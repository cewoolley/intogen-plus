"""
Synthetic cohorts with mutations, DNA methylation and RNA-seq data.

Planted signals (tumour samples T000...):

- SIL1: promoter hypermethylation in 30% of the tumours with loss of expression
  (epigenetic silencing). SIL1 is also mutated (truncating) in other tumours.
- HYP1: promoter hypermethylation in 30% of the tumours without expression changes
- MET1: promoter methylated in normal samples (not testable)
- OUT1: over-expression outliers in 20% of the tumours; mutated tumours overexpress it
- LOW1: under-expression outliers in 15% of the tumours
- NOEXP: not expressed
- TRUNC1: truncating mutations with lower expression (nonsense-mediated decay)
- XGENE / YGENE: genes in the sex chromosomes
- G0...: background genes
"""

import os

import numpy as np
import pandas as pd


PLANTED = ['SIL1', 'HYP1', 'MET1', 'OUT1', 'LOW1', 'NOEXP', 'TRUNC1', 'XGENE', 'YGENE']


def genes(n_background=200):
    return [f'G{i}' for i in range(n_background)] + PLANTED


def chromosome(gene):
    return {'XGENE': 'X', 'YGENE': 'Y'}.get(gene, '1')


def write_gene_annotation(path, gene_list):
    """cds_biomart.tsv-like annotation (MANE protein coding genes)"""
    rows = [[f'ENSG{i:011d}', g, f'ENSP{i:011d}', chromosome(g), 100, 200, 1, 100, 300, 1, f'ENST{i:011d}', 90, 210]
            for i, g in enumerate(gene_list)]
    pd.DataFrame(rows).to_csv(path, sep='\t', header=False, index=False)


class Cohort:

    def __init__(self, seed=0, n_tumors=120, n_normals=8, n_background=200):
        self.rng = np.random.default_rng(seed)
        self.genes = genes(n_background)
        self.gi = {g: i for i, g in enumerate(self.genes)}
        self.tumors = [f'T{i:03d}' for i in range(n_tumors)]
        self.normals = [f'N{i:03d}' for i in range(n_normals)]
        rng = self.rng
        n = n_tumors
        pick = lambda k: np.sort(rng.choice(n, k, replace=False))
        self.silenced = pick(int(0.3 * n))
        self.hyper_only = pick(int(0.3 * n))
        self.over = pick(int(0.2 * n))
        self.under = pick(int(0.15 * n))
        free = np.setdiff1d(np.arange(n), self.silenced)
        self.sil_mutated = np.sort(rng.choice(free, 6, replace=False))
        self.out_mutated = np.sort(rng.choice(self.over, 5, replace=False))
        self.trunc_mutated = pick(12)

    # --- methylation -----------------------------------------------------
    def beta(self):
        rng, gi = self.rng, self.gi
        n_samples = len(self.tumors) + len(self.normals)
        beta = rng.beta(2, 30, size=(len(self.genes), n_samples))
        # sporadic background hypermethylation (tumours only)
        n = len(self.tumors)
        hits = rng.uniform(size=(len(self.genes), n)) < 0.02
        beta[:, :n][hits] = rng.uniform(0.45, 0.9, hits.sum())
        beta[gi['SIL1'], self.silenced] = rng.uniform(0.5, 0.9, len(self.silenced))
        beta[gi['HYP1'], self.hyper_only] = rng.uniform(0.5, 0.9, len(self.hyper_only))
        beta[gi['MET1'], :] = rng.uniform(0.7, 0.9, n_samples)
        # X inactivation: females have ~0.5 promoter methylation
        female = rng.uniform(size=n_samples) < 0.5
        beta[gi['XGENE'], female] = rng.uniform(0.4, 0.6, female.sum())
        return beta

    def write_methylation(self, folder, name, probes_per_gene=3, values='beta'):
        """Probe-level matrix (with normals), promoter probes annotation and sample sheet"""
        beta = self.beta()
        rows, ids, symbols = [], [], []
        for g, i in self.gi.items():
            for k in range(probes_per_gene):
                rows.append(np.clip(beta[i] + self.rng.normal(0, 0.02, beta.shape[1]), 0.001, 0.999))
                ids.append(f'cg{i * 10 + k:08d}')
                symbols.append(g)
        m = pd.DataFrame(rows, index=ids, columns=self.tumors + self.normals)
        if values == 'm':
            m = np.log2(m / (1 - m))
        m.index.name = 'ID_REF'
        # probes not in promoters
        m.loc['cg99999999'] = 0.5
        m.iloc[5, 3] = np.nan
        paths = {
            'matrix': os.path.join(folder, f'{name}.beta.tsv.gz'),
            'probes': os.path.join(folder, 'promoter_probes.tsv.gz'),
        }
        m.to_csv(paths['matrix'], sep='\t')
        pd.DataFrame({'PROBE': ids, 'SYMBOL': symbols, 'CHROMOSOME': [chromosome(g) for g in symbols]}).to_csv(
            paths['probes'], sep='\t', index=False)
        return paths

    # --- expression ------------------------------------------------------
    def counts(self):
        rng, gi = self.rng, self.gi
        n = len(self.tumors)
        mu = rng.lognormal(4, 1.2, len(self.genes))
        lam = np.outer(mu, np.ones(n)) * rng.lognormal(0, 0.25, (len(self.genes), n))
        lam[gi['SIL1']] = 300 * rng.lognormal(0, 0.2, n)
        lam[gi['SIL1'], self.silenced] = 2
        lam[gi['OUT1']] = rng.lognormal(0, 0.3, n)
        lam[gi['OUT1'], self.over] = 500
        lam[gi['LOW1']] = 500 * rng.lognormal(0, 0.2, n)
        lam[gi['LOW1'], self.under] = 5
        lam[gi['NOEXP']] = 0.001
        lam[gi['TRUNC1']] = 400 * rng.lognormal(0, 0.2, n)
        lam[gi['TRUNC1'], self.trunc_mutated] = 60
        male = rng.uniform(size=n) < 0.5
        lam[gi['YGENE']] = np.where(male, 200, 0.001)
        # make the library size realistic (~20M reads for ~20,000 genes)
        return rng.poisson(lam * 1000)

    def write_expression(self, folder, name, ensembl=False):
        counts = self.counts()
        index = [f'ENSG{i:011d}.3' for i in range(len(self.genes))] if ensembl else self.genes
        df = pd.DataFrame(counts, index=index, columns=self.tumors)
        df.index.name = 'gene_id'
        path = os.path.join(folder, f'{name}.counts.tsv.gz')
        df.to_csv(path, sep='\t')
        return path

    # --- mutations -------------------------------------------------------
    def write_vep(self, path):
        """Mutations as in the processed VEP output (parse-vep)"""
        rng = self.rng
        rows = []

        def add(sample, gene, consequence):
            i = len(rows)
            rows.append([f'I{i:010d}__{sample}__C__T__{1000 + i}', f'1:{1000 + i}', 'T', 'ENSG', 'ENST',
                         'Transcript', consequence, gene, 'YES', 'NM_1'])

        for j, sample in enumerate(self.tumors):
            for g in rng.choice(self.genes[:200], 3, replace=False):
                add(sample, g, rng.choice(['missense_variant', 'synonymous_variant']))
        for j in self.sil_mutated:
            add(self.tumors[j], 'SIL1', 'stop_gained')
        for j in self.out_mutated:
            add(self.tumors[j], 'OUT1', 'missense_variant')
        for j in self.trunc_mutated:
            add(self.tumors[j], 'TRUNC1', 'frameshift_variant')
        df = pd.DataFrame(rows, columns=['#Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature',
                                         'Feature_type', 'Consequence', 'SYMBOL', 'CANONICAL', 'MANE_SELECT'])
        df.to_csv(path, sep='\t', index=False, compression='gzip' if path.endswith('.gz') else None)
        return path


def write_samplesheet(path, cohort):
    rows = [[n, n, 'normal'] for n in cohort.normals]
    pd.DataFrame(rows, columns=['ID', 'SAMPLE', 'TYPE']).to_csv(path, sep='\t', index=False)
    return path
