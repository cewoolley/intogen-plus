"""
Synthetic cohorts for the pathway analyses.

PathwayCohort (mutations and methylation):

- LONGTAIL: 40 genes with a modest excess of missense mutations each
  (selected through its long tail; no gene is individually significant)
- DRIVER_SET: contains DRV1, a strong driver, and neutral genes
- SILENCED: 30 genes with a mild excess of promoter hypermethylation each;
  their silencing co-occurs with DRV1 mutations
- NEUTRAL_*: random neutral gene sets

SequencingCohort (mutations only), with a selection network made of DRV2 and
the long tails of two gene sets:

- DRV1 and DRV2: drivers
- LONGTAIL: excess of missense mutations in random tumours (selected, not in the network)
- LOF_TAIL: excess of truncating mutations in the tumours with DRV2 mutations
- MIS_TAIL: excess of missense mutations in the tumours with DRV2 mutations
- DRIVER_SET: DRV1 and neutral genes
- NEUTRAL_*: random neutral gene sets
"""

import os

import numpy as np
import pandas as pd


def gene_names(n=600):
    return ['DRV1'] + [f'G{i}' for i in range(1, n)]


class PathwayCohort:

    def __init__(self, seed=0, n_genes=600, n_tumours=200):
        rng = self.rng = np.random.default_rng(seed)
        self.genes = np.array(gene_names(n_genes))
        self.tumours = np.array([f'T{i:03d}' for i in range(n_tumours)])
        self.sets = {
            'LONGTAIL': list(self.genes[1:41]),
            'DRIVER_SET': list(self.genes[[0] + list(range(41, 70))]),
            'SILENCED': list(self.genes[70:100]),
        }
        for i in range(20):
            self.sets[f'NEUTRAL_{i}'] = list(rng.choice(self.genes[100:], rng.integers(15, 40), replace=False))
        self.burden = rng.lognormal(0, 0.6, n_tumours)
        self.drv1_tumours = rng.choice(n_tumours, 60, replace=False)

    def mutations(self):
        """genemuts table and mutations per tumour (VEP-like)"""
        rng = self.rng
        n = len(self.genes)
        length = rng.lognormal(7.2, 0.5, n)
        pred = rng.lognormal(0, 0.3, n)
        rate = rng.gamma(8, pred / 8)
        base = 2.5e-3
        exp = {'syn': 0.25, 'mis': 0.65, 'non': 0.06, 'spl': 0.04}
        omega = np.ones(n)
        omega[1:41] = 2.2    # long tail
        counts = {k: rng.poisson(length * base * f * rate * (omega if k == 'mis' else 1)) for k, f in exp.items()}
        counts['mis'][0] += 45    # DRV1
        counts['non'][0] += 15
        genemuts = pd.DataFrame({'gene_name': self.genes, 'n_syn': counts['syn'], 'n_mis': counts['mis'],
                                 'n_non': counts['non'], 'n_spl': counts['spl']})
        for k, f in exp.items():
            genemuts[f'exp_{k}'] = length * base * f
        genemuts['exp_syn_cv'] = genemuts['exp_syn'] * pred

        consequence = {'syn': 'synonymous_variant', 'mis': 'missense_variant', 'non': 'stop_gained',
                       'spl': 'splice_donor_variant'}
        weights = self.burden / self.burden.sum()
        rows = []
        for gi, gene in enumerate(self.genes):
            for k in exp:
                for _ in range(counts[k][gi]):
                    if gi == 0 and k in ('mis', 'non'):
                        tumour = rng.choice(self.drv1_tumours)
                    else:
                        tumour = rng.choice(len(self.tumours), p=weights)
                    i = len(rows)
                    rows.append([f'I{i:010d}__{self.tumours[tumour]}__C__T__{i}', f'1:{i}', 'T', 'ENSG', 'ENST',
                                 'Transcript', consequence[k], gene, 'YES', 'NM'])
        vep = pd.DataFrame(rows, columns=['#Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature',
                                          'Feature_type', 'Consequence', 'SYMBOL', 'CANONICAL', 'MANE_SELECT'])
        return genemuts, vep

    def methylation(self):
        """methylation-analysis results, events and tumour list"""
        rng = self.rng
        n, t = len(self.genes), len(self.tumours)
        pi = np.clip(0.004 * self.burden, 0, 0.2)
        expected = np.full(n, pi.sum())
        excess = np.zeros(n)
        excess[70:100] = 2.5          # excess in the SILENCED genes (each gene still rarely altered)
        observed = rng.poisson(expected * np.where(excess > 0, excess, 1))
        events = []
        drv1 = set(self.drv1_tumours)
        for gi, gene in enumerate(self.genes):
            if 70 <= gi < 100:
                # the extra events happen in tumours with DRV1 mutations
                background = rng.choice(t, rng.poisson(expected[gi]), replace=False, p=pi / pi.sum())
                extra = rng.choice(self.drv1_tumours, max(0, observed[gi] - len(background)), replace=True)
                chosen = set(background) | set(extra)
            else:
                chosen = set(rng.choice(t, min(observed[gi], t), replace=False, p=pi / pi.sum()))
            observed[gi] = len(chosen)
            events += [(gene, self.tumours[s]) for s in chosen]
        results = pd.DataFrame({'SYMBOL': self.genes, 'STATUS': 'tested', 'HYPER_SAMPLES': observed,
                                'EXPECTED_HYPER': expected, 'FUNCTIONAL': 'True', 'Q_VALUE': 0.5})
        return results, pd.DataFrame(events, columns=['SYMBOL', 'SAMPLE']), list(self.tumours)

    def write(self, folder):
        os.makedirs(folder, exist_ok=True)
        p = lambda name: os.path.join(folder, name)
        genemuts, vep = self.mutations()
        genemuts.to_csv(p('C.dndscv_genemuts.tsv.gz'), sep='\t', index=False)
        vep.to_csv(p('C.tsv.gz'), sep='\t', index=False)
        pd.DataFrame({'SYMBOL': ['DRV1', 'G500'], 'TIER': [1, 4], 'FILTER': ['PASS', 'No driver'],
                      'QVALUE_COMBINATION': [1e-10, 0.4]}).to_csv(p('C.vet.tsv'), sep='\t', index=False)
        results, events, tumours = self.methylation()
        results.to_csv(p('C.methylation.tsv.gz'), sep='\t', index=False)
        events.to_csv(p('C.methylation_events.tsv.gz'), sep='\t', index=False)
        pd.DataFrame(columns=['SYMBOL'] + tumours).to_csv(p('C.promoter_methylation.tsv.gz'), sep='\t', index=False)
        rows = [(s, 'test', s, g) for s, genes in self.sets.items() for g in genes]
        pd.DataFrame(rows, columns=['SET', 'SOURCE', 'NAME', 'SYMBOL']).to_csv(p('gene_sets.tsv.gz'), sep='\t',
                                                                               index=False)
        return {
            'genemuts': p('C.dndscv_genemuts.tsv.gz'), 'mutations': p('C.tsv.gz'), 'vet': p('C.vet.tsv'),
            'gene_sets': p('gene_sets.tsv.gz'), 'methylation': p('C.methylation.tsv.gz'),
            'methylation_events': p('C.methylation_events.tsv.gz'),
            'methylation_matrix': p('C.promoter_methylation.tsv.gz'),
        }


class SequencingCohort:

    CONSEQUENCES = {'syn': 'synonymous_variant', 'mis': 'missense_variant', 'non': 'stop_gained',
                    'spl': 'splice_donor_variant'}
    FRACTIONS = {'syn': 0.25, 'mis': 0.65, 'non': 0.06, 'spl': 0.04}

    def __init__(self, seed=0, n_genes=800, n_tumours=250, lof_omega=4.0, mis_omega=2.2, network_mis_omega=2.5):
        rng = self.rng = np.random.default_rng(seed)
        self.genes = np.array(['DRV1', 'DRV2'] + [f'G{i}' for i in range(2, n_genes)])
        self.tumours = np.array([f'T{i:03d}' for i in range(n_tumours)])
        self.sets = {
            'LONGTAIL': list(self.genes[2:42]),
            'LOF_TAIL': list(self.genes[42:82]),
            'MIS_TAIL': list(self.genes[82:122]),
            'DRIVER_SET': list(self.genes[[0] + list(range(122, 151))]),
        }
        for i in range(30):
            self.sets[f'NEUTRAL_{i}'] = list(rng.choice(self.genes[160:], rng.integers(15, 60), replace=False))
        self.burden = rng.lognormal(0, 0.8, n_tumours)
        self.drv1 = rng.choice(n_tumours, 60, replace=False)
        self.drv2 = rng.choice(n_tumours, 50, replace=False)
        # excess of mutations: (genes, consequences, omega, tumours where they happen or None for any)
        self.excess = [(range(2, 42), ['mis'], mis_omega, None),
                       (range(42, 82), ['non', 'spl'], lof_omega, self.drv2),
                       (range(82, 122), ['mis'], network_mis_omega, self.drv2)]

    def mutations(self):
        rng = self.rng
        n = len(self.genes)
        length = rng.lognormal(7.2, 0.5, n)
        pred = rng.lognormal(0, 0.3, n)
        rate = rng.gamma(8, pred / 8)
        weights = self.burden / self.burden.sum()
        counts = {k: np.zeros(n, dtype=int) for k in self.FRACTIONS}
        rows = []

        def add(gi, k, tumour):
            i = len(rows)
            rows.append([f'I{i:010d}__{self.tumours[tumour]}__C__T__{i}', f'1:{i}', 'T', 'ENSG', 'ENST',
                         'Transcript', self.CONSEQUENCES[k], self.genes[gi], 'YES', 'NM'])
            counts[k][gi] += 1

        for gi in range(n):
            for k, f in self.FRACTIONS.items():
                mean = length[gi] * 2.5e-3 * f * rate[gi]
                for _ in range(rng.poisson(mean)):
                    add(gi, k, rng.choice(len(self.tumours), p=weights))
                for genes, consequences, omega, tumours in self.excess:
                    if gi in genes and k in consequences:
                        for _ in range(rng.poisson(mean * (omega - 1))):
                            add(gi, k, rng.choice(len(self.tumours), p=weights) if tumours is None
                                else rng.choice(tumours))
        for t in self.drv1:
            add(0, 'mis' if rng.uniform() < 0.8 else 'non', t)
        for t in self.drv2:
            add(1, 'non' if rng.uniform() < 0.6 else 'mis', t)

        genemuts = pd.DataFrame({'gene_name': self.genes, 'n_syn': counts['syn'], 'n_mis': counts['mis'],
                                 'n_non': counts['non'], 'n_spl': counts['spl']})
        for k, f in self.FRACTIONS.items():
            genemuts[f'exp_{k}'] = length * 2.5e-3 * f
        genemuts['exp_syn_cv'] = genemuts['exp_syn'] * pred
        vep = pd.DataFrame(rows, columns=['#Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature',
                                          'Feature_type', 'Consequence', 'SYMBOL', 'CANONICAL', 'MANE_SELECT'])
        return genemuts, vep

    def write(self, folder):
        os.makedirs(folder, exist_ok=True)
        p = lambda name: os.path.join(folder, name)
        genemuts, vep = self.mutations()
        genemuts.to_csv(p('S.dndscv_genemuts.tsv.gz'), sep='\t', index=False)
        vep.to_csv(p('S.tsv.gz'), sep='\t', index=False)
        pd.DataFrame({'SYMBOL': ['DRV1', 'DRV2'], 'TIER': [1, 1], 'FILTER': ['PASS', 'PASS'],
                      'QVALUE_COMBINATION': [1e-12, 1e-10]}).to_csv(p('S.vet.tsv'), sep='\t', index=False)
        rows = [(s, 'test', s, g) for s, genes in self.sets.items() for g in genes]
        pd.DataFrame(rows, columns=['SET', 'SOURCE', 'NAME', 'SYMBOL']).to_csv(p('gene_sets.tsv.gz'), sep='\t',
                                                                               index=False)
        return {'genemuts': p('S.dndscv_genemuts.tsv.gz'), 'mutations': p('S.tsv.gz'), 'vet': p('S.vet.tsv'),
                'gene_sets': p('gene_sets.tsv.gz')}
