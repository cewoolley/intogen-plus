"""Reactome and hallmark gene sets plus curated CRC pathways (canonical pathways of TCGA Network, Nature 2012)"""
import pandas as pd

CURATED = {
    'CRC_WNT': ['APC', 'CTNNB1', 'TCF7L2', 'SOX9', 'FBXW7', 'AMER1', 'AXIN1', 'AXIN2', 'RNF43', 'ZNRF3', 'DKK1', 'DKK2',
                'DKK3', 'DKK4', 'LRP5', 'FZD10', 'LEF1', 'TCF7', 'BCL9', 'BCL9L'],
    'CRC_TGFB_BMP': ['TGFBR1', 'TGFBR2', 'ACVR1B', 'ACVR2A', 'SMAD2', 'SMAD3', 'SMAD4', 'BMPR1A', 'BMPR2', 'ACVR1', 'ACVR2B'],
    'CRC_PI3K': ['PIK3CA', 'PIK3R1', 'PIK3CB', 'PTEN', 'IGF2', 'IRS2', 'AKT1', 'AKT2', 'AKT3', 'MTOR', 'TSC1', 'TSC2'],
    'CRC_RTK_RAS': ['KRAS', 'NRAS', 'HRAS', 'BRAF', 'ARAF', 'RAF1', 'MAP2K1', 'ERBB2', 'ERBB3', 'EGFR', 'MET', 'NF1',
                    'PTPN11', 'SOS1'],
    'CRC_TP53_DDR': ['TP53', 'ATM', 'ATR', 'CHEK2', 'MDM2', 'MDM4', 'CDKN2A', 'BRCA1', 'BRCA2'],
    'CRC_IMMUNE_ESCAPE': ['B2M', 'HLA-A', 'HLA-B', 'HLA-C', 'TAP1', 'TAP2', 'TAPBP', 'NLRC5', 'JAK1', 'JAK2', 'IRF1', 'STAT1',
                          'IFNGR1', 'IFNGR2', 'CIITA', 'CD58', 'CASP8'],
    'CRC_SWI_SNF': ['ARID1A', 'ARID1B', 'ARID2', 'SMARCA4', 'SMARCA2', 'SMARCB1', 'PBRM1', 'SMARCC1', 'SMARCC2', 'SMARCD1'],
}
base = pd.read_csv("data/gene_sets.tsv", sep="\t")
rows = [(s, 'CRC_curated', s, g) for s, genes in CURATED.items() for g in genes]
pd.concat([base, pd.DataFrame(rows, columns=['SET', 'SOURCE', 'NAME', 'SYMBOL'])]).to_csv(
    'gene_sets_crc.tsv.gz', sep='\t', index=False)
