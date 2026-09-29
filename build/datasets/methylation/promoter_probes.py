"""
Promoter probes of the Illumina methylation arrays (HM450, EPIC, EPICv2)
for the genes analysed by IntOGen.

A probe is a promoter probe of a gene if its CpG lies within a window
around the transcription start site (TSS) of the MANE transcript of the gene
(by default, from 1500 bp upstream to 500 bp downstream).
Genes whose MANE transcript is not present in the array annotation
fall back to the protein coding transcripts with the same gene name.
Probes flagged by the general quality mask (M_general) are discarded.

Uses the array annotations (hg38, GENCODE) from Zhou et al.
(https://zwdzwd.github.io/InfiniumAnnotation):
Zhou W, Laird PW, Shen H. Comprehensive characterization, annotation and
innovative use of Infinium DNA methylation BeadChip probes.
Nucleic Acids Res. 2017;45(4):e22.
"""

import csv
import gzip
from collections import defaultdict

import click


COLUMNS = ['PROBE', 'SYMBOL', 'CHROMOSOME', 'TRANSCRIPT', 'DIST_TO_TSS', 'CGI_POSITION', 'MATCH']


def load_genes(biomart):
    """MANE transcripts (cds_biomart.tsv) -> symbol and chromosome"""
    transcripts, chromosomes = {}, {}
    with open(biomart) as fd:
        for row in csv.reader(fd, delimiter='\t'):
            symbol, chromosome, transcript = row[1], row[3], row[10]
            transcripts[transcript] = symbol
            chromosomes[symbol] = chromosome
    return transcripts, chromosomes


def load_masked(masks):
    """Probes flagged by the general quality mask"""
    masked = set()
    for mask in masks:
        with gzip.open(mask, 'rt') as fd:
            for row in csv.DictReader(fd, delimiter='\t'):
                if row['M_general'].strip().upper() == 'TRUE':
                    masked.add(row['Probe_ID'])
    return masked


def promoter_probes(manifests, transcripts, symbols, masked, upstream, downstream):
    """
    Returns:
        dict symbol -> list of probe rows matched through the MANE transcript
        dict symbol -> list of probe rows matched through the gene name
    """
    by_transcript, by_name = defaultdict(dict), defaultdict(dict)
    for manifest in manifests:
        with gzip.open(manifest, 'rt') as fd:
            for row in csv.DictReader(fd, delimiter='\t'):
                probe = row['probeID']
                if probe in masked or row['distToTSS'] in ('', 'NA'):
                    continue
                cgi = row.get('CGIposition', 'NA') or 'NA'
                annotations = zip(row['geneNames'].split(';'), row['transcriptTypes'].split(';'),
                                  row['transcriptIDs'].split(';'), row['distToTSS'].split(';'))
                for name, type_, transcript, distance in annotations:
                    try:
                        distance = int(distance)
                    except ValueError:
                        continue
                    if not -upstream <= distance <= downstream:
                        continue
                    transcript = transcript.split('.')[0]
                    symbol = transcripts.get(transcript)
                    if symbol is not None:
                        target, match = by_transcript[symbol], 'MANE'
                    elif name in symbols and type_ == 'protein_coding':
                        symbol, target, match = name, by_name[name], 'gene_name'
                    else:
                        continue
                    previous = target.get(probe)
                    # keep the closest TSS for each probe and gene
                    if previous is None or abs(distance) < abs(previous[4]):
                        target[probe] = [probe, symbol, None, transcript, distance, cgi, match]
    return by_transcript, by_name


@click.command()
@click.option('--biomart', type=click.Path(exists=True), required=True, help='cds_biomart.tsv')
@click.option('--manifest', 'manifests', type=click.Path(exists=True), multiple=True, required=True,
              help='Array annotation with GENCODE transcripts (e.g. HM450.hg38.manifest.gencode.v41.tsv.gz)')
@click.option('--mask', 'masks', type=click.Path(exists=True), multiple=True, help='Array masks (e.g. HM450.hg38.mask.tsv.gz)')
@click.option('--upstream', type=int, default=1500, show_default=True, help='Promoter bp upstream of the TSS')
@click.option('--downstream', type=int, default=500, show_default=True, help='Promoter bp downstream of the TSS')
@click.option('-o', '--output', type=click.Path(), required=True)
def cli(biomart, manifests, masks, upstream, downstream, output):
    transcripts, chromosomes = load_genes(biomart)
    masked = load_masked(masks)
    by_transcript, by_name = promoter_probes(manifests, transcripts, set(chromosomes), masked, upstream, downstream)

    rows = []
    for symbol in chromosomes:
        # genes whose MANE transcript is not in the array annotation use the gene name
        probes = by_transcript.get(symbol) or by_name.get(symbol) or {}
        for row in probes.values():
            row[2] = chromosomes[symbol]
            rows.append(row)
    rows.sort(key=lambda r: (r[1], r[0]))

    with gzip.open(output, 'wt') as fd:
        writer = csv.writer(fd, delimiter='\t', lineterminator='\n')
        writer.writerow(COLUMNS)
        writer.writerows(rows)

    genes = len(set(r[1] for r in rows))
    by_mane = len(set(r[1] for r in rows if r[6] == 'MANE'))
    print(f'{len(rows)} promoter probes for {genes} of {len(chromosomes)} genes '
          f'({by_mane} through the MANE transcript). {len(masked)} masked probes discarded')


if __name__ == '__main__':
    cli()
