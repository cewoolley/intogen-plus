"""
Gene sets used by the pathway analyses.

- Reactome pathways (CC0): https://reactome.org
- MSigDB hallmark gene sets (CC BY 4.0): https://www.gsea-msigdb.org

Gene symbols are harmonised with the MANE genes analysed by IntOGen
(outdated symbols are updated; other genes are discarded).
Output: TSV with the columns SET, SOURCE, NAME and SYMBOL.
"""

import csv
import gzip
import io
import json
import zipfile

import click


def read_gmt(path):
    """Yield (name, description, genes) from a GMT file (optionally zipped)"""
    if path.endswith('.zip'):
        with zipfile.ZipFile(path) as z:
            member = [n for n in z.namelist() if n.endswith('.gmt')][0]
            fd = io.TextIOWrapper(z.open(member), encoding='utf-8')
            lines = fd.readlines()
    else:
        with open(path) as fd:
            lines = fd.readlines()
    for line in lines:
        fields = line.rstrip('\n').split('\t')
        if len(fields) >= 3:
            yield fields[0].strip(), fields[1].strip(), [g.strip() for g in fields[2:] if g.strip()]


@click.command()
@click.option('--reactome', type=click.Path(exists=True), help='ReactomePathways.gmt(.zip)')
@click.option('--hallmarks', type=click.Path(exists=True), help='MSigDB hallmarks GMT (symbols)')
@click.option('--biomart', type=click.Path(exists=True), required=True, help='cds_biomart.tsv')
@click.option('--symbols', type=click.Path(exists=True), help='Outdated to current HUGO symbols (JSON)')
@click.option('-o', '--output', type=click.Path(), required=True)
def cli(reactome, hallmarks, biomart, symbols, output):
    with open(biomart) as fd:
        genes = {row[1] for row in csv.reader(fd, delimiter='\t')}
    mapping = json.load(open(symbols)) if symbols else {}

    def harmonise(symbol):
        return symbol if symbol in genes else mapping.get(symbol) if mapping.get(symbol) in genes else None

    sets = []
    if reactome:
        sets += [(identifier, 'Reactome', name, members) for name, identifier, members in read_gmt(reactome)]
    if hallmarks:
        sets += [(name, 'MSigDB_Hallmark', name, members) for name, _, members in read_gmt(hallmarks)]

    rows, total, kept = [], 0, 0
    for identifier, source, name, members in sets:
        mapped = sorted({s for s in map(harmonise, members) if s is not None})
        total += len(set(members))
        kept += len(mapped)
        rows += [(identifier, source, name, s) for s in mapped]

    with gzip.open(output, 'wt') as fd:
        writer = csv.writer(fd, delimiter='\t', lineterminator='\n')
        writer.writerow(['SET', 'SOURCE', 'NAME', 'SYMBOL'])
        writer.writerows(rows)
    print(f'{len(sets)} gene sets; {kept} of {total} gene memberships mapped to MANE genes')


if __name__ == '__main__':
    cli()
