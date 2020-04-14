#!/usr/bin/env python

import pandas as pd
import click

def combine_taxa2(fcat, fdemovir, taxa2lineage):
    cat = pd.read_csv(fcat, sep='\t', names=['contig_id', 'taxid'])
    demovir = pd.read_csv(fdemovir, sep='\t', names=['contig_id', 'taxid'])
    df = pd.concat([cat, demovir], axis=0).drop_duplicates('contig_id')
    taxa2lineage = pd.read_csv(taxa2lineage, sep='\t')
    return df.merge(taxa2lineage, on='taxid', how='left').fillna('no_rank')


@click.command()
@click.option("--fcat", '-c', help='CAT_contig2taxid.tsv')
@click.option("--fdemovir", '-d', help='demovir_contig2taxid.tsv')
@click.option("--taxa2lineage", '-t', help='taxa2lineage.tsv')
@click.option("--assembler", '-a', help='[spades, megahit]')
@click.option("--out", '-o', help="taxa.tsv")
def main(fcat, fdemovir, taxa2lineage, assembler, out):
    df = combine_taxa2(fcat, fdemovir, taxa2lineage)
    
    if assembler == 'spades':
        df['length'] = df.apply(lambda x:int(x['contig_id'].split('_')[4]), axis=1)
    elif assembler == 'megahit':
        df['length'] = df.apply(lambda x:int(x['contig_id'].split('=')[3]), axis=1)
    else:
        print('assember should be spades or megahit!')
    df.to_csv(out, sep='\t', index=False)


if __name__ == '__main__':
    main()