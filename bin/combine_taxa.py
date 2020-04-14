#!/usr/bin/env python

import pandas as pd
import click

def read_taxa(fin, contig_id='contig_id'):
    df = pd.read_csv(fin, sep='\t')
    df.rename(columns={contig_id: 'contig_id'}, inplace=True)
    taxas = ['Superkingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    add_taxa = [x for x in taxas if x not in df.columns]
    for taxa in add_taxa:
        if taxa == 'Superkingdom':
            df[taxa] = 'Viruses'
        else:
            df[taxa] = 'noname'
    df.replace('Unassigned', 'noname', inplace=True)
    df.replace('noname', 'no_rank', inplace=True)
    df2 = df[['contig_id'] + taxas].copy()
#     df2['taxa_level'] = df2.apply(lambda x:list(x).index([i for i in x if i!='no_rank'][-1]), axis=1)
    return df2


@click.command()
@click.option("--fvcontact", '-v', default=False, help="vc_cluster.tsv")
@click.option("--fcat", '-c', help="cat_virus.tsv")
@click.option("--fdemovir", '-d', help="DemoVir_assignment.txt")
@click.option("--out", '-o', help="taxa.tsv")
def main(fcat, fdemovir, out, fvcontact=False):
    cat = read_taxa(fcat)
    demovir = read_taxa(fdemovir, contig_id='Sequence_ID')
    
    cat['taxa_priority'] = 2
    demovir['taxa_priority'] = 3
    
    if fvcontact:
        vcc = read_taxa(fvcontact)
        vcc['taxa_priority'] = 1
        taxa_tools = [vcc, cat, demovir]
    else:
        taxa_tools = [cat, demovir]

    df = pd.concat(taxa_tools, axis=0).sort_values('taxa_priority').drop_duplicates('contig_id').drop(['Phylum', 'Class', 'taxa_priority'], axis=1)
    df.replace('*', '', inplace=True)
    df['Species'] = df['Species'].str.replace('*', '')
    df['lowest_taxa'] = df.apply(lambda x:[i for i in x if i!='no_rank'][-1], axis=1)
    df.to_csv(out, sep='\t', index=False)

    
if __name__ == '__main__':
    main()
