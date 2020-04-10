#!/usr/bin/env python

import pandas as pd
import click
import re


def get_unclassified_contigs(fin):
    """
    input:
    -----
    sample.contig2classification_official.txt
    
    output:
    ------
    df:list of unclassified contigs
    """
    
    df = pd.read_csv(fin, sep='\t')
    df.rename(columns={'# contig':'contig_id'}, inplace=True)
    df_unclassified = df[df.classification=='unclassified']
    return df_unclassified[['contig_id']]


def get_classified_virus_lineage(fin):
    """
    input:
    -----
    sample.contig2classification_official.txt
    
    output:
    ------
    df: classified viral contigs
    col1: contig_id
    col2: lineage (lineage or unclassified)
    """
    
    df = pd.read_csv(fin, sep='\t')
    new_col = {
        '# contig': 'contig_id',
        'superkingdom': 'Superkingdom',
        'phylum': 'Phylum',
        'class': 'Class',
        'order': 'Order',
        'family': 'Family',
        'genus': 'Genus',
        'species': 'Species'
    }
    df.rename(columns=new_col, inplace=True)
    
    # get classified viruses
    df_classified = df[(df.classification=='classified') & ~(df['Superkingdom'].isna())]
    df_classified_virus = df_classified[df_classified['Superkingdom'].str.startswith('Viruses')].copy()
    df_classified_virus.replace(regex=r'\:.*$', value='', inplace=True)
    df_classified_virus.fillna('noname', inplace=True)
    df_classified_virus.replace('not classified', value='noname', inplace=True)
    # df_classified_virus['lineage'] = df_classified_virus.apply(lambda x:'k__{};p__{};c__{};o__{};f__{};g__{};s__{}'.format(x['superkingdom'], x['phylum'], x['class'], x['order'], x['family'], x['genus'], x['species']), axis=1)

    # return df_classified_virus[['contig_id', 'lineage']]
    return df_classified_virus[['contig_id', 'Superkingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']]


@click.command()
@click.option("--fin", '-i', help="sample.contig2classification_official.txt")
@click.option("--method", '-m', help="[virus | unclassified]")
@click.option("--out", '-o', help="out.tsv")
def main(fin, method, out):
    """
    method: [virus | unclassified]
    """
    if method == 'virus':
        df = get_classified_virus_lineage(fin)
    elif method == 'unclassified':
        df = get_unclassified_contigs(fin)
    df.to_csv(out, sep='\t', index=False)
    

if __name__ == '__main__':
    main()