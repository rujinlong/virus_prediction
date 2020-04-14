#!/usr/bin/env python

import pandas as pd
import click

def combine_multiple_predictor(dvf, cat, virsorter, vibrant):
    df_dvf = pd.read_csv('predictor_dvf.tsv', sep='\t', usecols=['name', 'qvalue'])
    df_dvf.rename(columns={'name': 'contig_id', 'qvalue':'dvf_qvalue'}, inplace=True)

    df_cat = pd.read_csv(cat, sep='\t')
    df_cat = df_cat[df_cat.Superkingdom!='noname'][['contig_id', 'Superkingdom']].copy()
    df_cat.rename(columns={'Superkingdom':'CAT_Superkingdom'}, inplace=True)

    df_vs = pd.read_csv(virsorter, sep='\t', usecols=['ctg_id', 'category_virsorter'])
    df_vs.rename(columns={'ctg_id':'contig_id', 'category_virsorter':'virsorter_category'}, inplace=True)
    
    df_vib = pd.read_csv(vibrant, sep='\t', usecols=['contig_id', 'phage_type'])
    df_vib.rename(columns={'phage_type':'VIBRANT_phage_type'}, inplace=True)
    
    df_all = df_dvf.merge(df_cat, on='contig_id', how='outer').merge(df_vs, on='contig_id', how='outer').merge(df_vib, on='contig_id', how='outer')
    return df_all


@click.command()
@click.option("--dvf", '-d', help="predictor_dvf.tsv")
@click.option("--cat", '-c', help="predictor_CAT.tsv")
@click.option("--virsorter", '-s', help="predictor_virsorter.tsv")
@click.option("--vibrant", '-b', help="predictor_VIBRANT.tsv")
@click.option("--out", '-o', help="predictor_all.tsv")
def main(dvf, cat, virsorter, vibrant, out):
    df = combine_multiple_predictor(dvf, cat, virsorter, vibrant)
    df.to_csv(out, sep='\t', index=False)
    
    
if __name__ == '__main__':
    main()