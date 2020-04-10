#!/usr/bin/env python

import pandas as pd
import click

def merge_virsorter_CATBAT_dvf(fvirsorter, fcat_virus, fcat_unclassified, fdvf):
    """
    ### Most confidence
    - virsorter category: 1, 2, 4, 5 (original 1, 2)
    - CAT viral

    ### Less confidence
    Union any of two prediction results
    - dvf predicted
    - virsorter category: 3, 6 (original 3)
    - CAT unclassified
    """
    # read data
    clms = ['contig_id', 'category', 'virus_type', 'circular']
    df_vsorter = pd.read_csv(fvirsorter, sep='\t')
    df_vsorter.columns = clms
    lcp_vsorter = df_vsorter[df_vsorter.category.isin([3,6])]
    cat_virus = pd.read_csv(fcat_virus, sep='\t')
    cat_unclassified = pd.read_csv(fcat_unclassified)
    dvf = pd.read_csv(fdvf, sep='\t')
    
    # lcp virus
    ctg_lcp_vsorter = set(lcp_vsorter.contig_id.tolist())
    ctg_cat_unclassified = set(cat_unclassified.contig_id.tolist())
    ctg_dvf = set(dvf.name.tolist())

    c1 = ctg_lcp_vsorter.intersection(ctg_dvf)
    c2 = ctg_lcp_vsorter.intersection(ctg_cat_unclassified)
    c3 = ctg_dvf.intersection(ctg_cat_unclassified)
    lcp_virus = c1.union(c2).union(c3)
    
    # merge
    df_mcp_cat = pd.DataFrame([[x, 7, 'CAT_virus', 7] for x in cat_virus.contig_id], columns=clms)
    df_lcp_virus = pd.DataFrame([[x, 8, 'lcp_virus', 8] for x in lcp_virus], columns=clms)
    df = pd.concat([df_vsorter, df_mcp_cat, df_lcp_virus], axis=0)
    df.sort_values(['contig_id', 'category'], inplace=True)
    
    df_dup = df[df.duplicated('contig_id', keep=False)].copy().drop_duplicates('contig_id')
    df_no36 = df.drop(df.category.isin([3,6])).copy()
    df_new = pd.concat([df_no36, df_dup], axis=0)
    df_final = df_new.sort_values('category', ascending=True).drop_duplicates('contig_id')
    return df_final


@click.command()
@click.option("--fvirsorter", '-s', help="out_virsorter.tsv")
@click.option("--fcat_virus", '-c', help="CAT_virus.tsv")
@click.option("--fcat_unclassified", '-u', help="CAT_unclassified.tsv")
@click.option("--fdvf", '-d', help="out_dvf.tsv")
@click.option("--out", '-o', help="out.tsv")
def main(fvirsorter, fcat_virus, fcat_unclassified, fdvf, out):
    df = merge_virsorter_CATBAT_dvf(fvirsorter, fcat_virus, fcat_unclassified, fdvf)
    df.to_csv(out, sep='\t', index=False)
    

if __name__ == '__main__':
    main()