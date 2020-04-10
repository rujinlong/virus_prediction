#!/usr/bin/env python

import pandas as pd
import click

def virsorter_to_df(fname):
    with open(fname, 'r') as fh:
        f = fh.read().splitlines()
    
    category_virsorter = 1
    phage_type = 'phage'
    entry_virsorter = []
    for row in f:
        if row.startswith('## 1'):
            category_virsorter = 1
            phage_type = 'phage'
        elif row.startswith('## 2'):
            category_virsorter = 2
            phage_type = 'phage'
        elif row.startswith('## 3'):
            category_virsorter = 3
            phage_type = 'phage'
        elif row.startswith('## 4'):
            category_virsorter = 4
            phage_type = 'prophage'
        elif row.startswith('## 5'):
            category_virsorter = 5
            phage_type = 'prophage'
        elif row.startswith('## 6'):
            category_virsorter = 6
            phage_type = 'prophage'
        elif row.startswith('VIRSorter_'):
            new_row = '{},{},{}'.format(row,category_virsorter,phage_type)
            entry_virsorter.append(new_row.split(','))
    
    clms = ['Contig_id',
            'Nb_genes_contigs',
            'Fragment',
            'Nb_genes',
            'raw_category',
            'Nb_phage_hallmark_genes',
            'Phage_gene_enrichment_sig',
            'Non-Caudovirales_phage_gene_enrichment_sig',
            'Pfam_depletion_sig',
            'Uncharacterized_enrichment_sig',
            'Strand_switch_depletion_sig',
            'Short_genes_enrichment_sig',
            'category_virsorter',
            'phage_type_virsorter']
    return pd.DataFrame(entry_virsorter, columns=clms)

@click.command()
@click.option("--fn_virsorter", '-i', help="virsorter output")
@click.option("--ctgid_mapping", '-d', help="virsorter contig id mapping")
@click.option("--out", '-o', help="virsorter.tsv")
def pretty_virsorter(fn_virsorter, ctgid_mapping, out):
    df = virsorter_to_df(fn_virsorter)
    
    def _detect_circular(x):
        if x.split('-')[-1] == 'circular':
            circular = 1
        else:
            circular = 0
        return circular
    
    def _modify_contig_id(x):
        if x['circular_virsorter'] == 1:
            vctg_id = x['Contig_id'].replace('-circular', '')
        else:
            vctg_id = x['Contig_id']
        return vctg_id
    
    df['circular_virsorter'] = df.apply(lambda x: _detect_circular(x['Contig_id']), axis=1)
    df['virsorter_ctgid'] = df.apply(lambda x: _modify_contig_id(x), axis=1)
 
    idmap = pd.read_csv(ctgid_mapping, sep='\t', names=['ctg_id', 'virsorter_ctgid'])
    final = df.merge(idmap, on='virsorter_ctgid', how='left')
    df_final = final[['ctg_id', 'category_virsorter', 'phage_type_virsorter', 'circular_virsorter']].copy()
    df_final.drop_duplicates(inplace=True)
    df_final.sort_values(['ctg_id','category_virsorter'], inplace=True, ascending=True)
    df_final.drop_duplicates('ctg_id', keep='first', inplace=True)
    df_final.to_csv(out, index=False, sep='\t')
    
if __name__ == '__main__':
    pretty_virsorter()