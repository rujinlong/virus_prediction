#!/usr/bin/env python

import pandas as pd
import click

@click.command()
@click.option("--virus_list", '-l', help="virus.list")
@click.option("--predicted_taxa", '-p', help="taxa_pred.tsv")
@click.option("--name2taxid", '-n', help="name2taxid.tsv")
@click.option("--assembler", '-a', help="spades or megahit")
@click.option("--taxid2lineage", '-t', help="taxid2lineage.tsv")
@click.option("--out", '-o', help="taxaonomy.tsv")
def main(virus_list, predicted_taxa, name2taxid, taxid2lineage, assembler, out):
    df_raw_taxa = pd.read_csv(predicted_taxa, sep='\t')
    df_name2taxid = pd.read_csv(name2taxid, sep='\t', names=['lowest_taxa', 'taxid'], dtype=str)
    df_taxid2lineage = pd.read_csv(taxid2lineage, sep='\t', dtype=str)
    df_name2lineage = df_name2taxid.merge(df_taxid2lineage, on='taxid', how='left')
    df_new_taxa = df_raw_taxa[['contig_id', 'lowest_taxa']].merge(df_name2lineage, on='lowest_taxa', how='left')
    df_new_taxa.drop(['Phylum', 'Class', 'lowest_taxa'], inplace=True, axis=1)
    
    taxa_by_taxonkit = df_new_taxa[~df_new_taxa.Superkingdom.isna()].copy()
    taxa_no_taxonkit = df_new_taxa[df_new_taxa.Superkingdom.isna()].copy()
    
    tbl = taxa_no_taxonkit[['contig_id']].merge(df_raw_taxa, on='contig_id', how='left')
    tbl['taxid'] = 'unknown'
    viral_taxa = pd.concat([taxa_by_taxonkit, tbl[taxa_by_taxonkit.columns]], axis=0)
    viral_taxa.fillna('no_rank', inplace=True)

    # concat with all virus list
    viral_all = pd.read_csv(virus_list, names=['contig_id'])
    viral = viral_all.merge(viral_taxa, on='contig_id', how='left')
    viral.taxid.fillna('unknown', inplace=True)
    viral.fillna('no_rank', inplace=True)
    viral.sort_values('Superkingdom', inplace=True)

    if assembler == 'spades':
        viral['length'] = viral.apply(lambda x:int(x['contig_id'].split('_')[4]), axis=1)
    elif assembler == 'megahit':
        viral['length'] = viral.apply(lambda x:int(x['contig_id'].split('=')[3]), axis=1)
    else:
        print('assember should be spades or megahit!')
    viral.to_csv(out, sep='\t', index=False)

if __name__ == '__main__':
    main()