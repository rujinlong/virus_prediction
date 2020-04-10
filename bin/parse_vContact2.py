#!/usr/bin/env python

import pandas as pd
import click


def vcontact_classified_virus(df_vc, assembler):
    if assembler == 'spades':
        assembler_id = '_NODE_'
    elif assembler == 'megahit':
        assembler_id = '_k141_'
    else:
        print("No assembler")
    df_vc_query = df_vc[df_vc.Genome.str.contains(assembler_id)][['Genome', 'VC', 'VC_clstr', 'VC Subcluster']].copy()
    df_vc_query.rename(columns={'Genome':'contig_id'}, inplace=True)
    df_vc_ref = df_vc[~df_vc.Genome.str.contains(assembler_id)][['VC', 'VC Subcluster', 'Genome', 'Order', 'Family', 'Genus']].copy()
    # by VC subcluster
    df_vc_query_subcluster = df_vc_query.merge(df_vc_ref.drop('VC', axis=1), on='VC Subcluster', how='inner').copy()
    # by VC cluster
    df_vc_query_cluster = df_vc_query.merge(df_vc_ref.drop('VC Subcluster', axis=1), on='VC', how='inner').copy()
    # subcluster at top, and then remove duplicates
    tbl = pd.concat([df_vc_query_subcluster, df_vc_query_cluster]).drop_duplicates('contig_id')
    tbl['VC_status'] = 'Clustered'
    tbl.reset_index(inplace=True, drop=True)
    tbl.rename(columns={'VC Subcluster':'VC_subcluster'}, inplace=True)
    tbl['taxa_by_overlapVC'] = 0
    if tbl.shape[0] > 0:
        tbl['Genome'] = tbl.apply(lambda x:x['Genome'].replace('~', ' '), axis=1)
    else:
        tbl['Genome'] = ''
    return tbl[['contig_id', 'taxa_by_overlapVC', 'VC_clstr', 'VC_subcluster', 'VC_status', 'Order', 'Family', 'Genus', 'Genome']].copy()


def vcontact_others(df_novc, assembler='spades'):
    if assembler=='spades':
        assembler_id = '_NODE_'
    df_others = df_novc[(~df_novc['VC Status'].str.startswith('Overlap')) & (df_novc['Genome'].str.contains(assembler_id))]
    tbl = df_others[['Genome', 'VC Status']].copy()
    tbl.rename(columns={'Genome':'contig_id', 'VC Status': 'VC_status'}, inplace=True)
    tbl.reset_index(inplace=True, drop=True)
    return tbl


def vcontact_overlap(df_vc, df_novc, ovl_id):
    df_novc_ovl = df_novc[df_novc['VC Status']==ovl_id][['Genome', 'Order', 'Family', 'Genus', 'VC Status', 'VC Subcluster']].copy()
    df_novc_ovl['VC_clstr'] = df_novc_ovl.apply(lambda x:x['VC Status'].replace('Overlap (', '').replace(')', ''), axis=1)
    split_col = df_novc_ovl['VC_clstr'].str.split('/').apply(pd.Series, 1).stack()
    split_col.index = split_col.index.droplevel(-1)
    split_col.name='VC_clstr'
    df_novc_ovl = df_novc_ovl.drop('VC_clstr', axis=1).join(split_col)
    df_vc_ovl = df_vc[df_vc['VC_clstr'].isin(df_novc_ovl['VC_clstr'].unique())][['Genome', 'Order', 'Family', 'Genus', 'VC Status', 'VC Subcluster', 'VC_clstr']].copy()
    
    df_ovl_all = pd.concat([df_vc_ovl, df_novc_ovl], axis=0)
    
    df_ref = df_ovl_all[~df_ovl_all.Genome.str.contains('_NODE_')]
    df_ref = df_ref[df_ref['Order'] != 'Unassigned'].copy()
    df_query = df_ovl_all[df_ovl_all.Genome.str.contains('_NODE_')].copy()
    df_query.rename(columns={'Genome': 'contig_id'}, inplace=True)
    
    if df_query.shape[0] != 0:
        taxa = find_lca(df_ref)
 
        # generate final df_query table
        df_query['Order'] = taxa[0]
        df_query['Family'] = taxa[1]
        df_query['Genus'] = taxa[2]
        df_query['Genome'] = taxa[3]
        df_query.rename(columns={'VC Status': 'VC_status', 'VC Subcluster': 'VC_subcluster'}, inplace=True)
        df_query['VC_clstr'] = df_query.apply(lambda x:x['VC_clstr'] if x['VC_status'].startswith('Clustered') else x['VC_status'], axis=1)
        df_query['VC_subcluster'] = df_query.apply(lambda x:x['VC_subcluster'] if x['VC_status'].startswith('Clustered') else x['VC_status'], axis=1)
        df_query.drop_duplicates(inplace=True)
        return df_query
    
def find_lca(df):
    taxa_levels = ['Order', 'Family', 'Genus', 'Genome']

    lca = 'none'
    order_name = [x for x in df['Order'].unique() if x!='noname']
    if len(order_name) == 0:
        lca = 'Order'
    elif len(order_name) == 1:
        lca = 'Order'
        for taxa_level in taxa_levels[1:]:
            taxa_names = [x for x in df[taxa_level].unique()]
            if len(taxa_names) == 1:
                lca = taxa_level
            else:
                lca = lca
                break
    else:
        lca = 'none'
    
    if lca == 'none':
        taxa = ['noname'] * 4
    else:
        if lca == 'Genome':
            taxa_sel = taxa_levels
        else:
            taxa_sel = taxa_levels[:taxa_levels.index(lca)+1]
        taxa = []
        for taxalevel in taxa_sel:
            taxaname = [x for x in df[taxalevel].unique()]
            if len(taxaname) > 1:
                if 'noname' in taxaname:
                    taxaname.remove('noname')
            taxa += taxaname
        # fill other taxonomy level with 'noname'
        for i in range(4):
            if len(taxa) > i:
                continue
            else:
                taxa += ['noname']
    return taxa


@click.command()
@click.option("--fin", '-i', help="genome_by_genome_overview.csv")
@click.option("--method", '-m', help="[cluster, overlap, others]")
@click.option("--assembler", '-a', default='spades', help="spades or megahit")
@click.option("--out", '-o', help="pid_orf.stats")
def main(fin, method, out, assembler='spades'):
    df = pd.read_csv(fin, index_col=0)
    df.replace('Unassigned', 'noname', inplace=True)
    df_vc = df[~df['VC'].isna()].copy()
    df_vc['VC_clstr'] = df_vc.apply(lambda x:'VC_'+x['VC'].split('_')[0], axis=1)
    df_novc = df[df['VC'].isna()].copy()
    
    ovl_ids = list(df_novc[df_novc['VC Status'].str.startswith('Overlap')]['VC Status'].unique())
    
    if method == 'cluster':
        tbl_taxa = vcontact_classified_virus(df_vc, assembler=assembler)
    elif method == 'overlap':
        dfs = []
        for ovl_id in ovl_ids:
            dfs.append(vcontact_overlap(df_vc, df_novc, ovl_id))
        dfall = pd.concat(dfs)
        dfall.drop_duplicates(inplace=True)
        # find LCA for conflict taxa at unique contigs
        dup = dfall[dfall.duplicated('contig_id', keep=False)].copy()
        non_dup = dfall.drop_duplicates('contig_id', keep=False).copy()
        dfs = []
        for ctg in dup['contig_id'].unique():
            df_sel = dup[dup['contig_id']==ctg].copy()
            taxa = find_lca(df_sel)
            df_sel['Order'] = taxa[0]
            df_sel['Family'] = taxa[1]
            df_sel['Genus'] = taxa[2]
            df_sel['Genome'] = taxa[3]
            df_sel.drop_duplicates(inplace=True)
            dfs.append(df_sel)
        dfs = [x for x in dfs if x is not None]
        if len(dfs) != 0:
            dup_dedup = pd.concat(dfs)
            tbl = pd.concat([non_dup, dup_dedup], axis=0)
        else:
            tbl = non_dup.copy()
        tbl['taxa_by_overlapVC'] = tbl.apply(lambda x:2 if x['Order']=='noname' else 1, axis=1)
        tbl = tbl[['contig_id', 'taxa_by_overlapVC', 'VC_clstr', 'VC_subcluster', 'VC_status', 'Order', 'Family', 'Genus', 'Genome']]
        tbl['Genome'] = tbl.apply(lambda x:x['Genome'].replace('~', ' '), axis=1)
        tbl_taxa = tbl.copy()
    elif method == 'others':
        tbl = vcontact_others(df_novc)
    if 'tbl_taxa' in locals():
        tbl = tbl_taxa.copy()
        tbl.rename(columns={'Genome': 'Species'}, inplace=True)
        tbl['Superkingdom'] = 'Viruses'
        tbl['Phylum'] = 'noname'
        tbl['Class'] = 'noname'
    tbl.sort_values('contig_id', inplace=True)
    tbl.reset_index(drop=True, inplace=True)
    tbl.to_csv(out, sep='\t', index=False)


if __name__ == '__main__':
    main()
