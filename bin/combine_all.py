#!/usr/bin/env python


import pandas as pd
import click


@click.command()
@click.option("--abuntaxa", '-a', help="abundance_taxonomy.tsv")
@click.option("--pred", '-p', help="predictor_all.tsv")
@click.option("--out", '-o', help="final.tsv") 
def main(abuntaxa, pred, out):
    abuntaxa = pd.read_csv(abuntaxa, sep='\t')
    pred = pd.read_csv(pred, sep='\t')
    df = pred.merge(abuntaxa, on='contig_id', how='outer')
    df.to_csv(out, sep='\t', index=False)
    
    
if __name__ == '__main__':
    main()