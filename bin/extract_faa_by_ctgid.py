#!/usr/bin/env python

from Bio import SeqIO
import click

@click.command()
@click.option("--faa", '-a', help="Protein FASTA file name.")
@click.option("--ctg", '-c', help="Contig list file name.")
@click.option("--header", '-t', help="Header in contig list file.")
@click.option("--out", '-o', help="Extracted protein FASTA file.")
def extract_fasta(faa, ctg, header, out):
    with open(ctg, 'r') as fh:
        ctg_list = fh.read().splitlines()
    
    if header in ctg_list:
        ctg_list.remove(header)
        
    records = []
    for record in SeqIO.parse(faa, "fasta"):
        ctg_id = '_'.join(record.id.split('_')[:-1])
        if ctg_id in ctg_list:
            records.append(record)
            
    SeqIO.write(records, out, "fasta")

if __name__ == '__main__':
    extract_fasta()