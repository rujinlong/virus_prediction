#!/usr/bin/env nextflow
// Usage: nextflow main.nf --fid project_id --contigs contigs.fna

//Creates working dir
fid=params.fid
workingpath = params.outdir + "/" + fid
workingdir = file(workingpath)
if( !workingdir.exists() ) {
    if( !workingdir.mkdirs() ) {
        exit 1, "Cannot create working directory: $workingpath"
    } 
}

Channel.value(file(params.contigs)).into { ctgs_to_dedup }


process contigs_dedup {
    cache "deep"
    publishDir "${workingpath}/p01_contigs"
    conda '/home/viro/jinlong.ru/conda3/envs/binfo2'

    input:
    file(contigs_all) from ctgs_to_dedup

    output:
    file("contigs_long_nr.fna") into ctgs_to_CAT
    file("contigs_long_nr.fna") into ctgs_to_dvf
    file("contigs_long_nr.fna") into ctgs_to_virsorter
    file("contigs_long_nr.fna") into ctgs_to_combine

    when:
    params.mode == 'all'

    """
    seqkit seq -m $params.min_contig $contigs_all > contigs_long.fna
    dedupe.sh in=contigs_long.fna out=contigs_long_nr.fna
    """
}


process run_CAT {
    cache "deep"
    publishDir "${workingpath}/p02_CAT"
    conda '/home/viro/jinlong.ru/conda3/envs/binfo2'

    input:
    file(contigs) from ctgs_to_CAT

    output:
    file("out.faa") into ctg_to_orfs
    file("out.faa") into faa_to_combine
    file("CAT_virus.tsv") into virus_cat
    file("CAT_unclassified.tsv") into unclassified_cat
    file("out.contig2classification_official.txt") into taxa_cat_to_taxonomy
    file("*")

    when:
    params.mode == 'all'

    """
    CAT contigs -c $contigs -d $params.db_cat_nr -t $params.db_cat_taxa -o out -n $params.ncpu_cat
    CAT add_names -i out.contig2classification.txt -o out.contig2classification_official.txt -t $params.db_cat_taxa --only_official
    CAT summarise -c $contigs -i out.contig2classification_official.txt -o out.contig2classification_official.sum
    sed 's/*//' out.predicted_proteins.faa |  grep -v '^\$' > out.faa
    parse_CAT.py -i out.contig2classification_official.txt -m virus -o CAT_virus.tsv
    parse_CAT.py -i out.contig2classification_official.txt -m unclassified -o CAT_unclassified.tsv
    """
}


process run_dvf {
    cache "deep"
    publishDir "${workingpath}/p02_dvf"
    conda '/home/viro/jinlong.ru/conda3/envs/binfo2'

    input:
    file(contigs) from ctgs_to_dvf

    output:
    file("dvf_virus.tsv") into virus_dvf

    when:
    params.mode == 'all'

    """
    ${params.dir_dvf}/dvf.py -i $contigs -o out_dvf -c $params.ncpu_dvf
    ln -s out_dvf/*_dvfpred.txt dvf_pred.tsv
    calc_qvalue.r dvf_pred.tsv params.qvalue dvf_virus.tsv
    """
}


process run_virsorter {
    cache "deep"
    publishDir "${workingpath}/p02_virsorter"
    conda '/home/viro/jinlong.ru/conda3/envs/iVirus'

    input:
    file(contigs) from ctgs_to_virsorter

    output:
    file("virsorter_virus.tsv") into virus_virsorter

    when:
    params.mode == 'all'

    """
    wrapper_phage_contigs_sorter_iPlant.pl -f $contigs --db 2 --wdir outdir --ncpu $params.ncpu_virsorter --data-dir $params.DB_VIRSORTER --no_c
    pretty_virsorter.py -i outdir/*.csv -d outdir/fasta/*sequences_id_translation.tsv -o virsorter_virus.tsv
    """
}


process viral_contigs_proteins {
    cache "deep"
    publishDir "${workingpath}/p03_viral_ctgs_prots"
    conda '/home/viro/jinlong.ru/conda3/envs/binfo2'

    input:
    file(vs_virus) from virus_virsorter
    file(cat_virus) from virus_cat
    file(cat_unclassified) from unclassified_cat
    file(dvf_virus) from virus_dvf
    file(ctgs_nr) from ctgs_to_combine
    file(proteins) from faa_to_combine
    

    output:
    file("viral_contigs_nr.fna") into ctgs_to_abundance
    file("viral_contigs_nr.fna") into virus_nr_to_taxa
    file("viral_proteins.faa") into faa_to_VIBRANT
    file("viral_proteins.faa") into faa_to_vcontact
    file("viral_proteins.faa") into faa_to_demovir
    file("viral_proteins.faa") into faa_to_taxa
    file("viral_contigs_nr.list") into virus_nr_list

    when:
    params.mode == 'all'

    """
    parse_viral_predictor.py -s $vs_virus -c $cat_virus -u $cat_unclassified -d $dvf_virus -o viral_contigs.tsv
    cut -f1 viral_contigs.tsv | sed '1d' > viral_contigs.list
    seqkit grep -f viral_contigs.list $ctgs_nr > viral_contigs.fna
    cd-hit-est -i viral_contigs.fna -o viral_contigs_nr.fna -c 0.95 -n 9 -aS 0.8 -d 0 -T params.ncpu_cdhit -M 0 1>cdhit.log 2>&1
    grep '^>' viral_contigs_nr.fna | sed 's/^>//' > viral_contigs_nr.list
    extract_faa_by_ctgid.py -a out.faa -c viral_contigs_nr.list -o viral_proteins.faa
    """
}


process run_VIBRANT {
    cache "deep"
    publishDir "${workingpath}/p04_VIBRANT"
    conda '/home/viro/jinlong.ru/conda3/envs/binfo2'

    input:
    file(proteins) from faa_to_VIBRANT


    output:
    file("VIBRANT_viral_proteins") into func_anno

    when:
    params.mode == 'all'

    """
    VIBRANT_run.py -i $proteins -f prot -t $params.ncpu_VIBRANT -o 4 -virome
    """
}


process run_vcontact {
    cache "deep"
    publishDir "${workingpath}/p04_vcontact"
    conda '/home/viro/jinlong.ru/conda3/envs/binfo2'

    input:
    file(proteins) from faa_to_vcontact


    output:
    file("out_vContact2/genome_by_genome_overview.csv") into out_vcontact
    file("*")

    when:
    params.mode == 'all'

    """
    gene_to_genome.py -a $proteins -o vcontact_gene2genome.tsv
    vcontact --raw-proteins $proteins --rel-mode 'Diamond' --proteins-fp vcontact_gene2genome.tsv --db 'ProkaryoticViralRefSeq85-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin $params.jarfile --output-dir out_vContact2 -t $params.ncpu_vcontact --pc-inflation 1.5 --vc-inflation 1.5
    """
}


process run_demovir {
    cache "deep"
    publishDir "${workingpath}/p04_demovir"
    conda '/home/viro/jinlong.ru/conda3/envs/binfo2'

    input:
    file(proteins) from faa_to_demovir


    output:
    file("DemoVir_assignments.txt") into out_demovir

    when:
    params.mode == 'all'

    """
    ${params.dir_demovir}/demovir.sh $proteins $params.ncpu_demovir
    """
}


process taxonomy {
    cache "deep"
    publishDir "${workingpath}/p05_taxonomy"
    conda '/home/viro/jinlong.ru/conda3/envs/binfo2'

    input:
    file(out_vc) from out_vcontact
    file(taxa_cat) from taxa_cat_to_taxonomy
    file(taxa_demovir) from out_demovir
    file(ctgs_nr) from virus_nr_list
    file(proteins) from faa_to_taxa


    output:
    file("taxonomy.tsv")

    when:
    params.mode == 'all'

    """
    parse_vContact2.py -i $out_vc -a $params.assembler -m cluster -o taxa_pred_vc.tsv
    parse_CAT.py -i $taxa_cat -m virus -o taxa_pred_CAT.tsv
    combine_taxa.py -v taxa_pred_vc.tsv -c taxa_pred_CAT.tsv -d $taxa_demovir -o taxa_pred.tsv
    pretty_taxa.sh taxa_pred.tsv $ctgs_nr
    pretty_taxa.py -l $ctgs_nr -p taxa_pred.tsv -n nf_tmp_name2taxid.tsv -t nf_tmp_taxaid2lineage.tsv -o taxonomy.tsv -a $params.assembler
    """
}


workflow.onComplete { 
    println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
