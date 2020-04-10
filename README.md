# virus_prediction
Identify viral sequences from metagenomic assembled contigs

# Usage

```sh
nextflow main.nf --fid project_id --contigs contigs.fna
```

Final output file will be saved to `results/p05_taxonomy/taxonomy.tsv`.

# Test
The test may take more than 1 hour.

```bash
cd test
./test d01_data_raw
```

If succeed, you will find `taxonomy.tsv` in `test/result/test/p05_taxonomy` folder.

# TODO

1. Create a new conda environment to install all required softwares and packages, and export new 'requirements' file.
2. Remove redundant code. 
