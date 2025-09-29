# A Novel Algorithm for Allele Reconstruction of Pharmacogenes from Oxford Nanopore Technologies’ Amplicon-Derived Long Reads

Pre-print: [INSERT_BIORXIV]

## Install Required Software

- ont-guppy_6.5.7 (for basecalling)
- conda v22.9.0
- snakemake v6.3.0
- R v4.1.3

## Download Required References

- CYP2D6 reference sequence (data/NG_008376.fasta)
- CYP2D6 star allele definition table (data/CYP2D6_allele_definition_table_250121_cleaned250121.xls) (Original downloaded from https://www.pharmgkb.org/page/cyp2d6RefMaterials on January 21, 2025)
- Human reference genome (hg38)

## Set Up Config File

- An example config file is provided.

```yaml
## Reference files
NG_008376_REFSEQ: "../data/NG_008376.fasta"
HUMAN_REF: "/path/to/ref/hg38/genome/hg38.analysisSet.fa"
STAR_ALLELE_DEFINITIONS: "../data/CYP2D6_allele_definition_table_250121_cleaned250121.xls"

## Source directories
BASECALLED_DIR: "/path/to/basecalled/data/"

## Working directories
WORKING_DIR: "/path/to/output/directory/"
SCRATCH_DIR: "/path/to/scratch/space/"
```

- Set paths for your system.

## Basecalling raw .fast5 data

- Basecalling was run using GPU-enabled Guppy v6.5.7 using default settings for our flowcell version. An example `slurm` script is provided.

## Running Snakemake Pipeline

- The snakemake pipeline is broken into 3 separate files, to avoid working with Snakemake checkpoints.
- Run the files sequentially, as follows:

```sh
$ snakemake -p --cores 32 -s allele_reconstruction_pipeline_step1.smk --configfile example_config.yaml --resources single_concurrent=1 --use-conda
$ snakemake -p --cores 32 -s allele_reconstruction_pipeline_step2.smk --configfile example_config.yaml --resources single_concurrent=1 --use-conda
$ snakemake -p --cores 32 -s allele_reconstruction_pipeline_step3.smk --configfile example_config.yaml --resources single_concurrent=1 --use-conda
```

- Upon the first run, Snakemake will build the required conda environments for rules that require them, based on the environment definition files in `scripts/env`, and so may take a long time. Subsequent runs will be faster.

## Running R markdown

- Once the snakemake pipeline has finished running, the output files can be inspected directly, or can be parsed and summarized using the provided R markdown script (CYP2D6-ONT_allele_sequence_reconstruction_report.Rmd). Knitting this `Rmd` will produce a summary report file, and will write two `.tsv` files with results in tabular form.