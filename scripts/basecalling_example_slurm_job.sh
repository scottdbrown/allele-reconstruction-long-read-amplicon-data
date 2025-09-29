#!/bin/bash
#SBATCH --partition=gpu3090
#SBATCH -n 4
#SBATCH --gres=gpu:1
#SBATCH --mem=20G
#SBATCH --job-name=r10.4.1flongle
#SBATCH -e /path/to/logs/NA02016/00_fastq/%j.%N.err
#SBATCH -o /path/to/logs/NA02016/00_fastq/%j.%N.out

source /home/sbrown/.bashrc;
/path/to/bin/ont-guppy_6.5.7/bin/guppy_basecaller --input_path /path/to/flowcell/fast5/ --save_path /path/to/output/NA02016/00_fastq/basecall --config dna_r10.4.1_e8.2_400bps_sup.cfg --device "cuda:0"