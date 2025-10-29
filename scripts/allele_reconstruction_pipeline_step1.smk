## Allele Reconstruction Pipeline
## Step 1

'''
1. Merge fastqs
2. Run Porechop to remove any remaining adapters, and to split possibly chimeric reads
3. Get read lengths
4. Align reads to human genome, extract regions of high coverage.
'''


from Bio import SeqIO
import os

## Load Config
NG_008376_REFSEQ = config["NG_008376_REFSEQ"]
HUMAN_REF = config["HUMAN_REF"]

BASECALLED_DIR = config["BASECALLED_DIR"]

WORKING_DIR = config["WORKING_DIR"]
SCRATCH_DIR = config["SCRATCH_DIR"]



SAMPLES = os.listdir(BASECALLED_DIR)

## get specific basecall root dir for each 
BD = {}
for samp in SAMPLES:
    BD[samp] = BASECALLED_DIR



rule all:
    input:
        expand(os.path.join(SCRATCH_DIR, "{sample}", "01_alignment", "alignment_genome_01.inferred_amplicon_seqs.fasta"), sample = SAMPLES),
        expand(os.path.join(WORKING_DIR, "{sample}", "00_fastq", "readlens.txt"), sample = SAMPLES),
        expand(os.path.join(WORKING_DIR, "{sample}", "00_fastq", "all.fastq.numreads"), sample = SAMPLES),
   




#######################
##  00  Merge fastq  ##
#######################

rule merge_fastq:
    '''Merge fastq files'''
    params:
        fastq_dir = lambda wildcards: os.path.join(BD[wildcards.sample], wildcards.sample, "basecalls", "*")
    output:
        fastq_file = os.path.join(SCRATCH_DIR, "{sample}", "00_fastq", "all.fastq")
    shell:
        "cat {params.fastq_dir}/*.fastq > {output.fastq_file}"

rule porechop:
    '''Run Porechop to cleave possibly chimeric reads'''
    input:
        fastq_file = rules.merge_fastq.output.fastq_file
    output:
        fastq_file = os.path.join(SCRATCH_DIR, "{sample}", "00_fastq", "porechopped.fastq")
    conda: "env/porechop_env.yaml"
    shell:
        """porechop -i {input.fastq_file} -o {output.fastq_file}"""

rule get_read_lengths:
    '''Get file of read lengths of all reads for summary stats'''
    input:
        fastq = rules.porechop.output.fastq_file
    output:
        readlens = os.path.join(WORKING_DIR, "{sample}", "00_fastq", "readlens.txt")
    shell:
        "awk '{{if(NR%4==0){{print length($1)}}}}' {input.fastq} > {output.readlens}"

rule get_number_reads_00:
    '''Get the number of merged reads'''
    input:
        fastq = rules.porechop.output.fastq_file
    output:
        txt = os.path.join(WORKING_DIR, "{sample}", "00_fastq", "all.fastq.numreads")
    shell:
        "echo $(cat {input.fastq} | wc -l)/4|bc > {output.txt}"


#######################################################
##  01  Align to Human and get High Coverage Regions ##
#######################################################

rule alignment_sam_01:
    '''Generate the .sam alignment file for extracting CYP2D6'''
    input:
        fastq = rules.porechop.output.fastq_file,
        refseq = HUMAN_REF
    output:
        sam = os.path.join(SCRATCH_DIR, "{sample}", "01_alignment", "alignment_genome_01.sam")
    threads: 4
    conda: "env/conda_env.yaml"
    shell:
        "minimap2 -x map-ont --secondary=no -t {threads} -a {input.refseq} {input.fastq} > {output.sam}"

rule sort_and_index:
    '''Sort and Index the sam file'''
    input:
        sam = rules.alignment_sam_01.output.sam
    output:
        bam = os.path.join(SCRATCH_DIR, "{sample}", "01_alignment", "alignment_genome_01.bam")
    conda: "env/conda_env.yaml"
    shell:
        """samtools sort -o {output.bam} {input.sam}; samtools index {output.bam}"""

rule depth_genome:
    '''Run samtools depth on genome alignment'''
    input:
        bam = rules.sort_and_index.output.bam
    output:
        depthfile = os.path.join(SCRATCH_DIR, "{sample}", "01_alignment", "alignment_genome_01.depth")
    conda: "env/conda_env.yaml"
    shell:
        """samtools depth {input.bam} > {output.depthfile}"""

rule get_high_coverage:
    '''Get the high coverage regions'''
    input:
        depthfile = rules.depth_genome.output.depthfile
    output:
        high_coverage_bed = os.path.join(SCRATCH_DIR, "{sample}", "01_alignment", "alignment_genome_01.highcov.bed")
    params:
        top_percent = 0.50
    shell:
        """awk 'BEGIN {{max = 0}} {{if($3 > max) max=$3}} END {{close(FILENAME); while ((getline < FILENAME) > 0) if($3 >= max*{params.top_percent}){{print $1"\t"$2"\t"($2+1)"\t"$3}}}}' {input.depthfile} > {output.high_coverage_bed}"""

rule merge_regions:
    '''Merge contiguous regions'''
    input:
        high_coverage_bed = rules.get_high_coverage.output.high_coverage_bed
    output:
        inferred_amplicon = os.path.join(SCRATCH_DIR, "{sample}", "01_alignment", "alignment_genome_01.amplicons.bed")
    conda: "env/conda_env.yaml"
    shell:
        """bedtools merge -d 100 -i {input.high_coverage_bed} > {output.inferred_amplicon}"""

rule extract_ref_sequence:
    '''Get the genomic sequences for the inferred amplicon regions'''
    input:
        inferred_amplicon = rules.merge_regions.output.inferred_amplicon
    output:
        inferred_amplicon_seqs = os.path.join(SCRATCH_DIR, "{sample}", "01_alignment", "alignment_genome_01.inferred_amplicon_seqs.fasta")
    params:
        refseq = HUMAN_REF
    conda: "env/conda_env.yaml"
    shell:
        """bedtools getfasta -fi {params.refseq} -bed {input.inferred_amplicon} -fo {output.inferred_amplicon_seqs}"""

