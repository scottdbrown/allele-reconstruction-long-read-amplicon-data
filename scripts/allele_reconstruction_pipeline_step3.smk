## Allele Reconstruction Pipeline
## Step 3

'''
6. For each length cluster with abundance greater than minimum size, align reads to their inferred amplicon reference.
7. do iterative consensus workflow for all aligned read groups
8. Will be left with n number of consensus sequences, each derived from m reads. These should all be 'clean' consensus sequences
9. do variant calling on the consensus sequences as normal.
10. Score each set of variants for best match ot known CYP2D6 allele.
'''


from Bio import SeqIO
import os
import glob
import re

WORKFLOW_DIR = os.path.dirname(os.path.abspath(workflow.snakefile))

## Load Config
NG_008376_REFSEQ = config["NG_008376_REFSEQ"]
HUMAN_REF = config["HUMAN_REF"]
STAR_ALLELE_DEFINITIONS = config["STAR_ALLELE_DEFINITIONS"]

BASECALLED_DIR = config["BASECALLED_DIR"]

WORKING_DIR = config["WORKING_DIR"]
SCRATCH_DIR = config["SCRATCH_DIR"]





SAMPLES = os.listdir(BASECALLED_DIR)

## get specific basecall root dir for each 
BD = {}
for samp in SAMPLES:
    BD[samp] = BASECALLED_DIR


MIN_FRACTION_LENGTH_CLUSTER_SIZE = 0.25

## get length clusters to use, and check if inferred amplicon is of the gene of interest.
LENGTH_CLUSTERS = {}
for samp in SAMPLES:
    if samp not in LENGTH_CLUSTERS:
        LENGTH_CLUSTERS[samp] = {}
    for ampid in [re.search("cluster_stats_amp(\d+)\.tsv", f).group(1) for f in glob.glob(os.path.join(WORKING_DIR, f"{samp}", "02_clusters", 'cluster_stats_amp*.tsv'))]:
        ## now have ampids
        LENGTH_CLUSTERS[samp][ampid] = [[],""]
        cluster_stats_file = os.path.join(WORKING_DIR, f"{samp}", "02_clusters", f'cluster_stats_amp{ampid}.tsv')
        max_num_reads = 0
        HEADER = True
        for line in open(cluster_stats_file, "r"):
            if HEADER:
                HEADER = False
            else:
                if int(line.split("\t")[1]) > max_num_reads:
                    max_num_reads = int(line.split("\t")[1])
        ## now that we have max, determine min threshold
        min_num_reads = max_num_reads * MIN_FRACTION_LENGTH_CLUSTER_SIZE
        ## now loop through again and save clusters that meet this min size
        HEADER = True
        for line in open(cluster_stats_file, "r"):
            if HEADER:
                HEADER = False
            else:
                if int(line.split("\t")[1]) >= min_num_reads and int(line.split("\t")[0]) != -1:
                    LENGTH_CLUSTERS[samp][ampid][0].append(line.split("\t")[0])

        ## also pick the cluster closest in length to the inferred amplicon
        ampseqlen = len(open(os.path.join(WORKING_DIR, f"{samp}", "01_alignment", f"alignment_genome_01.inferred_amplicon_seq_{ampid}.fasta")).readlines()[-1].strip())
        closest_cluster = [-1,999999]
        HEADER = True
        for line in open(cluster_stats_file, "r"):
            if HEADER:
                HEADER = False
            else:
                if abs(float(line.split("\t")[2]) - ampseqlen) < closest_cluster[1] and int(line.split("\t")[0]) != -1:
                    closest_cluster[1] = abs(float(line.split("\t")[2]) - ampseqlen)
                    closest_cluster[0] = line.split("\t")[0]
        if closest_cluster[0] not in LENGTH_CLUSTERS[samp][ampid][0]:
            LENGTH_CLUSTERS[samp][ampid][0].append(closest_cluster[0])

        ## check which reference to use
        amp_coords = open(os.path.join(WORKING_DIR, f"{samp}", "01_alignment", f"alignment_genome_01.inferred_amplicon_seq_{ampid}.fasta")).readlines()[0].strip()
        amp_coords = amp_coords.split(">")[1]
        chrom, coords = amp_coords.split(":")
        start_coord, end_coord = coords.split("-")
        LENGTH_CLUSTERS[samp][ampid][1] = os.path.join(WORKING_DIR, f"{samp}", "01_alignment", f"alignment_genome_01.inferred_amplicon_seq_{ampid}.fasta")



rule all:
    input:
        [os.path.join(WORKING_DIR, f"{samp}", "05_consensus", f"{samp}_amp{ampid}_lenclust{cluster_id}_alleles.tsv") for samp in SAMPLES for ampid in LENGTH_CLUSTERS[samp] for cluster_id in LENGTH_CLUSTERS[samp][ampid][0]],

        




###########################
##  04  Alignment (paf)  ##
###########################

rule alignment_for_consensus_paf:
    '''Generate the .paf alignment file for iterative consensus generation'''
    input:
        fastq = os.path.join(SCRATCH_DIR, "{sample}", "02_clusters", "amp{ampid}_cluster{cluster_id}.fastq"),
        refseq = lambda wildcards: LENGTH_CLUSTERS[wildcards.sample][wildcards.ampid][1]
    output:
        paf = os.path.join(SCRATCH_DIR, "{sample}", "04_alignment", "amp{ampid}_lenclust{cluster_id}_alignment.paf")
    conda: "env/conda_env.yaml"
    shell:
        "minimap2 -x map-ont --cs {input.refseq} {input.fastq} > {output.paf}"



##############################
##  05  Consensus Building  ##
##############################

rule consensus_05:
    '''Pileup base called at each position and iteratively generate consensus sequences based on making "pure" consensus'''
    input:
        refseq = lambda wildcards: LENGTH_CLUSTERS[wildcards.sample][wildcards.ampid][1],
        raw_reads = os.path.join(SCRATCH_DIR, "{sample}", "02_clusters", "amp{ampid}_cluster{cluster_id}.fastq"),
        paf = rules.alignment_for_consensus_paf.output.paf
    output:
        consensus = os.path.join(WORKING_DIR, "{sample}", "05_consensus", "{sample}_amp{ampid}_lenclust{cluster_id}_consensus.fasta"),
        chromat = os.path.join(WORKING_DIR, "{sample}", "05_consensus", "{sample}_amp{ampid}_lenclust{cluster_id}_chromatogram-data.tsv"),
    params:
        script = os.path.join(WORKFLOW_DIR,"recursive_allele_reconstruction.py"),
        min_depth_factor = 0.25,
        signal_to_noise = 2,
        global_min_depth_thresh = 0.05,
        min_base_qual_percentile = 33,
        min_base_qual_global = 10,
        name = lambda wildcards: f"amp{wildcards.ampid}_lenclust{wildcards.cluster_id}",
    conda: "env/conda_env.yaml"
    shell:
        """python {params.script} \
        -v \
        --ref {input.refseq} \
        --read_seqs {input.raw_reads} \
        --paf {input.paf} \
        --name {params.name} \
        --consensus {output.consensus} \
        --chromat {output.chromat} \
        --min_depth_factor {params.min_depth_factor} \
        --global_min_depth_thresh {params.global_min_depth_thresh} \
        --min_base_qual_percentile {params.min_base_qual_percentile} \
        --min_base_qual_global {params.min_base_qual_global} \
        --sig_to_noise_thresh {params.signal_to_noise}"""


rule annotate_variants_05:
    '''From the variants, select best matching star alleles'''
    input:
        refseq = NG_008376_REFSEQ,
        consensus = rules.consensus_05.output.consensus
    output:
        variants = os.path.join(WORKING_DIR, "{sample}", "05_consensus", "{sample}_amp{ampid}_lenclust{cluster_id}_variants.txt")
    params:
        script = os.path.join(WORKFLOW_DIR,"CYP2D6_variant_enumeration.py"),
        GENE_START_POS = 5001,
        GENE_END_POS = 9312
    conda: "env/conda_env.yaml"
    shell:
        "python {params.script} --reference {input.refseq} --consensus {input.consensus} --gene_start {params.GENE_START_POS} --gene_end {params.GENE_END_POS} --output {output.variants}"


rule genotype_variants_05:
    '''From the conensus sequence(s), annotate variants'''
    input:
        refseq = NG_008376_REFSEQ,
        hapdef = STAR_ALLELE_DEFINITIONS,
        variants = rules.annotate_variants_05.output.variants
    output:
        genotypes = os.path.join(WORKING_DIR, "{sample}", "05_consensus", "{sample}_amp{ampid}_lenclust{cluster_id}_alleles.tsv")
    params:
        script = os.path.join(WORKFLOW_DIR,"CYP2D6_allele_matcher.py")
    conda: "env/conda_env.yaml"
    shell:
        "python {params.script} --ref {input.refseq} --hapdef {input.hapdef} --variants {input.variants} --output {output.genotypes}"

