## Allele Reconstruction Pipeline
## Step 2

'''
5. Cluster reads by length.
'''


from Bio import SeqIO
import os

WORKFLOW_DIR = os.path.dirname(os.path.abspath(workflow.snakefile))

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




## get inferred amplicons to use
INFERRED_AMPS = {}
for samp in SAMPLES:
    if not os.path.exists(os.path.join(WORKING_DIR, f"{samp}", "01_alignment")):
        os.makedirs(os.path.join(WORKING_DIR, f"{samp}", "01_alignment"))

    INFERRED_AMPS[samp] = [[],[],[]]   ## list of ids, regions, seqs
    inferred_amplicon_regions_file = os.path.join(SCRATCH_DIR, f"{samp}", "01_alignment", "alignment_genome_01.amplicons.bed")
    inferred_amplicon_seqs_file = os.path.join(SCRATCH_DIR, f"{samp}", "01_alignment", "alignment_genome_01.inferred_amplicon_seqs.fasta")
    
    inferred_regions = []
    for line in open(inferred_amplicon_regions_file, "r"):
        inferred_regions.append(line)
    ## now inferred_seqs is like ["chr4    7587067 7589193\n", "chr22   42126034        42131139\n", ...]
    
    inferred_seqs = []
    for line in open(inferred_amplicon_seqs_file, "r"):
        if line.startswith(">"):
            inferred_seqs.append(line)
        else:
            inferred_seqs[-1] += line.upper()
    ## now inferred_seqs is like [">seq1\nACTGCTACA\n", ">seq2\nTCGTCGTCGC\n", ...]

    for infseqid in range(len(inferred_regions)):
        INFERRED_AMPS[samp][0].append(infseqid)
        
        INFERRED_AMPS[samp][1].append(os.path.join(WORKING_DIR, f"{samp}", "01_alignment", f"alignment_genome_01.amplicon_{infseqid}.bed"))
        out = open(INFERRED_AMPS[samp][1][-1], "w")
        out.write(inferred_regions[infseqid])
        out.close()

        INFERRED_AMPS[samp][2].append(os.path.join(WORKING_DIR, f"{samp}", "01_alignment", f"alignment_genome_01.inferred_amplicon_seq_{infseqid}.fasta"))
        out = open(INFERRED_AMPS[samp][2][-1], "w")
        out.write(inferred_seqs[infseqid])
        out.close()





rule all:
    input:
        [os.path.join(WORKING_DIR, f"{samp}", "02_clusters", f"cluster_stats_amp{ampid}.tsv") for samp in SAMPLES for ampid in INFERRED_AMPS[samp][0]],
        [os.path.join(WORKING_DIR, f"{samp}", "01_alignment", f"candidate_reads_amp{ampid}.numreads") for samp in SAMPLES for ampid in INFERRED_AMPS[samp][0]]


##########################
##  01  Align to Human  ##
##########################

rule get_candidate_readnames:
    '''extract of names of reads that align to the high coverage region'''
    input:
        bam = os.path.join(SCRATCH_DIR, "{sample}", "01_alignment", "alignment_genome_01.bam"),
        amplicon_region = os.path.join(WORKING_DIR, "{sample}", "01_alignment", "alignment_genome_01.amplicon_{ampid}.bed")
    output:
        txt = os.path.join(SCRATCH_DIR, "{sample}", "01_alignment", "candidate_readnames_amp{ampid}.txt")
    resources:
        single_concurrent=1
    conda: "env/conda_env.yaml"
    shell:
        """samtools view -L {input.amplicon_region} {input.bam} | awk '{{print $1}}' | sort | uniq > {output.txt}"""

rule get_candidate_reads:
    '''Subset fastq to identified reads'''
    input:
        fastq = os.path.join(SCRATCH_DIR, "{sample}", "00_fastq", "porechopped.fastq"),
        candidate_readnames = rules.get_candidate_readnames.output.txt
    output:
        fastq = os.path.join(SCRATCH_DIR, "{sample}", "01_alignment", "candidate_reads_amp{ampid}.fastq")
    run:
        readnames = set()
        for line in open(input.candidate_readnames, "r"):
            readnames.add(line.rstrip())
        
        records = []
        for record in SeqIO.parse(input.fastq, "fastq"):
            if record.id in readnames:
                records.append(record)

        with open(output.fastq, "w") as out:
            SeqIO.write(records, out, "fastq")

rule get_number_candidate_reads_01:
    '''Get the number of candidate reads'''
    input:
        txt = rules.get_candidate_readnames.output.txt,
    output:
        txt = os.path.join(WORKING_DIR, "{sample}", "01_alignment", "candidate_reads_amp{ampid}.numreads")
    shell:
        "wc -l {input.txt} | awk '{{print $1}}' > {output.txt}"


##########################
##  02  Length Cluster  ##
##########################

rule length_based_clustering:
    '''Cluster reads by length'''
    input:
        FASTQ = rules.get_candidate_reads.output.fastq,
    output:
        LENGTH_CLUSTERING_FLAG = os.path.join(WORKING_DIR, "{sample}", "02_clusters", "02_clusters_amp{ampid}.done"),
        cluster_stats = os.path.join(WORKING_DIR, "{sample}", "02_clusters", "cluster_stats_amp{ampid}.tsv")
    params:
        allowed_neighbour_length_variability = 10,
        min_core_point_neighbours = lambda wildcards: int((sum(1 for line in open(os.path.join(SCRATCH_DIR, f"{wildcards.sample}", "01_alignment", f"candidate_reads_amp{wildcards.ampid}.fastq")))/4) * 0.001),     ## 0.1 % of reads.
        output_dir = lambda wildcards: os.path.join(SCRATCH_DIR, f"{wildcards.sample}", "02_clusters"),
        script = os.path.join(WORKFLOW_DIR,"read_length_clustering.py"),
    conda: "env/conda_env.yaml"
    shell:
        "python {params.script} --fastq {input.FASTQ} --anlv {params.allowed_neighbour_length_variability} --mcpc {params.min_core_point_neighbours} --outdir {params.output_dir} --cluster_stats {output.cluster_stats} --lcf {output.LENGTH_CLUSTERING_FLAG} --ampid {wildcards.ampid}"
    # run:
    #     import numpy as np
    #     from sklearn.cluster import DBSCAN
    #     from collections import defaultdict

    #     def read_fastq(file_path):
    #         sequences = []
    #         headers = []
    #         qualities = []
    #         with open(file_path, "r") as f:
    #             while True:
    #                 header = f.readline().strip()
    #                 if not header:
    #                     break
    #                 sequence = f.readline().strip()
    #                 plus = f.readline().strip()
    #                 quality = f.readline().strip()
    #                 headers.append(header)
    #                 sequences.append(sequence)
    #                 qualities.append(quality)
    #         return headers, sequences, qualities

    #     def write_fastq(file_path, headers, sequences, qualities):
    #         with open(file_path, 'w') as f:
    #             for header, sequence, quality in zip(headers, sequences, qualities):
    #                 f.write(f"{header}\n{sequence}\n+\n{quality}\n")

    #     if not os.path.isdir(params.output_dir):
    #         os.mkdir(params.output_dir)

    #     headers, sequences, qualities = read_fastq(input.FASTQ)
    #     lengths = np.array([len(seq) for seq in sequences]).reshape(-1, 1)

    #     # Perform DBSCAN clustering
    #     dbscan = DBSCAN(eps=params.allowed_neighbour_length_variability, min_samples=params.min_core_point_neighbours).fit(lengths)
    #     clusters = dbscan.labels_

    #     # Create subsets based on clusters
    #     clustered_data = defaultdict(list)
    #     for cluster_id, header, sequence, quality, length in zip(clusters, headers, sequences, qualities, lengths):
    #         clustered_data[cluster_id].append((header, sequence, quality, length))
    #     out = open(output.cluster_stats, "w")
    #     out.write("cluster\tnum_reads\tmean_length\n")
    #     # Write subsets to files
    #     for cluster_id, data in clustered_data.items():
    #         subset_file = f'{params.output_dir}/amp{wildcards.ampid}_cluster{cluster_id}.fastq'
    #         subset_headers, subset_sequences, subset_qualities, subset_lengths = zip(*data)
    #         out.write("{}\t{}\t{}\n".format(cluster_id, len(subset_headers), (sum(subset_lengths)/len(subset_lengths))[0]))
    #         write_fastq(subset_file, subset_headers, subset_sequences, subset_qualities)
    #     out.close()

    #     # write flag.
    #     out = open(output.LENGTH_CLUSTERING_FLAG, "w")
    #     out.write("{}\n".format(params.min_core_point_neighbours))
    #     out.close()

