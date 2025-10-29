TITLE = "Cluster Reads by Length"
DESC = "Given read lengths, cluster reads."


import argparse
import os
import sys
import time
import numpy as np
from sklearn.cluster import DBSCAN
from collections import defaultdict

def read_fastq(file_path):
    sequences = []
    headers = []
    qualities = []
    with open(file_path, "r") as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            sequence = f.readline().strip()
            plus = f.readline().strip()
            quality = f.readline().strip()
            headers.append(header)
            sequences.append(sequence)
            qualities.append(quality)
    return headers, sequences, qualities

def write_fastq(file_path, headers, sequences, qualities):
    with open(file_path, 'w') as f:
        for header, sequence, quality in zip(headers, sequences, qualities):
            f.write(f"{header}\n{sequence}\n+\n{quality}\n")

def statprint(msg, msg_type = "STATUS"):
    print("{message_type} [{datetime}]: {message}".format(message_type = msg_type,
             datetime = time.strftime("%Y/%m/%d %T"), 
             message = msg), flush=True)

## Main

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = TITLE)
    parser.add_argument("--fastq", dest = "FASTQ", help = "Fastq of reads", type = str)
    parser.add_argument("--anlv", dest = "allowed_neighbour_length_variability", help = "allowed_neighbour_length_variability", type = int)
    parser.add_argument("--mcpc", dest = "min_core_point_neighbours", help = "min_core_point_neighbours", type = int)
    parser.add_argument("--outdir", dest = "output_dir", help = "output directory", type = str)
    parser.add_argument("--cluster_stats", dest = "cluster_stats", help = "cluster stats file", type = str)
    parser.add_argument("--lcf", dest = "LENGTH_CLUSTERING_FLAG", help = "LENGTH_CLUSTERING_FLAG", type = str)
    parser.add_argument("--ampid", dest = "ampid", help = "Amplicon ID", type = str)
    parser.add_argument("-d", "--debug", action = "store_true", dest = "DEBUG", help = "Flag for setting debug/test state.")
    parser.add_argument("-v", "--verbose", action = "store_true", dest = "VERB", help = "Flag for setting verbose output.")
    args = parser.parse_args()

    ## Set Global Vars
    DEBUG = args.DEBUG
    VERB = args.VERB

    print("=======================================================")
    print("Python version: {}".format(sys.version))
    print("Python environment: {}".format(sys.prefix))
    print("Server: {}".format(os.uname()[1]))
    print("Current directory: {}".format(os.getcwd()))
    print("Command: {}".format(" ".join(sys.argv)))
    print("Time: {}".format(time.strftime("%Y/%m/%d %T")))
    print("=======================================================\n")

if not os.path.isdir(args.output_dir):
    os.mkdir(args.output_dir)

headers, sequences, qualities = read_fastq(args.FASTQ)
lengths = np.array([len(seq) for seq in sequences]).reshape(-1, 1)

# Perform DBSCAN clustering
dbscan = DBSCAN(eps=args.allowed_neighbour_length_variability, min_samples=args.min_core_point_neighbours).fit(lengths)
clusters = dbscan.labels_

# Create subsets based on clusters
clustered_data = defaultdict(list)
for cluster_id, header, sequence, quality, length in zip(clusters, headers, sequences, qualities, lengths):
    clustered_data[cluster_id].append((header, sequence, quality, length))
out = open(args.cluster_stats, "w")
out.write("cluster\tnum_reads\tmean_length\n")
# Write subsets to files
for cluster_id, data in clustered_data.items():
    subset_file = f'{args.output_dir}/amp{args.ampid}_cluster{cluster_id}.fastq'
    subset_headers, subset_sequences, subset_qualities, subset_lengths = zip(*data)
    out.write("{}\t{}\t{}\n".format(cluster_id, len(subset_headers), (sum(subset_lengths)/len(subset_lengths))[0]))
    write_fastq(subset_file, subset_headers, subset_sequences, subset_qualities)
out.close()

# write flag.
out = open(args.LENGTH_CLUSTERING_FLAG, "w")
out.write("{}\n".format(args.min_core_point_neighbours))
out.close()