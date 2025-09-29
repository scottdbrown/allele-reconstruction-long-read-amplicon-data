TITLE = "CYP2D6 Variant Caller"
DESC = "Given the CYP2D6 reference sequence and a consensus sequence, report all SNV and indel variants."


## Import Libraries

import sys
import argparse
import os
import time
from Bio import pairwise2
from Bio import SeqIO
from Bio.Seq import Seq



## Declare global variables

DEBUG = False
VERB = False


## Classes and functions

class bcolors:
    CYAN = '\033[1;36;40m'
    BLUE = '\033[1;34;40m'
    GREEN = '\033[1;32;40m'
    YELLOW = '\033[1;33;40m'
    RED = '\033[1;31;40m'
    BOLDWHITE = '\033[1;37;40m'
    DARK = '\033[1;30;40m'
    PURPLE = '\033[1;35;40m'
    ENDC = '\033[0m'

def statprint(msg, msg_type = "STATUS"):
    typeColour = ""
    if msg_type == "ERROR":
        typeColour = bcolors.RED
    elif msg_type == "WARNING":
        typeColour = bcolors.YELLOW
    elif msg_type == "DEBUG":
        typeColour = bcolors.GREEN
    elif msg_type == "SUBPROCESS":
        typeColour = bcolors.GREEN
        msg_type = "     " + msg_type
    else:
        typeColour = bcolors.BOLDWHITE

    print("{message_color}{message_type}{end_color} {time_color}[{datetime}]{end_color}: {message}".format(message_color = typeColour, 
             message_type = msg_type,
             end_color = bcolors.ENDC, 
             time_color = bcolors.BLUE, 
             datetime = time.strftime("%Y/%m/%d %T"), 
             message = msg))





## Main

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = TITLE)
    parser.add_argument("--reference", help = "Reference sequence of plasmid (FASTA)", type = str)
    parser.add_argument("--consensus", help = "Consensus sequence.", type = str)
    parser.add_argument("--gene_start", help = "Start of gene in reference.", type = int)
    parser.add_argument("--gene_end", help = "End of gene in reference.", type = int)
    parser.add_argument("--output_file", help = "Prefix name to give output files", type = str)
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


    # read in reference
    statprint("Loading Reference...")
    ref = ""
    for line in open(args.reference, "r"):
        if not line.startswith(">"):
            ref += line.rstrip().upper()

    out = open(args.output_file, "w")
    out.write("consensus_id\tvariant\n")
    
    # read in consensus
    statprint("Loading Consensus...")
    for record in SeqIO.parse(args.consensus, "fasta"):
        con_name = record.id
        con_seq = record.seq
    
        # align them
        statprint("Aligning sequences...")
        alignment = pairwise2.align.localxs(ref[::-1], con_seq[::-1], -2, 0, one_alignment_only=True)   # penalize gaps with -1 to open, 0 to extend.

        ref_align = alignment[0].seqA[::-1]
        con_align = alignment[0].seqB[::-1]

        if DEBUG: print("Ref align string: {}".format(ref_align))
        if DEBUG: print("Consensus align string: {}".format(con_align))

        # parse alignment to get SNVs, indels.
        statprint("Detecting variants...")

        variants = {}
        # this will hold variants. Key will be start position and mutation type (for ins/del) or ref (snv), value will be mut bases.

        i = 0
        
        for r, c in zip(ref_align, con_align):
            ## r is ref base
            ## c is consensus base
            
            ## increment position i if there is not a gap in reference.
            if r != "-":
                i += 1

            if i >= args.gene_start and i <= args.gene_end:
                ## check if bases don't match
                if r != c:
                    ## check if SNV
                    if r != "-" and c != "-":
                        variants["{}{}".format(r, i)] = c
                    elif r == "-":
                        ## insertion
                        if "{}_{}ins".format(i, i+1) not in variants:
                            variants["{}_{}ins".format(i, i+1)] = ""
                        ## add to it if already exists
                        variants["{}_{}ins".format(i, i+1)] += c
                    elif c == "-":
                        ## deletion
                        ## it is complicated to check for a continuous deletion, so will just track each position as separate deletions
                        variants["{}del".format(i)] = r

        if DEBUG: statprint("Variants Found: {}".format(variants))

        ## do reverse complement (if using unguided consensus generation, do not know which strand consensus was made from)
        con_seq_rc = str(Seq(con_seq).reverse_complement())

        # align them
        statprint("Aligning reverse complement sequences...")
        alignment = pairwise2.align.localxs(ref[::-1], con_seq_rc[::-1], -2, 0, one_alignment_only=True)   # penalize gaps with -1 to open, 0 to extend.

        ref_align = alignment[0].seqA[::-1]
        con_align = alignment[0].seqB[::-1]

        if DEBUG: print("Ref align string: {}".format(ref_align))
        if DEBUG: print("Consensus align string: {}".format(con_align))

        # parse alignment to get SNVs, indels.
        statprint("Detecting variants...")

        variants_rc = {}
        # this will hold variants. Key will be start position and mutation type (for ins/del) or ref (snv), value will be mut bases.

        i = 0

        for r, c in zip(ref_align, con_align):
            ## r is ref base
            ## c is consensus base
            
            ## increment position i if there is not a gap in reference.
            if r != "-":
                i += 1

            if i >= args.gene_start and i <= args.gene_end:
                ## check if bases don't match
                if r != c:
                    ## check if SNV
                    if r != "-" and c != "-":
                        variants_rc["{}{}".format(r, i)] = c
                    elif r == "-":
                        ## insertion
                        if "{}_{}ins".format(i, i+1) not in variants_rc:
                            variants_rc["{}_{}ins".format(i, i+1)] = ""
                        ## add to it if already exists
                        variants_rc["{}_{}ins".format(i, i+1)] += c
                    elif c == "-":
                        ## deletion
                        ## it is complicated to check for a continuous deletion, so will just track each position as separate deletions
                        variants_rc["{}del".format(i)] = r

        if DEBUG: statprint("Variants Found: {}".format(variants_rc))

        ## pick best of variants and variants_rc
        if len(variants_rc) < len(variants):
            statprint("Reverse complement was a better match, using these results.")
            variants = variants_rc

        # report SNVs and indels. All positions relative to reference sequence.
        if len(variants) > 0:
            statprint("Variants detected.", "WARNING")

            for var in variants:
                out.write("{}\t{}{}\n".format(con_name, var, variants[var]))
        else:
            statprint("No variants detected.")

    out.close()
    statprint("Done.")
    