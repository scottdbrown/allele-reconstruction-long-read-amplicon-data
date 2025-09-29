TITLE = "Recursive Consensus Generator"
DESC = "Parse the minimap2 mapped reads, and generate a set of pure consensus sequences"

## Import Libraries

import sys
import argparse
import os
import time
from Bio import SeqIO
import copy
import numpy as np


## Declare global variables

DEBUG = False
VERB = False

## Compliment bases
BASE_COMPLIMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", "Q": "Q"}

## Classes and functions

def statprint(msg, msg_type = "STATUS"):
    print("{message_type} [{datetime}]: {message}".format(message_type = msg_type,
             datetime = time.strftime("%Y/%m/%d %T"), 
             message = msg), flush=True)

def subsetDictionaryOnReadset(d, rs):
    '''
    Given a dictionary of the form
    {'A': {"+": set("r1", "r3"), "-": set("r5", "r8")}, 'T': {"+": set("r2"), "-": set("r9")}, 'C': {"+": set(), "-": set()}, 'G': {"+": set(), "-": set()}}
    and a set of readnames,
    return the same dictionary but with only data from those readnames
    '''

    # Create a new dictionary to store filtered values
    new_dict = {}

    # Iterate through the original dictionary
    for key, value_set in d.items():
        for strand, read_set in value_set.items():
            # Filter the values in the set based on whether they are present in the given set
            filtered_values = read_set.intersection(rs)
            # Add the key-value pair to the new dictionary if there are filtered values
            if key not in new_dict:
                new_dict[key] = {"+": set(), "-": set()}
            new_dict[key][strand] = filtered_values

    return(new_dict)


def processBaseString_leftIndel(obsarr, qualarr, obsarr_i, this_qual, baseseqstr, read_name, strand):
    ## since I want to right-align all the base strings, but i don't know what the longest string will be,
    ## I am inserting new positions to the front of the list.
    ## 
    ## Read 1:  TAAT
    ## Read 2:  AT
    ## 
    ## Right-aligned: TAAT
    ##                  AT
    ##
    ## Stored (tracking read names that support each base):
    ##   obsarr[i][0] = {'A': {'+': set(), '-': set()}, 'T': {'+': set(r1), '-': set()}, 'C': {'+': set(), '-': set()}, 'G': {'+': set(), '-': set()}, 'Q': {'+': set(), '-': set()}}
    ##   obsarr[i][1] = {'A': {'+': set(r1), '-': set()}, 'T': {'+': set(), '-': set()}, 'C': {'+': set(), '-': set()}, 'G': {'+': set(), '-': set()}, 'Q': {'+': set(), '-': set()}}
    ##   obsarr[i][2] = {'A': {'+': set(r1, r2), '-': set()}, 'T': {'+': set(), '-': set()}, 'C': {'+': set(), '-': set()}, 'G': {'+': set(), '-': set()}, 'Q': {'+': set(), '-': set()}}
    ##   obsarr[i][3] = {'A': {'+': set(), '-': set()}, 'T': {'+': set(r1, r2), '-': set()}, 'C': {'+': set(), '-': set()}, 'G': {'+': set(), '-': set()}, 'Q': {'+': set(), '-': set()}}
    
    ## we iterate through the base string to add, in reverse (right-aligned string, left-aligned indel)
    ## bi will index from 0:len(baseseqstr), and we use that to iterate through bases of baseseqstr from right to left.
    for bi in range(len(baseseqstr)):  ## iterate through positions of bases string
        b = baseseqstr[len(baseseqstr) - bi - 1]  ## get base starting from 3' end
        if len(obsarr[obsarr_i]) <= bi:
            # insert a position at the front of the list - not ideal for efficiency
            obsarr[obsarr_i].insert(0, {'A': {'+': set(), '-': set()}, 'T': {'+': set(), '-': set()}, 'C': {'+': set(), '-': set()}, 'G': {'+': set(), '-': set()}, 'Q': {'+': set(), '-': set()}})
            qualarr[obsarr_i].insert(0, {'A': {'+': {}, '-': {}}, 'T': {'+': {}, '-': {}}, 'C': {'+': {}, '-': {}}, 'G': {'+': {}, '-': {}}, 'Q': {'+': {}, '-': {}}})
        ## update the bi position (from end of list)
        obsarr[obsarr_i][-(1+bi)][b][strand].add(read_name)
        qualarr[obsarr_i][-(1+bi)][b][strand][read_name] = this_qual[bi]

    return (obsarr, qualarr)

def processBaseString_rightIndel(obsarr, qualarr, obsarr_i, this_qual, baseseqstr, read_name, strand):
    ## we iterate through the base string to add (left-aligned string, right-aligned indel)
    ## bi will index through obsarr[obsarr_i] left to right
    for bi in range(len(baseseqstr)):  ## iterate through positions of bases string
        b = baseseqstr[bi]  ## get base starting from 5' end
        if len(obsarr[obsarr_i]) <= bi:
            obsarr[obsarr_i].append({'A': {'+': set(), '-': set()}, 'T': {'+': set(), '-': set()}, 'C': {'+': set(), '-': set()}, 'G': {'+': set(), '-': set()}, 'Q': {'+': set(), '-': set()}})
            qualarr[obsarr_i].append({'A': {'+': {}, '-': {}}, 'T': {'+': {}, '-': {}}, 'C': {'+': {}, '-': {}}, 'G': {'+': {}, '-': {}}, 'Q': {'+': {}, '-': {}}})
        obsarr[obsarr_i][bi][b][strand].add(read_name)
        qualarr[obsarr_i][bi][b][strand][read_name] = this_qual[bi]

    return (obsarr, qualarr)

def processOperation(obsarr, qualarr, i, operator, operand, refarr, read_name, qual, read_i, strand):
    ## process cstag from .paf alignment.
    if operator == ":":
        ## bases match
        for x in range(int(operand)):
            ## x is used just to iterate through the number of matched bases, but i is the index we care about.
            if (args.MIN_BASE_QUAL_GLOBAL and qual[read_i] >= args.MIN_BASE_QUAL_GLOBAL) or (args.MIN_BASE_QUAL_PERCENTILE and qual[read_i] >= POS_QUAL_THRESH[strand][(2*i)+1]):
                ## base quality is sufficient, add the support for this base.
                this_qual = [qual[read_i]]
                (obsarr, qualarr) = processBaseString_leftIndel(obsarr, qualarr, (2*i)+1, this_qual, refarr[(2*i)+1], read_name, strand)
            else:
                ## quality insufficient, track this. This is tracked as a "Q" base.
                if len(obsarr[(2*i)+1]) == 0:
                    obsarr[(2*i)+1].insert(0, {'A': {'+': set(), '-': set()}, 'T': {'+': set(), '-': set()}, 'C': {'+': set(), '-': set()}, 'G': {'+': set(), '-': set()}, 'Q': {'+': set(), '-': set()}})
                    qualarr[(2*i)+1].insert(0, {'A': {'+': {}, '-': {}}, 'T': {'+': {}, '-': {}}, 'C': {'+': {}, '-': {}}, 'G': {'+': {}, '-': {}}, 'Q': {'+': {}, '-': {}}})
                obsarr[(2*i)+1][0]["Q"][strand].add(read_name)
                qualarr[(2*i)+1][0]["Q"][strand][read_name] = qual[read_i]
            i += 1
            read_i += 1

    elif operator == "+":
        ## insertion of bases
        inserted = ""
        this_qual = []
        for k in operand:
            ## k == each base that is inserted
            if (args.MIN_BASE_QUAL_GLOBAL and qual[read_i] >= args.MIN_BASE_QUAL_GLOBAL) or (args.MIN_BASE_QUAL_PERCENTILE and qual[read_i] >= POS_QUAL_THRESH[strand][(2*i)]):
                inserted += k
                this_qual.append(qual[read_i])
            ## we don't track poor quality insertions since likely noise, and not going to lead to heterozygous position
            read_i += 1
        if inserted != "":  ## (there was insertion with sufficient quality)
            (obsarr, qualarr) = processBaseString_leftIndel(obsarr, qualarr, 2*i, this_qual, inserted.upper(), read_name, strand)

    elif operator == "-":
        ## deletion of bases
        ## don't write anything, does not contribute to depth.
        for x in range(len(operand)):
            i += 1

    elif operator == "*":
        ## substitution of bases
        ## write second base of operand
        if (args.MIN_BASE_QUAL_GLOBAL and qual[read_i] >= args.MIN_BASE_QUAL_GLOBAL) or (args.MIN_BASE_QUAL_PERCENTILE and qual[read_i] >= POS_QUAL_THRESH[strand][(2*i)+1]):
            this_qual = [qual[read_i]]
            (obsarr, qualarr) = processBaseString_leftIndel(obsarr, qualarr, (2*i)+1, this_qual, operand[-1].upper(), read_name, strand)
        else:
            ## quality insufficient, track this.
            if len(obsarr[(2*i)+1]) == 0:
                obsarr[(2*i)+1].insert(0, {'A': {'+': set(), '-': set()}, 'T': {'+': set(), '-': set()}, 'C': {'+': set(), '-': set()}, 'G': {'+': set(), '-': set()}, 'Q': {'+': set(), '-': set()}})
                qualarr[(2*i)+1].insert(0, {'A': {'+': {}, '-': {}}, 'T': {'+': {}, '-': {}}, 'C': {'+': {}, '-': {}}, 'G': {'+': {}, '-': {}}, 'Q': {'+': {}, '-': {}}})
            obsarr[(2*i)+1][0]["Q"][strand].add(read_name)
            qualarr[(2*i)+1][0]["Q"][strand][read_name] = qual[read_i]
        i += 1
        read_i += 1

    elif operator == "Z":
        pass

    else:
        ## unknown operator, error.
        sys.exit("Unknown operator: {}".format(operator))
    
    return (obsarr, qualarr, i, read_i)


## Main

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = TITLE)
    parser.add_argument("--ref", dest = "REF", help = "Reference fasta file of sequences", type = str)
    parser.add_argument("--read_seqs", dest = "READSEQ", help = "Raw reads fastq file", type = str)
    parser.add_argument("--paf", dest = "PAF", help = "Mapped reads .paf file", type = str)
    parser.add_argument("--consensus", help = "Consensus output file", type = str)
    parser.add_argument("--chromat", help = "Chromatogram data output file", type = str)
    parser.add_argument("--name", help = "Name of sequences, add to consensus names.", type = str)
    parser.add_argument("--min_depth_factor", dest = "MIN_DEPTH_FACTOR", help = "Minimum number of reads needed to call a base is set to max_depth*MIN_DEPTH_FACTOR", type = float)
    parser.add_argument("--global_min_depth_thresh", dest = "GLOBAL_MIN_DEPTH_THRESH", help = "Fraction of total reads to be the minimum number used for a consensus.", type = float)
    parser.add_argument("--sig_to_noise_thresh", dest = "SIG_TO_NOISE_THRESH", help = "Value of minimum ratio of most frequent to second most frequent base needed to consider a position homozygous.", type = float)
    parser.add_argument("--min_base_qual_global", dest = "MIN_BASE_QUAL_GLOBAL", help = "Minimum base quality, globally, to use.", type = int)
    parser.add_argument("--min_base_qual_percentile", dest = "MIN_BASE_QUAL_PERCENTILE", help = "Minimum base quality percentile to use, per position.", type = int)
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
    
    ## General Overview:
    ## 1.  Initialize an array of the reference sequence
    ## 2.  Read in .paf file and extract relevant fields for each alignment
    ## 3.  Read in raw reads to extract seq before and after alignments
    ## 4.  Parse through each alignment from .paf file and iterate through the cs tag
    ##     Determine the observed nucleotide frequency at each position of the reference
    ## 5.  Determine the consensus sequence based on the most observed, high quality, nucleotide at each position
    ## 6.  Scan consensus (chromat data) for heterogeneous positions. If found, use bases at that position to split reads into new groups for new consensuses.
    ## 7.  Output the consensus(es) and frequency at each position for plotting.


    #q_string = """!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"""
    #qual_conversion = {k: v for k, v in zip([q for q in q_string], range(len(q_string)))}


    ########################################
    ## Step 1: Initialize Reference Array ##
    ########################################

    ## Position n of reference sequence will be at position 2n+1 of reference array
    ## Position 2n of reference array will hold additional bases in the case of insertions between position n-1 and n of reference sequence
    '''
    Reference:  A C C A T G
                0 1 2 3 4 5
    Array:     ['', A, '',  C, '',  C, '',  A, '',  T, '',  G, '']
                0   1   2   3   4   5   6   7   8   9  10  11  12
    '''


    statprint("Initializing reference array...")
    
    refarr = None
    for record in SeqIO.parse(args.REF, "fasta"):

        if refarr:
            ## already have read in a reference sequence, throw error
            sys.exit("Multiple references present in ref file - only accepts one reference.")

        read_name = record.id
        refseq = record.seq.upper()
    
        refarr = ["" for x in range((2*len(refseq))+1)] ## to hold the reference sequence
        obsarr = [[] for x in range((2*len(refseq))+1)] ## to hold the actual observed bases
        qualarr = [[] for x in range((2*len(refseq))+1)] ## to hold the actual base qualities for calculating average consensus base quality

        ## ultimately obsarr will be obsarr[array_position][base_position_from_5'_end]{dict of stranded read sets supporting each base}
        ## this allows for one array_position to hold more than a single base for each read (ie an insertion of 2 bases)
        ## In cases where there is more than 1 base, we right-align* all the possible bases across the reads (left-aligning indel):
        ## * except for read sequence downstream of an alignment end.
        ##    base pos:    0 1 2
        ##     Read 1:       A T
        ##     Read 2:         T
        ##     Read 3:     C A A
        ## Then, when checking each position, we step through base_positions at array_position

        ## initialize the reference sequence
        ## This is used to determine the base when an alignment matches the reference, because the reference base is not provided in the .paf file
        for i in range(len(refseq)):
            refarr[(2*i)+1] = refseq[i]

    statprint("Ref array length: {}".format(len(refarr)))

    DEPTH_THRESHOLDS = {}
    consensuses = {}
    chromats = {}

    INDICES_TO_USE = [0]
    COMPLETE_INDICES = set()

    
    ################################
    ## Step 2: Read in alignments ##
    ################################

    ## read in paf file
    statprint("Reading in aligments...")
    num_mapped_reads = 0    ## this actually counts alignments, not mapped reads, per se.
    paf = {}
    readnames = [set()]

    for line in open(args.PAF, "r"):
        line = line.rstrip().split("\t")
        ## PAF format:
        # Ind   Type        Description
        # 0     string      Query sequence name
        # 1     int         Query sequence length
        # 2     int         Query start coordinate (0-based)
        # 3     int         Query end coordinate (0-based)
        # 4     char        ‘+’ if query/target on the same strand; ‘-’ if opposite
        # 5     string      Target sequence name
        # 6     int         Target sequence length
        # 7     int         Target start coordinate on the original strand
        # 8     int         Target end coordinate on the original strand
        # 9     int         Number of matching bases in the mapping
        # 10    int         Number bases, including gaps, in the mapping
        # 11    int         Mapping quality (0-255 with 255 for missing)

        read_name = line[0]

        if VERB and num_mapped_reads % 1000 == 0:
            statprint("{} alignments have been read in.".format(num_mapped_reads))

        
        read_length = int(line[1])
        align_start_pos = int(line[2])
        align_end_pos = int(line[3])    ## note this is the index you would use in python subsetting [align_start_pos:align_end_pos]
        ref_start_pos = int(line[7])
        strand = line[4]

        num_mapped_reads += 1

        ## reverse indices if reverse alignment
        if strand == "-":
            align_start_pos = read_length - align_end_pos    # 0 based index of reveresed fastq record.
            align_end_pos = read_length - int(line[2])

        ## the position of the cs tag is variable, so step through and find it
        cstag = [e for e in line if e.startswith("cs:")][0]
        
        # cleave off 'cs:'
        cstag = cstag[3:]

        ## save record
        if read_name not in paf:
            ## we only take first alignment from a read to the appropriate reference
            paf[read_name] = {"align_start_pos": align_start_pos,
                            "align_end_pos": align_end_pos,
                            "ref_start_pos": ref_start_pos,
                            "strand": strand,
                            "cstag": cstag}
            
            readnames[0].add(read_name)

    GLOBAL_MIN_DEPTH = num_mapped_reads * args.GLOBAL_MIN_DEPTH_THRESH


    #####################################
    ## Step 3: Read in raw reads fastq ##
    #####################################

    ## this is used to get the read sequence before and after the alignment tracked in the .paf file
    statprint("Reading in fasta file of reads and saving upstream and downstream of alignments...")

    read_name = ""
    seq = ""
    qual = []
    for record in SeqIO.parse(args.READSEQ, "fastq"):
        read_name = record.id
        seq = record.seq.upper()
        qual = record.letter_annotations["phred_quality"]

        if read_name != "" and read_name in paf:
            if paf[read_name]["strand"] == "-":
                ## revcomp seq
                seq = "".join([BASE_COMPLIMENT[x.upper()] for x in seq[::-1]])
                qual = qual[::-1]
            paf[read_name]["upstream_seq"] = seq[:paf[read_name]["align_start_pos"]]
            paf[read_name]["downstream_seq"] = seq[paf[read_name]["align_end_pos"]:]
            paf[read_name]["quality"] = qual
    


    ############################################
    ## Step 4: Parse alignments for qualities ##
    ############################################
    
    statprint("Parsing aligned qualities...")
    
    POS_QUAL = {"+": {}, "-": {}}
    
    SPECIAL_CHARS = [":","Z","+","-","*"]

    for read_name in paf:

        operand = ""    ## specifics of variant, or number of matched bases
        operator = ""   ## one of the SPECIAL_CHARS that tells you what action is happening
        i = paf[read_name]["ref_start_pos"]

        ## get start position of read that is aligned.
        read_i = paf[read_name]["align_start_pos"]

        ## get strand
        strand = paf[read_name]["strand"]

        ## step through and parse
        for c in paf[read_name]["cstag"]:
            if c in SPECIAL_CHARS:
                ## check if existing operator to do
                if operand != "":
                    ## do the previous operator ##
                    if operator == ":":
                        ## bases match
                        for x in range(int(operand)):
                            ## x is used just to iterate through the number of matched bases, but i and read_i are the indexes we care about.
                            ## save qual
                            if (2*i)+1 not in POS_QUAL[strand]:
                                POS_QUAL[strand][(2*i)+1] = []
                            POS_QUAL[strand][(2*i)+1].append(paf[read_name]["quality"][read_i])
                            
                            i += 1
                            read_i += 1

                    elif operator == "+":
                        ## insertion of bases
                        for k in operand:
                            if (2*i) not in POS_QUAL[strand]:
                                POS_QUAL[strand][(2*i)] = []
                            POS_QUAL[strand][(2*i)].append(paf[read_name]["quality"][read_i])
                            read_i += 1
                    elif operator == "-":
                        ## deletion of bases
                        ## don't write anything, no quality associated with absense of a base.
                        for x in range(len(operand)):
                            i += 1
                    elif operator == "*":
                        ## substitution of bases
                        if (2*i)+1 not in POS_QUAL[strand]:
                            POS_QUAL[strand][(2*i)+1] = []
                        POS_QUAL[strand][(2*i)+1].append(paf[read_name]["quality"][read_i])
                        i += 1
                        read_i += 1
                    elif operator == "Z":
                        pass

                ## reset the operator
                operator = c
                operand = ""
            else:
                operand += c ## add character to the operator, since it can be multiple characters

        ## end of cstag, do the last operator
        if operator == ":":
            ## bases match
            for x in range(int(operand)):
                ## x is used just to iterate through the number of matched bases, but i and read_i are the indexes we care about.
                ## save qual
                if (2*i)+1 not in POS_QUAL[strand]:
                    POS_QUAL[strand][(2*i)+1] = []
                POS_QUAL[strand][(2*i)+1].append(paf[read_name]["quality"][read_i])
                
                i += 1
                read_i += 1

        elif operator == "+":
            ## insertion of bases
            for k in operand:
                if (2*i) not in POS_QUAL[strand]:
                    POS_QUAL[strand][(2*i)] = []
                POS_QUAL[strand][(2*i)].append(paf[read_name]["quality"][read_i])
                read_i += 1
        elif operator == "-":
            ## deletion of bases
            ## don't write anything, no quality associated with absense of a base.
            for x in range(len(operand)):
                i += 1
        elif operator == "*":
            ## substitution of bases
            if (2*i)+1 not in POS_QUAL[strand]:
                POS_QUAL[strand][(2*i)+1] = []
            POS_QUAL[strand][(2*i)+1].append(paf[read_name]["quality"][read_i])
            i += 1
            read_i += 1
        elif operator == "Z":
            pass


        ## end of alignment

    if args.MIN_BASE_QUAL_PERCENTILE:
        ## Calculate per-position thresholds, once, so I am not repeatedly calculating thresholds.
        statprint("Calculating position-specific quality thresholds...")
        POS_QUAL_THRESH = {"+": {}, "-": {}}
        for s in POS_QUAL:
            for pqi in POS_QUAL[s]:
                POS_QUAL_THRESH[s][pqi] = np.percentile(POS_QUAL[s][pqi], args.MIN_BASE_QUAL_PERCENTILE)



    ########################################
    ## Step 5: Parse alignments for bases ##
    ########################################

    ## step through alignments and parse the cs tag to determine the matching and mismatching bases
    statprint("Parsing alignments...")

    ## cs tag like:
    ## cs:Z::7-t:112+tt:71-c:49+a:34-ga:3*ct:16*ga:43-t:29+a:1-tta:52+gg:23-a:17+c:2*ag:1-g:198*ca*tc:59
    ## but "cs:" already sliced off the start.
    SPECIAL_CHARS = [":","Z","+","-","*"]
    num_parsed = 0
    for read_name in paf:
        num_parsed += 1

        if VERB and num_parsed % 1000 == 0:
            statprint("{} alignments have been parsed and processed.".format(num_parsed))

        operand = ""    ## specifics of variant, or number of matched bases
        operator = ""   ## one of the SPECIAL_CHARS that tells you what action is happening
        i = paf[read_name]["ref_start_pos"]

        ## get start position of read that is upstream.
        read_i_up = paf[read_name]["align_start_pos"] - len(paf[read_name]["upstream_seq"])

        ## get start position of read that is aligned.
        read_i = paf[read_name]["align_start_pos"]

        ## process prefix of read before it aligns
        (obsarr, qualarr) = processBaseString_leftIndel(obsarr, qualarr, 2*i, paf[read_name]["quality"][read_i_up:read_i], paf[read_name]["upstream_seq"], read_name, paf[read_name]["strand"])
        

        ## step through and parse
        for c in paf[read_name]["cstag"]:
            if c in SPECIAL_CHARS:
                ## check if existing operator to do
                if operand != "":
                    ## do the previous operator ##
                    (obsarr, qualarr, i, read_i) = processOperation(obsarr, qualarr, i, operator, operand, refarr, read_name, paf[read_name]["quality"], read_i, paf[read_name]["strand"])
                
                ## reset the operator
                operator = c
                operand = ""
            else:
                operand += c ## add character to the operator, since it can be multiple characters

        ## end of cstag, do the last operator
        (obsarr, qualarr, i, read_i) = processOperation(obsarr, qualarr, i, operator, operand, refarr, read_name, paf[read_name]["quality"], read_i, paf[read_name]["strand"])

        ## process suffix of read after it aligns
        (obsarr, qualarr) = processBaseString_rightIndel(obsarr, qualarr, 2*i, paf[read_name]["quality"][read_i:], paf[read_name]["downstream_seq"], read_name, paf[read_name]["strand"])

        ## end of alignment

    ################################
    ## Step 6: Consensus Building ##
    ################################
    
    CONSENSUS_DONE = False

    while not CONSENSUS_DONE:
        break_out = False
        for cn in INDICES_TO_USE:

            if cn not in COMPLETE_INDICES:  ## this consensus is not yet pure

                statprint("Attempting consensus id{} of [{}]...".format(cn, ",".join([str(j) for j in INDICES_TO_USE])))

                #####################################
                ## Step 6a: Calculate max coverage ##
                #####################################

                statprint("Calculating max depth across sequence...")
                max_depth = 0
                for i in range(len(obsarr)):
                    if len(obsarr[i]) > 0:
                        ## take last base position of each array position as this will have the max (left-aligned indels)
                        filtered = subsetDictionaryOnReadset(obsarr[i][-1], readnames[cn])
                        this_depth = sum(len(rs["+"]) for rs in filtered.values()) + sum(len(rs["-"]) for rs in filtered.values())       ## NOTE depth including reads with low qual at this position, and is for both strands combined
                        if this_depth > max_depth: 
                            max_depth = this_depth
                
                DEPTH_THRESHOLDS[cn] = max_depth*args.MIN_DEPTH_FACTOR

                statprint("Max depth is {}.".format(max_depth))
                statprint("DEPTH_THRESHOLD is {}.".format(DEPTH_THRESHOLDS[cn]))


                ################################################################
                ## Step 6b: Calculate consensus and chromatogram trace values ##
                ################################################################

                statprint("Calculating consensus...")
                
                consensuses[cn] = ""
                chromats[cn] = [{"+": [["",0,0], ["",0,0], ["",0,0], ["",0,0], ["",0,0]], "-": [["",0,0], ["",0,0], ["",0,0], ["",0,0], ["",0,0]]} for x in range((2*len(refseq))+1)]
                
                for x in range(len(obsarr)):    ## x is "position", with odds being 2n+1 of pos of ref, evens being in between.

                    if VERB and x % 1000 == 0:
                        statprint("{} positions have been checked...".format(x))

                    ## step through possible multiple bases from 5' to 3' (end of array to beginning)
                    for base_position in range(len(obsarr[x])):

                        ## convert to list of tuples with no zeros
                        filtered = subsetDictionaryOnReadset(obsarr[x][base_position], readnames[cn])
                        list_of_tuples_pos = [(k, len(v["+"])) for k,v in filtered.items() if len(v["+"]) > 0]
                        list_of_tuples_neg = [(k, len(v["-"])) for k,v in filtered.items() if len(v["-"]) > 0]
                        list_of_tuples_all = [(k, len(v["+"]) + len(v["-"])) for k,v in filtered.items() if len(v["+"]) > 0 or len(v["-"]) > 0]

                        ## sort this on the values
                        summary_x_pos = sorted(list_of_tuples_pos, key=lambda tup: tup[1])[::-1]
                        summary_x_neg = sorted(list_of_tuples_neg, key=lambda tup: tup[1])[::-1]
                        summary_x_all = sorted(list_of_tuples_all, key=lambda tup: tup[1])[::-1]
                        ## now summary_x is like [('C', 60), ('T', 14), ('G', 6), ('A', 5)]

                        ## remove "Q" bases (low quality)
                        summary_x_pos = [tup for tup in summary_x_pos if tup[0] != "Q"]
                        summary_x_neg = [tup for tup in summary_x_neg if tup[0] != "Q"]
                        summary_x_all = [tup for tup in summary_x_all if tup[0] != "Q"]

                        bases_pos = [summary_x_pos[tb][0] for tb in range(len(summary_x_pos))]
                        counts_pos = [summary_x_pos[tb][1] for tb in range(len(summary_x_pos))]
                        bases_neg = [summary_x_neg[tb][0] for tb in range(len(summary_x_neg))]
                        counts_neg = [summary_x_neg[tb][1] for tb in range(len(summary_x_neg))]
                        bases_all = [summary_x_all[tb][0] for tb in range(len(summary_x_all))]
                        counts_all = [summary_x_all[tb][1] for tb in range(len(summary_x_all))]

                        ## check if insufficient signal to noise, and splitting of reads needed. 
                        if ( (len(bases_all) >= 2 and counts_all[0] > DEPTH_THRESHOLDS[cn] and (counts_all[0]/counts_all[1]) < args.SIG_TO_NOISE_THRESH) and
                                (len(bases_pos) >= 2 and counts_pos[0]/counts_pos[1] < args.SIG_TO_NOISE_THRESH*1.5) and
                                (len(bases_neg) >= 2 and counts_neg[0]/counts_neg[1] < args.SIG_TO_NOISE_THRESH*1.5) ):
                            ## scrap this consensus generation, split reads, and try again.
                            statprint("Consensus id{} had a heterogeneous position ({}) - splitting reads and retrying.".format(cn, x))
                            if DEBUG: statprint(summary_x_pos, "DEBUG")
                            if DEBUG: statprint(summary_x_neg, "DEBUG")
                            ## first check number of bins - use both strands
                            num_bins = 2
                            if len(bases_all) == 3 and counts_all[0]/counts_all[2] < args.SIG_TO_NOISE_THRESH:
                                ## third base also insufficient signal to noise
                                num_bins = 3
                            elif len(bases_all) == 4 and counts_all[0]/counts_all[3] < args.SIG_TO_NOISE_THRESH:
                                ## fourth base also insufficient signal to noise
                                num_bins = 4
                            
                            statprint("Splitting reads into {} bins.".format(num_bins))
                            ## create new readsets.
                            
                            #cn is current index that we are splitting, to be deleted.
                            INDICES_TO_USE.remove(cn)
                            newi = len(readnames) ## this will be new first index
                            for newreadseti in range(num_bins):
                                this_bin_readnames = filtered[bases_all[newreadseti]]["+"] | filtered[bases_all[newreadseti]]["-"]
                                if len(this_bin_readnames) > GLOBAL_MIN_DEPTH: ## this consensus has enough depth to be output. Number of reads minus weighting for ambiguous reads
                                    readnames.append(this_bin_readnames)  ## append to readnames list, set of readnames supporting this base
                                    INDICES_TO_USE.append(newi+newreadseti)
                                else:
                                    readnames.append(set()) ## in case the above if is not met, need to keep the size of readnames and the INDICES_TO_USE in sync.
                            statprint("New consensus ids: [{}]".format(",".join([str(j) for j in INDICES_TO_USE])))
                            statprint("Read set sizes: {}".format(",".join([str(len(readnames[j])) for j in INDICES_TO_USE])))

                            ## break out of for loops and start trying to build consensuses again
                            break_out = True
                            break
                        
                        #else:
                            ## no issue at this position, moving on.

                    
                        ## MINIMUM COVERAGE CHECK
                        ## check that meets minimum depth threshold
                        if len(bases_all) > 0:
                            if counts_all[0] > DEPTH_THRESHOLDS[cn]:
                                consensuses[cn] += bases_all[0]

                                ## save raw chromat bases
                                for cb in range(len(bases_pos)):
                                    chromats[cn][x]["+"][cb][0] += bases_pos[cb]
                                    chromats[cn][x]["+"][cb][1] += counts_pos[cb]
                                    if DEBUG: print(f"x: {x}, base_position: {base_position}, cb: {cb}")
                                    chromats[cn][x]["+"][cb][2] += sum([qualarr[x][base_position][bases_pos[cb]]["+"][rn] for rn in readnames[cn] if rn in qualarr[x][base_position][bases_pos[cb]]["+"]])
                                for cb in range(len(bases_neg)):
                                    chromats[cn][x]["-"][cb][0] += bases_neg[cb]
                                    chromats[cn][x]["-"][cb][1] += counts_neg[cb]
                                    if DEBUG: print(f"x: {x}, base_position: {base_position}, cb: {cb}")
                                    chromats[cn][x]["-"][cb][2] += sum([qualarr[x][base_position][bases_neg[cb]]["-"][rn] for rn in readnames[cn] if rn in qualarr[x][base_position][bases_neg[cb]]["-"]])
                            else:
                                ## not enough support for a base here, but save chromat details.
                                for cb in range(len(bases_pos)):
                                    chromats[cn][x]["+"][cb][0] += bases_pos[cb]
                                    chromats[cn][x]["+"][cb][1] += counts_pos[cb]
                                    if DEBUG: print(f"x: {x}, base_position: {base_position}, cb: {cb}")
                                    chromats[cn][x]["+"][cb][2] += sum([qualarr[x][base_position][bases_pos[cb]]["+"][rn] for rn in readnames[cn] if rn in qualarr[x][base_position][bases_pos[cb]]["+"]])
                                for cb in range(len(bases_neg)):
                                    chromats[cn][x]["-"][cb][0] += bases_neg[cb]
                                    chromats[cn][x]["-"][cb][1] += counts_neg[cb]
                                    if DEBUG: print(f"x: {x}, base_position: {base_position}, cb: {cb}")
                                    chromats[cn][x]["-"][cb][2] += sum([qualarr[x][base_position][bases_neg[cb]]["-"][rn] for rn in readnames[cn] if rn in qualarr[x][base_position][bases_neg[cb]]["-"]])
                        #else: (implicit)
                        ## nothing at this position
                            

                    if break_out:
                        break   
            
                if break_out:
                    break
            
                ## if we have gotten here, this consensus sequence is pure and complete, do not need to process it again.
                COMPLETE_INDICES.add(cn)
        
            if break_out:
                break
        
        if not break_out:
            CONSENSUS_DONE = True
    
    
    ##########################
    ## Step 7: Write output ##
    ##########################
    
    
    con_out = open(args.consensus, "w")
    chromat_out = open(args.chromat, "w")
    chromat_out.write("consensus_id\tpos\tbase\tcount\ttop\tstrand\tmean_quality\n")

    for ind in INDICES_TO_USE:
        ## check that there is actually a consensus (if only supported by one strand, nothing written)
        if len(consensuses[ind]) > 0:
            con_out.write(">{}_{}_{}reads\n{}\n".format(args.name, ind, len(readnames[ind]), consensuses[ind]))
            
            for x in range(len(chromats[ind])):
                if x%2==1:
                    # position of reference
                    top = True
                    for ind2 in range(len(chromats[ind][x]["+"])):
                        chromat_out.write("{}_{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(args.name, ind,
                                                                        ((x-1)/2)+1,
                                                                        chromats[ind][x]["+"][ind2][0],
                                                                        chromats[ind][x]["+"][ind2][1],
                                                                        top,
                                                                        "+",
                                                                        chromats[ind][x]["+"][ind2][2]/chromats[ind][x]["+"][ind2][1] if chromats[ind][x]["+"][ind2][1] > 0 else 0))
                        top = False
                    top = True
                    for ind2 in range(len(chromats[ind][x]["-"])):
                        chromat_out.write("{}_{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(args.name, ind,
                                                                        ((x-1)/2)+1,
                                                                        chromats[ind][x]["-"][ind2][0],
                                                                        chromats[ind][x]["-"][ind2][1],
                                                                        top,
                                                                        "-",
                                                                        chromats[ind][x]["-"][ind2][2]/chromats[ind][x]["-"][ind2][1] if chromats[ind][x]["-"][ind2][1] > 0 else 0))
                        top = False

                else:
                    # insertion position
                    top = True
                    for ind2 in range(len(chromats[ind][x]["+"])):
                        chromat_out.write("{}_{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(args.name, ind,
                                                                        ((x-1)/2)+1,
                                                                        chromats[ind][x]["+"][ind2][0],
                                                                        chromats[ind][x]["+"][ind2][1],
                                                                        top,
                                                                        "+",
                                                                        chromats[ind][x]["+"][ind2][2]/chromats[ind][x]["+"][ind2][1] if chromats[ind][x]["+"][ind2][1] > 0 else 0))
                        top = False

                    top = True
                    for ind2 in range(len(chromats[ind][x]["-"])):
                        chromat_out.write("{}_{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(args.name, ind,
                                                                        ((x-1)/2)+1,
                                                                        chromats[ind][x]["-"][ind2][0],
                                                                        chromats[ind][x]["-"][ind2][1],
                                                                        top,
                                                                        "-",
                                                                        chromats[ind][x]["-"][ind2][2]/chromats[ind][x]["-"][ind2][1] if chromats[ind][x]["-"][ind2][1] > 0 else 0))
                        top = False

    con_out.close()
    chromat_out.close()

    statprint("Done.")