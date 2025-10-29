TITLE = "CYP2D6 allele from variants"
DESC = "Parse the haplotype definition file, and determine the closest match for each set of variants"


## Import Libraries

import sys
import argparse
import os
import re
import time
from Bio import SeqIO # type: ignore
import pandas as pd # type: ignore



## Declare global variables

DEBUG = False
VERB = False

## Compliment bases
BASE_COMPLIMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

## IUPAC
IUPAC = {"A": ["A"],
         "C": ["C"],
         "G": ["G"],
         "T": ["T"],
         "R": ["A","G"],
         "Y": ["C","T"],
         "S": ["G","C"],
         "W": ["A","T"],
         "K": ["G","T"],
         "M": ["A","C"],
         "B": ["C","G","T"],
         "D": ["A","G","T"],
         "H": ["A","C","T"],
         "V": ["A","C","G"]}
## all possible IUPAC codes for a given nucleotide
IUPAC_INV = {"A": ["A","R","W","M","D","H","V","N"],
             "C": ["C","Y","S","M","B","H","V","N"],
             "T": ["T","Y","W","K","B","D","H","N"],
             "G": ["G","R","S","K","B","D","V","N"]}


## Classes and functions

def statprint(msg, msg_type = "STATUS"):
    print("{message_type} [{datetime}]: {message}".format(message_type = msg_type,
             datetime = time.strftime("%Y/%m/%d %T"), 
             message = msg), flush=True)


## Main

if __name__ == "__main__":

    ## Deal with command line arguments
    parser = argparse.ArgumentParser(description = TITLE)
    parser.add_argument("--ref", dest = "REF", help = "Reference sequence", type = str)
    parser.add_argument("--hapdef", dest = "HAPDEF", help = "CYP2D6 Haplotype definition file", type = str)
    parser.add_argument("--variants", dest = "VARIANTS", help = "Consensus variant list", type = str)
    parser.add_argument("--output", help = "haplotypes output file", type = str)
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


    #######################################
    ## Step 1: Initialize Reference Array ##
    ########################################

    ## Position n of reference sequence will be at position 2n+1 of reference array
    ## Position 2n of reference array will hold additional bases in the case of insertions between position n-1 and n of reference sequence
    '''
    Reference:  A C C A T G
                0 1 2 3 4 5
    Array:     ['',  A, '',  C, '',  C, '',  A, '',  T, '',  G, '']
                 0   1   2   3   4   5   6   7   8   9  10  11  12
    
    array_index (0-based) = 2*ref_ind + 1 (0-based)
    array_index (0-based) = 2*(ref_ind - 1) + 1 (1-based)
    '''

    statprint("Initializing reference array(s)...")

    ref = {}
    REFNAME = ""
    num_ref = 0
    for record in SeqIO.parse(args.REF, "fasta"):
        num_ref += 1
        REFNAME = record.id
        ref = {"refseq": record.seq.upper(),
               "refarr": ["" for x in range((2*len(record.seq.upper()))+1)] ## to hold the reference sequence
        }

        blank_refarr = ["" for x in range((2*len(record.seq.upper()))+1)]


        ## initialize the reference sequence
        ## This is used to determine the base when an alignment matches the reference.
        for i in range(len(ref["refseq"])):
            ref["refarr"][(2*i)+1] = ref["refseq"][i]

    if num_ref > 1:
        statprint("More than 1 reference sequence received - only the last one is being used!", "ERROR")



    ## define exons and splices
    exons_and_splice = set(x for x in range(5001,5199+1))
    exons_and_splice |= set(x for x in range(5902,6073+1))
    exons_and_splice |= set(x for x in range(6626,6778+1))
    exons_and_splice |= set(x for x in range(6867,7027+1))
    exons_and_splice |= set(x for x in range(7461,7637+1))
    exons_and_splice |= set(x for x in range(7828,7969+1))
    exons_and_splice |= set(x for x in range(8177,8364+1))
    exons_and_splice |= set(x for x in range(8819,8960+1))
    exons_and_splice |= set(x for x in range(9059,9237+1))

    exons_and_splice |= set([5901, 6866, 7460, 7959, 7970, 8008])   ## splice variants


    ################################################
    ## Step 2: Load in Haplotype Definition Table ##
    ################################################

    ## NOTE These variants are given in 1-based positions, so to get the array position:
    ## array_index (0-based) = 2*(ref_ind - 1) + 1 (1-based)

    ## NOTE logic based on haplotype definition file format as of 10/31/2024 release. Downloaded 01/21/2025. Minor edits made to make variant reporting consistent.

    df = pd.read_excel(args.HAPDEF, sheet_name="Alleles")
    headers = df.iloc[0] # create headers from NG_008376.4 ATG coordinates - note need to add 5019 to get full NG_008376.4 coords
    COORD_CORRECTION = 5019
    df = pd.DataFrame(df.values[6:], columns=headers)
    df.rename(columns={'Position at NG_008376.4 (CYP2D6 RefSeqGene; reverse relative to chromosome)': 'haplotype'}, inplace=True)

    # iterate over dataframe to construct CPIC haplotypes including iupac nucleotides
    
    positions = set()
    cpic = {}
    ref_row = None

    for idx, row in df.iterrows():

        allele = 'CYP2D6' + row[0]
        statprint("Processing {}".format(allele))
        
        if allele not in ['CYP2D6*5', 'CYP2D6*13', 'CYP2D6*61', 'CYP2D6*63', 'CYP2D6*68']:
            # *5 is a complete deletion, *13,*61,*63,*68 are 2D6-2D7 hybrids and no consensus sequence exists 

            if allele == "CYP2D6*1":
                ## save this row.
                ref_row = row

            cpic[allele] = [set() for x in range((2*len(record.seq.upper()))+1)]

            ## loop through each known variant position (columns)
            for i in range(1,len(row)-1):   ## last column is structural vars, which we are not tracking here.
                

                if "ins" in row.index[i]:
                    ## insertion position
                    # 5156_5157insT
                    # or 
                    # 6907_6908insTA
                    pos = int(row.index[i].split(";")[0].split("_")[0]) + COORD_CORRECTION + 1   ## rightmost position.
                    
                    if isinstance(row[i], str) and allele != "CYP2D6*1":
                        ## there is a defined base at this position for this variant allele
                        ## insertions and deletions are not using IUPAC degeneracy. So, only one possiblity for each.
                        ## if this is a "true" insertion, ref will be "del", allele will be "insB"
                        ## if this is a duplication, ref will be the base(s) that is duplicated, allele will be "BB" or "BS(2)"
                        ## --> this requires that we cleave off the first occurrence of the wildtype from the allele string to determine what the inserted bases are
                        if ref_row[i] == "del":
                            ## allele value will be like "insTA"
                            base_string = [BASE_COMPLIMENT[x.strip()] for x in row[i][3:][::-1]]    ## reverse compliments after cleaving off leading "ins"
                        else:
                            if "(" in row[i]:
                                copies_in_this_allele = int(row[i].split("(")[1][:-1])
                                allele_base_string = row[i].split("(")[0] * copies_in_this_allele
                            else:
                                allele_base_string = row[i]
                            base_string = [BASE_COMPLIMENT[x.strip()] for x in allele_base_string[len(ref_row[i].strip()):][::-1]]  ## what is actually inserted

                    else:
                        ## no variant for this allele
                        base_string = ""

                    ## build haplotype
                    cpic[allele][int((2*(pos-0.5-1))+1)].add("".join(base_string))
                    positions.add(pos-0.5)  #subtracting 0.5 to make this position appear as the insertion.

                elif "del" in row.index[i]:
                    ## deletion position
                    # 6727delT
                    # or
                    # 5083delTCGA

                    ## get first deleted position
                    pos = int(row.index[i].split(";")[0].split("del")[0]) + COORD_CORRECTION

                    if isinstance(row[i], str) and allele != "CYP2D6*1":
                        ## there is a defined base at this position for this variant allele
                        ## note there are some that are deleting more than one base, eg. g.7559delAACT
                        # in this case, the positions are 7559,7560,7561,7562 all deleted.
                        ## but there are some variants that are not mutually exclusive:
                        ## g.7947delGATCCTACATCCGGATGTG	g.7955A>C	g.7959G>A
                        ## If the deletion is present, the next two SNVs can not exist (since they will be deleted).
                        ## I handle this by just updating the array.
                        ## i also use the row.index string as these are the bases that are deleted...don't have to deal with cleaving off "del" or rev comp or deal with 

                        base_string = row.index[i].split(";")[0].split("del")[1]
                        for j in range(len(base_string)):
                            ## build haplotype
                            cpic[allele][int((2*(pos-1))+1)].add("-")
                            
                            positions.add(pos)
                            pos += 1
                    else:
                        base_string = row.index[i].split(";")[0].split("del")[1]
                        for i in range(len(base_string)):
                            ## build haplotype(s)
                            if "-" not in cpic[allele][int((2*(pos-1))+1)]:  # there is not an upstream deletion that has already deleted this position.
                                cpic[allele][int((2*(pos-1))+1)].add(base_string[i])
                            pos += 1

                    

                else:
                    ## substitution position
                    # 5033C>T
                    pos = int(row.index[i].split(";")[0][:-3]) + COORD_CORRECTION   ## if there are multiple possible (sep ;), the position is the same so just look at the first.

                    if isinstance(row[i], str) and allele != "CYP2D6*1":
                        ## there is a defined base at this position for this variant allele
                        bases = [BASE_COMPLIMENT[x] for x in IUPAC[row[i].strip()]]     ## use IUPAC dict to get all possible bases, then rev comp because table is genomic, while gene is on - strand.

                        ## build haplotype
                        ## note if there was an upstream long deletion variant that exists but is wildtype in this allele, then this position will already have the wildtype base added as an acceptable base, which is not quite true.
                        ## need to account for this by removing that base if it exists.
                        if len(cpic[allele][int((2*(pos-1))+1)]) > 0 and BASE_COMPLIMENT[IUPAC[ref_row[i].strip()][0]] in cpic[allele][int((2*(pos-1))+1)]:
                            cpic[allele][int((2*(pos-1))+1)].remove(BASE_COMPLIMENT[IUPAC[ref_row[i].strip()][0]])
                        ## now add variant base for this allele as usual:
                        for b in bases:
                            ## make sure not deleted by upstream variant
                            if "-" not in cpic[allele][int((2*(pos-1))+1)]:
                                cpic[allele][int((2*(pos-1))+1)].add(b)
                       
                    else:
                        bases = [BASE_COMPLIMENT[x] for x in IUPAC[ref_row[i].strip()]]
                        for b in bases:
                            ## make sure not deleted by upstream variant
                            if "-" not in cpic[allele][int((2*(pos-1))+1)]:
                                cpic[allele][int((2*(pos-1))+1)].add(b)                    

                    positions.add(pos)
       
    
    ## now cpic[allele] contains an array of all possible variant position values for each haplotype.
    statprint("Finished reading in CPIC vars.")
    if DEBUG: 
        for allele in cpic:
            if allele == "CYP2D6*59":
                print(allele)
                #print("".join(["".join(cpic[allele][i]) for i in range(len(cpic[allele])) if ((i-1)/2)+1 in positions]))
                print([cpic[allele][i] for i in range(len(cpic[allele])) if ((i-1)/2)+1 in positions])
                print(cpic[allele][int((2*(7959-1))+1)])
    if DEBUG: sys.exit()


    ######################################################
    ## Step 3: Read in variant lists for each consensus ##
    ##         Split by coding and non-coding.          ##
    ######################################################

    ## file format:
    '''
    consensus_id    variant
    c1_7_32652reads C5119T
    c1_7_32652reads G5329T
    c1_7_32652reads T5861G
    c1_7_32652reads C6057T
    c1_7_32652reads G6681C
    c1_7_32652reads A7117G
    c1_7_32652reads A8404C
    c1_7_32652reads A8602G
    c1_7_32652reads G9200C
    c1_11_18601reads        G5233C
    c1_11_18601reads        C5240A
    c1_11_18601reads        C5242G
    c1_11_18601reads        T5246C
    c1_11_18601reads        G5251C
    c1_11_18601reads        A5252C
    c1_11_18601reads        A5264G
    c1_11_18601reads        G5329T
    c1_11_18601reads        C5764G
    c1_11_18601reads        T5861G
    c1_11_18601reads        G6681C
    c1_11_18601reads        C7870T
    c1_11_18601reads        G8008A
    c1_11_18601reads        A8404C
    c1_11_18601reads        G8604A
    c1_11_18601reads        C8810T
    c1_11_18601reads        G9200C
    '''

    convars = {}

    header = True
    for line in open(args.VARIANTS, "r"):
        if header:
            header = False
        else:
            cid, var = line.rstrip().split("\t")
            if cid not in convars:
                convars[cid] = []
            convars[cid].append(var)
    

    ## generate list of variant bases to compare to cpic

    consensus_haplotypes = {}
    
    for cid in convars:
        consensus_haplotypes[cid] = blank_refarr.copy()
        
        for var in convars[cid]:

            if "ins" in var:
                ## insertion position
                # 5157_5158insT
                pos = int(var.split("_")[0]) + 1   ## rightmost position.
                
                base_string = var.split("ins")[1]

                if pos in exons_and_splice or pos-1 in exons_and_splice:
                    ## build haplotype
                    consensus_haplotypes[cid][int((2*(pos-0.5-1))+1)] = base_string
                ## else ignore because this is not a coding mutation

            elif "del" in var:
                ## deletion position
                # 6727delT
                pos = float(var.split("del")[0])

                ## note there are some that are deleting more than one base, eg. g.7559delAACT
                # in this case, the positions are 7559,7560,7561,7562 all deleted.
                ## but there are some variants that are not mutually exclusive:
                ## g.7947delGATCCTACATCCGGATGTG	g.7955A>C	g.7959G>A
                ## If the deletion is present, the next two SNVs can not exist (since they will be deleted).
                ## I handle this by just updating the array.
                base_string = var.split("del")[1]
                for i in range(len(base_string)):
                    ## build haplotype
                    if pos in exons_and_splice:
                        consensus_haplotypes[cid][int((2*(pos-1))+1)] = "-"
                    pos += 1

            else:
                ## substitution position
                # C5033T
                pos = float(var[1:-1])
                bases = var[-1]

                ## build haplotype
                if pos in exons_and_splice:
                    consensus_haplotypes[cid][int((2*(pos-1))+1)] = bases


        ## now we fill in the array with ref bases at variant positions so that it matches what we are comparing to for the cpic array.
        for pos in positions:
            if consensus_haplotypes[cid][int((2*(pos-1))+1)] == "":
                if len(cpic["CYP2D6*1"][int((2*(pos-1))+1)]) > 0:
                    consensus_haplotypes[cid][int((2*(pos-1))+1)] = next(iter(cpic["CYP2D6*1"][int((2*(pos-1))+1)])) ## gets a value from the set of possible values





    #####################################
    ## Step 4: Find Top-Scoring Allele ##
    #####################################

    ## for each consensus, find top scoring allele, and enumerate mismatches.
    ## Ignore non-coding

    statprint("Calculating top-scoring allele...")

    best_matches = {}
    best_scores = {}

    for cid in consensus_haplotypes:
        statprint("Assessing cid {}".format(cid))
        best_matches[cid] = []
        best_scores[cid] = 9999
        allele_scores = {}

        for allele in cpic:
            if allele not in allele_scores:
                allele_scores[allele] = 9999

            this_score = 0

            for i in range(len(consensus_haplotypes[cid])):
                if (len(cpic[allele][i]) > 0 and consensus_haplotypes[cid][i] not in cpic[allele][i]) or (len(cpic[allele][i]) == 0 and consensus_haplotypes[cid][i] != ""):
                    this_score += 1
            
            if this_score < allele_scores[allele]:
                allele_scores[allele] = this_score
                

        ## get best scoring allele, enumerate any differences.
        best_score = min(value for value in allele_scores.values())
        best_scores[cid] = best_score
        statprint("Best scoring allele has a score of {}".format(best_score))

        alleles_with_best_score = set("{}".format(allele) for allele, value in allele_scores.items() if value == best_score)

        statprint("Alelles with this best score: {}".format(alleles_with_best_score))

        ## now I have allele(s) with best score. If multiple equally-scoring hits, enumerate all.
        while len(alleles_with_best_score) > 0:
            allele = alleles_with_best_score.pop()

            # step through variants again and enumerate any discrepancies for this top hit.
            this_res = {"allele": allele, "missing":[], "extra":[]}
            for i in range(len(consensus_haplotypes[cid])):
                if (consensus_haplotypes[cid][i] in cpic[allele][i]) or (consensus_haplotypes[cid][i] == "" and len(cpic[allele][i]) == 0):
                    ## they match, no issue
                    pass
                else:
                    ## mismatch
                    pos = float(((i-1)/2)+1)

                    if consensus_haplotypes[cid][i] == "" or len(cpic[allele][i]) > 0:
                        ## missing variant, either no variant in consensus or it doesnt match expected variant
                        if ref["refarr"][i] == "":
                            # insertion
                            this_res["missing"].append("{}_{}ins{}".format(int(pos-0.5), int(pos+0.5), "".join(cpic[allele][i])))
                        elif "-" in cpic[allele][i]:
                            # deletion
                            this_res["missing"].append("{}del{}".format(int(pos), ref["refarr"][i]))
                        else:
                            # SNV
                            ## only if ref and alt are differnet. Otherwise missing var is also recorded as extra var, ie missing G9200C shows as extra G9200G.
                            if ref["refarr"][i] not in cpic[allele][i]:
                                this_res["missing"].append("{}{}[{}]".format(ref["refarr"][i], int(pos), "".join(cpic[allele][i])))
                    
                    if len(cpic[allele][i]) == 0 or consensus_haplotypes[cid][i] != "":
                        ## extra variant, either no variant in cpic or it doesnt match expected variant.
                        #if pos in positions and ref["refarr"][i] == "":
                        if ref["refarr"][i] == "":
                            # insertion
                            this_res["extra"].append("{}_{}ins{}".format(int(pos-0.5), int(pos+0.5), consensus_haplotypes[cid][i]))
                        elif consensus_haplotypes[cid][i] == "-":
                            # deletion
                            this_res["extra"].append("{}del{}".format(int(pos), ref["refarr"][i]))
                        else:
                            # snv
                            if ref["refarr"][i] != consensus_haplotypes[cid][i]:
                                this_res["extra"].append("{}{}{}".format(ref["refarr"][i], int(pos), consensus_haplotypes[cid][i]))
                    
            
            best_matches[cid].append(this_res)
        


    ##########################
    ## Step 5: Write output ##
    ##########################

    out = open(args.output, "w")
    out.write("consensus_id\ttop_allele\tscore\tmissing_var\textra_var\n")
    for cid in best_matches:
        for result in best_matches[cid]:
            out.write(f"{cid}\t{result['allele']}\t{best_scores[cid]}\t{','.join(result['missing'])}\t{','.join(result['extra'])}\n")
    out.close()
