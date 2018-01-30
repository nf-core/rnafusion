#!/usr/bin/env python 
import os
import re
import argparse
import sys

import pdb

def read_files_store_data(input_files,output_file):
    fusion_dict={}
    for input_file in input_files:
        if input_file.endswith("star-fusion.fusion_candidates.final.abridged"):
           #We have a star fusion file
            with open(input_file, 'r') as f:
                for line in f:
                    if line.startswith("#"):
                        #Skip header
                            continue
                    else:
                        fusion=line.split("\t")[0]
                        # If we want to store to metadata then that can be inserted here
                        if fusion in fusion_dict.keys():
                            fusion_dict[fusion]='Both'
                        else:
                            fusion_dict[fusion]='STAR'

        elif input_file.endswith("summary_candidate_fusions.txt"):
            #We have a Fusion catcher file
            pdb.set_trace()
            with open(input_file, 'r') as f:
                for line in f:
                    if line.startswith("  * "):
                        fusion=line.split(" ")[3]
                        if fusion in fusion_dict.keys():
                            fusion_dict[fusion]='Both'
                        else:
                            fusion_dict[fusion]='FusionCatcher'

        else:
           print"Found file with incorect file ending, omitting file {}".format(input_file)
    make_report(fusion_dict, output_file)


def group_files(input_files, outputfile):
    sample_dict = {}
    # Look through the input files and find sample names.
    for input_file in input_files:
        #Check for Star-fusion
        if input_file.endswith('star-fusion.fusion_candidates.final.abridged'):
            key = input_file.rstrip('star-fusion.fusion_candidates.final.abridged')    
            try:
                #We have already encountered the fusioncatcher mate
                sample_dict[key].append(input_file)
            except KeyError:
                sample_dict[key]=[input_file]
        #We have fusioncatcher
        elif input_file.endswith("summary_candidate_fusions.txt"):    
            try:
                key = input_file.rstrip('summary_candidate_fusions.txt')   
                try:
                    #We have already encountered the star-fusion mate
                    sample_dict[key].append(input_file)
                except KeyError:
                    sample_dict[key]=[input_file]
            except KeyError:
                continue

    outfile="{}.fusion_comparison.txt".format(sample_dict.keys()[0])
    read_files_store_data(sample_dict.values()[0],outfile)  



def make_report(fusion_dict, output_file):
    content=str()
    gene_in_both=[]
    gene_star_only=[]
    gene_fc_only=[]
    
    len_fc=0
    len_star=0

    for fusion_gene in fusion_dict.keys():
        if fusion_dict[fusion_gene] == 'Both':
            gene_in_both.append(fusion_gene)
            len_fc+=1
            len_star+=1
        elif fusion_dict[fusion_gene] == 'STAR':
            gene_star_only.append(fusion_gene)
            len_star+=1
        elif fusion_dict[fusion_gene] == 'FusionCatcher':
            gene_fc_only.append(fusion_gene)
            len_fc+=1
    
    content+="## Number of Fusion genes detected with STAR-fusion: {} \n".format(len_star)
    content+="## Number of Fusion genes detected with FusionCatcher: {} \n".format(len_fc)
    content +="##FUSIONCATCHER\tSTAR-FUSION\tBOTH\n"
    ##cleanup
    
    
    gene_in_both=[item.rstrip() for item in gene_in_both]
    gene_star_only=[item.rstrip() for item in gene_star_only]
    gene_fc_only=[item.rstrip() for item in gene_fc_only]
    
    pdb.set_trace()
    maxlen = max([len(l) for l in [gene_in_both,gene_star_only,gene_fc_only]])
    for idx in range(0, maxlen):
        both_str = gene_in_both[idx] if len(gene_in_both) > idx else ''
        star_str = gene_star_only[idx] if len(gene_star_only) > idx else ''
        fc_str = gene_fc_only[idx] if len(gene_fc_only) > idx else ''
        content += "{}\t{}\t{}\n".format(fc_str, star_str, both_str)    
 
    with open(output_file, 'w') as f:
        f.write(content)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Compare two list of fusion genes and give which fusions are found in both """)
    parser.add_argument("input_files", metavar='Input file', nargs='+', default='.',
                                   help="Input files from STAR fusion and Fusion catcher ")
    parser.add_argument("output_file", metavar='Output file', nargs='?', default='fusion_comparison.txt',
                                   help="File to save output to. ")
    args = parser.parse_args() 
    group_files(args.input_files,args.output_file)
