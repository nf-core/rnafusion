#!/usr/bin/env python 
import os
import re
import argparse
import sys


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


def group_NGI_files(input_files,outputfile):
    sample_pattern=re.compile("^(P[0-9]+_[0-9]+)")
    matches=[]
    for input_file in input_files:
        try:
            match=sample_pattern.search(os.path.basename(input_file)).group(1)
            if match:
                matches.append(match)
        except AttributeError:
            continue
    NGI_names=matches    
    for NGI_name in NGI_names:
        sample_files=[]
        for fusion_file in input_files:
            if os.path.basename(fusion_file).startswith(NGI_name):
                sample_files.append(fusion_file)
        outfile="{}.fusion_comparison.txt".format(NGI_name)
        read_files_store_data(sample_files,outfile)


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
    
    maxlen = max([len(l) for l in [gene_in_both,gene_star_only,gene_fc_only]])
    for idx in range(0, maxlen-1):
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
    group_NGI_files(args.input_files,args.output_file)
    read_files_store_data(args.input_files,args.output_file)
