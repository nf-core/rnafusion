#!/usr/bin/env python
from __future__ import print_function
from __future__ import with_statement
from yaml import dump
import argparse
import os.path
import sys

TOOLS = ['star_fusion' ,'fusioncatcher']
SUMMARY = 'summary.yaml'
OUTPUT = 'fusions.txt'

def fi_format(p_gene1, p_gene2):
    return '{}--{}\n'.format(p_gene1, p_gene2)

def star_fusion():
    exit('Not Implemented!')

def fusioncatcher(p_file):
    fusions = ''
    next(p_file)    # skip header
    for line in p_file:
        tmp = line.split('\t')
        fusions += fi_format(tmp[0], tmp[1])
    
    return fusions

def transform(p_input, p_tool, p_output):
    if not os.path.exists(p_input):
        sys.exit('Defined {} doesn\'t exist'.format(p_input))

    if p_tool not in TOOLS:
       sys.exit('Defined {} not in the supported list of transformations!'.format(p_tool))

    try:
        with open(p_input, 'r') as in_file, open(p_output, 'a') as out_file, open(SUMMARY, 'a') as summary:
            func = getattr(sys.modules[__name__], p_tool)   # get function from parameter
            fusions = func(in_file).rstrip()  # call function
            if len(fusions) > 0:
                out_file.write(fusions + '\n')
                # append to summary
                summary_data = [x.split('--') for x in fusions.split('\n')]
                summary.write(dump({p_tool : dict((k,v) for k,v in summary_data)}, default_flow_style=False, allow_unicode=True))
            
            # closing files
            in_file.close()
            out_file.close()
            summary.close()
    except IOError as error:
        sys.exit(error)
    except Exception as error:
        sys.exit(error)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Utility for transforming the results from apps to FusionInspector format""")
    parser.add_argument('-i', '--input', nargs='?', help='Input file', type=str, required=True)
    parser.add_argument('-t', '--tool', nargs='?', metavar='|'.join(TOOLS), type=str, required=True)
    args = parser.parse_args()
    transform(args.input, args.tool, OUTPUT)
