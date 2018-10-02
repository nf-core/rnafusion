#!/usr/bin/env python
from __future__ import print_function
from __future__ import with_statement
import argparse
import os.path
import sys

TOOLS = ['star_fusion' ,'fusioncatcher']

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
        with open(p_input, 'r') as in_file, open(p_output, 'w') as out_file:
            func = getattr(sys.modules[__name__], p_tool)   # get function from parameter
            fusions = func(in_file)  # call function
            out_file.write(fusions)
            in_file.close()
            out_file.close()
    except IOError as error:
        sys.exit(error)
    except Exception as error:
        sys.exit(error)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Utility for transforming the results from apps to FusionInspector format""")
    parser.add_argument('-i', '--input', nargs='?', help='Input file', type=str, required=True)
    parser.add_argument('-t', '--tool', nargs='?', metavar='|'.join(TOOLS), type=str, required=True)
    parser.add_argument("-o", '--output', nargs='?', help='Output file', type=str, required=True)
    args = parser.parse_args()
    transform(args.input, args.tool, args.output)