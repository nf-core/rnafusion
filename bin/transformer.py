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

def save(p_output, p_fusions):
    old = p_output.read().split('\n')
    new = p_fusions.split('\n')
    unique_fusions = set().union(old, new)
    unique_fusions = '\n'.join(unique_fusions).lstrip()
    # clear file
    p_output.seek(0)
    p_output.truncate()
    # write fusions
    p_output.write(unique_fusions)

def fi_format(p_gene1, p_gene2):
    return '{}--{}\n'.format(p_gene1, p_gene2)

def star_fusion(p_file):
    fusions = ''
    next(p_file)    #skip header
    for line in p_file:
        tmp = line.split('\t')[0].split('--')
        fusions += fi_format(tmp[0], tmp[1])

    return fusions

def fusioncatcher(p_file):
    fusions = ''
    next(p_file)    # skip header
    for line in p_file:
        tmp = line.split('\t')
        fusions += fi_format(tmp[0], tmp[1])
    
    return fusions

def transform(p_input, p_tool, p_output):
    if not os.path.exists(p_input):
        print('Defined {} doesn\'t exist'.format(p_input))
    
    if not os.path.exists(OUTPUT):
        with open(OUTPUT, 'w'): pass

    if p_tool not in TOOLS:
       print('Defined {} not in the supported list of transformations!'.format(p_tool))

    try:
        fusions = ''
        with open(p_input, 'r') as in_file:
            func = getattr(sys.modules[__name__], p_tool)   # get function from parameter
            fusions = func(in_file).rstrip()  # call function
            in_file.close()    

        with open(p_output, 'r+') as out_file, open(SUMMARY, 'a') as summary:
            if len(fusions) > 0:
                save(out_file, fusions)
                
                summary_data = [x.split('--') for x in fusions.split('\n')]
                summary.write(dump(
                    {
                        p_tool : dict((left_gene,right_gene) for left_gene,right_gene in summary_data)
                    }, 
                    default_flow_style=False, allow_unicode=True
                ))
            else:
                summary.write(dump({p_tool: None}, default_flow_style=False, allow_unicode=True))
            
            out_file.close()
            summary.close()
    except IOError as error:
        print(error)
    except Exception as error:
        print(error)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Utility for transforming the results from apps to FusionInspector format""")
    parser.add_argument('-i', '--input', nargs='?', help='Input file', type=str, required=True)
    parser.add_argument('-t', '--tool', nargs='?', metavar='|'.join(TOOLS), type=str, required=True)
    args = parser.parse_args()
    transform(args.input, args.tool, OUTPUT)
