#!/usr/bin/env python3
import argparse
import os.path
import sys
import json
from yaml import dump

TOOLS = ['star_fusion', 'fusioncatcher', 'ericscript', 'pizzly', 'squid']
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

def squid(p_file):
    fusions = ''
    next(p_file)    # skip header
    for line in p_file:
        tmp = line.rstrip('\n').split('\t')[11].split(':')
        if len(tmp) == 2:
            fusions += fi_format(tmp[0], tmp[1])

    return fusions

def pizzly(p_file):
    fusions = ''
    data = json.load(p_file)
    for item in data['genes']:
        fusions += fi_format(item['geneA']['name'], item['geneB']['name'])

    return fusions

def ericscript(p_file):
    fusions = ''
    next(p_file)    # skip header
    for line in p_file:
        tmp = line.split('\t')
        fusions += fi_format(tmp[0], tmp[1])

    return fusions

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
        exit('Defined {} doesn\'t exist'.format(p_input))

    if not os.path.exists(OUTPUT):
        with open(OUTPUT, 'w'):
            pass

    if p_tool not in TOOLS:
        exit('Defined {} not in the supported list of transformations!'.format(p_tool))

    try:
        fusions = ''
        with open(p_input, 'r') as in_file:
            func = getattr(sys.modules[__name__], p_tool)   # get function from parameter
            fusions = func(in_file).rstrip()  # call function
            in_file.close()

        with open(p_output, 'r+') as out_file, open(SUMMARY, 'a') as summary:
            if fusions:
                save(out_file, fusions)
                summary.write(
                    dump(
                        {p_tool : list(set(fusions.split('\n')))},
                        default_flow_style=False,
                        allow_unicode=True)
                )
            else:
                summary.write(
                    dump(
                        {p_tool: None},
                        default_flow_style=False,
                        allow_unicode=True
                    )
                )

            out_file.close()
            summary.close()
    except IOError as error:
        print(error)

def main():
    parser = argparse.ArgumentParser(
        description='Tool for extracting fusion-genes from various fusion detection tools'
    )
    parser.add_argument(
        '-i', '--input',
        help='Output file of a fusion tool',
        type=str,
        required=True
    )
    parser.add_argument(
        '-t', '--tool',
        help='',
        choices='|'.join(TOOLS),
        type=str,
        required=True
    )
    args = parser.parse_args()
    transform(args.input, args.tool, OUTPUT)

if __name__ == "__main__":
    main()
