#!/usr/bin/env python
from __future__ import print_function
from __future__ import with_statement
from collections import OrderedDict
from yaml import dump
import argparse
import yaml
import sys
import os

OUTPUT = 'fusion_genes_mqc.yaml'
TEMPLATE = OrderedDict([
    ('id', 'fusion_genes'),
    ('s', 'Fusion genes'),
    ('description', 'Number of fusion genes found by various tools'),
    ('plot_type', 'bargraph'),
    ('pconfig', {
        'id': 'barplot_config_only',
        'title': 'Detected fusion genes',
        'ylab': 'Number of detected fusion genes'
    })
])

def findings(p_yaml, p_sample):
    data = TEMPLATE
    fusions = {}

    # Counts per tool
    tools = p_yaml.keys()
    for tool in tools:
        if p_yaml[tool] == None:
            fusions[tool] = 0
        else:
            fusions[tool] = len(p_yaml[tool].keys())

    # If only one tool was found, there is no need to make intercept
    if len(tools) == 1:
        data['data'] = { p_sample: fusions}
        return OrderedDict(data)

    # Intersect
    intersect_genes = {}
    if p_yaml[tools[0]] != None:
        for (fusion_left, fusion_right) in p_yaml[tools[0]].items():
            for tool in tools[1:]:
                if p_yaml[tool] != None:
                    # check if the fusion is not swapped
                    if (fusion_left in p_yaml[tool] and p_yaml[tool][fusion_left] == fusion_right) or (fusion_right in p_yaml[tool] and p_yaml[tool][fusion_right] == fusion_left):
                        intersect_genes[fusion_left] = fusion_right
    fusions['together'] = len(intersect_genes)
    
    # Group results
    data['data'] = { p_sample: fusions}
    return OrderedDict(data)

def summary(p_input, p_sample):
    if not os.path.exists(p_input):
        sys.exit('Defined {} doesn\'t exist'.format(p_input))
    try:
        with open(p_input, 'r') as stream, open(OUTPUT, 'w') as out_file:
            yaml_data = yaml.safe_load(stream)
            yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
            out_file.write(dump(findings(yaml_data, p_sample), default_flow_style=False, allow_unicode=True))
            stream.close()
            out_file.close()
    except IOError as error:
        sys.exit(error)
    except yaml.YAMLError as error:
        sys.exit(error)
    except Exception as error:
        sys.exit(error)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Utility for generating data structure for MultiQC""")
    parser.add_argument('-i', '--input', nargs='?', help='Input file', type=str, required=True)
    parser.add_argument('-s', '--sample', nargs='?', help='Sample name', type=str, required=True)
    args = parser.parse_args()
    summary(args.input, args.sample)