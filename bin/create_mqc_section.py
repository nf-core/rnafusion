#!/usr/bin/env python3
from collections import OrderedDict
from yaml import dump
import argparse
import yaml
import sys
import os

OUTPUT = 'fusion_genes_config_mqc.yaml'
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

def findings(p_yaml, p_sample_name):
    template = TEMPLATE
    result = {}

    if p_yaml is None:
        return

    # Counts per tool
    for tool, fusions in p_yaml.items():
        result[tool] = len(fusions) if fusions is not None else 0

    # If only one tool was found, there is no need to make intercept
    if len(result) == 1:
        template['data'] = { p_sample_name: result }
        return OrderedDict(template)

    # Intersect
    result['together'] = len(set.intersection(*map(set, [fusions for _, fusions in p_yaml.items()])))
    
    # Group results
    template['data'] = { p_sample_name: result }
    return OrderedDict(template)

def summary(p_input, p_sample_name):
    if not os.path.exists(p_input):
        sys.exit('Defined {} doesn\'t exist'.format(p_input))
    try:
        with open(p_input, 'r') as stream, open(OUTPUT, 'w') as out_file:
            yaml_data = yaml.safe_load(stream)
            # Conversion to nice yaml file
            yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))
            # Find and store
            out_file.write(dump(findings(yaml_data, p_sample_name), default_flow_style=False, allow_unicode=True))
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