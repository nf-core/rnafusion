#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re
import os

regexes = {
    'nf-core/rnafusion': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'STAR-Fusion': ['v_star_fusion.txt', r"=(\S+)"],
    'FusionCatcher': ['v_fusioncatcher.txt', r"\"(.*?)\""],
    'Fusion-Inspector': ['v_fusion_inspector.txt', r"=(\S+)"],
    'EricScript': ['v_ericscript.txt', r"=(\S+)"],
    'Pizzly': ['v_pizzly.txt', r"=(\S+)"],
    'Squid': ['v_squid.txt', r"=(\S+)"]
}
results = OrderedDict()
results['nf-core/rnafusion'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'
results['STAR-Fusion'] = '<span style="color:#999999;\">N/A</span>'
results['FusionCatcher'] = '<span style="color:#999999;\">N/A</span>'
results['Fusion-Inspector'] = '<span style="color:#999999;\">N/A</span>'
results['Pizzly'] = '<span style="color:#999999;\">N/A</span>'
results['Squid'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    if os.path.exists(v[0]):
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'nf-core/rnafusion-software-versions'
section_name: 'nf-core/rnafusion Software Versions'
section_href: 'https://github.com/nf-core/rnafusion'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")
