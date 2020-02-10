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
    'Arriba': ['v_arriba.txt', r"arriba=(\S+)"],
    'EricScript': ['v_ericscript.txt', r"ericscript=(\S+)"],
    'FusionCatcher': ['v_fusioncatcher.txt', r"fusioncatcher=(\S+)"],
    'Fusion-Inspector': ['v_fusion_inspector.txt', r"fusion-inspector=(\S+)"],
    'fusion-report': ['v_fusion_report.txt', r"fusion-report=(\S+)"],
    'Pizzly': ['v_pizzly.txt', r"pizzly=(\S+)"],
    'STAR-Fusion': ['v_star_fusion.txt', r"star-fusion=(\S+)"],
    'Squid': ['v_squid.txt', r"squid=(\S+)"]
}
results = OrderedDict()
results['nf-core/rnafusion'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'
results['Arriba'] = '<span style="color:#999999;\">N/A</span>'
results['EricScript'] = '<span style="color:#999999;\">N/A</span>'
results['FusionCatcher'] = '<span style="color:#999999;\">N/A</span>'
results['Fusion-Inspector'] = '<span style="color:#999999;\">N/A</span>'
results['fusion-report'] = '<span style="color:#999999;\">N/A</span>'
results['Pizzly'] = '<span style="color:#999999;\">N/A</span>'
results['STAR-Fusion'] = '<span style="color:#999999;\">N/A</span>'
results['Squid'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del(results[k])

# Remove software set to false in results
for k in results:
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'nf-core/rnafusion Software Versions'
section_href: 'https://github.com/nf-core/rnafusion'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k,v in results.items():
        f.write("{}\t{}\n".format(k,v))
