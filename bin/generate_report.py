#!/usr/bin/env python3

import argparse
import yaml
from summary.db import Db
from summary.page import Page
from summary.report import Report

def ppi_network_graph(p_data):
    graph_data = []
    if len(p_data) < 1:
        return graph_data

    if 'h_gene' in p_data[0].keys() and 't_gene' not in p_data[0].keys():
        graph_data.append({
            'id': 1,
            'name': p_data[0]['h_gene'],
        })
        for idx, value in enumerate(p_data[0]['h_gene_interactions'].split(',')):
            graph_data.append({
                'id': idx + 1,
                'name': value
            })

    if 'h_gene' not in p_data[0].keys() and 't_gene' in p_data[0].keys():
        graph_data.append({
            'id': 1,
            'name': p_data[0]['t_gene'],
        })
        for idx, value in enumerate(p_data[0]['t_gene_interactions'].split(',')):
            graph_data.append({
                'id': idx + 1,
                'name': value
            })

    if 'h_gene' in p_data[0].keys() and 't_gene' in p_data[0].keys():
        graph_data.append({
            'id': 1,
            'name': 'fusion',
        })

        graph_data.append({
            'id': 2,
            'name': p_data[0]['h_gene'],
            'parent': 1
        })
        graph_data.append({
            'id': 3,
            'name': p_data[0]['t_gene'],
            'parent': 1
        })
        counter = 4
        left_fusion = set(map(str.strip, p_data[0]['h_gene_interactions'].split(',')))
        right_fusion = set(map(str.strip, p_data[0]['t_gene_interactions'].split(',')))
        intersect = left_fusion & right_fusion
        left_fusion -=  intersect
        right_fusion -= intersect

        for gene in left_fusion:
            graph_data.append({
                'id': counter,
                'name': gene,
                'parent': 2
            })
            counter += 1

        for gene in right_fusion:
            graph_data.append({
                'id': counter,
                'name': gene,
                'parent': 3
            })
            counter += 1

        for gene in list(intersect):
            graph_data.append({
                'id': counter,
                'name': gene,
                'parent': 1
            })
            counter += 1

    return graph_data

def parse_summary(p_summary):
    try:
        with open(p_summary, 'r') as in_file:
            return yaml.safe_load(in_file.read())
    except IOError:
        exit('File ' + p_summary + ' was not found!')

def build_fusion_summary(p_summary, p_known_fusions, p_unknown_fusions):
    fusions = {}
    all_fusions = p_known_fusions + p_unknown_fusions
    for fusion in all_fusions:
        if fusion not in fusions.keys():
            tools = {}
            tmp_split = fusion.split('--')
            for tool in p_summary.keys():
                tools[tool] = True if tmp_split[0] in p_summary[tool] and p_summary[tool][tmp_split[0]] == tmp_split[1] else False
            fusions[fusion] = {
                'left_gene': tmp_split[0],
                'right_gene': tmp_split[1],
                'known': fusion in p_known_fusions,
                'tools': tools,
                'tools_total': sum(tools.values())
            }

    return { 'fusions': fusions, 'tools': p_summary.keys() }

def generate(p_fusion_list, p_summary, p_config, p_output):
    if p_fusion_list is None or p_summary is None:
        exit('Fusion list or summary was not provided')
    
    db = Db()
    known_fusions = []
    unknown_fusions = []
    summary_file = parse_summary(p_summary)
    report = Report(p_config, p_output)

    # Get all fusions from DB
    db_fusions = db.select('SELECT DISTINCT (h_gene || "--" || t_gene) as fusion_pair FROM TCGA_ChiTaRS_combined_fusion_information_on_hg19')
    db_fusions = [x['fusion_pair'] for x in db_fusions]

    # Create page per fusion
    with open(p_fusion_list, 'r') as fusions:
        for fusion_pair in fusions:
            fusion_pair = fusion_pair.rstrip()

            if fusion_pair not in db_fusions:
                unknown_fusions.append(fusion_pair)
                continue

            known_fusions.append(fusion_pair)
            fusion_page = Page(db, fusion_pair, 'fusion')
            fusion = fusion_pair.split('--')
            fusion_page.add_section(
                'variations',
                'Fusion gene variations',
                'SELECT * FROM TCGA_ChiTaRS_combined_fusion_information_on_hg19 WHERE h_gene = ? AND t_gene = ?',
                fusion
            )
            fusion_page.add_section(
                'transcripts',
                'Ensembl transcripts',
                'SELECT * FROM TCGA_ChiTaRS_combined_fusion_ORF_analyzed_gencode_h19v19 WHERE h_gene = ? AND t_gene = ?',
                fusion
            )
            fusion_page.add_section(
                'ppi',
                'Protein-Protein interaction',
                'SELECT DISTINCT h_gene, h_gene_interactions, t_gene, t_gene_interactions FROM fusion_ppi WHERE h_gene = ? AND t_gene = ?',
                fusion
            )
            graph = ppi_network_graph(fusion_page.get_section('ppi')['data'])
            fusion_page.add_graph('ppi', 'ppi_graph', graph)
            fusion_page.add_section(
                'targeting_drugs',
                'Targeting drugs',
                '''
                SELECT gene_symbol, drug_status, drug_bank_id, drug_name, drug_action, fusion_uniprot_related_drugs.uniprot_acc FROM fusion_uniprot_related_drugs
                INNER JOIN uniprot_gsymbol ON fusion_uniprot_related_drugs.uniprot_acc = uniprot_gsymbol.uniprot_acc
                WHERE gene_symbol = ? OR gene_symbol = ?
                ''',
                fusion
            )
            fusion_page.add_section(
                'related_diseases',
                'Related diseases',
                'SELECT * FROM fgene_disease_associations WHERE (gene = ? OR gene = ?) AND disease_prob > 0.2001 ORDER BY disease_prob DESC',
                fusion
            )
            report.add_page(fusion_page)

    # Create index page
    fusions_summary = build_fusion_summary(summary_file, known_fusions, unknown_fusions)
    index_page = Page(db, 'index', 'index')
    index_page.add_custom_section(
        'intro',
        'Dashboard fusion summary',
        {}
    )
    index_page.add_graph('intro', 'tool_detection_chart', [[key, len(value)] for key, value in summary_file.items()])
    index_page.add_graph('intro', 'known_unknown_chart', [['known', len(known_fusions)], ['unknown', len(unknown_fusions)]])
    index_page.add_custom_section(
        'summary_fusions',
        'List of all found fusions',
        fusions_summary
    )
    
    report.add_page(index_page)
    # Render final report
    report.render()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Utility for generating HTML summary report""")
    parser.add_argument('fusions', help='Input fusion file', type=str)
    parser.add_argument('summary', help='Input fusion file', type=str)
    parser.add_argument('-c', '--config', nargs='?', help='Input config file', type=str, required=False)
    parser.add_argument('-o', '--output', nargs='?', help='Output directory', type=str, required=True)
    args = parser.parse_args()

    generate(args.fusions, args.summary, args.config, args.output)