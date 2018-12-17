#!/usr/bin/env python3

import argparse
import yaml
from summary.db import Db
from summary.page import Page
from summary.report import Report
from summary.section import Section
from summary.graph import Graph

def parse_summary(p_summary):
    try:
        with open(p_summary, 'r') as in_file:
            return yaml.safe_load(in_file.read())
    except IOError:
        exit('File ' + p_summary + ' was not found!')

def create_tool_detection_chart(p_summary):
    result = []
    all_fusions = []
    for tool, fusions in p_summary.items():
        all_fusions.append(fusions)
        result.append([tool, len(fusions)])
    result.append(['all tools', len(set.intersection(*map(set, all_fusions)))])

    return result 

def create_ppi_graph(p_data):
    graph_data = []
    if len(p_data) < 1:
        return graph_data

    graph_data = [
        {'data': {'id': 'fusion' }, 'classes': 'core'},
        {'data': {'id': p_data[0]['h_gene']}, 'classes': 'core'},
        {'data': {'id': p_data[0]['t_gene']}, 'classes': 'core'},
        {'data': {
            'id': 'fusion' + p_data[0]['h_gene'],
            'source': 'fusion',
            'target': p_data[0]['h_gene']
            },
            'classes': 'core-connection'
        },
        {'data': {
            'id': 'fusion' + p_data[0]['t_gene'],
            'source': 'fusion',
            'target': p_data[0]['t_gene']
            },
            'classes': 'core-connection'
        },
    ]

    left_fusion = set(map(str.strip, p_data[0]['h_gene_interactions'].split(',')))
    right_fusion = set(map(str.strip, p_data[0]['t_gene_interactions'].split(',')))
    intersect = left_fusion & right_fusion
    left_fusion -=  intersect
    right_fusion -= intersect

    for gene in left_fusion:
        graph_data.append({'data': {'id': gene }})
        graph_data.append({
            'data': {
               'id': gene + '--' + p_data[0]['h_gene'],
               'source': p_data[0]['h_gene'],
               'target': gene
            }
        })

    for gene in right_fusion:
        graph_data.append({'data': {'id': gene}})
        graph_data.append({
            'data': {
               'id': gene + '--' + p_data[0]['t_gene'],
               'source': p_data[0]['t_gene'],
               'target': gene
            }
        })

    for gene in list(intersect):
        graph_data.append({'data': {'id': gene}})
        graph_data.append({
            'data': {
               'id': 'fusion' + '--' + gene,
               'source': 'fusion',
               'target': gene
            }
        })

    return graph_data

def create_distribution_chart(p_summary):
    graph_data = [set() for i in range(len(p_summary.keys()))]
    all_fusions = [fusions for _, fusions in p_summary.items()]
    all_fusions = sum(all_fusions, [])
    for fusion in all_fusions:
        index = all_fusions.count(fusion)
        graph_data[index - 1].add(fusion)

    return [[str(index + 1) + ' tool/s', len(value)] for index, value in enumerate(graph_data)]

def create_fusions_table(p_summary, p_known_fusions):
    fusions = {}
    all_fusions = [fusions for _, fusions in p_summary.items()]
    unique_fusions = set(sum(all_fusions, []))

    for fusion in unique_fusions:
        tools = [fusion in x for x in all_fusions]
        fusions[fusion] = {
            'known': fusion in p_known_fusions,
            'tools': tools,
            'tools_total': sum(tools)
        }

    return { 'fusions': fusions, 'tools': p_summary.keys() }

def generate(p_args):
    if p_args.fusions is None or p_args.summary is None:
        exit('Fusion list or summary was not provided')
    
    db = Db()
    known_fusions = []
    unknown_fusions = []
    summary_file = parse_summary(p_args.summary)
    report = Report(p_args.config, p_args.output)

    # Get all fusions from DB
    db_fusions = db.select('SELECT DISTINCT (h_gene || "--" || t_gene) as fusion_pair FROM TCGA_ChiTaRS_combined_fusion_information_on_hg19')
    db_fusions = [x['fusion_pair'] for x in db_fusions]

    # Create page per fusion
    with open(p_args.fusions, 'r') as fusions:
        for fusion_pair in fusions:
            fusion_pair = fusion_pair.rstrip()

            if fusion_pair not in db_fusions:
                unknown_fusions.append(fusion_pair)
                continue

            known_fusions.append(fusion_pair)
            fusion_page_variables = {
                'sample': p_args.sample
            }
            fusion_page = Page(fusion_pair, fusion_page_variables, 'fusion')
            fusion = fusion_pair.split('--')

            variations_section = Section()
            variations_section.id = 'variations'
            variations_section.title = 'Fusion gene variations'
            variations_section.content = '''
            Fusion gene information taken from three different sources ChiTars (NAR, 2018), tumorfusions (NAR, 2018) and Gao et al. (Cell, 2018).
            Genome coordinates are lifted-over GRCh37/hg19 version. <br>Note: LD (Li Ding group, RV: Roel Verhaak group, ChiTaRs fusion database).
            '''
            variations_section.set_query(
                db,
                'SELECT * FROM TCGA_ChiTaRS_combined_fusion_information_on_hg19 WHERE h_gene = ? AND t_gene = ?',
                fusion
            )
            fusion_page.add_section(variations_section)
            
            transcripts_section = Section()
            transcripts_section.id = 'transcripts'
            transcripts_section.title = 'Ensembl transcripts'
            transcripts_section.content = '''Open reading frame (ORF) analsis of fusion genes based on Ensembl gene isoform structure.'''
            transcripts_section.set_query(
                db,
                'SELECT * FROM TCGA_ChiTaRS_combined_fusion_ORF_analyzed_gencode_h19v19 WHERE h_gene = ? AND t_gene = ?',
                fusion
            )
            fusion_page.add_section(transcripts_section)

            ppi_section = Section()
            ppi_section.id = 'ppi'
            ppi_section.title = 'Chimeric Protein-Protein interactions'
            ppi_section.content = '''
            Protein-protein interactors with each fusion partner protein in wild-type.
            Data are taken from <a href="http://chippi.md.biu.ac.il/index.html">here</a>
            '''
            ppi_section.set_query(
                db,
                'SELECT DISTINCT h_gene, h_gene_interactions, t_gene, t_gene_interactions FROM fusion_ppi WHERE h_gene = ? AND t_gene = ?',
                fusion
            )
            ppi_graph = Graph(
                'ppi_graph',
                'Network graph of gene interactions',
                '',
                create_ppi_graph(ppi_section.data))
            ppi_section.add_graph(ppi_graph)
            fusion_page.add_section(ppi_section)

            drugs_section = Section()
            drugs_section.id = 'targeting_drugs'
            drugs_section.title = 'Targeting drugs'
            drugs_section.content = '''Drugs targeting genes involved in this fusion gene (DrugBank Version 5.1.0 2018-04-02).'''
            drugs_section.set_query(
                db,
                '''
                SELECT gene_symbol, drug_status, drug_bank_id, drug_name, drug_action, fusion_uniprot_related_drugs.uniprot_acc FROM fusion_uniprot_related_drugs
                INNER JOIN uniprot_gsymbol ON fusion_uniprot_related_drugs.uniprot_acc = uniprot_gsymbol.uniprot_acc
                WHERE gene_symbol = ? OR gene_symbol = ?
                ''',
                fusion
            )
            fusion_page.add_section(drugs_section)

            diseases_section = Section()
            diseases_section.id = 'related_diseases'
            diseases_section.title = 'Related diseases'
            diseases_section.content = '''Diseases associated with fusion partners (DisGeNet 4.0).'''
            diseases_section.set_query(
                db,
                'SELECT * FROM fgene_disease_associations WHERE (gene = ? OR gene = ?) AND disease_prob > 0.2001 ORDER BY disease_prob DESC',
                fusion
            )
            fusion_page.add_section(diseases_section)
            report.add_page(fusion_page)

    index_page_variables = {
        'sample': p_args.sample,
        'total_fusions': len(unknown_fusions) + len(known_fusions),
        'known_fusions': len(known_fusions),
        'tools': summary_file.keys()
    }
    index_page = Page('index', index_page_variables, 'index')

    dashboard_section = Section()
    dashboard_section.id = 'dashboard'
    dashboard_section.title = 'Dashboard fusion summary'
    dashboard_graph1 = Graph(
        'tool_detection_chart',
        'Tool detection',
        'Displays number of found fusions per tool.',
        create_tool_detection_chart(summary_file)
    )
    dashboard_graph2 = Graph(
        'known_unknown_chart',
        'Known Vs Unknown',
        'Shows the ration between found and unknown missing fusions in the local database.',
        [['known', len(known_fusions)], ['unknown', len(unknown_fusions)]])
    dashboard_graph3 = Graph(
        'distribution_chart',
        'Tool detection distribution',
        'Sum of counts detected by different tools per fusion.',
        create_distribution_chart(summary_file))
    dashboard_section.add_graph(dashboard_graph1)
    dashboard_section.add_graph(dashboard_graph2)
    dashboard_section.add_graph(dashboard_graph3)
    index_page.add_section(dashboard_section)
    
    fusion_list_section = Section()
    fusion_list_section.id = 'fusion_list'
    fusion_list_section.title = 'List of all detected fusions'
    fusion_list_section.data = create_fusions_table(summary_file, known_fusions)
    index_page.add_section(fusion_list_section)
    report.add_page(index_page)
    
    # Render final report
    report.render()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Utility for generating HTML summary report""")
    parser.add_argument('fusions', help='Input fusion file', type=str)
    parser.add_argument('summary', help='Input fusion file', type=str)
    parser.add_argument('-c', '--config', nargs='?', help='Input config file', type=str, required=False)
    parser.add_argument('-s', '--sample', nargs='?', help='Sample name', type=str, required=True)
    parser.add_argument('-o', '--output', nargs='?', help='Output directory', type=str, required=True)
    args = parser.parse_args()

    generate(args)