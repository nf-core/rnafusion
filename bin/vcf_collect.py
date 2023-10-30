#!/usr/bin/env python

import argparse
import logging
import sys
from pathlib import Path
import pandas as pd
import ast
from gtfparse import read_gtf

logger = logging.getLogger()


def vcf_collect(fusioninspector_in_file, fusionreport_in_file, gtf, hgnc, sample, out):
    """
    Process FusionInspector and FusionReport data,
    merge with GTF from FusionInspector and HGNC database,
    and write a VCF file.

    Args:
        fusioninspector_in_file (str): Path to FusionInspector input file.
        fusionreport_in_file (str): Path to FusionReport input file.
        sample (str): Sample name for the header.
        hgnc (str): Path to HGNC file.
        gtf (str): Path to GTF file.
        out (str): Output VCF file path.

    Adapted from: https://github.com/J35P312/MegaFusion
    """
    merged_df = build_fusioninspector_dataframe(fusioninspector_in_file).join(
        read_build_fusionreport(fusionreport_in_file), how="outer", on='FUSION').reset_index()

    df = build_hgnc_dataframe(hgnc).merge(merged_df, how='right', left_on='ensembl_gene_id',
                                          right_on='Left_ensembl_gene_id')
    df = df.rename(columns={"hgnc_id": "Left_hgnc_id"})
    df = build_hgnc_dataframe(hgnc).merge(df, how='right', left_on='ensembl_gene_id', right_on='Right_ensembl_gene_id')
    df = df.rename(columns={"hgnc_id": "Right_hgnc_id"})
    gtf_df = build_gtf_dataframe(gtf)
    all_df = df.merge(gtf_df, how='left', left_on='CDS_LEFT_ID', right_on='Transcript_id')
    all_df = all_df[(all_df['PosA'] >= all_df['orig_start']) & (all_df['PosA'] <= all_df['orig_end'])]
    all_df = all_df.rename(columns={"transcript_version": "Left_transcript_version"})
    all_df = all_df.rename(columns={"exon_number": "Left_exon_number"})
    all_df = all_df[
        ['FUSION', 'GeneA', 'GeneB', 'PosA', 'PosB', 'ChromosomeA', 'ChromosomeB', 'TOOLS_HITS', 'SCORE', 'FOUND_DB',
         'FOUND_IN', 'JunctionReadCount', 'SpanningFragCount', 'FFPM', 'PROT_FUSION_TYPE', 'CDS_LEFT_ID',
         'CDS_RIGHT_ID', 'Left_transcript_version', 'Left_exon_number', 'Left_hgnc_id', 'Right_hgnc_id', 'Strand1',
         'Strand2', 'annots']].drop_duplicates()
    all_df = all_df.merge(gtf_df, how='left', left_on='CDS_RIGHT_ID', right_on='Transcript_id')
    all_df = all_df[(all_df['PosB'] >= all_df['orig_start']) & (all_df['PosB'] <= all_df['orig_end'])]
    all_df = all_df.rename(columns={"transcript_version": "Right_transcript_version"})
    all_df = all_df.rename(columns={"exon_number": "Right_exon_number"})
    all_df = all_df[
        ['FUSION', 'GeneA', 'GeneB', 'PosA', 'PosB', 'ChromosomeA', 'ChromosomeB', 'TOOLS_HITS', 'SCORE', 'FOUND_DB',
         'FOUND_IN', 'JunctionReadCount', 'SpanningFragCount', 'FFPM', 'PROT_FUSION_TYPE', 'CDS_LEFT_ID',
         'CDS_RIGHT_ID', 'Left_transcript_version', 'Left_exon_number', 'Left_hgnc_id', 'Right_transcript_version',
         'Right_exon_number', 'Right_hgnc_id', 'Strand1', 'Strand2', 'annots']].drop_duplicates()

    return write_vcf(column_manipulation(all_df), header_def(sample), out)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "--fusioninspector",
        metavar="FUSIONINSPECTOR",
        type=Path,
        help="FusionInspector output in TSV format.",
    )
    parser.add_argument(
        "--fusionreport",
        metavar="FUSIONREPORT",
        type=Path,
        help="Fusionreport output in TSV format.",
    )
    parser.add_argument(
        "--fusioninspector_gtf",
        metavar="GTF",
        type=Path,
        help="FusionInspector GTF output.",
    )
    parser.add_argument(
        "--hgnc",
        metavar="HGNC",
        type=Path,
        help="HGNC database.",
    )
    parser.add_argument("--sample", metavar="SAMPLE", type=Path, help="Sample name.", default="Sample")
    parser.add_argument(
        "--out",
        metavar="OUT",
        type=Path,
        help="VCF output path.",
    )
    return parser.parse_args(argv)


def header_def(sample):
    """
    Define the header of the VCF file
    """
    return '##fileformat=VCFv4.1\n\
##ALT=<ID=BND,Description="Break end">\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n\
##INFO=<ID=CHRA,Number=1,Type=String,Description="Chromosome A">\n\
##INFO=<ID=CHRB,Number=1,Type=String,Description="Chromosome B">\n\
##INFO=<ID=GENEA,Number=.,Type=String,Description="Gene A">\n\
##INFO=<ID=GENEB,Number=.,Type=String,Description="Gene B">\n\
##INFO=<ID=POSA,Number=.,Type=String,Description="Breakpoint position A">\n\
##INFO=<ID=POSB,Number=.,Type=String,Description="Breakpoint position B">\n\
##INFO=<ID=ORIENTATION,Number=.,Type=String,Description="Strand1 and strand2 directions">\n\
##INFO=<ID=FOUND_DB,Number=.,Type=String,Description="Databases in which the fusion has been found">\n\
##INFO=<ID=FOUND_IN,Number=.,Type=String,Description="Callers that have found the fusion">\n\
##INFO=<ID=TOOL_HITS,Number=.,Type=Integer,Description="Number of tools that found the fusion">\n\
##INFO=<ID=SCORE,Number=.,Type=Float,Description="Score from fusionreport">\n\
##INFO=<ID=FRAME_STATUS,Number=.,Type=Float,Description="Frame status of the fusion">\n\
##INFO=<ID=TRANSCRIPT_ID_A,Number=.,Type=Float,Description="Transcript id A ">\n\
##INFO=<ID=TRANSCRIPT_ID_B,Number=.,Type=Float,Description="Transcript id B">\n\
##INFO=<ID=TRANSCRIPT_VERSION_A,Number=.,Type=Float,Description="Transcript version A">\n\
##INFO=<ID=TRANSCRIPT_VERSION_B,Number=.,Type=Float,Description="Transcript version B">\n\
##INFO=<ID=HGNC_ID_A,Number=.,Type=Float,Description="HGNC id A">\n\
##INFO=<ID=HGNC_ID_B,Number=.,Type=Float,Description="HGNC id A">\n\
##INFO=<ID=EXON_NUMBER_A,Number=.,Type=Float,Description="Exon number A">\n\
##INFO=<ID=EXON_NUMBER_B,Number=.,Type=Float,Description="Exon number B">\n\
##INFO=<ID=ANNOTATIONS,Number=.,Type=Float,Description="Annotations from FusionInspector">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n\
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of paired-ends that support the event">\n\
##FORMAT=<ID=RV,Number=1,Type=Integer,Description="Number of split reads that support the event">\n\
##FORMAT=<ID=FFPM,Number=1,Type=Float,Description="Fusion fragments per million total RNA-seq fragments">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}'.format(
        sample
    )


def build_fusioninspector_dataframe(file):
    """
    Read FusionInspector output from a CSV file, preprocess the data, and set 'FUSION' as the index.
    """
    df = pd.read_csv(file, sep="\t")
    df = df.rename(columns={"#FusionName": "FUSION"})
    df[['ChromosomeA', 'PosA', 'Strand1']] = df['LeftBreakpoint'].str.split(':', expand=True)
    df[['ChromosomeB', 'PosB', 'Strand2']] = df['RightBreakpoint'].str.split(':', expand=True)
    df[['GeneA', 'GeneB']] = df['FUSION'].str.split('--', expand=True)
    df[['LeftGeneName', 'Left_ensembl_gene_id']] = df['LeftGene'].str.split('^', expand=True)
    df[['RightGeneName', 'Right_ensembl_gene_id']] = df['RightGene'].str.split('^', expand=True)
    return df.set_index(['FUSION'])


def replace_value_with_column_name(row, value_to_replace, column_name):
    """
    Replace a specific value in a row with the corresponding column name.
    """
    new_values = ''
    for col_name, value in row.items():
        if col_name == column_name:
            if value == value_to_replace:
                new_values = col_name
            else:
                new_values = ''
    return new_values


def concatenate_columns(row):
    """
    Concatenate non-empty values in a row into a single string separated by commas.
    """
    non_empty_values = [str(value) for value in row if value != '']
    return ','.join(non_empty_values)


def read_build_fusionreport(fusionreport_file):
    """
    Read and preprocess fusion-report data from a file, including handling missing tool columns,
    getting the columns with each tool and create a new FOUND_IN column with all the tool hits.
    Convert the list of databases in FOUND_DB into a joined string with a comma separator.
    Make all column headers uppercase.
    """
    with open(fusionreport_file) as f:
        from_html = [line.split('rows": [')[1] for line in f if 'name="fusion_list' in line]
        expression = from_html[0].split('], "tool')[0]
    fusion_report = pd.DataFrame.from_dict(ast.literal_eval(expression))
    if "arriba" not in fusion_report.columns:
        fusion_report["arriba"] = ""
    if "fusioncatcher" not in fusion_report.columns:
        fusion_report["fusioncatcher"] = ""
    if "starfusion" not in fusion_report.columns:
        fusion_report["starfusion"] = ""
    fusion_report['arriba'] = fusion_report[['arriba']].apply(replace_value_with_column_name,
                                                              args=('true', 'arriba'), axis=1)
    fusion_report['fusioncatcher'] = fusion_report[['fusioncatcher']].apply(replace_value_with_column_name,
                                                                            args=('true', 'fusioncatcher'), axis=1)
    fusion_report['starfusion'] = fusion_report[['starfusion']].apply(replace_value_with_column_name,
                                                                      args=('true', 'starfusion'), axis=1)
    fusion_report['FOUND_IN'] = fusion_report[['arriba', 'starfusion',
                                               'fusioncatcher']].apply(concatenate_columns, axis=1)
    fusion_report.columns = fusion_report.columns.str.upper()
    fusion_report['FOUND_DB'] = fusion_report['FOUND_DB'].apply(lambda x: ', '.join(x))
    return fusion_report[['FUSION', 'TOOLS_HITS', 'SCORE', 'FOUND_DB', 'FOUND_IN']].set_index(['FUSION'])


def column_manipulation(df):
    """
    Manipulate and prepare DataFrame for VCF file creation.
    """
    df["ALT"] = ""
    df = df.reset_index()
    df["FORMAT"] = "GT:DV:RV:FFPM"
    df["ID"] = "."
    df["QUAL"] = "."
    df["FILTER"] = "PASS"
    df["REF"] = "N"
    df["INFO"] = ""
    df["Sample"] = ""

    for index, row in df.iterrows():
        # ALT
        if not row["Strand1"] in ["+", "-"] or not row["Strand2"] in ["+", "-"]:
            df.loc[index, "ALT"] = "N[{}:{}[".format(df["ChromosomeB"], row["PosB"])
        elif row["Strand1"] == "-" and row["Strand2"] == "-":
            df.loc[index, "ALT"] = "[{}:{}[N".format(row["ChromosomeB"], row["PosB"])
        elif row["Strand1"] == "+" and row["Strand2"] == "-":
            df.loc[index, "ALT"] = "N]{}:{}]".format(row["ChromosomeB"], row["PosB"])
        elif row["Strand1"] == "-" and row["Strand2"] == "+":
            df.loc[index, "ALT"] = "N]{}:{}]".format(row["ChromosomeB"], row["PosB"])
        else:
            df.loc[index, "ALT"] = "N[{}:{}[".format(row["ChromosomeB"], row["PosB"])
        # INFO
        df.loc[index, "INFO"] = (
            "SVTYPE=BND;CHRA={};CHRB={};GENEA={};GENEB={};POSA={};POSB={};ORIENTATION={},{};FOUND_DB={};"
            "FOUND_IN={};;TOOL_HITS={};SCORE={};FRAME_STATUS={};TRANSCRIPT_ID_A={};TRANSCRIPT_ID_B={};"
            "TRANSCRIPT_VERSION_A={};TRANSCRIPT_VERSION_B={};HGNC_ID_A={};HGNC_ID_B={};EXON_NUMBER_A={};"
            "EXON_NUMBER_B={};ANNOTATIONS={}".format(
                row["ChromosomeA"],
                row["ChromosomeB"],
                row["GeneA"],
                row["GeneB"],
                row['PosA'],
                row['PosB'],
                row["Strand1"],
                row["Strand2"],
                row["FOUND_DB"],
                row["FOUND_IN"],
                row["TOOLS_HITS"],
                row["SCORE"],
                row["PROT_FUSION_TYPE"],
                row["CDS_LEFT_ID"],
                row["CDS_RIGHT_ID"],
                row["Left_transcript_version"],
                row["Right_transcript_version"],
                row["Left_hgnc_id"],
                row["Right_hgnc_id"],
                row["Left_exon_number"],
                row["Right_exon_number"],
                row["annots"],
            )
        )
        df.loc[index, "Sample"] = "./1:{}:{}:{}".format(row["JunctionReadCount"], row["SpanningFragCount"], row["FFPM"])
    return df


def write_vcf(df_to_print, header, out_file):
    """
    Write a VCF file with a specified DataFrame, header, and output file path.
    """
    df_to_print[
        [
            "ChromosomeA",
            "PosA",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "Sample",
        ]
    ].to_csv(
        path_or_buf=out_file,
        sep="\t",
        header=None,
        index=False,
    )

    with open(out_file, "r+") as f:
        content = f.read()
        f.seek(0, 0)
        f.write(header.rstrip("\r\n") + "\n" + content)


def build_hgnc_dataframe(file):
    """
    Build a DataFrame from HGNC input file, extracting 'hgnc_id' and 'ensembl_gene_id' columns.
    """
    df = pd.read_csv(file, sep="\t", low_memory=False)
    return df[['hgnc_id', 'ensembl_gene_id']].dropna()


def build_gtf_dataframe(file):
    """
    Build a DataFrame from GTF file, extracting relevant columns.
    """
    df = read_gtf(file)
    df[['fusion_dump', 'Transcript_id']] = df['transcript_id'].str.split('^', expand=True)
    df[['orig_chromosome', 'orig_start', 'orig_end', 'orig_dir']] = df['orig_coord_info'].str.split(',', expand=True)
#     return df
    return df[['Transcript_id', 'transcript_version', 'exon_number', 'exon_id', 'orig_start', 'orig_end']]


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    if not args.fusioninspector.is_file() or not args.fusionreport.is_file():
        logger.error(f"The given input file {args.fusioninspector} or {args.fusionreport} was not found!")
        sys.exit(2)
    vcf_collect(args.fusioninspector, args.fusionreport, args.fusioninspector_gtf, args.hgnc, args.sample, args.out)


if __name__ == "__main__":
    sys.exit(main())
