#!/usr/bin/env python3

import argparse
import logging
import sys
from pathlib import Path
import pandas as pd
import ast
import numpy as np
import csv

logger = logging.getLogger()


def vcf_collect(
    fusioninspector_in_file: str,
    fusionreport_in_file: str,
    gtf: str,
    fusionreport_csv: str,
    hgnc: str,
    sample: str,
    out_file,
) -> None:
    """
    Process FusionInspector and FusionReport data,
    merge with GTF from FusionInspector and HGNC database,
    and write a VCF file.

    Args:
        fusioninspector_in_file (str): Path to FusionInspector input file.
        fusionreport_in_file (str): Path to Fusion-report input file.
        sample (str): Sample name for the header.
        hgnc (str): Path to HGNC file.
        gtf (str): Path to output GTF file from FusionInspector in TSV format.
        fusionreport_csv (str): Path to Fusion-report CSV output file.
        out (str): Output VCF file path.

    Adapted from: https://github.com/J35P312/MegaFusion
    """
    merged_df = (
        build_fusioninspector_dataframe(fusioninspector_in_file)
        .join(read_build_fusionreport(fusionreport_in_file), how="outer", on="FUSION")
        .reset_index()
    )
    hgnc_df = build_hgnc_dataframe(hgnc)
    df_symbol = merged_df[merged_df["Left_ensembl_gene_id"].isna()]
    df_not_symbol = merged_df[merged_df["Left_ensembl_gene_id"].notna()]

    df_not_symbol = hgnc_df.merge(
        df_not_symbol,
        how="right",
        left_on="ensembl_gene_id",
        right_on="Left_ensembl_gene_id",
    )
    df_symbol = hgnc_df.merge(
        df_symbol, how="right", left_on="symbol", right_on="GeneA"
    )
    df = pd.concat([df_not_symbol, df_symbol])
    df = df.rename(columns={"hgnc_id": "Left_hgnc_id"})

    df_symbol = df[df["Right_ensembl_gene_id"].isna()]
    df_not_symbol = df[df["Right_ensembl_gene_id"].notna()]

    df_not_symbol = hgnc_df.merge(
        df_not_symbol,
        how="right",
        left_on="ensembl_gene_id",
        right_on="Right_ensembl_gene_id",
    )
    df_symbol = hgnc_df.merge(
        df_symbol, how="right", left_on="symbol", right_on="GeneB"
    )
    df = pd.concat([df_not_symbol, df_symbol])
    df = df.rename(columns={"hgnc_id": "Right_hgnc_id"})

    gtf_df = build_gtf_dataframe(gtf)
    all_df = df.merge(
        gtf_df, how="left", left_on="CDS_LEFT_ID", right_on="Transcript_id"
    )
    all_df[["PosA", "orig_start", "orig_end"]] = (
        all_df[["PosA", "orig_start", "orig_end"]].fillna(0).astype(int)
    )

    all_df = all_df[
        (
            (all_df["PosA"] >= all_df["orig_start"])
            & (all_df["PosA"] <= all_df["orig_end"])
        )
        | ((all_df["orig_start"] == 0) & (all_df["orig_end"] == 0))
    ]

    all_df.replace("", np.nan, inplace=True)
    all_df = all_df.drop_duplicates()

    all_df[["exon_number", "transcript_version"]] = all_df[
        ["exon_number", "transcript_version"]
    ].replace(0, np.nan)
    # Fill non-empty values within each group for 'exon_number' and 'transcript_version'
    all_df["exon_number"] = all_df.groupby("PosA")["exon_number"].transform(
        lambda x: x.fillna(method="ffill").fillna(method="bfill")
    )
    all_df["transcript_version"] = all_df.groupby("PosA")[
        "transcript_version"
    ].transform(lambda x: x.fillna(method="ffill").fillna(method="bfill"))

    all_df = all_df.rename(columns={"transcript_version": "Left_transcript_version"})
    all_df = all_df.rename(columns={"exon_number": "Left_exon_number"})
    all_df = all_df[
        [
            "FUSION",
            "GeneA",
            "GeneB",
            "PosA",
            "PosB",
            "ChromosomeA",
            "ChromosomeB",
            "TOOLS_HITS",
            "SCORE",
            "FOUND_DB",
            "FOUND_IN",
            "JunctionReadCount",
            "SpanningFragCount",
            "FFPM",
            "PROT_FUSION_TYPE",
            "CDS_LEFT_ID",
            "CDS_RIGHT_ID",
            "Left_transcript_version",
            "Left_exon_number",
            "Left_hgnc_id",
            "Right_hgnc_id",
            "Strand1",
            "Strand2",
            "annots",
        ]
    ].drop_duplicates()
    all_df["CDS_RIGHT_ID"] = all_df["CDS_RIGHT_ID"].astype("str")
    all_df = all_df.merge(
        gtf_df, how="left", left_on="CDS_RIGHT_ID", right_on="Transcript_id"
    )
    all_df[["PosB", "orig_start", "orig_end"]] = all_df[
        ["PosB", "orig_start", "orig_end"]
    ].fillna(0)
    all_df[["PosB", "orig_start", "orig_end"]] = all_df[
        ["PosB", "orig_start", "orig_end"]
    ].astype(int)
    all_df = all_df[
        (
            (all_df["PosB"] >= all_df["orig_start"])
            & (all_df["PosB"] <= all_df["orig_end"])
        )
        | ((all_df["orig_start"] == 0) & (all_df["orig_end"] == 0))
    ]

    all_df[["PosA", "PosB"]] = all_df[["PosA", "PosB"]].replace(0, np.nan)
    all_df = all_df.replace("", np.nan)

    all_df[["exon_number", "transcript_version"]] = all_df[
        ["exon_number", "transcript_version"]
    ].replace(0, np.nan)
    # Fill non-empty values within each group for 'exon_number' and 'transcript_version'
    all_df["exon_number"] = all_df.groupby("PosB")["exon_number"].transform(
        lambda x: x.fillna(method="ffill").fillna(method="bfill")
    )
    all_df["transcript_version"] = all_df.groupby("PosB")[
        "transcript_version"
    ].transform(lambda x: x.fillna(method="ffill").fillna(method="bfill"))

    all_df = all_df.rename(columns={"transcript_version": "Right_transcript_version"})
    all_df = all_df.rename(columns={"exon_number": "Right_exon_number"})

    all_df = all_df[
        [
            "FUSION",
            "GeneA",
            "GeneB",
            "PosA",
            "PosB",
            "ChromosomeA",
            "ChromosomeB",
            "TOOLS_HITS",
            "SCORE",
            "FOUND_DB",
            "FOUND_IN",
            "JunctionReadCount",
            "SpanningFragCount",
            "FFPM",
            "PROT_FUSION_TYPE",
            "CDS_LEFT_ID",
            "CDS_RIGHT_ID",
            "Left_transcript_version",
            "Left_exon_number",
            "Left_hgnc_id",
            "Right_transcript_version",
            "Right_exon_number",
            "Right_hgnc_id",
            "Strand1",
            "Strand2",
            "annots",
        ]
    ].drop_duplicates()
    all_df = all_df.rename(columns={"FUSION": "Fusion"})
    all_df = all_df.set_index("Fusion")

    all_df = all_df.combine_first(read_fusionreport_csv(fusionreport_csv))

    return write_vcf(column_manipulation(all_df), header_def(sample), out_file)


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
        help="Fusionreport output in index/html format.",
    )
    parser.add_argument(
        "--fusionreport_csv",
        metavar="FUSIONREPORT_CSV",
        type=Path,
        help="Fusionreport output in CSV format.",
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
    parser.add_argument(
        "--sample", metavar="SAMPLE", type=Path, help="Sample name.", default="Sample"
    )
    parser.add_argument(
        "--out",
        metavar="OUT",
        type=Path,
        help="VCF output path.",
    )
    return parser.parse_args(argv)


def header_def(sample: str) -> str:
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
##INFO=<ID=FRAME_STATUS,Number=.,Type=String,Description="Frame status of the fusion">\n\
##INFO=<ID=TRANSCRIPT_ID_A,Number=.,Type=String,Description="Transcript id A ">\n\
##INFO=<ID=TRANSCRIPT_ID_B,Number=.,Type=String,Description="Transcript id B">\n\
##INFO=<ID=TRANSCRIPT_VERSION_A,Number=.,Type=Float,Description="Transcript version A">\n\
##INFO=<ID=TRANSCRIPT_VERSION_B,Number=.,Type=Float,Description="Transcript version B">\n\
##INFO=<ID=HGNC_ID_A,Number=.,Type=Float,Description="HGNC id A">\n\
##INFO=<ID=HGNC_ID_B,Number=.,Type=Float,Description="HGNC id A">\n\
##INFO=<ID=EXON_NUMBER_A,Number=.,Type=Float,Description="Exon number A">\n\
##INFO=<ID=EXON_NUMBER_B,Number=.,Type=Float,Description="Exon number B">\n\
##INFO=<ID=ANNOTATIONS,Number=.,Type=String,Description="Annotations from FusionInspector">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n\
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of paired-ends that support the event">\n\
##FORMAT=<ID=RV,Number=1,Type=Integer,Description="Number of split reads that support the event">\n\
##FORMAT=<ID=FFPM,Number=1,Type=Float,Description="Fusion fragments per million total RNA-seq fragments">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}'.format(
        sample
    )


def convert_to_list(annots_str: str) -> list:
    try:
        return ast.literal_eval(annots_str)
    except (SyntaxError, ValueError):
        return np.nan


def build_fusioninspector_dataframe(file: str) -> pd.DataFrame:
    """
    Read FusionInspector output from a CSV file, preprocess the data, and set 'FUSION' as the index.
    """
    df = pd.read_csv(file, sep="\t")
    df = df.rename(columns={"#FusionName": "FUSION"})
    if not (df.empty):
        df[["ChromosomeA", "PosA", "Strand1"]] = df["LeftBreakpoint"].str.split(
            ":", expand=True
        )
        df[["ChromosomeB", "PosB", "Strand2"]] = df["RightBreakpoint"].str.split(
            ":", expand=True
        )
        df[["LeftGeneName", "Left_ensembl_gene_id"]] = df["LeftGene"].str.split(
            "^", expand=True
        )
        df[["RightGeneName", "Right_ensembl_gene_id"]] = df["RightGene"].str.split(
            "^", expand=True
        )
        df["annots"] = (
            df["annots"]
            .apply(convert_to_list)
            .apply(
                lambda x: (
                    ",".join(map(str, x))
                    if isinstance(x, list)
                    else str(x) if pd.notna(x) else ""
                )
            )
        )
    else:
        for i in [
            "ChromosomeA",
            "Strand1",
            "ChromosomeB",
            "Strand2",
            "LeftGeneName",
            "Left_ensembl_gene_id",
            "RightGeneName",
            "Right_ensembl_gene_id",
            "annots",
        ]:
            df[i] = ""
        for j in [
            "PosA",
            "PosB",
        ]:
            df[j] = np.nan

    return df.set_index(["FUSION"])


def replace_value_with_column_name(
    row: pd.Series, value_to_replace: str, column_name: str
) -> str:
    """
    Replace a specific value in a row with the corresponding column name.
    """
    new_values = ""
    for col_name, value in row.items():
        if col_name == column_name:
            if value == value_to_replace:
                new_values = col_name
            else:
                new_values = ""
    return new_values


def concatenate_columns(row: pd.Series) -> str:
    """
    Concatenate non-empty values in a row into a single string separated by commas.
    """
    non_empty_values = [str(value) for value in row if value != ""]
    return ",".join(non_empty_values)


def read_build_fusionreport(fusionreport_file: str) -> pd.DataFrame:
    """
    Read and preprocess fusion-report data from a file, including handling missing tool columns,
    getting the columns with each tool and create a new FOUND_IN column with all the tool hits.
    Convert the list of databases in FOUND_DB into a joined string with a comma separator.
    Make all column headers uppercase.
    """
    with open(fusionreport_file) as f:
        from_html = [
            line.split('rows": ')[1] for line in f if 'name="fusion_list' in line
        ]
        tmp = str(from_html)[2:]
        tmp2 = tmp.split(', "tools": ')[0]
        fusion_report = pd.DataFrame(ast.literal_eval(tmp2))
    if not "arriba" in fusion_report.columns:
        fusion_report["arriba"] = ""
    if not "fusioncatcher" in fusion_report.columns:
        fusion_report["fusioncatcher"] = ""
    if not "starfusion" in fusion_report.columns:
        fusion_report["starfusion"] = ""
    fusion_report["arriba"] = fusion_report[["arriba"]].apply(
        replace_value_with_column_name, args=("true", "arriba"), axis=1
    )
    fusion_report["fusioncatcher"] = fusion_report[["fusioncatcher"]].apply(
        replace_value_with_column_name, args=("true", "fusioncatcher"), axis=1
    )
    fusion_report["starfusion"] = fusion_report[["starfusion"]].apply(
        replace_value_with_column_name, args=("true", "starfusion"), axis=1
    )
    fusion_report["FOUND_IN"] = fusion_report[
        ["arriba", "starfusion", "fusioncatcher"]
    ].apply(concatenate_columns, axis=1)
    fusion_report.columns = fusion_report.columns.str.upper()
    fusion_report["FOUND_DB"] = fusion_report["FOUND_DB"].apply(
        lambda x: ",".join(x) if len(x) > 0 else ""
    )
    fusion_report[["GeneA", "GeneB"]] = fusion_report["FUSION"].str.split(
        "--", expand=True
    )

    return fusion_report[
        ["FUSION", "GeneA", "GeneB", "TOOLS_HITS", "SCORE", "FOUND_DB", "FOUND_IN"]
    ].set_index(["FUSION"])


def read_fusionreport_csv(file: str) -> pd.DataFrame:
    df = pd.read_csv(file)
    columns_to_iterate = ["starfusion", "arriba", "fusioncatcher"]
    for column in columns_to_iterate:
        if column not in df.columns:
            df[column] = ""
    df[["starfusion", "arriba", "fusioncatcher"]] = df[
        ["starfusion", "arriba", "fusioncatcher"]
    ].astype("str")
    for index, row in df.iterrows():
        for column in columns_to_iterate:
            cell_value = row[column]

            if "#" in cell_value:
                df.at[index, column] = df.at[index, column].split(",")[0]
                df.at[index, column] = df.at[index, column].replace("position: ", "")
                df.at[index, "A"] = df.at[index, column].split("#")[0]
                df.at[index, "B"] = df.at[index, column].split("#")[1]
                df.at[index, "ChromosomeA"] = df.at[index, "A"].split(":")[0]
                df.at[index, "PosA"] = df.at[index, "A"].split(":")[1]
                if "+" in df.at[index, "A"] or "-" in df.at[index, "A"]:
                    df.at[index, "StrandA"] = df.at[index, "A"].split(":")[2]
                else:
                    df.at[index, "StrandA"] = ""

                df.at[index, "ChromosomeB"] = df.at[index, "B"].split(":")[0]
                df.at[index, "PosB"] = df.at[index, "B"].split(":")[1]
                if "+" in df.at[index, "B"] or "-" in df.at[index, "B"]:
                    df.at[index, "StrandB"] = df.at[index, "B"].split(":")[2]
                else:
                    df.at[index, "StrandB"] = ""

                break
    df[["GeneA", "GeneB"]] = df["Fusion"].str.split("--", expand=True)
    df = df.set_index("Fusion")
    df.to_csv("tmp.csv")
    return df[
        [
            "GeneA",
            "GeneB",
            "ChromosomeA",
            "PosA",
            "StrandA",
            "ChromosomeB",
            "PosB",
            "StrandB",
        ]
    ]


def column_manipulation(df: pd.DataFrame) -> pd.DataFrame:
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
    df["Strand1"] = df["Strand1"].astype(str)
    df["JunctionReadCount"] = df["JunctionReadCount"].fillna(0).astype(int).astype(str)
    df["SpanningFragCount"] = df["SpanningFragCount"].fillna(0).astype(int).astype(str)
    df["FFPM"] = df["FFPM"].fillna(0).astype(float).astype(str)
    df["ChromosomeA"] = df["ChromosomeA"].fillna(0).astype(str)
    df["ChromosomeB"] = df["ChromosomeB"].fillna(0).astype(str)
    df["Left_hgnc_id"] = df["Left_hgnc_id"].fillna(0).astype(int).astype(str)
    df["Right_hgnc_id"] = df["Right_hgnc_id"].fillna(0).astype(int).astype(str)
    df["Left_exon_number"] = df["Left_exon_number"].fillna(0).astype(int).astype(str)
    df["Right_exon_number"] = df["Right_exon_number"].fillna(0).astype(int).astype(str)
    df["Left_transcript_version"] = (
        df["Left_transcript_version"].fillna(0).astype(int).astype(str)
    )
    df["Right_transcript_version"] = (
        df["Right_transcript_version"].fillna(0).astype(int).astype(str)
    )
    df["PosA"] = df["PosA"].fillna(0).astype(int).astype(str)
    df["PosB"] = df["PosB"].fillna(0).astype(int).astype(str)
    df["PROT_FUSION_TYPE"] = df["PROT_FUSION_TYPE"].replace(".", "nan")
    df["CDS_LEFT_ID"] = df["CDS_LEFT_ID"].replace(".", "nan")
    df["CDS_RIGHT_ID"] = df["CDS_RIGHT_ID"].replace(".", "nan")

    for index, row in df.iterrows():
        if row["Strand1"] == "-" and row["Strand2"] == "-":
            df.loc[index, "ALT"] = f'[{row["ChromosomeB"]}:{row["PosB"]}[N'
        elif row["Strand1"] == "+" and row["Strand2"] == "-":
            df.loc[index, "ALT"] = f'N]{row["ChromosomeB"]}:{row["PosB"]}]'
        elif row["Strand1"] == "-" and row["Strand2"] == "+":
            df.loc[index, "ALT"] = f'N]{row["ChromosomeB"]}:{row["PosB"]}]'
        else:
            df.loc[index, "ALT"] = f'N[{row["ChromosomeB"]}:{row["PosB"]}['

        df.loc[index, "INFO"] = (
            f"SVTYPE=BND;CHRA={row['ChromosomeA']};CHRB={row['ChromosomeB']};GENEA={row['GeneA']};GENEB={row['GeneB']};"
            f"POSA={row['PosA']};POSB={row['PosB']};ORIENTATION={row['Strand1']},{row['Strand2']};FOUND_DB={row['FOUND_DB']};"
            f"FOUND_IN={row['FOUND_IN']};TOOL_HITS={row['TOOLS_HITS']};SCORE={row['SCORE']};FRAME_STATUS={row['PROT_FUSION_TYPE']};"
            f"TRANSCRIPT_ID_A={row['CDS_LEFT_ID']};TRANSCRIPT_ID_B={row['CDS_RIGHT_ID']};"
            f"TRANSCRIPT_VERSION_A={row['Left_transcript_version']};TRANSCRIPT_VERSION_B={row['Right_transcript_version']};"
            f"HGNC_ID_A={row['Left_hgnc_id']};HGNC_ID_B={row['Right_hgnc_id']};"
            f"EXON_NUMBER_A={row['Left_exon_number']};EXON_NUMBER_B={row['Right_exon_number']};"
            f"ANNOTATIONS={row['annots']}"
        )
        df.loc[index, "Sample"] = (
            f"./1:{row['JunctionReadCount']}:{row['SpanningFragCount']}:{row['FFPM']}"
        )

    return df


def write_vcf(df_to_print: pd.DataFrame, header: str, out_file: str) -> None:
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
        path_or_buf=out_file, sep="\t", header=None, index=False, quoting=csv.QUOTE_NONE
    )

    with open(out_file, "r+") as f:
        content = f.read()
        f.seek(0, 0)
        f.write(header.rstrip("\r\n") + "\n" + content)


def build_hgnc_dataframe(file: str) -> pd.DataFrame:
    """
    Build a DataFrame from HGNC input file, extracting 'hgnc_id' and 'ensembl_gene_id' columns.
    """
    df = pd.read_csv(file, sep="\t", low_memory=False)
    df["hgnc_id"] = df["hgnc_id"].str.replace("HGNC:", "")
    return df[["hgnc_id", "ensembl_gene_id", "symbol"]].dropna()


def build_gtf_dataframe(file: str) -> pd.DataFrame:
    """
    Build a DataFrame from GTF file converted in TSV, extracting relevant columns.
    """
    df = pd.read_csv(file, sep="\t")
    df[["fusion_dump", "Transcript_id"]] = df["transcript_id"].str.split(
        "^", expand=True
    )
    df[["orig_chromosome", "orig_start", "orig_end", "orig_dir"]] = df[
        "orig_coord_info"
    ].str.split(",", expand=True)
    return df[
        ["Transcript_id", "transcript_version", "exon_number", "orig_start", "orig_end"]
    ]


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    if (
        not args.fusioninspector.is_file()
        or not args.fusionreport.is_file()
        or not args.fusioninspector_gtf
        or not args.fusionreport_csv
        or not args.hgnc
    ):
        logger.error(
            f"The given input file {args.fusioninspector} or {args.fusionreport} was not found!"
        )
        sys.exit(2)
    vcf_collect(
        args.fusioninspector,
        args.fusionreport,
        args.fusioninspector_gtf,
        args.fusionreport_csv,
        args.hgnc,
        args.sample,
        args.out,
    )


if __name__ == "__main__":
    sys.exit(main())
