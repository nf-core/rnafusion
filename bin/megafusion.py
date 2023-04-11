#!/usr/bin/env python

import argparse
import csv
import logging
import os
import sys
from collections import Counter
from pathlib import Path
import pandas as pd
import ast
import re

logger = logging.getLogger()

FUSIONINSPECTOR_MAP = {
    "fusion": {"column": 0, "delimiter": "\t", "element": 0},
    "chromosomeA": {"column": 7, "delimiter": ":", "element": 0},
    "chromosomeB": {"column": 10, "delimiter": ":", "element": 0},
    "posA": {"column": 7, "delimiter": ":", "element": 1},
    "posB": {"column": 10, "delimiter": ":", "element": 1},
    "strand1": {"column": 7, "delimiter": ":", "element": 2},
    "strand2": {"column": 10, "delimiter": ":", "element": 2},
    "geneA": {"column": 0, "delimiter": "--", "element": 0},
    "geneB": {"column": 0, "delimiter": "--", "element": 1},
    "split_reads": {"column": 1, "delimiter": "\t", "element": 0},
    "discordant_pairs": {"column": 2, "delimiter": "\t", "element": 0},
    "ffpm": {"column": 25, "delimiter": "\t", "element": 0},
}


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
    parser.add_argument("--sample", metavar="SAMPLE", type=Path, help="Sample name.", default="Sample")
    parser.add_argument(
        "--out",
        metavar="OUT",
        type=Path,
        help="Output path.",
    )
    return parser.parse_args(argv)


def header_def(sample):
    return '##fileformat=VCFv4.1\n\
##ALT=<ID=BND,Description="Break end">\n\
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n\
##INFO=<ID=CHRA,Number=1,Type=String,Description="Chromosome A">\n\
##INFO=<ID=CHRB,Number=1,Type=String,Description="Chromosome B">\n\
##INFO=<ID=GENEA,Number=.,Type=String,Description="Gene A">\n\
##INFO=<ID=GENEB,Number=.,Type=String,Description="Gene B">\n\
##INFO=<ID=ORIENTATION,Number=.,Type=String,Description="Strand1 and strand2 directions">\n\
##INFO=<ID=FOUND_DB,Number=.,Type=String,Description="Databases in which the fusion has been found">\n\
##INFO=<ID=ARRIBA,Number=.,Type=String,Description="Found by arriba">\n\
##INFO=<ID=FUSIONCATCHER,Number=.,Type=String,Description="Found by fusioncatcher">\n\
##INFO=<ID=PIZZLY,Number=.,Type=String,Description="Found by pizzly">\n\
##INFO=<ID=SQUID,Number=.,Type=String,Description="Found by squid">\n\
##INFO=<ID=STARFUSION,Number=.,Type=String,Description="Found by starfusion">\n\
##INFO=<ID=TOOL_HITS,Number=.,Type=String,Description="Number of tools that found the fusion">\n\
##INFO=<ID=SCORE,Number=.,Type=String,Description="Score from fusionreport">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n\
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of paired-ends that support the event">\n\
##FORMAT=<ID=RV,Number=1,Type=Integer,Description="Number of split reads that support the event">\n\
##FORMAT=<ID=FFPM,Number=1,Type=Integer,Description="Fusion fragments per million total RNA-seq fragments">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}'.format(
        sample
    )


def read_fusioninspector(fusioninspector_file, col_num, delimiter, element):
    with open(fusioninspector_file) as fusioninspector:
        return [line.split()[col_num].split(delimiter)[element] for line in fusioninspector if not line.startswith("#")]


def build_fusioninspector_dataframe(file, map):
    new_dict = {}
    for key in FUSIONINSPECTOR_MAP:
        new_dict[key] = read_fusioninspector(
            file,
            map[key]["column"],
            map[key]["delimiter"],
            map[key]["element"],
        )
    return pd.DataFrame.from_dict(new_dict).set_index("fusion")


def read_build_fusionreport(fusionreport_file):
    with open(fusionreport_file) as fusionreport:
        from_html = [line.split('rows": [')[1] for line in fusionreport if 'name="fusion_list' in line]
        expression = from_html[0].split('], "tool')[0]
        return pd.DataFrame.from_dict(ast.literal_eval(expression)).set_index("fusion")


def column_manipulation(df):
    df["ALT"] = ""
    df = df.reset_index()
    df["FORMAT"] = "GT:DV:RV:FFPM"
    df["ID"] = "."
    df["QUAL"] = "."
    df["FILTER"] = "PASS"
    df["REF"] = "N"

    for index, row in df.iterrows():
        # ALT
        if not row["strand1"] in ["+", "-"] or not row["strand2"] in ["+", "-"]:
            df.loc[index, "ALT"] = "N[{}:{}[".format(df["chromosomeB"], row["posB"])
        elif row["strand1"] == "-" and row["strand2"] == "-":
            df.loc[index, "ALT"] = "[{}:{}[N".format(row["chromosomeB"], row["posB"])
        elif row["strand1"] == "+" and row["strand2"] == "-":
            df.loc[index, "ALT"] = "N]{}:{}]".format(row["chromosomeB"], row["posB"])
        elif row["strand1"] == "-" and row["strand2"] == "+":
            df.loc[index, "ALT"] = "N]{}:{}]".format(row["chromosomeB"], row["posB"])
        else:
            df.loc[index, "ALT"] = "N[{}:{}[".format(row["chromosomeB"], row["posB"])
        # INFO
        df.loc[index, "INFO"] = (
            "SVTYPE=BND;CHRA={};CHRB={};GENEA={};GENEB={};ORIENTATION={},{};FOUND_DB={};"
            "ARRIBA={};FUSIONCATCHER={};PIZZLY={};SQUID={};STARFUSION={};TOOL_HITS={};SCORE={}".format(
                row["chromosomeA"],
                row["chromosomeB"],
                row["geneA"],
                row["geneB"],
                row["strand1"],
                row["strand2"],
                row["found_db"],
                row["arriba"],
                row["fusioncatcher"],
                row["pizzly"],
                row["squid"],
                row["starfusion"],
                row["tools_hits"],
                row["score"],
            )
        )
        # FORMAT
        df.loc[index, "Sample"] = "./1:{}:{}:{}".format(row["split_reads"], row["discordant_pairs"], row["ffpm"])
    return df


def write_vcf(df_to_print, header, out_file):
    df_to_print[["chromosomeA", "posA", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample",]].to_csv(
        path_or_buf=out_file,
        sep="\t",
        header=None,
        index=False,
    )

    with open(out_file, "r+") as f:
        content = f.read()
        f.seek(0, 0)
        f.write(header.rstrip("\r\n") + "\n" + content)


def megafusion(fusioninspector_in_file, fusionreport_in_file, sample, out):
    merged_df = build_fusioninspector_dataframe(fusioninspector_in_file, FUSIONINSPECTOR_MAP).join(
        read_build_fusionreport(fusionreport_in_file), how="left"
    )
    write_vcf(column_manipulation(merged_df), header_def(sample), out)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    if not args.fusioninspector.is_file() or not args.fusionreport.is_file():
        logger.error(f"The given input file {args.fusioninspector} or {args.fusionreport} was not found!")
        sys.exit(2)
    megafusion(args.fusioninspector, args.fusionreport, args.sample, args.out)


if __name__ == "__main__":
    sys.exit(main())
