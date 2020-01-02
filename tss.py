#!/usr/bin/env python3
import argparse
import re
import sys
from collections import namedtuple
from enum import Enum

# GTF feature based on:
# https://useast.ensembl.org/info/website/upload/gff.html
Feature = namedtuple(
    typename = "Feature",
    field_names = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute"
    ]
)


ATTR_REG = re.compile(r"(\w+) \"(.*?)\";")


def _line_to_feature(line):
    fields = line.strip().split("\t")
    if len(fields) < 9:
        raise ValueError(f"Expected 9 fields, received {len(fields)}")
    attributes = dict(ATTR_REG.findall(fields[8]))
    return Feature(
        "chr" + fields[0],
        fields[1],
        fields[2],
        int(fields[3]),
        int(fields[4]),
        fields[5],
        fields[6],
        fields[7],
        attributes
    )


Feature.from_line = _line_to_feature

def feature_to_bed(f, start_tss=None, end_tss=None):
    start = f.start
    end = f.end
    if start_tss is not None:
        start = f.start + start_tss
    if end_tss is not None:
        end = f.start + end_tss
    if "transcript_id" in f.attribute:
        name = f.attribute["transcript_id"]
    elif "gene_id" in f.attribute:
        name = f.attribute["gene_id"]
    else:
        name = "."
    if "gene_name" in f.attribute:
        name = f'{f.attribute["gene_name"]}-{name}'
    return f"{f.seqname}\t{start}\t{end}\t{name}\t{f.strand}\t{f.score}"

class FType(Enum):
    """enum representing the supported Feature types"""
    GENE = 0
    TRANSCRIPT = 1

GENE_ID_REG = re.compile(r"\t.*?\tgene.*?gene_id \"(.*?)\"")
TRANSCRIPT_ID_REG = re.compile(r"\t.*?\ttranscript.*?transcript_id \"(.*?)\"")


def find_ids(gtf_file, ids, feature_type):
    """
    returns  a list of lines with matching IDs
    gtf_file: an iterable of strings
    ids: a collection that supports "contains"
    feature_type: FType.GENE or FType.TRANSCRIPT (see FType Enum)
    """
    results = []
    if feature_type == FType.GENE:
        pattern = GENE_ID_REG
    elif feature_type == FType.TRANSCRIPT:
        pattern = TRANSCRIPT_ID_REG
    else:
        raise ValueError("Expected GENE or TRANSCRIPT for feature_type")
    for index, line in enumerate(gtf_file):
        match = pattern.search(line)
        if match:
            feature_id = match.group(1)
            if feature_id in ids:
                results.append(line)
    return results

gtf = open("annotated/Homo_sapiens.GRCh38.head.gtf")

parser = argparse.ArgumentParser(description = "Find regions around the TSS of different transcripts.")
parser.add_argument("--gtf", required=True, help="A GTF file containing genes / transcripts")
parser.add_argument("--ids", nargs="+", required=True,
                    help="(whitespace-separated) list of Ensembl IDs")
#parser.add_argument("--feature-type", nargs=1, choices=["transcript", "gene"], default="TRANSCRIPT", help="Type of features [GENE or TRANSCRIPT]")

if __name__ == "__main__":
    args = parser.parse_args()
    ftype = FType.TRANSCRIPT
    with open(args.gtf) as gtf_file:
        features = find_ids(gtf_file, set(args.ids), ftype)
        for f in features:
            f = Feature.from_line(f)
            print(feature_to_bed(f, start_tss=-500, end_tss=500))