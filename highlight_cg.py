#!/usr/bin/env python
from __future__ import print_function
import sys
from Bio import SeqIO
from Bio import SeqFeature as SF
import re

USAGE='''highlight_cg.py
Usage: highlight_cg.py in.gb sequence [out.gb]
Add annotations to all sequences matching [seq]
Output a GB file
'''

CPG = re.compile("CG")
GPC = re.compile("GC")

def remove_gc_features(features):
    return [feature for feature in features if "CpG" not in feature.qualifiers and "GpC" not in feature.qualifiers]

def match_to_feature(match, name=None, color=None):
    feat = SF.SeqFeature(SF.FeatureLocation(SF.ExactPosition(match.start()), 
                                            SF.ExactPosition(match.end())), 
                                            type="misc_feature")
    if name:
        feat.qualifiers["label"] = name
    if color:
        feat.qualifiers["ApEinfo_revcolor"] = color
        feat.qualifiers["ApEinfo_fwdcolor"] = color
    return feat

def get_all_gc(sequence):
    cpgs = [ match_to_feature(pos, "CpG", "#cc0a1d") for pos in CPG.finditer(sequence)]
    gpcs = [ match_to_feature(pos, "GpC", "#1dcc0a") for pos in GPC.finditer(sequence)]
    return (cpgs, gpcs)


def read_genbank(filename):
    return SeqIO.read(filename, "genbank")


def write_genbank(record, handle):
    output = record.format('gb')
    handle.write(output)

def main():
    if len(sys.argv) < 2:
        print(USAGE, file=sys.stderr)
        print("Error: Expected 1-2 arguments, received %s" 
              % (len(sys.argv)-1), file=sys.stderr)
        exit(1)
    filename = sys.argv[1]
    output = sys.stdout
    if len(sys.argv) > 2:
        output = open(sys.argv[2], 'w')
    record = read_genbank(filename)
    cpgs, gpcs = get_all_gc(str(record.seq))
    print(cpgs, gpcs)
    print(record.features)
    record.features = remove_gc_features(record.features)
    print(record.features)
    record.features += cpgs + gpcs
    write_genbank(record, output)
    output.close()

if __name__ == "__main__":
	main()