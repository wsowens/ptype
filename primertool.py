#!/usr/bin/env python
from __future__ import print_function
import sys
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp
from Bio import Alphabet

usage ='''primertool [accepting better names]
get primer related information
'''

def bisulfite_convert(seq, p_type):
    '''Get a bisulfite converted sequence based on Seq [seq] and a primer type
    This function is intended to produce strands complimentary to the p_type
    '''
    bases = str(seq)
    if p_type == "a1" or p_type == "b1":
        bases = bases.replace("CG", "yG").replace("GC", "Gy").replace("C", "t")
    elif p_type == "a2" or p_type == "b2":
        bases = bases.replace("CG", "Cr").replace("GC", "rC").replace("G", "a")
    else:
        raise ValueError("Invalid primer type \'%s\'" % p_type)
    return Seq(bases, Alphabet.generic_dna)

def get_stats(bases, p_type):
    bases = bases.upper()
    stats = {"primer type" : p_type}
    if p_type == "b1" or "a2":
        stats["original"] = Seq(bases, Alphabet.generic_dna)
        stats["template"] = stats["original"].reverse_complement()
    else:
        stats["template"] = Seq(bases, Alphabet.generic_dna)
    stats["converted"] = bisulfite_convert(bases, p_type)
    stats["primer"] = stats["converted"].reverse_complement()
    stats["tm"] = MeltingTemp.Tm_NN(stats["primer"])
    stats["length"] = len(bases)
    return stats

def exhaustive_search(region, p_type, minlength=20, maxlength=30):
    results = []
    for length in range(minlength, maxlength):
        for start in range(0, len(region) - length + 1):
            results.append(get_stats(region[start:start+length], p_type))
    results.sort(key=(lambda x: x["tm"]), reverse=True)
    return results

def print_stats(stats):
    #find the max length of a value
    max_length = 0
    for key in stats:
        length = len(key)
        if length > max_length:
            max_length = length
    
    for key in ["primer type", "original", "template", "converted", "primer", "tm", "length"]:
        if key not in stats:
            continue
        spaces = " " * (max_length - len(key) + 1)
        print("%s:" % key, "%s" % stats[key], sep=spaces)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(usage)
        exit(1)
    for arg in sys.argv[1::]:
        mode = "default"
        try:
            p_type, bases = arg.split(":")
        except ValueError:
            try: 
                p_type, bases, mode = arg.split(":")
            except ValueError:
                print("Invalid argument: %s" % arg, file=sys.stderr)
                continue
        if p_type not in ["a1", "a2", "b1", "b2"]:
            print("Primer type not recognized: %s" % p_type, file=sys.stderr)
            continue
        if mode not in ["default", "best"]:
            print("Search mode not recognized: %s" % mode, file=sys.stderr)
            continue
        if mode == "best":
            print_stats(exhaustive_search(bases, p_type)[0])
        else:
            print_stats(get_stats(bases, p_type))