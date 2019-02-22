#!/usr/bin/env python
from __future__ import print_function
import sys
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp
from Bio import Alphabet

usage ='''primertool [accepting better names]
get primer related information
'''

def bisulfite_convert(seq):
    '''Get a bisulfite converted sequence based on Seq [seq] and a primer type
    This function is intended to produce strands complimentary to the p_type
    '''
    bases = str(seq)
    bases = bases.replace("CG", "yG").replace("GC", "Gy").replace("C", "t")
    return Seq(bases, Alphabet.generic_dna)

def get_stats(bases, p_type):
    bases = bases.upper()
    stats = {"primer type" : p_type}
    # note that we automatically orientation for type b primers
    stats["unconverted"] = Seq(bases)
    if p_type == "b1" or p_type == "b2":
        stats["unconverted"] = stats["unconverted"].reverse_complement()
    # convert the strand
    stats["converted"] = bisulfite_convert(stats["unconverted"])
    if p_type == "a2" or p_type == "b2":
        stats["primer"] = stats["converted"]
    elif p_type == "a1" or p_type == "b1":
        stats["primer"] = stats["converted"].reverse_complement()
    else:
        raise ValueError("Invalid p_type %s" % (p_type))
    stats["template"] = stats["primer"].reverse_complement()
    stats["tm"] = MeltingTemp.Tm_NN(stats["primer"], dnac1=250, Na=100)
    stats["length"] = len(bases)
    return stats

# 

def has_bad_cg(primer):
    '''Check that if farthest cg/gc in a primer
    is too close too the 3' end
    '''
    return (primer.rfind("r") > (len(primer) // 2) or
            primer.rfind("y") > (len(primer) // 2))

def exhaustive_search(region, p_type, minlength=20, maxlength=30):
    results = []
    for length in range(minlength, maxlength):
        for start in range(0, len(region) - length + 1):
            result = get_stats(region[start:start+length], p_type)
            cg_count = result["primer"].count("y") + result["primer"].count("r")
            if has_bad_cg(result["primer"]):
                continue
            if cg_count > 2:
                continue
            if len(result["primer"]) > 28:
                continue
            results.append(result)
    results.sort(key=(lambda x: x["tm"]), reverse=True)
    #sprint(list(map(lambda x: x["tm"], results)))
    return results

def print_stats(stats):
    #find the max length of a value
    max_length = 0
    for key in stats:
        length = len(key)
        if length > max_length:
            max_length = length
    
    for key in ["primer type", "unconverted", "converted", "template", "primer", "tm", "length"]:
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
            p_types, bases = arg.split(":")
        except ValueError:
            try: 
                p_types, bases, mode = arg.split(":")
            except ValueError:
                print("Invalid argument: %s" % arg, file=sys.stderr)
                continue
        p_types = p_types.split(",")
        if mode not in ["default", "best", "top3"]:
            print("Search mode not recognized: %s" % mode, file=sys.stderr)
            continue
        for p_type in p_types:
            if p_type not in ["a1", "a2", "b1", "b2"]:
                print("Primer type not recognized: %s" % p_type, file=sys.stderr)
                continue
            if mode == "best":
                print_stats(exhaustive_search(bases, p_type)[0])
            elif mode == "top3":
                for result in exhaustive_search(bases, p_type)[:3]:
                    print_stats(result)
            else:
                print_stats(get_stats(bases, p_type))