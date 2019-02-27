#!/usr/bin/env python
from __future__ import print_function
import sys
import math
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp
from Bio import Alphabet
from AligoOnalyzer import analyze

usage ='''primertool [accepting better names]
get primer related information'''


'''
Basing these concentrations on this methods paper:
https://www.ncbi.nlm.nih.gov/pubmed/25827879
    Use 2-4 uL deaminated DNA as template in the following 50 uL PCR reaction:
    dH2O to volume, 1x Qiagen HotStar Coral PCR buffer, 2.25 mM MgCl2,
    0.2 mM dNTPs, 0.8 uM b1 (or a1) primer, 0.8 uM b2 (or a2) primer,
    and 1.25 U Qiagen HotStar Plus Taq DNA polymerase (see Note 21 ).
[Tris] and [K] are not provided by Qiagen (as far as I can see) so
estimates were made.
'''
# all values in units of mM, except dnac1, which is in nM
PCR_MIX = {
    "dnac1": 800, 
    "Mg" : 2,  
    "dNTPs" : 0.2,
    "Na" : 0,     
    "K" : 50,     
    "Tris": 20    
}

def bisulfite_convert(seq):
    '''Get a bisulfite converted sequence based on Seq [seq] and a primer type
    This function is intended to produce strands complimentary to the p_type
    '''
    bases = str(seq)
    bases = bases.replace("CG", "yG").replace("GC", "Gy").replace("C", "t")
    return Seq(bases, Alphabet.generic_dna)

def get_stats(bases, p_type):
    '''get stats for a set of bases, if designing for p_type'''
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
    stats["tm"] = MeltingTemp.Tm_NN(stats["primer"], **PCR_MIX)
    stats["length"] = len(bases)
    return stats


def format_stats(stats):
    #find the max length of a value
    max_length = 0
    for key in stats:
        length = len(key)
        if length > max_length:
            max_length = length
    lines = []
    for key in ["primer type", "unconverted", "converted", "template", "primer", "tm", "length"]:
        if key not in stats:
            continue
        spaces = " " * (max_length - len(key) + 1)
        lines.append("%s:" % key + spaces + "%s" % stats[key])
    return "\n".join(lines)

def analyze_stats(stats):
    return analyze(stats["primer"], "DNA", **PCR_MIX)

def has_bad_cg(primer):
    '''Check that if farthest cg/gc in a primer
    is too close too the 3' end
    '''
    return (primer.rfind("r") > (len(primer) // 2) or
            primer.rfind("y") > (len(primer) // 2))


def exhaustive_search(region, p_type, minlength=18, maxlength=25, maxtemp=59, max_cg=0, mintemp=54):
    results = []
    for length in range(minlength, maxlength):
        for start in range(0, len(region) - length + 1):
            result = get_stats(region[start:start+length], p_type)
            cg_count = result["primer"].count("y") + result["primer"].count("r")
            if has_bad_cg(result["primer"]):
                continue
            if cg_count > max_cg:
                continue
            if result["tm"] > maxtemp or result["tm"] < mintemp:
                continue
            results.append(result)
    # sort by tm
    results.sort(key=(lambda x: x["tm"]), reverse=True)
    return results

def pair_diff(p1, p2):
    return abs(p1["tm"] - p2["tm"])


def pair_sort(p1, p2, maxtemp, mintemp, maxdiff):
    return merge_num((maxdiff - pair_diff(p1, p2) )/ maxdiff * 100, float(int(max(p1["tm"], p2["tm"]) )))


def merge_num(num1, num2):
    return float(int(num1)) + num2 * 10 ** (-(1+int(math.log10(num2))))

def format_pair(p1, p2):
    p1 = format_stats(p1).split("\n")
    p2 = format_stats(p2).split("\n")
    max_length = 0
    for line in p1:
        length = len(line)
        if length > max_length:
            max_length = length
    output = []
    for line1, line2 in zip(p1,p2):
        spaces = " " * (max_length - len(line1) + 1)
        output.append(line1 + spaces + line2)
    return "\n".join(output)


def exhaustive_pairing(first_region, second_region, strand, maxdiff=2, search_params={}):
    if strand == "a":
        results1 = exhaustive_search(first_region, "a2", **search_params)
        results2 = exhaustive_search(second_region, "a1", **search_params)
    elif strand == "b":
        results1 = exhaustive_search(first_region, "b1", **search_params)
        results2 = exhaustive_search(second_region, "b2", **search_params)
    else:
        raise ValueError("'%s' is not recognized as a stand. (Valid inputs: 'a' or 'b'.)" % strand)
    # find matches
    # this can be done in at least nlogn time
    pairings = []
    for result in results1:
        diff = float('inf')
        pairing = ()
        for other in results2:
            temp = (result, other)
            temp_diff = pair_diff(*temp)
            if temp_diff < diff and temp_diff < maxdiff:
                diff = pair_diff(*temp)
                pairing = temp
            if temp_diff > diff:
                break
        if pairing:
            pairings.append(pairing)
    pairings.sort(key=(lambda x: pair_sort(x[0], x[1], 54, 59, maxdiff)), reverse=True)
    return pairings
     
def sort_pairings(pair_list, maxtemp, mintemp, maxdiff):
    ranked_pairs = {}
    maxtemp = lambda pair : max(pair[0]["tm"], pair[1]["tm"])
    for first, second in pair_list:
        tier = int((maxdiff - pair_diff(first, second)) / maxdiff * 10)
        if tier not in ranked_pairs or maxtemp(ranked_pairs[tier]) < maxtemp((first, second)):
            ranked_pairs[tier] = (first, second)
    return ranked_pairs.items() 
          

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
        if mode not in ["default", "best", "top3", "all", "pair", "ao"]:
            print("Search mode not recognized: %s" % mode, file=sys.stderr)
            continue
        for p_type in p_types:
            if not (mode == "pair" and p_type in ["a", "b"]) and p_type not in ["a1", "a2", "b1", "b2"]:
                print("Primer type not recognized: %s" % p_type, file=sys.stderr)
                continue
            if mode == "best":
                print(format_stats(exhaustive_search(bases, p_type)[0]))
            elif mode == "top3":
                for result in exhaustive_search(bases, p_type)[:3]:
                    print(format_stats(result))
            elif mode == "all":
                for result in exhaustive_search(bases, p_type):
                    print(format_stats(result))
            elif mode == "pair":
                bases = bases.split(",")
                if len(bases) < 2:
                    print("Expected input with the following format: [ab]:(region1),(region2):pair.", file=sys.stderr)
                    continue
                for tier, (first, second) in sort_pairings(exhaustive_pairing(bases[0], bases[1], p_type), 54, 59, 2):
                    print(format_pair(first, second))
                    print(tier)
                    print(pair_sort(first, second, 54, 59, 2))
            elif mode == "ao":
                print(analyze_stats(get_stats(bases, p_type)))
            else:
                print(format_stats((get_stats(bases, p_type))))
