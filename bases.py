from __future__ import print_function
import exrex
import itertools

DEGEN_CODE_UPPER = {
    "W" : "[AT]",
    "S" : "[CG]",
    "M" : "[AC]",
    "K" : "[GT]",
    "R" : "[AG]",
    "Y" : "[CT]",
    "B" : "[CGT]",
    "D" : "[AGT]",
    "H" : "[ACT]",
    "V" : "[ACG]",
    "N" : "[ACTG]"
}

DEGEN_CODE_LOWER = {
    "w" : "[at]",
    "s" : "[cg]",
    "m" : "[ac]",
    "k" : "[gt]",
    "r" : "[ag]",
    "y" : "[ct]",
    "b" : "[cgt]",
    "d" : "[agt]",
    "h" : "[act]",
    "v" : "[acg]",
    "n" : "[actg]"
}

DEGEN_CODE_SENSITIVE = dict(list(DEGEN_CODE_LOWER.items()) + list(DEGEN_CODE_UPPER.items()))

print(DEGEN_CODE_SENSITIVE)

def make_regex(bases, table=DEGEN_CODE_SENSITIVE):
    output = []
    for base in bases:
        if base in table:
            output += table[base]
        else:
            output.append(base)
    return "".join(output)

def generate_all(bases):
    regex = make_regex(bases)
    return list(exrex.generate(regex))

bases = "AAAAAAAAACCCCCCCCCTGYY"
print(bases)
print(make_regex(bases))
print(generate_all(bases))
print(itertools.product(generate_all(), generate_all()))