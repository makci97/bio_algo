import sys
import numpy as np
from math import factorial
from functools import lru_cache
from collections import Counter, defaultdict

mapping = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}

compl_nucs = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
}

mass_of = {
    'G': 57,
    'A': 71,
    'S': 87,
    'P': 97,
    'V': 99,
    'T': 101,
    'C': 103,
    'L': 113,
    'N': 114,
    'D': 115,
    'Q': 128,
    'E': 129,
    'M': 131,
    'H': 137,
    'F': 147,
    'R': 156,
    'Y': 163,
    'W': 186,

}


@lru_cache(maxsize=100)
def frac(n):
    return factorial(n)


def num_combinations(peptide):
    return frac(len(peptide)) / np.prod([frac(n) for n in Counter(peptide).values()])


def pep_num_by_mass(amins, mass, peptide):
    if mass < 0 or (len(amins) == 0 and mass > 0):
        return 0
    elif mass == 0:
        return num_combinations(peptide)
    else:
        res = pep_num_by_mass(amins, mass - amins[0][1], peptide + amins[0][0])
        res += pep_num_by_mass(amins[1:], mass, peptide)
        return res


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        parent_mass = int(file.readline().strip())

    # write
    with open(output_filename, 'w') as file:
        file.write(str(int(pep_num_by_mass(list(mass_of.items()), parent_mass, ''))))
