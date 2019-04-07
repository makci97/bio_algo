import sys
import numpy as np
from math import factorial
from functools import lru_cache
from collections import Counter, defaultdict

import time

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


def expand(peptides):
    return [peptide + (mass,) for peptide in peptides for mass in mass_of.values()]


def peptide_mass(peptide):
    return sum(peptide)


def is_equal_spectrum(peptide, spectrum):
    return sorted(gen_cyclo_spectrum(peptide)) == sorted(spectrum)


def is_consistent_spectrum(peptide, spectrum):
    return len(set(gen_lin_spectrum(peptide)).difference(set(spectrum))) == 0


def gen_cyclo_spectrum(peptide):
    yield 0
    for i in range(len(peptide)):
        # for j in range(i + 1, len(peptide) + 1):
        for j in range(1, len(peptide) + 1):
            if i < j:
                yield sum(peptide[i:j])
            elif i > j:
                yield sum(peptide[i:]) + sum(peptide[:j])


def gen_lin_spectrum(peptide):
    yield 0
    for i in range(len(peptide)):
        for j in range(i + 1, len(peptide) + 1):
            yield sum(peptide[i:j])


def cyclopeptide_sequencing(spectrum):
    result = []
    peptides = [()]
    parent_mass = max(spectrum)
    while len(peptides) > 0:
        peptides = expand(peptides)
        i = 0
        while i < len(peptides):
            if peptide_mass(peptides[i]) == parent_mass:
                if is_equal_spectrum(peptides[i], spectrum):
                    result.append(peptides[i])
                peptides.pop(i)
                i -= 1
            elif not is_consistent_spectrum(peptides[i], spectrum):
                peptides.pop(i)
                i -= 1
            i += 1

    return result


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        spectrum = [int(mass) for mass in file.readline().split()]

    result = cyclopeptide_sequencing(spectrum)
    # write
    with open(output_filename, 'w') as file:
        file.write(' '.join(map(lambda peptide: '-'.join(map(str, peptide)), result)))
