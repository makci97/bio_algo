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


def expand(leaderboard, spectrum):
    peptides = [pair[0] + (mass,) for pair in leaderboard for mass in mass_of.values()]
    return list(map(lambda peptide: (peptide, score(peptide, spectrum)), peptides))


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


def score(peptide, spectrum):
    score_res = 0
    pep_spectrum = sorted(gen_cyclo_spectrum(peptide))

    i = 0
    j = 0
    while i < len(pep_spectrum) and j < len(spectrum):
        if pep_spectrum[i] == spectrum[j]:
            score_res += 1
            i += 1
            j += 1
        elif pep_spectrum[i] < spectrum[j]:
            i += 1
        else:
            j += 1

    return score_res


def cut(leaderboard, n):
    return sorted(leaderboard, key=lambda pair: -pair[1])[:n]


def leaderboard_cyclopeptide_sequencing(spectrum, n):
    best_peptide = ()
    best_score = 0

    leaderboard = [((), 0)]

    parent_mass = max(spectrum)
    while len(leaderboard) > 0:
        leaderboard = expand(leaderboard, spectrum)
        i = 0
        while i < len(leaderboard):
            if peptide_mass(leaderboard[i][0]) == parent_mass:
                if leaderboard[i][1] > best_score:
                    best_peptide = leaderboard[i][0]
                    best_score = leaderboard[i][1]
            elif peptide_mass(leaderboard[i][0]) > parent_mass:
                leaderboard.pop(i)
                i -= 1
            i += 1
        leaderboard = cut(leaderboard, n)

    return best_peptide


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        n = int(file.readline())
        spectrum = [int(mass) for mass in file.readline().split()]

    peptide = leaderboard_cyclopeptide_sequencing(spectrum, n)
    # write
    with open(output_filename, 'w') as file:
        file.write('-'.join(map(str, peptide)))
