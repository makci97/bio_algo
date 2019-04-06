import sys
import numpy as np
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


def substring_encoding(dna, peptide):
    substrings = []
    for i in range(0, 3):
        k = 0
        for j in range(i, len(dna), 3):
            if peptide[k] == mapping.get(dna[j:j + 3]):
                k += 1
            else:
                k = 0
            if k == len(peptide):
                substrings.append(dna[j - 3 * (len(peptide) - 1):j + 3])
                k = 0
    return substrings


def get_reverse_complement_dna(dna):
    return ''.join(map(lambda nuc: compl_nucs[nuc], dna))[::-1]


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        dna = str(file.readline().strip())
        peptide = str(file.readline().strip())

    substrings = list(map(
        get_reverse_complement_dna,
        substring_encoding(get_reverse_complement_dna(dna), peptide)
    ))
    substrings.extend(substring_encoding(dna, peptide))

    # write
    with open(output_filename, 'w') as file:
        file.write('\n'.join(map(lambda seq: ''.join(map(str, seq)), sorted(substrings))))


"""
AAAGAAGTTTTCGAACCACATTATTAC
AAAGAAGTTTTCGAGCCGCACTACTAC
AAAGAGGTGTTTGAACCTCATTACTAT
AAGGAAGTATTCGAACCACATTACTAT
AAGGAAGTATTTGAGCCTCATTATTAC
AAGGAAGTGTTTGAACCTCACTATTAT
AAGGAGGTATTTGAACCCCACTATTAC
ATAATAATGCGGCTCGAATACTTCCTT
ATAATAATGTGGCTCGAACACTTCTTT
ATAATAGTGAGGCTCAAAAACTTCCTT
ATAGTAATGAGGTTCGAAAACCTCTTT
ATAGTAATGGGGTTCGAAGACTTCCTT
ATAGTAGTGAGGTTCGAAGACTTCCTT
ATAGTAGTGTGGTTCAAATACCTCCTT
GTAATAGTGCGGTTCAAAAACTTCCTT
GTAGTAATGAGGTTCAAAAACCTCCTT
GTAGTAATGGGGCTCAAACACCTCTTT
GTAGTAATGGGGCTCGAAAACCTCCTT
GTAGTAATGGGGTTCGAAGACTTCCTT
GTAGTAGTGCGGCTCAAAAACTTCCTT

"""