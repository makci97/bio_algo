import sys
import numpy as np
from collections import Counter


def dist(profile, motif):
    return sum(p != m for p, m in zip(profile, motif))


def dist_btw_pattern_and_dna(pattern, dna):
    k = len(pattern)
    res_dist = 0
    for line in dna:
        h_dist = np.inf
        for i in range(len(line) - k + 1):
            cur_dist = dist(pattern, line[i:i+k])
            if h_dist > cur_dist:
                h_dist = cur_dist
        res_dist += h_dist
    return res_dist


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        pattern = [nuc for nuc in str(file.readline().strip())]
        dna = [
            [nuc for nuc in line]
            for line in file.readline().strip().split(' ')
        ]

    # write
    with open(output_filename, 'w') as file:
        file.write(str(dist_btw_pattern_and_dna(pattern, dna)))
