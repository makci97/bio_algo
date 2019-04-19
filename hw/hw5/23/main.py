import sys
import numpy as np

from enum import Enum
from functools import lru_cache
from collections import Counter, defaultdict


class Move(Enum):
    none = 0
    down = 1
    right = 2
    diag = 3


def longest_path(word1, word2, weights, sigma):
    n, m = len(word1), len(word2)
    max_way = []
    for i in range(n + 1):
        row = []
        for j in range(m + 1):
            if i == 0 and j == 0:
                row.append([0, Move.none])
                continue
            candidates = []
            if i > 0:
                candidates.append([max_way[i - 1][j][0] + sigma, Move.down])
            if j > 0:
                candidates.append([row[j - 1][0] + sigma, Move.right])
            if i > 0 and j > 0:
                candidates.append([
                    max_way[i - 1][j - 1][0] + weights[word1[i - 1]][word2[j - 1]],
                    Move.diag
                ])
            row.append(sorted(candidates, key=lambda pair: -pair[0])[0])
        max_way.append(row)
    return max_way


def unroll_max_way(max_way, word1, word2):
    i = len(max_way) - 1
    j = len(max_way[i]) - 1

    first_line = []
    second_line = []
    while i > 0 or j > 0:
        if max_way[i][j][1] == Move.down:
            first_line.append(word1[i - 1])
            second_line.append('-')
            i -= 1
        elif max_way[i][j][1] == Move.right:
            first_line.append('-')
            second_line.append(word2[j - 1])
            j -= 1
        elif max_way[i][j][1] == Move.diag:
            first_line.append(word1[i - 1])
            second_line.append(word2[j - 1])
            i -= 1
            j -= 1
    return str(max_way[-1][-1][0]), ''.join(first_line[::-1]), ''.join(second_line[::-1])


def read_weights(filename):
    weights = defaultdict(dict)
    with open(filename, 'r') as file:
        col_names = file.readline().strip().split()
        for row in file.readlines():
            from_pep = row.strip().split()[0]
            for i, w in enumerate(row.strip().split()[1:]):
                weights[from_pep][col_names[i]] = int(w)
    return weights


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'
    weights_filename = 'weights.txt'

    # read
    with open(input_filename, 'r') as file:
        word1 = file.readline().strip()
        word2 = file.readline().strip()

    weights = read_weights(weights_filename)
    sigma = -5

    # write
    with open(output_filename, 'w') as file:
        file.write('\n'.join(
            unroll_max_way(longest_path(word1, word2, weights, sigma), word1, word2)
        ))
