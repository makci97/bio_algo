import sys
import numpy as np
from functools import lru_cache
from collections import Counter, defaultdict


def breakpoints_num(perm):
    perm = [0] + perm + [len(perm) + 1]
    n_breakponits = 0
    for i in range(len(perm) - 1):
        if perm[i + 1] - perm[i] != 1:
            n_breakponits += 1
    return n_breakponits


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    graph = defaultdict(Counter)
    with open(input_filename, 'r') as file:
        perm = list(map(int, file.readline().strip()[1:-1].split()))

    # write
    with open(output_filename, 'w') as file:
        file.write(str(breakpoints_num(perm)))
