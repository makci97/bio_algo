import sys
import numpy as np
from functools import lru_cache
from collections import Counter, defaultdict


def greedy_sort(perm):
    seq = []
    for i in range(1, len(perm) + 1):
        if abs(perm[i - 1]) != i:
            j = i
            while j < len(perm) and abs(perm[j]) != i:
                j += 1
            perm[i - 1: j + 1] = -perm[j: i - 2 - len(perm): -1]
            seq.append(perm.copy())
        if perm[i - 1] < 0:
            perm[i - 1] = abs(perm[i - 1])
            seq.append(perm.copy())

    return seq


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    graph = defaultdict(Counter)
    with open(input_filename, 'r') as file:
        perm = np.array(list(map(int, file.readline().strip()[1:-1].split())))

    seq = greedy_sort(perm)
    res = '\n'.join(map(
        lambda p: '(' + ' '.join(map(
            lambda num: str(num) if num < 0 else '+'+str(num),
            p
        )) + ')',
        seq
    ))

    # write
    with open(output_filename, 'w') as file:
        file.write(res)
