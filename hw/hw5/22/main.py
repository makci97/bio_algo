import sys
import numpy as np
from functools import lru_cache
from collections import Counter, defaultdict


def length_of_longest_path(graph, n, m):
    max_way_weight = np.zeros((n + 1, m + 1))
    for i in range(n + 1):
        for j in range(m + 1):
            if i == 0 and j == 0:
                continue
            candidates = []
            if i > 0:
                candidates.append(max_way_weight[i - 1, j] + graph[(i - 1, j)][(i, j)])
            if j > 0:
                candidates.append(max_way_weight[i, j - 1] + graph[(i, j - 1)][(i, j)])
            max_way_weight[i][j] = max(candidates)
    return max_way_weight[n][m]


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    graph = defaultdict(Counter)
    with open(input_filename, 'r') as file:
        n, m = map(int, file.readline().strip().split())
        for i in range(n):
            for j, w in enumerate(map(int, file.readline().strip().split())):
                graph[(i, j)][(i + 1, j)] = w
        file.readline()
        for i in range(n + 1):
            for j, w in enumerate(map(int, file.readline().strip().split())):
                graph[(i, j)][(i, j + 1)] = w

    # write
    with open(output_filename, 'w') as file:
        file.write(str(int(length_of_longest_path(graph, n, m))))
