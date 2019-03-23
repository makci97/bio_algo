import sys
import numpy as np
from collections import Counter, defaultdict


def create_de_bruijn_graph(pairs):
    graph = defaultdict(list)
    d_from, d_to = Counter(), Counter()
    for pair in pairs:
        v_from = (pair[0][:-1], pair[1][:-1])
        v_to = (pair[0][1:], pair[1][1:])
        graph[v_from].append(v_to)
        d_from[v_from] += 1
        d_to[v_to] += 1

        d_from[v_to] *= 1
        d_to[v_from] *= 1

    first = None
    for v in graph.keys():
        if d_from.get(v, 0) > d_to.get(v, 0):
            first = v
            break

    return graph, first


def get_euler_path(pairs):
    graph, first = create_de_bruijn_graph(pairs)

    euler_path = []
    vertices = [first]
    while len(vertices) > 0:
        v = vertices[-1]
        if len(graph[v]) == 0:
            euler_path.append(v)
            vertices.pop()
        else:
            u = graph[v].pop()
            vertices.append(u)

    return euler_path[::-1]


def euler_path2genome(euler_path, k, d):
    genome_beg = ''.join(
        [euler_path[0][0]] + [pair[0][-1] for pair in euler_path[1:]]
    )
    genome_end = ''.join(
        [euler_path[0][1]] + [pair[1][-1] for pair in euler_path[1:]]
    )
    genome = genome_beg + genome_end[-(k + d):]
    return genome


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        k, d = map(int, file.readline().strip().split())
        pairs = [tuple(pair.strip().split('|')) for pair in file.readlines()]

    # write
    with open(output_filename, 'w') as file:
        file.write(euler_path2genome(get_euler_path(pairs), k, d))
