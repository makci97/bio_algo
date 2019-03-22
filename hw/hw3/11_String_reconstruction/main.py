import sys
import numpy as np
from collections import Counter, defaultdict


def create_de_bruijn_graph(mers):
    graph = defaultdict(list)
    d_in, d_out = Counter(), Counter()
    for mer in mers:
        graph[mer[:-1]].append(mer[1:])
        d_in[mer[1:]] += 1
        d_out[mer[:-1]] += 1

    first = None
    for v in d_out.keys():
        if d_in.get(v, 0) != d_out.get(v, 0):
            first = v
            break

    return graph, first


def get_euler_path(mers):
    graph, first = create_de_bruijn_graph(mers)

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


def euler_path2genome(euler_path):
    genome = [euler_path[0]] + [mer[-1] for mer in euler_path[1:]]
    return ''.join(genome)


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        k = map(int, file.readline().strip().split())
        mers = [mer.strip() for mer in file.readlines()]

    # write
    with open(output_filename, 'w') as file:
        file.write(euler_path2genome(get_euler_path(mers)))
