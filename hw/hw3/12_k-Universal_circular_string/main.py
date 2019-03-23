import sys
import numpy as np
from collections import Counter, defaultdict


def gen_binary_seq(n, pref=None):
    if pref is None:
        pref = []

    if n == 0:
        yield ''.join(map(str, pref))
    else:
        pref.append(0)
        yield from gen_binary_seq(n - 1, pref)
        pref[-1] = 1
        yield from gen_binary_seq(n - 1, pref)
        pref.pop()


def create_de_bruijn_graph(mers):
    graph = defaultdict(list)
    d_in, d_out = Counter(), Counter()
    for mer in mers:
        graph[mer[:-1]].append(mer[1:])
        d_in[mer[1:]] += 1
        d_out[mer[:-1]] += 1

    first = ''.join(map(str, np.zeros(len(mers[0]) - 1, dtype=int)))

    return graph, first


def get_euler_path(mers):
    graph, first = create_de_bruijn_graph(mers)

    # print(graph, first)
    euler_path = []
    vertices = [first]
    while len(vertices) > 0:
        v = vertices[-1]
        if len(graph[v]) == 0:
            euler_path.append(v)
            vertices.pop()
        else:
            u = graph[v].pop(0)
            vertices.append(u)

    return euler_path[:0:-1]


def euler_path2genome(euler_path, k):
    # print(euler_path)
    genome = [euler_path[0]] + [mer[-1] for mer in euler_path[1:-(k-2)]]
    return ''.join(genome)


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        k = list(map(int, file.readline().strip().split()))[0]

    mers = [mer for mer in gen_binary_seq(k)]
    # print(mers)

    # # write
    with open(output_filename, 'w') as file:
        file.write(euler_path2genome(get_euler_path(mers), k))
