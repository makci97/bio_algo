import sys
import numpy as np
from enum import Enum
from functools import lru_cache
from collections import Counter, defaultdict


def make_graph(gen):
    graph = defaultdict(list)
    add_gen_to_graph(graph, gen)
    return graph


def add_gen_to_graph(graph, gen):
    for subgen in gen:
        for i in range(len(subgen)):
            next_i = (i + 1) % len(subgen)

            from_type = 2 if subgen[i] > 0 else 1
            to_type = 1 if subgen[next_i] > 0 else 2
            from_v = (abs(subgen[i]) - 1) * 2 + from_type
            to_v = (abs(subgen[next_i]) - 1) * 2 + to_type

            graph[from_v].append(to_v)
            graph[to_v].append(from_v)


def graph2gen(graph, old_gen):
    gen = []
    used_nums = defaultdict(lambda: False)
    max_num = max(graph.keys()) // 2
    flatten_old_gen = [num for subgen in old_gen for num in subgen]

    for num in range(1, max_num + 1):
        if used_nums[num]:
            continue
        # init subgen collecting process
        subgen = []
        sign = 1
        cur_num = num
        subgen.append(cur_num * sign)
        used_nums[cur_num] = True
        cur_num, sign = next_signed_num(graph, cur_num, sign)

        # collect whole subgen
        while cur_num != num:
            subgen.append(cur_num * sign)
            used_nums[cur_num] = True
            cur_num, sign = next_signed_num(graph, cur_num, sign)

        # check sign of subgen
        if intersec_len(flatten_old_gen, subgen) < intersec_len(flatten_old_gen, invert_sign(subgen)):
            subgen = invert_sign(subgen)

        gen.append(subgen)
    return gen


def next_signed_num(graph, num, sign):
    from_type = 2 if sign > 0 else 1
    v = (num - 1) * 2 + from_type
    u = graph[v][0]
    next_num = (u + 1) // 2
    next_sign = 1 if u % 2 == 1 else -1
    return next_num, next_sign


def intersec_len(flatten_old_gen, subgen):
    return len(set(flatten_old_gen).intersection(set(subgen)))


def invert_sign(subgen):
    return [-num for num in subgen]


def two_break_on_gen_graph(graph, i_1, i_2, j_1, j_2):
    # remove
    graph[i_1].remove(i_2)
    graph[i_2].remove(i_1)
    graph[j_1].remove(j_2)
    graph[j_2].remove(j_1)

    # add
    graph[i_1].append(j_1)
    graph[j_1].append(i_1)
    graph[i_2].append(j_2)
    graph[j_2].append(i_2)


def gen2str(gen):
    str_subgens = []
    for subgen in gen:
        str_subgen = '(' + ' '.join([
            '+' + str(num) if num > 0 else str(num)
            for num in subgen
        ]) + ')'
        str_subgens.append(str_subgen)
    return ' '.join(str_subgens)


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        gen = list(map(
            lambda block: list(map(int, block[1:].split())),
            file.readline().strip().split(')')[:-1]
        ))
        i_1, i_2, j_1, j_2 = map(int, file.readline().strip().split(', '))

    # print(gen)
    # print(i_1, i_2, j_1, j_2)
    graph = make_graph(gen)
    # print(graph.items())
    two_break_on_gen_graph(graph, i_1, i_2, j_1, j_2)
    # print(graph.items())
    # print(graph2gen(graph, gen))
    # print(gen2str(graph2gen(graph, gen)))

    # write
    with open(output_filename, 'w') as file:
        file.write(gen2str(graph2gen(graph, gen)))
