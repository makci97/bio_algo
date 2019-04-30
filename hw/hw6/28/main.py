import sys
import numpy as np
from enum import Enum
from functools import lru_cache
from collections import Counter, defaultdict


class VColor(Enum):
    white = 0
    black = 1


def make_graph(gen_1, gen_2):
    graph = defaultdict(list)
    add_gen_to_graph(graph, gen_1)
    add_gen_to_graph(graph, gen_2)
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


def find_n_cycles(graph):
    n_cycles = 0
    v_colors = defaultdict(lambda: VColor.white)
    for v in list(graph.keys()):
        if v_colors[v] == VColor.white:
            bfs(v, graph, v_colors)
            n_cycles += 1
    return n_cycles


def bfs(v, graph, v_colors):
    v_list = [v]
    while len(v_list) > 0:
        cur_v = v_list.pop(0)
        v_colors[cur_v] = VColor.black
        v_list.extend([u for u in graph[cur_v] if v_colors[u] == VColor.white])


def dfs(v, graph, v_colors):
    v_colors[v] = VColor.black
    v_stack = [v]
    while len(v_stack) > 0:
        cur_v = v_stack[-1]
        is_proccesed = True
        while len(graph[cur_v]) > 0:
            u = graph[cur_v].pop()
            graph[u].remove(cur_v)
            if v_colors[u] == VColor.white:
                v_colors[u] = VColor.black
                v_stack.append(u)
                is_proccesed = False
                break
        if is_proccesed:
            v_stack.pop()


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        gen_1 = list(map(
            lambda block: list(map(int, block[1:].split())),
            file.readline().strip().split(')')[:-1]
        ))
        gen_2 = list(map(
            lambda block: list(map(int, block[1:].split())),
            file.readline().strip().split(')')[:-1]
        ))

    graph = make_graph(gen_1, gen_2)

    # write
    with open(output_filename, 'w') as file:
        file.write(str(len(graph) // 2 - find_n_cycles(graph)))
