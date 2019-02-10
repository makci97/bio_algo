import sys
import numpy as np


# def get_skew(text):
#     cum_skew = 0
#     skew = []
#     for nuc in text:
#         skew.append(cum_skew)
#         if nuc == 'G':
#             cum_skew += 1
#         elif nuc == 'C':
#             cum_skew -= 1
#     return skew

def get_skew(text):
    skew = np.where(list(map(lambda t: t == 'G', text)), 1, 0)
    skew += np.where(list(map(lambda t: t == 'C', text)), -1, 0)
    return np.cumsum(skew)


def args_min(array):
    cond_min = (array == np.min(array))
    return np.arange(array.shape[0])[cond_min]


def get_args_min_skew(text):
    return args_min(get_skew(text)) + 1


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        text = str(file.readline().strip())

    # write
    with open(output_filename, 'w') as file:
        file.write(' '.join(map(str, get_args_min_skew(text))))
