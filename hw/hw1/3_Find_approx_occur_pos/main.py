import sys


def approx_compare(pattern, sub_text, d):
    return sum(pattern[i] != sub_text[i] for i in range(len(pattern))) <= d


def get_approx_occur_pos(pattern, text, d):
    possitions = []
    pattern_size = len(pattern)
    for i in range(len(text) - pattern_size):
        if approx_compare(pattern, text[i:i + pattern_size], d):
            possitions.append(i)
    return possitions


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        pattern = str(file.readline().strip())
        text = str(file.readline().strip())
        d = int(file.readline().strip())

    # write
    with open(output_filename, 'w') as file:
        file.write(' '.join(map(str, get_approx_occur_pos(pattern, text, d))))
