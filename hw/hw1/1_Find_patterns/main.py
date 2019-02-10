import sys
from collections import Counter


def create_freq_dict(text, k):
    return Counter([text[i:i + k] for i in range(len(text) - k)])


def update_freq_patterns(freq_patterns, freq_dict, t):
    freq_patterns.update({
        pattern for pattern, count in freq_dict.items() if count >= t
    })


def find_clumps(genome, k, t, l):
    freq_patterns = set()
    freq_dict = create_freq_dict(genome[:l], k)
    update_freq_patterns(freq_patterns, freq_dict, t)

    for i in range(1, len(genome) - l):
        old_pattern = genome[(i - 1):(i - 1 + k)]
        new_pattern = genome[(i + l - k):(i + l)]

        freq_dict[old_pattern] -= 1
        freq_dict[new_pattern] += 1
        if freq_dict[new_pattern] >= t:
            freq_patterns.add(new_pattern)

    return sorted(freq_patterns)


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        genome = str(file.readline().strip())
        k, l, t = map(int, file.readline().strip().split())

    # write
    with open(output_filename, 'w') as file:
        file.write(' '.join(map(lambda seq: ''.join(map(str, seq)), find_clumps(genome, k, t, l))))
