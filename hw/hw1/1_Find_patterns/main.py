import sys


def pattern_count(text, pattern):
    count = 0
    pattern_size = len(pattern)
    for i in range(len(text) - pattern_size):
        if text[i:i + pattern_size] == pattern:
            count += 1
    return count


def t_frequent_words(text, k, t):
    counts_cache = dict()
    # counts = []
    for i in range(len(text) - k):
        pattern = text[i:i + k]
        if pattern not in counts_cache:
            counts_cache[pattern] = pattern_count(text, pattern)
        # counts.append(counts_cache[pattern])
    # print(*zip(text, counts), sep='\n')

    return sorted([pattern for pattern, count in counts_cache.items() if count >= t])


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        text = str(file.readline().strip())
        k, l, t = map(int, file.readline().strip().split())

    # write
    with open(output_filename, 'w') as file:
        file.write(' '.join(map(lambda seq: ''.join(map(str, seq)), t_frequent_words(text, k, t))))
