import sys

NUCS = ['A', 'C', 'G', 'T']
COMPL_NUCS = dict(
    A='T',
    C='G',
    G='C',
    T='A',
)


def get_compl_pair(mer):
    return ''.join([COMPL_NUCS[nuc] for nuc in mer])


def get_reverse_compl_pair(mer):
    return get_compl_pair(mer[::-1])


def approx_compare(pattern, sub_text, d):
    return sum(pattern[i] != sub_text[i] for i in range(len(pattern))) <= d


def gen_seq(k, prefix=None):
    if prefix is None:
        prefix = []

    if k <= 0:
        yield prefix
    else:
        for nuc in NUCS:
            prefix.append(nuc)
            yield from gen_seq(k - 1, prefix)
            prefix.pop()


def get_freq_approx_mers(text, k, d):
    occurrences = {''.join(mer): 0 for mer in gen_seq(k)}
    for i in range(len(text) - k):
        for mer in occurrences.keys():
            if approx_compare(mer, text[i:i + k], d):
                occurrences[mer] += 1

    max_occur = max([
        occur + occurrences[get_reverse_compl_pair(mer)]
        for mer, occur in occurrences.items()
    ])
    freq_mers = [
        mer
        for mer, occur in occurrences.items()
        if (occur + occurrences[get_reverse_compl_pair(mer)]) == max_occur
    ]
    return freq_mers


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        text = str(file.readline().strip())
        k, d = map(int, file.readline().strip().split())

    # write
    with open(output_filename, 'w') as file:
        file.write(' '.join(map(str, get_freq_approx_mers(text, k, d))))
