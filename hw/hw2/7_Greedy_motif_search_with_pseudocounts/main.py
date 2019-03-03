import sys
import numpy as np
from collections import Counter


def make_profile(motifs, k):
    t_motifs = np.array([
        nuc for motif in motifs for nuc in motif
    ]).reshape((-1, k)).T
    profile = []
    for col in t_motifs:
        val, counts = np.unique(col, return_counts=True)
        profile.append(val[np.argmax(counts)])
    return profile


def dist(profile, motif):
    return sum(p != m for p, m in zip(profile, motif))


def get_most_likely_motif(motifs, line, k):
    t_motifs = np.array([
        nuc for motif in motifs for nuc in motif
    ]).reshape((-1, k)).T
    pseudo_counter = Counter(['A', 'C', 'G', 'T'])
    distr = [Counter(col) + pseudo_counter for col in t_motifs]
    probs = [
        np.prod([distr[j][nuc] for j, nuc in enumerate(line[i:i+k])])
        for i in range(line.shape[0] - k + 1)
    ]
    i = np.argmax(probs)
    likely_motif = line[i:i+k]
    return likely_motif


def score(motifs, k):
    profile = make_profile(motifs, k)
    return sum(dist(profile, motif) for motif in motifs)


def greedy_motif_search(dna, k, t):
    best_motifs = [line[:k] for line in dna]

    for i in range(dna.shape[1] - k + 1):
        motifs = [dna[0, i:i+k]]
        for line in dna[1:]:
            motif = get_most_likely_motif(motifs, line, k)
            motifs.append(motif)
        if score(motifs, k) < score(best_motifs, k):
            best_motifs = motifs

    return best_motifs


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        k, t = map(int, file.readline().strip().split())
        dna = []
        for _ in range(t):
            dna.extend([nuc for nuc in str(file.readline().strip())])
        dna = np.array(dna).reshape((t, -1))

    # write
    with open(output_filename, 'w') as file:
        file.write('\n'.join(map(lambda seq: ''.join(map(str, seq)), greedy_motif_search(dna, k, t))))
