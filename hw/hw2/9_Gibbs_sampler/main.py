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


def get_closest_motif(profile, line, k):
    i = sorted([
        (i, dist(profile, line[i:i+k]))
        for i in range(line.shape[0] - k + 1)
    ], key=lambda pair: pair[1])[0][0]
    return line[i:i+k]


def score(motifs, k):
    profile = make_profile(motifs, k)
    return sum(dist(profile, motif) for motif in motifs)


    # GIBBSSAMPLER(Dna, k, t, N)
    #     randomly select k-mers Motifs = (Motif1, …, Motift) in each string
    #         from Dna
    #     BestMotifs ← Motifs
    #     for j ← 1 to N
    #         i ← Random(t)
    #         Profile ← profile matrix constructed from all strings in Motifs
    #                    except for Motifi
    #         Motifi ← Profile-randomly generated k-mer in the i-th sequence
    #         if Score(Motifs) < Score(BestMotifs)
    #             BestMotifs ← Motifs
    #     return BestMotifs

def get_random_motif(motifs, line, k):
    t_motifs = np.array([
        nuc for motif in motifs for nuc in motif
    ]).reshape((-1, k)).T
    pseudo_counter = Counter(['A', 'C', 'G', 'T'])
    distr = [Counter(col) + pseudo_counter for col in t_motifs]
    # distr = [Counter(col) for col in t_motifs]
    probs = [
        np.prod([distr[j][nuc] for j, nuc in enumerate(line[i:i+k])])
        for i in range(line.shape[0] - k + 1)
    ]
    i = np.random.choice(range(line.shape[0] - k + 1), p=probs/np.sum(probs))
    random_motif = line[i:i+k]
    return random_motif


def gibbs_sampler(dna, k, t, n):
    motifs = [
        dna[i][j:j+k]
        for i, j in enumerate(
            np.random.randint(0, dna.shape[1] - k, dna.shape[0])
        )
    ]
    best_motifs = motifs
    for _ in range(n):
        i = np.random.randint(0, dna.shape[0])
        motifs[i] = get_random_motif(motifs[:i] + motifs[i+1:], dna[i], k)

        if score(motifs, k) < score(best_motifs, k):
            best_motifs = motifs
    return best_motifs


if __name__ == "__main__":
    input_filename = sys.argv[1]
    output_filename = 'output.txt'

    # read
    with open(input_filename, 'r') as file:
        k, t, n = map(int, file.readline().strip().split())
        dna = []
        for _ in range(t):
            dna.extend([nuc for nuc in str(file.readline().strip())])
        dna = np.array(dna).reshape((t, -1))

    # random recalling
    best_motifs = gibbs_sampler(dna, k, t, n)
    for _ in range(15):
        motifs = gibbs_sampler(dna, k, t, n)
        if score(motifs, k) < score(best_motifs, k):
            best_motifs = motifs
    # write
    with open(output_filename, 'w') as file:
        file.write('\n'.join(map(lambda seq: ''.join(map(str, seq)), best_motifs)))
