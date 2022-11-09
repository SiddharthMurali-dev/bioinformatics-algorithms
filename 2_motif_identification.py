from bioinformatics_1_replication_origins import hamming_distance, neighbours
"""
Given a collection of strings Dna and an integer d, a k-mer is a (k,d)-motif if it appears in every string from Dna with
at most d mismatches.

Implanted Motif Problem: Find all (k, d)-motifs in a collection of strings.
    Input: A collection of strings 'Dna', and integers k and d.
    Output: All (k, d)-motifs in 'Dna'.
    ---Brute force algorithm for finding motifs---
"""


def motif_enumeration(k, d, dna):
    """
    :param k: length of pattern
    :param d: number of mismatches tolerated
    :param dna: space-separated collection of strings
    :return: All (k, d) motifs
    """
    patterns = []
    dna_list = dna.split(" ")
    first_string = dna_list[0]
    for i in range(len(first_string) - k + 1):
        k_mer = first_string[i: i + k]
        for neighbour in neighbours(k_mer, d):
            dna_list_len = len(dna_list)
            for string in dna_list:
                for j in range(len(string) - k + 1):
                    k_mer_2 = string[j: j + k]
                    if hamming_distance(neighbour, k_mer_2) <= d:
                        dna_list_len -= 1
                        break
                if dna_list_len == 0 and neighbour not in patterns:
                    patterns.append(neighbour)
    return ' '.join(patterns)


patterns_string = motif_enumeration(5, 2, "AGTTAGTTATCAACGACCCATCCGC ACTTGTATTCGCAATGACCTGCTCA TGCGGAGAGTGCTATACACCTAGCG AAACGGAAATTACCGCTTACACCTA AATGGGTCGTGCATGGGGATGCGAT AGATCGGTTGGGTATCCCGTTCCCA")
print(patterns_string)
