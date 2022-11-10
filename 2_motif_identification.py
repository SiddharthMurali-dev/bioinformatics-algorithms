import math
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


# patterns_string = motif_enumeration(5, 2, "")
# print(patterns_string)

# --Scoring Motifs--

def freq_nt_col(col):
    """
    :param col: list of nts in a column of the motif matrix
    :return: list of frequencies--f(A,C,G,T)
    """
    len_col = len(col)
    nt_list = ["A", "C", "G", "T"]
    # print(col)
    freqs = [col.count(nt)/len_col for nt in nt_list]
    # print(freqs)
    return freqs


# print(freq_nt_col(['A', 'G', 'C', 'A']))


def entropy(matrix):
    """
    :param matrix: a list of rows of the motif matrix
    :return: entropy of the matrix
    """
    entropy_matrix = 0
    for i in range(len(matrix[0])):
        freqs_i = freq_nt_col([row[i] for row in matrix])
        entropy_i = sum(math.log(f, 2)*f for f in freqs_i if f != 0)
        # print(entropy_i)
        entropy_matrix += entropy_i

    return entropy_matrix * -1


# print(entropy([ "TCGGGGGTTTTT", "CCGGTGACTTAC", "ACGGGGATTTTC", "TTGGGGACTTTT", "AAGGGGACTTCC", "TTGGGGACTTCC",
#                 "TCGGGGATTCAT", "TCGGGGATTCCT", "TAGGGGAACTAC", "TCGGGTATAACC"]))
