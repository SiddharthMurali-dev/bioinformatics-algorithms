def skew_i(text):
    """
    Input dna string
    Output: Returns dictionary of values at each i of string
    """
    skew_i = {}
    # print(len(text))
    for i in range(0, len(text) + 1):
        nuc = text[i - 1]
        if i == 0:
            skew_i[i] = 0
        elif nuc == 'C':
            skew_i[i] = skew_i[i - 1] - 1
        elif nuc == 'G':
            skew_i[i] = skew_i[i - 1] + 1
        else:
            skew_i[i] = skew_i[i - 1]
    return skew_i

print(skew_i('GCATACACTTCCCAGTAGGTACTG'))


#Minimum Skew Problem
def min_skew(dna_string):
    skew_i = {}
    # print(len(dna_string))
    for i in range(0, len(dna_string) + 1):
        nuc = dna_string[i - 1]
        if i == 0:
            skew_i[i] = 0
            min_skew_i = [0]
        elif nuc == 'C':
            skew_i[i] = skew_i[i - 1] - 1
        elif nuc == 'G':
            skew_i[i] = skew_i[i - 1] + 1
        else:
            skew_i[i] = skew_i[i - 1]
        value = skew_i[i]
        if i != 0 and value < skew_i[min_skew_i[0]]:
            min_skew_i = [i]
        elif i != 0 and value == skew_i[min_skew_i[0]]:
            min_skew_i.append(i)
        # print(min_skew_i)
    # print(skew_i)
    return min_skew_i


print(min_skew('TAAAGACTGCCGAGAGGCCAACACGAGTGCTAGAACGAGGGGCGTAAACGCGGGTCCGAT'))
