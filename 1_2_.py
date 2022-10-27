#The reverse complement problem: Find the reverse complement of a DNA string
#Input: A DNA string 'pattern' ; Output: Pattern_rc, the reverse complement of 'pattern'

def reverse_complement(dna_string):
	"""
	Input: a DNA string
	Returns a string that is its reverse complement
	"""
	list_complementary_nucleotides = {'A':'T','T':'A','G':'C','C':'G'}
	reverse_dna_string = dna_string[::-1]
	reverse_complement = ''
	for nucleotide in reverse_dna_string:
		nucleotide_complement = list_complementary_nucleotides[nucleotide]
		reverse_complement = reverse_complement + nucleotide_complement
	return reverse_complement

# print(reverse_complement(''))

