genome_file = open('C:/Users/Siddharth Murali/Downloads/E_coli.txt', 'r')


# text_ori = input("give text sequence ").upper()
# pattern = input("give pattern ")

# def pattern_count(text,pattern):
# 	k = len(pattern)
# 	text_length = len(text)
# 	count = 0
# 	for i in range(0,text_length-k+1):
# 		# print(text[i:i+k])
# 		if text[i:i+k] == pattern :
# 			count += 1
# 	# print(count)
# 	return count	

# # pattern_count(text,pattern)										#---RUN line---	

# def most_freq_kmer(text, k):
# 	count_threshold = 0
# 	checked_patterns = []
# 	freq_patterns = []
# 	for i in range(0,len(text)-k+1) :
# 		pattern = text[i:i+k]
# 		if pattern not in checked_patterns :
# 			checked_patterns.append(pattern)
# 			count = pattern_count(text,pattern)
# 			# print(pattern)
# 			if count > count_threshold :
# 				count_threshold = count
# 				freq_patterns = [pattern]	
# 			elif count == count_threshold :
# 				freq_patterns.append(pattern)
# 	print(freq_patterns)
# 	return freq_patterns

# def most_freq_kmer_dic(text, k): #same as most_freq_kmer but using python dictionaries
# 	counts_patterns_dic = {}
# 	freq_patterns = []
# 	max_count = 0
# 	for i in range(0,len(text)-k+1):
# 		pattern = text[i:i+k]
# 		if pattern not in counts_patterns_dic :
# 			counts_patterns_dic[pattern] = pattern_count(text,pattern)
# 			if counts_patterns_dic[pattern] == max_count :
# 				print(pattern)
# 				freq_patterns.append(pattern)
# 			elif counts_patterns_dic[pattern] > max_count :
# 				freq_patterns = [pattern]
# 			max_count = max(counts_patterns_dic.values())
# 	# print(freq_patterns)
# 	# print(counts_patterns_dic)
# 	return freq_patterns

# most_freq_kmer_dic(text,3)											#---RUN line---	


def frequency_table(text, k):
	freq_map = {}
	for i in range(0, len(text) - k + 1):
		pattern = text[i:i + k]
		if pattern not in freq_map:
			freq_map[pattern] = 1
		elif pattern in freq_map:
			freq_map[pattern] += 1
	# print(freq_map)
	return freq_map


# print(frequency_table(text,9))																		#---RUN line---

# def frequency_table_rc(text, k):
# 	"""
# 	Input: a string and natural number k
# 	Returns frequency map of the frequencies of each pattern(wither as itself or its reverse compliment)
# 	"""
# 	freq_map = {}
# 	for i in range(0,len(text)-k+1):
# 		pattern = text[i:i+k]
# 		# print(pattern)
# 		rc_pattern = reverse_complement(pattern)
# 		# print(pattern + ' ' + rc_pattern)
# 		if pattern not in freq_map and rc_pattern not in freq_map:
# 			freq_map[pattern] = 1
# 		elif pattern in freq_map:
# 			freq_map[pattern] += 1
# 		elif rc_pattern in freq_map:
# 			# print('aye')
# 			freq_map[rc_pattern + ' & ' + pattern] = freq_map[rc_pattern] + 1
# 	# print(freq_map)
# 	return freq_map


# print(frequency_table_rc(text, 9))														#---RUN line---	 

def better_freq_patterns(text, k):
	freq_patterns = []
	freq_map = frequency_table(text, k)
	max_count = max(freq_map.values())
	for pattern in freq_map.keys():
		if freq_map[pattern] == max_count:
			freq_patterns.append([pattern, max_count])
	# print(freq_patterns)
	return freq_patterns


def better_freq_patterns_mincount(text, k, min_count):
	freq_patterns = {}
	freq_map = frequency_table(text, k)
	for pattern in freq_map.keys():
		frequency = freq_map[pattern]
		if frequency >= min_count:
			freq_patterns[pattern] = frequency
	# print(freq_patterns)
	return freq_patterns


# print(better_freq_patterns_mincount('', 9, 3))

# print(better_freq_patterns(text_ori, 9))										#---RUN line---


def reverse_complement(dna_string):
	"""
	Input: a DNA string
	Returns a string that is its reverse complement
	"""
	list_complementary_nucleotides = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
	reverse_dna_string = dna_string[::-1]
	# print(reverse_dna_string)
	bases = list(reverse_dna_string)
	bases = [list_complementary_nucleotides[base] for base in bases]  # list comprehension
	return ''.join(bases)


# reverse_complement = ''
# for nucleotide in reverse_dna_string:
# 	# print(nucleotide)
# 	nucleotide_complement = list_complementary_nucleotides[nucleotide]
# 	reverse_complement = reverse_complement + nucleotide_complement
# return reverse_complement

# print(reverse_complement('GCTAGCT'))													 #---RUN line---

def better_freq_patterns_rc(text, k, min_count):
	'''
	input: dna string, k-mer and minimum count/repeats of k-mer to look for
	output: freq map with the most frequent patterns(inclusive of their reverse complements)
	'''
	freq_patterns = []
	freq_map = frequency_table(text, k)
	# print(freq_map.keys())
	used_patterns = []
	for pattern in freq_map.keys():
		used_patterns.append(pattern)
		rc_pattern = reverse_complement(pattern)
		# if pattern == 'ATGATCAAG' :
		# 	print(rc_pattern)
		if rc_pattern in freq_map.keys() and rc_pattern not in used_patterns:
			updated_count = freq_map[pattern] + freq_map[rc_pattern]
			if updated_count >= min_count:
				freq_patterns.append([pattern + '*', updated_count])
				freq_patterns.append([rc_pattern + '**', updated_count])
		elif freq_map[pattern] >= min_count:
			freq_patterns.append([pattern, freq_map[pattern]])
	return freq_patterns


# print(better_freq_patterns_rc(text, 9, 3))

# def better_freq_patterns(text, k):
# 	freq_patterns = []
# 	freq_map = frequency_table(text,k)
# 	for pattern in freq_map.keys():
# 		reverse_complement = reverse_complement(pattern)
# 		if freq_map[reverse_complement] in freq_map:
# 			freq_map[pattern] += freq_map[reverse_complement]
# 			freq_map[reverse_complement] = freq_map[pattern]
# 	max_count = max(freq_map.values())

# 	# print(freq_patterns)
# 	return freq_patterns
###


# The reverse complement problem: Find the reverse complement of a DNA string
# Input: A DNA string 'pattern' ; Output: Pattern_rc, the reverse complement of 'pattern'


# def freq_patterns_with_rc(text, k):
# 	"""
# 	Input: a string and integer k
# 	Returns most frequent k-mers(after clubbing them with their reverse complements if they exist)
#	Checks only the max count values from freq_map
# 	"""
# 	text = text.upper()
# 	# print(text)
# 	freq_patterns = []
# 	freq_map = frequency_table(text,k)
# 	max_count = max(freq_map.values())
# 	for pattern in freq_map.keys():
# 		# print(pattern)
# 		if freq_map[pattern] >= 3 :
# 			freq_patterns.append(pattern)
# 	rc_freq_pattern = {}
# 	for key in freq_patterns:
# 		print(key)
# 		reverse_complement_key = reverse_complement(key)
# 		# rc_freq_pattern = {key + ',' + reverse_complement_key : freq_map[key] + freq_map[reverse_complement_key] for reverse_complement_key in freq_patterns}		
# 		if reverse_complement_key in freq_patterns:
# 			print('aye')
# 			freq_patterns.remove(reverse_complement_key)
# 			rc_freq_pattern[key + ',' + reverse_complement_key] = freq_map[key] + freq_map[reverse_complement_key]
# 		# else:
# 			# return freq_patterns
# 	return rc_freq_pattern


# print(freq_patterns_with_rc(text, 9))										#---RUN line---

# pattern matching problem

def find_occ_pattern(pattern, genome):
	"""
	Input: Strings Pattern and Genome.
	Output: All starting positions in Genome where Pattern appears as a substring.
	"""
	# genome = genome.read()
	gen_length = len(genome)
	pattern_length = len(pattern)
	pattern_positions = []
	j = 0
	for i in range(0, gen_length - pattern_length + 1):
		if pattern == genome[i:i + pattern_length]:
			pattern_positions.append(i)
			j += 1
	# print(j)
	return pattern_positions


# print(find_occ_pattern('CGC', 'ATGACTTCGCTGTTACGCGC'))


# ----------------------------------------------------------------------end of lesson 1.3------------------------------------------------------------------

# Clump finding problem

def clumping_patterns(genome, k, L, t):
	'''
	Input: a string of genome, 3 integers k, L, t; k is the length of pattern/k-mer, L is the length of region where we wish to find a clump, t is the minimum number of times pattern/kmer needs to occur
	Output: Returns all distinct k-mers forming (L,t) clumps
	'''
	genome = genome_file.read()
	gen_length = len(genome)
	# print(gen_length)
	kmer_clumps = []

	# freq_patterns = better_freq_patterns_mincount(genome, k, t)
	# for pattern in freq_patterns.keys():
	# 	pattern_positions = find_occ_pattern(pattern, genome)
	# 	for i in range(0,len(pattern_positions) - t):
	# 			# print(i)
	# 			# print(pattern_positions[i+t-1])
	# 			# print(pattern_positions[i])
	# 			if pattern_positions[i + t] - pattern_positions[i] <= L:
	# 				kmer_clumps.append(pattern)
	# 				break

	# freq_map = frequency_table(genome, k)
	# # print(len(freq_map))
	# for pattern in freq_map:
	# 	if freq_map[pattern] >= t:
	# 		pattern_positions = find_occ_pattern(pattern, genome)
	# 		# print(pattern_positions)
	# 		# print(pattern)
	# 		# print(len(pattern_positions))
	# 		for i in range(0,len(pattern_positions) - t):
	# 			# print(i)
	# 			# print(pattern_positions[i+t-1])
	# 			# print(pattern_positions[i])
	# 			if pattern_positions[i + t] - pattern_positions[i] <= L:
	# 				kmer_clumps.append(pattern)
	# 				break

	# for i in range(0, gen_length - L + 1):
	# 	region = text[i:i+L]
	# 	for pattern in freq_map.keys():
	# 		freq = freq_map[pattern]
	# 		if freq_map[pattern] >= t and pattern not in kmer_clumps:
	# 			kmer_clumps.append(pattern)

	for i in range(0, gen_length - L + 1):
		region = genome[i:i + L]
		freq_map = frequency_table(region, k)
		for pattern in freq_map.keys():
			frequency = freq_map[pattern]
			if frequency >= t and pattern not in kmer_clumps:
				kmer_clumps.append(pattern)
	# print("number of kmers = " + str(len(kmer_clumps)))
	return kmer_clumps


# with open("D:\Productive\Courses\Bioinformatics Algorithms\output.txt", 'a') as file:
# 	for item in clumping_patterns(genome_file, 9, 500, 3):
# 		file.write(item + ", ")
# frequency_table(genome, 9)

# print(clumping_patterns('',10,24,4))

def skew_i(text):
	"""
	Input DNA string 
	Output: Returns dictionary of values at each i of string
	"""
	skew_i = {}
	for i in range(0, len(text)):
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


# print(skew_i(genome_file))                                        					#---RUN LINE---

# print(find_occ_pattern('CTTGATCAT',genome_file))

def hamming_distance(x, y):
	"""
	Input: two strings x and y of equal length
	Output: integer, hamming distance between the two strings
	"""
	h_d = 0
	for i in range(len(x)):
		if x[i] != y[i]:
			h_d += 1
	return h_d


# print(hamming_distance('CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG','CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA'))


def approx_pattern_match(pattern, text, d):
	'''
	input strings pattern and text, and integer d
	output: starting positions of approx matches (hamming distance <= d) as a list of integers
	'''
	len_text = len(text)
	len_pattern = len(pattern)
	occurrences = []
	for i in range(len(text) - len(pattern) + 1):
		boxed_string = text[i:i + len(pattern)]
		if hamming_distance(pattern, boxed_string) <= d:
			occurrences.append(i)
	# print(len(occurrences))
	return occurrences


# print(*approx_pattern_match('AA', 'TACGCATTACAAAGCACA', 1))


def neighbours(pattern, d):
	"""
	Input: a string pattern and a positive integer d
	Output: the list of all strings that are at most a hamming distance d away from pattern
	---Recursive function---
	"""
	if d == 0:
		return [pattern]
	len_pattern = len(pattern)
	if len_pattern == 1:
		return ['A', 'T', 'G', 'C']
	neighbourhood = []
	suffix_pattern = pattern[1:len_pattern]
	suffix_neighbourhood = neighbours(suffix_pattern, d)
	for string in suffix_neighbourhood:
		h_d = hamming_distance(suffix_pattern, string)
		if h_d < d:
			for letter in ['A', 'T', 'G', 'C']:
				neighbourhood.append(letter + string)
		elif h_d == d:
			neighbourhood.append(pattern[0] + string)
	return neighbourhood


# list_str_neighbours = '\n'.join(map(str, neighbours('TGCAT', 2)))
# print(len(neighbours('TGCAT', 2)))


def frequent_words_with_mismatches(text, k, d):
	'''
	Input: takes in a string text, integer k (k-mer) and d for maximum hamming distance
	Output: most freq k-mers with matches upto d mismatches in text, a list
	'''
	patterns = []
	freq_map = {}
	len_text = len(text)
	for i in range(len_text - k + 1):
		pattern = text[i:i + k]
		neighbourhood = neighbours(pattern, d)
		# print(pattern, neighbourhood)
		for neighbour in neighbourhood:
			rc_neighbour = reverse_complement(neighbour)
			if neighbour not in freq_map.keys() and rc_neighbour not in freq_map.keys():
				freq_map[neighbour] = 1
				freq_map[rc_neighbour] = 1
			else:
				freq_map[neighbour] += 1
				freq_map[rc_neighbour] += 1
	m = max(freq_map.values())
	for pattern in freq_map.keys():
		if freq_map[pattern] == m:
			patterns.append(pattern)
	return patterns

# print(*frequent_words_with_mismatches('AATGATGATGACGTCAAAAGGATCCGGATAAAACATGGTGATTGCCTCGCATAACGCGGTATGAAAATGGATTGAAGCCCGGGCCGTGGATTCTACTCAACTTTGTCGGCTTGAGAAAGACCTGGGATCCTGGGTATTAAAAAGAAGATCTATTTATTTAGAGATCTGTTCTATTGTGATCTCTTATTAGGATCGCACTGCCCTGTGGATAACAAGGATCCGGCTTTTAAGATCAACAACCTGGAAAGGATCATTAACTGTGAATGATCGGTGATCCTGGACCGTATAAGCTGGGATCAGAATGAGGGGTTATACACAACTCAAAAACTGAACAACAGTTGTTCTTTGGATAACTACCGGTTGATCCAAGCTTCCTGACAGAGTTATCCACAGTAGATCGCACGATCTGTATACTTATTTGAGTAAATTAACCCACGATCCCAGCCATTCTTCTGCCGGATCTTCCGGAATGTCGTGATCAAGAATGTTGATCTTCAGTG', 9, 1))
