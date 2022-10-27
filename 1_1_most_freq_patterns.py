# text = input("give text sequence ")
# pattern = input("give pattern ")
def pattern_count(text,pattern):
	k = len(pattern)
	text_length = len(text)
	count = 0
	for i in range(0,text_length-k+1):
		# print(text[i:i+k])
		if text[i:i+k] == pattern :
			count += 1
	# print(count)
	return count	

# print(pattern_count('GACCATCAAAACTGATAAACTACTTAAAAATCAGT', 'AAA'))

def most_freq_kmer(text, k):
	count_threshold = 0
	checked_patterns = []
	freq_patterns = []
	for i in range(0,len(text)-k+1) :
		pattern = text[i:i+k]
		if pattern not in checked_patterns :
			checked_patterns.append(pattern)
			count = pattern_count(text,pattern)
			# print(pattern)
			if count > count_threshold :
				count_threshold = count
				freq_patterns = [pattern]	
			elif count == count_threshold :
				freq_patterns.append(pattern)
	print(freq_patterns)
	return freq_patterns

def most_freq_kmer_dic(text, k): #same as most_freq_kmer but using python dictionaries
	counts_patterns_dic = {}
	freq_patterns = []
	max_count = 0
	for i in range(0,len(text)-k+1):
		pattern = text[i:i+k]
		if pattern not in counts_patterns_dic :
			counts_patterns_dic[pattern] = pattern_count(text,pattern)
			if counts_patterns_dic[pattern] == max_count :
				# print(pattern)
				freq_patterns.append(pattern)
			elif counts_patterns_dic[pattern] > max_count :
				freq_patterns = [pattern]
			max_count = max(counts_patterns_dic.values())
	# print(freq_patterns)
	# print(counts_patterns_dic)
	return freq_patterns

# print(most_freq_kmer_dic('TAAACGTGAGAGAAACGTGCTGATTACACTTGTTCGTGTGGTAT?',3))

# frequency_table and better_freq_patterns are the modified versions of the above code that have a lower complexity. The following code runs independent of the above functions.
def frequency_table(text, k):
	freq_map = {}
	for i in range(0,len(text)-k+1):
		pattern = text[i:i+k]
		if pattern not in freq_map:
			freq_map[pattern] = 1
		elif pattern in freq_map:
			freq_map[pattern] += 1
	# print(freq_map)
	return freq_map

def better_freq_patterns(text, k):
	freq_patterns = []
	freq_map = frequency_table(text,k)
	max_count = max(freq_map.values())
	for pattern in freq_map.keys():
		if freq_map[pattern] == max_count :
			freq_patterns.append([pattern,max_count])
	# print(freq_patterns)
	return freq_patterns

# print(better_freq_patterns(text, 9))
# print(most_freq_kmer_dic(text,12))