import pandas
import time
import re

start = time.time()

# input_file = 'example_data_2_masked.txt'
# masked = pandas.read_csv(input_file, delim_whitespace=True, header=None) # This is a pandas dataframe

# num_snps = masked.shape[0]
# # Impute genotypes based on most likely

# # Make a series of each row in str form for easier/faster searching
# snp_strs = pandas.Series([''.join(masked.iloc[i,:].apply(str).values) for i in range(num_snps)])

# # Create a dictionary mapping individual (row) to SNP (col) which is location of *
# missing_indices = masked[masked=='*'].stack().index.tolist()
# snpdict = {key: list() for key in range(num_snps)}
# for index in missing_indices:
# 	snpdict[index[0]].append(index[1])

# # Search each individual and replace * with whichever more commonly appears after the 5 that precede it
# for row in snpdict.keys():
# 	for col in snpdict[row]:
# 		five = snp_strs[row][col-5:col]
# 		if five == '':
# 			five = snp_strs[row][col+1:col+6]
# 		five_zero = five + '0'
# 		five_two = five + '2'
# 		count_0 = snp_strs.str.count(re.escape(five_zero)).sum()
# 		count_2 = snp_strs.str.count(re.escape(five_two)).sum()
# 		imputed = '0' if count_0 >= count_2 else '2'
# 		snp_strs[row] = snp_strs[row][:col] + imputed + snp_strs[row][col+1:]
# 		masked.iloc[row, col] = imputed


# masked.to_csv('imputed_test1.csv', index=False, header=False, sep=' ')

# imp_time = time.time() - start
# print(imp_time)

# could probably combine methods by imputing columns? 
input_file = 'example_data_1.txt'
masked = pandas.read_csv(input_file, delim_whitespace=True, header=None) # This is a pandas dataframe
num_snps = masked.shape[0]
# Clark's for haplotype phasing with chunks of 4
chunk = 8
encountered = dict()
# encountered should be a dict of genotype: (hap1, hap2)
used = set()
sol_cols = [str(i//2) + '_' + str(i % 2 + 1) for i in range(masked.shape[1])]
solution = pandas.DataFrame(index = range(num_snps), columns = sol_cols)

genotype_strs = pandas.Series([''.join(masked.iloc[:, i].apply(str).values) for i in range(masked.shape[1])])
# check understanding: genotype_strs should have 50 members

def get_haplotypes(string, haplo):
	option = ''
	for index in range(len(string)):
		if string[index] == '0':
			option += '0'
		else:
			option += '1'
			if string[index] == '1':
				get_haplotypes(string[:index] + '0' + string[index+1:], haplo)
	haplo.add(option)
	return

def run_clarks(geno):
	if geno in encountered:
		return encountered[geno][0], encountered[geno][1]
	else:
		haplotypes = set()
		get_haplotypes(geno, haplotypes)
		once = []
		found = False
		for h in haplotypes:
			h_complement = ''.join([str(int(gchar) - int(hchar)) for (gchar, hchar) in zip(geno, h)])
			if found == False:
				if h in used:
					if h_complement in used:
						encountered[geno] = (h, h_complement)
						used.add(h)
						used.add(h_complement)
						found = True
						return h, h_complement
					else:
						once.append((h, h_complement))
		if found == False:
			for h in haplotypes:
				if found == False:
					if once != []:
						encountered[geno] = once[0]
						used.add(once[0][0])
						used.add(once[0][1])
						found = True
						return once[0][0], once[0][1]
					else:
						encountered[geno] = (h, h_complement)
						used.add(h)
						used.add(h_complement)
						found = True
						return h, h_complement
	return

skipped = []
count1 = 0
countd = 0
for col in range(len(genotype_strs)):
	for row in range(0, num_snps, chunk):
		current_gchunk = genotype_strs[col][row:row+chunk]
		if '1' in current_gchunk:
			skipped.append((row, col))
			count1+=1
			continue
		else:
			countd+=1
			h_list = [str(int(char)//2) for char in current_gchunk]
			h_string = ''.join(h_list)
			encountered[current_gchunk] = (h_string, h_string)
			used.add(h_string)
			solution.loc[[row+i for i in range(chunk)], str(col)+'_1'] = h_list
			solution.loc[[row+i for i in range(chunk)], str(col)+'_2'] = h_list

for index in skipped:
	row, col = index
	current_gchunk = genotype_strs[col][row:row+chunk]
	h, hc = run_clarks(current_gchunk)
	solution.loc[[row+i for i in range(chunk)], str(col)+'_1'] = list(h)
	solution.loc[[row+i for i in range(chunk)], str(col)+'_2'] = list(hc)

solution.to_csv('mysol_skip_test_genofix.csv', index=False, header=False, sep=' ')



#haplo.add(''.join([str(int(i)//2) for i in even]))



