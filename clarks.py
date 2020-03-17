import pandas
import time
import re

start = time.time()

input_file = 'example_data_2_masked.txt'
masked = pandas.read_csv(input_file, delim_whitespace=True, header=None) # This is a pandas dataframe

num_snps = masked.shape[0]
# Impute genotypes based on most likely

# Make a series of each row in str form for easier/faster searching
snp_strs = pandas.Series([''.join(masked.iloc[i,:].apply(str).values) for i in range(num_snps)])

# Create a dictionary mapping individual (row) to SNP (col) which is location of *
missing_indices = masked[masked=='*'].stack().index.tolist()
snpdict = {key: list() for key in range(num_snps)}
for index in missing_indices:
	snpdict[index[0]].append(index[1])

# Search each individual and replace * with whichever more commonly appears after the 5 that precede it
for row in snpdict.keys():
	for col in snpdict[row]:
		five = snp_strs[row][col-5:col]
		if five == '':
			five = snp_strs[row][col+1:col+6]
		five_zero = five + '0'
		five_two = five + '2'
		count_0 = snp_strs.str.count(re.escape(five_zero)).sum()
		count_2 = snp_strs.str.count(re.escape(five_two)).sum()
		imputed = '0' if count_0 >= count_2 else '2'
		snp_strs[row] = snp_strs[row][:col] + imputed + snp_strs[row][col+1:]
		masked.iloc[row, col] = imputed


masked.to_csv('imputed_test2.csv', index=False, header=False, sep=' ')

imp_time = time.time() - start
print(imp_time)

# could probably combine methods by imputing columns? 
# input_file = 'example_data_1.txt'
# masked = pandas.read_csv(input_file, delim_whitespace=True, header=None) # This is a pandas dataframe
# num_snps = masked.shape[0]
# # Clark's for haplotype phasing with chunks of 4
# chunk = 4
# encountered = dict()
# # encountered should be a dict of genotype: (hap1, hap2)
# used = set()
# sol_cols = [str(i//2) + '_' + str(i % 2 + 1) for i in range(masked.shape[1])]
# solution = pandas.DataFrame()

# genotype_strs = pandas.Series([''.join(masked.iloc[:, i].apply(str).values) for i in range(masked.shape[1])])
# # check understanding: genotype_strs should have 50 members

# def get_haplotypes(string, haplo):
# 	option = ''
# 	for index in range(len(string)):
# 		if string[index] == '0':
# 			option += '0'
# 		else:
# 			option += '1'
# 			if string[index] == '1':
# 				get_haplotypes(string[:index] + '0' + string[index+1:], haplo)
# 	haplo.add(option)
# 	return

# def run_clarks(geno, list1, list2)
# 	if current_gchunk in encountered:
# 			list1.extend(list(encountered[current_gchunk][0]))
# 			list2.extend(list(encountered[current_gchunk][1]))
# 		else:
# 			haplotypes = set()
# 			get_haplotypes(current_gchunk, haplotypes)
# 			once = []
# 			found = False
# 			for h in haplotypes:
# 				h_complement = ''.join([str(int(gchar) - int(hchar)) for (gchar, hchar) in zip(current_gchunk, h)])
# 				if found == False:
# 					if h in used:
# 						if h_complement in used:
# 							encountered[current_gchunk] = (h, h_complement)
# 							list1.extend(list(h))
# 							list2.extend(list(h_complement))
# 							used.add(h)
# 							used.add(h_complement)
# 							found = True
# 						else:
# 							once.append((h, h_complement))
# 			if found == False:
# 				for h in haplotypes:
# 					if found == False:
# 						if once != []:
# 							encountered[current_gchunk] = once[0]
# 							list1.extend(list(once[0][0]))
# 							list2.extend(list(once[0][1]))
# 							used.add(once[0][0])
# 							used.add(once[0][1])
# 							found = True
# 						else:
# 							encountered[current_gchunk] = (h, h_complement)
# 							list1.extend(list(h))
# 							list2.extend(list(h_complement))
# 							used.add(h)
# 							used.add(h_complement)
# 							found = True
# 	return

# for col in range(len(genotype_strs)):
# 	individual1 = []
# 	individual2 = []
# 	skipped = []
# 	for row in range(0, num_snps, chunk):
# 		current_gchunk = genotype_strs[col][row:row+chunk]
# 		if '1' in current_gchunk:
# 			skipped.append((row, col))
# 			continue
# 		individual1, individual2 = run_clarks(current_gchunk)
# 	solution[str(col)+'_1'] = individual1
# 	solution[str(col)+'_2'] = individual2

# for index in skipped:


# solution.to_csv('mysol_test1.csv', index=False, header=False, sep=' ')



# #haplo.add(''.join([str(int(i)//2) for i in even]))



