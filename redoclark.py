import pandas
import time
import re

input_file = 'first_imputed_test1.csv'
masked = pandas.read_csv(input_file, delim_whitespace=True, header=None) # This is a pandas dataframe
num_snps = masked.shape[0]


# Make 2 chunks dictionaries of (row, col): chunk

def get_chunks(chunk_size, genotypes):
	chunks_ambig = dict()
	chunks_unambig = dict()
	for g in range(len(genotypes)):
		for index in range(0, num_snps, chunk_size):
			substr = genotypes[g][index:index+chunk_size]
			if '1' in substr:
				chunks_ambig[(index,g)] = substr
			else:
				chunks_unambig[(index, g)] = substr
	return chunks_ambig, chunks_unambig

def add_to_sol(to_add, index, chunk_size, ind, output_df):
	output_df.loc[[index[0]+i for i in range(chunk_size)], str(index[1])+ind] = to_add	


def clarks_init(unambig_dict, chunk_size, output_df):
	encountered_genotypes = dict()
	used_haplotypes = set()
	for index in unambig_dict:
		geno = unambig_dict[index]
		haplo = [str(int(char)//2) for char in geno]
		haplo_str = ''.join(haplo)
		encountered_genotypes[geno] = (haplo_str, haplo_str)
		used_haplotypes.add(haplo_str)
		add_to_sol(haplo, index, chunk_size, '_1', output_df)
		add_to_sol(haplo, index, chunk_size, '_2', output_df)
	return encountered_genotypes, used_haplotypes

def get_haplotypes(string, haplo_set):
	option = ''
	for index in range(len(string)):
		if string[index] == '0':
			option += '0'
		else:
			option += '1'
			if string[index] == '1':
				get_haplotypes(string[:index] + '0' + string[index+1:], haplo_set)
	haplo_set.add(option)
	return

def check_membership(genotype, haplotype, complement, encountered, used):
	if genotype in encountered:
		return encountered[genotype]
	elif haplotype in used and complement in used:
		encountered[genotype] = (haplotype, complement)
		return (haplotype, complement)
	return False

def clarks_iterative(encountered, used, ambig_dict, chunk_size, output_df):
	for index in ambig_dict:
		geno = ambig_dict[index]
		haplotypes = set()
		try_later = dict()
		get_haplotypes(geno, haplotypes)
		for haplo in haplotypes:
			complement = ''.join([str(int(gchar) - int(hchar)) for (gchar, hchar) in zip(geno, haplo)])
			check = check_membership(geno, haplo, complement, encountered, used)
			if check != False:
				add_to_sol(list(haplo), index, chunk_size, '_1', output_df)
				add_to_sol(list(complement), index, chunk_size, '_2', output_df)
				break
			else:
				if haplo in used or complement in used:
					add_to_sol(list(haplo), index, chunk_size, '_1', output_df)
					add_to_sol(list(complement), index, chunk_size, '_2', output_df)
					encountered[geno] = (haplo, complement)
					if haplo in used:
						used.add(complement)
					elif complement in used:
						used.add(haplo)
				else:
					try_later[index] = (geno, haplo, complement)
	return try_later

def clarks_repeat(encountered, used, to_try, chunk_size, output_df):
	for index in to_try:
		haps = to_try[index]
		geno = haps[0]
		haplo = haps[1]
		complement = haps[2]
		check = check_membership(geno, haplo, complement, encountered, used)
		if check != False:
			add_to_sol(list(haplo), index, chunk_size, '_1', output_df)
			add_to_sol(list(complement), index, chunk_size, '_2', output_df)
			to_try.pop(index)
			break
		else:
			if haplo in used or complement in used:
				add_to_sol(list(haplo), index, chunk_size, '_1', output_df)
				add_to_sol(list(complement), index, chunk_size, '_2', output_df)
				encountered[geno] = (haplo, complement)
				to_try.pop(index)
				if haplo in used:
					used.add(complement)
				elif complement in used:
					used.add(haplo)
	return to_try

def run(genotype_strs):
	ambiguous, unambiguous = get_chunks(4, genotype_strs)
	sol_cols = [str(i//2) + '_' + str(i % 2 + 1) for i in range(masked.shape[1])]
	solution = pandas.DataFrame(index = range(num_snps), columns = sol_cols)
	encountered, used = clarks_init(unambiguous, 4, solution)
	unmatched = clarks_iterative(encountered, used, ambiguous, 4, solution)
	if len(unmatched) != 0:
		prev_unmatched = None
		while len(prev_unmatched) != len(unmatched):
			prev_unmatched = unmatched
			unmatched = clarks_repeat(encountered, used, prev_unmatched, 4, solution)
			if len(unmatched) == 0:
				break
		if len(unmatched) != 0:
			for key in unmatched:
				encountered[unmatched[key][0]] = (unmatched[key][1], unmatched[key][2])
				used.add(unmatched[key][1])
				used.add(unmatched[key][2])
				add_to_sol(list(unmatched[key][1]), key, 4, '_1', solution)
				add_to_sol(list(unmatched[key][2]), key, 4, '_2', solution)
				break
		print("Sanity Check: Exited while loop.")
	solution.to_csv('clark_myimp.csv', index=False, header=False, sep=' ')

genotype_strs = pandas.Series([''.join(masked.iloc[:, i].apply(str).values) for i in range(masked.shape[1])])
run(genotype_strs)