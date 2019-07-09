#pip3 install biopython
#pip3 install statsmodels
#pip3 install scipy

import sys
import Bio
import statsmodels
import math
import numpy as np
import pandas as pd
from operator import itemgetter
from statsmodels.compat.python import range
from statsmodels.compat.collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import defaultdict
from itertools import islice, chain, product, combinations_with_replacement
from time import process_time


########## main function ##############################################################################

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    nullomers_file = sys.argv[2]
    threshold = float(sys.argv[3])
    level = sys.argv[4]
    correction_method = sys.argv[5]
else:
    print("\n**An error occured**\n")
    raise SystemExit() 

if (correction_method != "FDR") and (correction_method != "ADJ-BONF") and (correction_method != "BONF"):
    print("\n**The correction method you declared is incorrect. Value should be either FDR or BONF or ADJ-BONF**\n")
    raise SystemExit()


if (level == "DNA"):
    alphabet = "TCGA"
elif (level == "PROT"):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
else:
    print("\n**The level you declared is incorrect. Value should be either DNA or PROT**\n")
    raise SystemExit()         
                                                 
#################################################################################################################


######### secondary functions ###################################################################################

def count_occurrences(subword, word):
    return word.count(subword)

def return_indices(array, value):
    return np.where(array==value)

def first_order_probs(sequences):   
    distinct = set(chain.from_iterable(sequences))
    n = len(distinct)
    coding = {j:i for i, j in enumerate(distinct)}
    counts = np.zeros((n, n))
    for seq in sequences:
        coded_seq = np.fromiter((coding[i] for i in seq), dtype=np.int64)
        pairs = coded_seq[:-1] + n * coded_seq[1:]
        counts += np.bincount(pairs, minlength=n*n).reshape(n, n)
    totals = counts.sum(axis=0)
    totals[totals == 0] = 1     # avoid division by zero
    probs = counts / totals
    return {a+b: p for a, b in product(distinct, repeat=2) for p in (probs[coding[b], coding[a]],) if p }

def second_order_probs(sequences):
    num = {}   
    for seq in sequences:
        for ind, current_letter in enumerate(str(seq)):
            if (ind == 0):
                two_previous_letters = current_letter
            elif (ind == 1):
                one_previous_letter = current_letter
            else:
                key = str(two_previous_letters)+str(one_previous_letter)+str(current_letter)
                if key in num:
                    num[key] += 1
                else:
                    num[key] = 1
                two_previous_letters = one_previous_letter
                one_previous_letter = current_letter
    denom = {}
    for key, value in (num.items()):
        if key[0:2] in denom:
           denom[key[0:2]] += value
        else:
           denom[key[0:2]] = value      
    second_order_dict = {key: value/value2 for key, value in num.items() for key2, value2 in denom.items() if key[0:2] == key2 }
    return {key: value for key, value in (second_order_dict.items()) }

def third_order_probs(sequences):
    num = {}   
    for seq in sequences:
        for ind, current_letter in enumerate(str(seq)):
            if (ind == 0):
                three_previous_letters = current_letter
            elif (ind == 1):
                two_previous_letters = current_letter
            elif (ind == 2):
                one_previous_letter = current_letter
            else:
                key = str(three_previous_letters)+str(two_previous_letters)+str(one_previous_letter)+str(current_letter)
                if key in num:
                    num[key] += 1
                else:
                    num[key] = 1
                three_previous_letters = two_previous_letters
                two_previous_letters = one_previous_letter
                one_previous_letter = current_letter
    denom = {}
    for key, value in (num.items()):
        if key[0:3] in denom:
           denom[key[0:3]] += value
        else:
           denom[key[0:3]] = value      
    third_order_dict = {key: value/value2 for key, value in num.items() for key2, value2 in denom.items() if key[0:3] == key2 }
    return {key: value for key, value in (third_order_dict.items()) }

def _ecdf(x):
    nobs = len(x)
    return np.arange(1,nobs+1)/float(nobs)

def fdrcorrection(pvals, thresh, is_sorted=False):
    ###FDR correction -- more info at: http://www.statsmodels.org/devel/_modules/statsmodels/stats/multitest.html#multipletests
    pvals = np.asarray(pvals)

    if not is_sorted:
        pvals_sortind = np.argsort(pvals)
        pvals_sorted = np.take(pvals, pvals_sortind)
    else:
        pvals_sorted = pvals  # alias

    ecdffactor = _ecdf(pvals_sorted)
    reject = pvals_sorted <= ecdffactor*thresh
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True

    pvals_corrected_raw = pvals_sorted / ecdffactor
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected>1] = 1
    if not is_sorted:
        pvals_corrected_ = np.empty_like(pvals_corrected)
        pvals_corrected_[pvals_sortind] = pvals_corrected
        del pvals_corrected
        reject_ = np.empty_like(reject)
        reject_[pvals_sortind] = reject
        return reject_, pvals_corrected_
    else:
        return reject, pvals_corrected


def calc_comb(substring, n):
    res = set()
    for combination in itertools.product(alphabet, repeat=n-len(substring)):
        for idx in range(len(combination)+1):
            res.add(''.join((*combination[:idx], substring, *combination[idx:])))
    return list(res)


def create_comb(sequence_length):
    sequence_list = []
    for sequence in product(alphabet, repeat=sequence_length):
        sequence_string = ''.join(sequence)
        sequence_list.append(sequence_string)
    return list(sequence_list)


######### static arguments ##############
codon_table = {
    'A': ('GCT', 'GCC', 'GCA', 'GCG'),
    'C': ('TGT', 'TGC'),
    'D': ('GAT', 'GAC'),
    'E': ('GAA', 'GAG'),
    'F': ('TTT', 'TTC'),
    'G': ('GGT', 'GGC', 'GGA', 'GGG'),
    'H': ('CAT', 'CAC'),
    'I': ('ATT', 'ATC', 'ATA'),
    'K': ('AAA', 'AAG'),
    'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
    'M': ('ATG',),
    'N': ('AAT', 'AAC'),
    'P': ('CCT', 'CCC', 'CCA', 'CCG'),
    'Q': ('CAA', 'CAG'),
    'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
    'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
    'T': ('ACT', 'ACC', 'ACA', 'ACG'),
    'V': ('GTT', 'GTC', 'GTA', 'GTG'),
    'W': ('TGG',),
    'Y': ('TAT', 'TAC'),
}

charge_table = { # '*' = neutral, '-' = negative charge, '+' = positive charge
    'A': ('*'),
    'C': ('*'),
    'D': ('-'),
    'E': ('-'),
    'F': ('*'),
    'G': ('*'),
    'H': ('-'),
    'I': ('*'),
    'K': ('+'),
    'L': ('*'),
    'M': ('*'),
    'N': ('*'),
    'P': ('*'),
    'Q': ('*'),
    'R': ('+'),
    'S': ('*'),
    'T': ('*'),
    'V': ('*'),
    'W': ('*'),
    'Y': ('*'),
}

polarity_table = { # '-' = non polar, '+' = polar
    'A': ('-'),
    'C': ('+'),
    'D': ('+'),
    'E': ('+'),
    'F': ('-'),
    'G': ('-'),
    'H': ('+'),
    'I': ('-'),
    'K': ('+'),
    'L': ('-'),
    'M': ('-'),
    'N': ('+'),
    'P': ('-'),
    'Q': ('+'),
    'R': ('+'),
    'S': ('+'),
    'T': ('+'),
    'V': ('-'),
    'W': ('-'),
    'Y': ('+'),
}

title = "***** Statistical evaluation of significance of Nullomers *****\n"
print(title)


########## calculations ######################
start_time = process_time()

if (level == 'PROT'):
    print("- The reference proteome is now processed")
elif (level == 'DNA'):
    print("- The reference genome is now processed")

full_seq_list = []
for rec in SeqIO.parse(fasta_file, "fasta"):
   x = ProteinAnalysis(str(rec.seq))
   full_seq_list.append(str(x.sequence))

print("- The calculation of transition matrices of probabilities is in progress")

aminoacid_percentage = ProteinAnalysis(''.join(full_seq_list)).get_amino_acids_percent()
aminoacid_counter = ProteinAnalysis(''.join(full_seq_list)).count_amino_acids()
total_sequence_length = sum(aminoacid_counter.values())
print("- The percentage of residues was calculated successfully")

transition_dictionary_first_order = {}
for p in list(first_order_probs(full_seq_list).items()):
    transition_dictionary_first_order[p[0][0] + p[0][1]] = p[1]
print("- The 1st order transition matrix was calculated successfully")
    
transition_dictionary_second_order = {}
for p in list(second_order_probs(full_seq_list).items()):
    transition_dictionary_second_order[p[0][0] + p[0][1] + p[0][2]] = p[1]    
print("- The 2nd order transition matrix was calculated successfully")

transition_dictionary_third_order = {}
for p in list(third_order_probs(full_seq_list).items()):
    transition_dictionary_third_order[p[0][0] + p[0][1] + p[0][2] + p[0][3]] = p[1]    
print("- The 3rd order transition matrix was calculated successfully")


##empty full_seq_list to free up memory
full_seq_list.clear()


print("- Absent words are currently processed")

line_length = {}
exclusions = {}
nullomer_list = []
exp_num_occur_zero_order = []
exp_num_occur_first_order = []
exp_num_occur_second_order = []
exp_num_occur_third_order = []
prob_corr_zero_occur_list_zero_order = []
prob_corr_zero_occur_list_first_order = []
prob_corr_zero_occur_list_second_order = []
prob_corr_zero_occur_list_third_order = []


with open(nullomers_file, encoding='utf8') as f:

    lines = f.read().splitlines()
    if (level == 'PROT'):
        lines = [ x for x in lines if x and (">" not in x) ] ##exclude blank lines and lines contain the '>' symbol
    elif (level == 'DNA'):
        lines = [ x for x in lines if x and (">" not in x) and ("N" not in x) ] ##exclude blank lines, lines that include 'N's and lines contain the '>' symbol
            
    print("- Absent words were parsed successfully\n")
    print("- The assesment of absent words has started")
    print("** Please be patient this might be a time-consuming process **")
   
    max_len = len(max(lines, key=len))
    min_len = len(min(lines, key=len))

    if (level == 'PROT'):
        print("- The shortest peptide in the nullomer's file is a " + str(min_len) + "-mer while the longest is a " + str(max_len) + "-mer")
    elif (level == 'DNA'):
        print("- The shortest nucleotide sequence in the nullomer's file is a " + str(min_len) + "-mer while the longest is a " + str(max_len) + "-mer")


    if (level == 'PROT' and max_len > 6):
        print("\n**Nullomers should be up to 6 amino acids in length. Please try again with shorter sequences**\n")
        raise SystemExit()
    elif (level == 'DNA' and correction_method == 'FDR' and max_len > 18):
        print("\n**Nullomers should be up to 18 nucleotides in length using the FDR method. Please try again with shorter sequences**\n")
    elif (level == 'DNA' and correction_method == 'BONF' and max_len > 18):
        print("\n**Nullomers should be up to 18 nucleotides in length using the BONF method. Please try again with shorter sequences**\n")
    elif (level == 'DNA' and correction_method == 'ADJ-BONF' and max_len > 14):
        print("\n**Nullomers should be up to 14 nucleotides in length using the ADJ-BONF method. Please try again with shorter sequences**\n")
        raise SystemExit()



    if (correction_method == 'FDR'):
        print("- The selected correction method is: " + str(correction_method) + "")
        print("- The p-value threshold is: " + str(threshold) + "\n")

        probability_zero_occurrence_zero_order = []
        probability_zero_occurrence_first_order = []
        probability_zero_occurrence_second_order = []
        probability_zero_occurrence_third_order = []
        pvalues_zero_order = np.ones(len(lines), dtype=int)

        for index in range(len(lines)):
            if (index % 1000000) == 0 and index != 0:
                print(str(index) + ' rows have been parsed')
            
            motif_length = len(lines[index])

            probability_one_occurrence_zero_order = 1
            probability_one_occurrence_first_order = 1
            probability_one_occurrence_second_order = 1
            probability_one_occurrence_third_order = 1
            
            for ind, current_letter in enumerate(str(lines[index])):
                
                if (ind == 0):
                    probability_one_occurrence_zero_order = probability_one_occurrence_zero_order * aminoacid_percentage[str(current_letter)]
                    probability_one_occurrence_first_order = probability_one_occurrence_first_order * aminoacid_percentage[str(current_letter)]
                    probability_one_occurrence_second_order = probability_one_occurrence_second_order * aminoacid_percentage[str(current_letter)]
                    probability_one_occurrence_third_order = probability_one_occurrence_third_order * aminoacid_percentage[str(current_letter)]
                    one_previous_letter = current_letter
                if (ind==1):
                    probability_one_occurrence_zero_order = probability_one_occurrence_zero_order * aminoacid_percentage[str(current_letter)]
                    probability_one_occurrence_first_order = probability_one_occurrence_first_order * transition_dictionary_first_order.get(str(one_previous_letter)+str(current_letter))
                    probability_one_occurrence_second_order = probability_one_occurrence_second_order * transition_dictionary_first_order.get(str(one_previous_letter)+str(current_letter))
                    probability_one_occurrence_third_order = probability_one_occurrence_third_order * transition_dictionary_first_order.get(str(one_previous_letter)+str(current_letter))
                    two_previous_letters = one_previous_letter
                    one_previous_letter = current_letter
                if (ind==2):
                    probability_one_occurrence_zero_order = probability_one_occurrence_zero_order * aminoacid_percentage[str(current_letter)]
                    probability_one_occurrence_first_order = probability_one_occurrence_first_order * transition_dictionary_first_order.get(str(one_previous_letter)+str(current_letter))
                    if transition_dictionary_second_order.get(str(two_previous_letters) + str(one_previous_letter) + str(current_letter)) is None:
                        probability_one_occurrence_second_order = probability_one_occurrence_second_order * 1
                        probability_one_occurrence_third_order = probability_one_occurrence_third_order * 1
                    else:
                        probability_one_occurrence_second_order = probability_one_occurrence_second_order * transition_dictionary_second_order.get(str(two_previous_letters)+str(one_previous_letter)+str(current_letter))
                        probability_one_occurrence_third_order = probability_one_occurrence_third_order * transition_dictionary_second_order.get(str(two_previous_letters)+str(one_previous_letter)+str(current_letter))
                    three_previous_letters = two_previous_letters
                    two_previous_letters = one_previous_letter
                    one_previous_letter = current_letter
                if (ind>=3):
                    probability_one_occurrence_zero_order = probability_one_occurrence_zero_order * aminoacid_percentage[str(current_letter)]
                    probability_one_occurrence_first_order = probability_one_occurrence_first_order * transition_dictionary_first_order.get(str(one_previous_letter)+str(current_letter))                    
                    if transition_dictionary_second_order.get(str(two_previous_letters) + str(one_previous_letter) + str(current_letter)) is None:
                        probability_one_occurrence_second_order = probability_one_occurrence_second_order * 1
                    else:
                        probability_one_occurrence_second_order = probability_one_occurrence_second_order * transition_dictionary_second_order.get(str(two_previous_letters)+str(one_previous_letter)+str(current_letter))
                    if transition_dictionary_third_order.get(str(three_previous_letters) + str(two_previous_letters) + str(one_previous_letter) + str(current_letter)) is None:
                        probability_one_occurrence_third_order = probability_one_occurrence_third_order * 1
                    else:
                        probability_one_occurrence_third_order = probability_one_occurrence_third_order * transition_dictionary_third_order.get(str(three_previous_letters) + str(two_previous_letters) + str(one_previous_letter) + str(current_letter)) 
                    three_previous_letters = two_previous_letters
                    two_previous_letters = one_previous_letter
                    one_previous_letter = current_letter

            expected_number_occurrences_zero_order = probability_one_occurrence_zero_order * (total_sequence_length - motif_length + 1)
            expected_number_occurrences_first_order = probability_one_occurrence_first_order * (total_sequence_length - motif_length + 1)
            expected_number_occurrences_second_order = probability_one_occurrence_second_order * (total_sequence_length - motif_length + 1)
            expected_number_occurrences_third_order = probability_one_occurrence_third_order * (total_sequence_length - motif_length + 1)

            probability_zero_occurrence_zero_order.append(math.exp(-expected_number_occurrences_zero_order))
            probability_zero_occurrence_first_order.append(math.exp(-expected_number_occurrences_first_order))
            probability_zero_occurrence_second_order.append(math.exp(-expected_number_occurrences_second_order))
            probability_zero_occurrence_third_order.append(math.exp(-expected_number_occurrences_third_order))
        
        prob_corr_zero_occur_list_zero_order = fdrcorrection(probability_zero_occurrence_zero_order,threshold)
        prob_corr_zero_occur_list_first_order = fdrcorrection(probability_zero_occurrence_first_order,threshold)
        prob_corr_zero_occur_list_second_order = fdrcorrection(probability_zero_occurrence_second_order,threshold)
        prob_corr_zero_occur_list_third_order = fdrcorrection(probability_zero_occurrence_third_order,threshold)

        idx = return_indices(prob_corr_zero_occur_list_zero_order[0], True)
        idx1 = return_indices(prob_corr_zero_occur_list_first_order[0], True)
        idx2 = return_indices(prob_corr_zero_occur_list_second_order[0], True)
        idx3 = return_indices(prob_corr_zero_occur_list_third_order[0], True)
        
        #for i in idx:
            #print(str(i))
        #for i1 in idx1:
            #print(str(i1))
        #for i2 in idx2:
            #print(str(i2))
        #for i3 in idx3:
            #print(str(i3))
        
        ids_in_common = np.intersect1d(idx, idx1) 
        #print(ids_in_common)
        ids_in_common = np.intersect1d(ids_in_common, idx2)
        #print(ids_in_common)
        ids_in_common = np.intersect1d(ids_in_common, idx3)
        #print(ids_in_common)

      

        if ids_in_common.size:
            print("\n** Results **\n")
            for index in ids_in_common:
                print(str(lines[index]) + ' - ' + str(max(prob_corr_zero_occur_list_zero_order[1][index], prob_corr_zero_occur_list_first_order[1][index], prob_corr_zero_occur_list_second_order[1][index], prob_corr_zero_occur_list_third_order[1][index])))
        else:        
            print("\nNo significant results found")
            
######################

    elif (correction_method == 'BONF'): ##bonferroni correction
        print("- The selected correction method is: " + str(correction_method) + "")
        print("- The p-value threshold is: " + str(threshold) + "\n")
        
        for index in range(len(lines)):
            if (index % 1000000) == 0 and index != 0:
                print(str(index) + ' rows have been parsed')             
            if not ">" in str(lines[index]) and len(str(lines[index]))!=0:
                motif_length = len(lines[index])

                probability_one_occurrence_zero_order = 1
                probability_one_occurrence_first_order = 1
                probability_one_occurrence_second_order = 1
                probability_one_occurrence_third_order = 1
            
                for ind, current_letter in enumerate(str(lines[index])):
                    if (ind == 0):
                        probability_one_occurrence_zero_order = probability_one_occurrence_zero_order * aminoacid_percentage[str(current_letter)]
                        probability_one_occurrence_first_order = probability_one_occurrence_first_order * aminoacid_percentage[str(current_letter)]
                        probability_one_occurrence_second_order = probability_one_occurrence_second_order * aminoacid_percentage[str(current_letter)]
                        probability_one_occurrence_third_order = probability_one_occurrence_third_order * aminoacid_percentage[str(current_letter)]
                        one_previous_letter = current_letter
                    if (ind==1):
                        probability_one_occurrence_zero_order = probability_one_occurrence_zero_order * aminoacid_percentage[str(current_letter)]
                        probability_one_occurrence_first_order = probability_one_occurrence_first_order * transition_dictionary_first_order.get(str(one_previous_letter)+str(current_letter))
                        probability_one_occurrence_second_order = probability_one_occurrence_second_order * transition_dictionary_first_order.get(str(one_previous_letter)+str(current_letter))
                        probability_one_occurrence_third_order = probability_one_occurrence_third_order * transition_dictionary_first_order.get(str(one_previous_letter)+str(current_letter))
                        two_previous_letters = one_previous_letter
                        one_previous_letter = current_letter
                    if (ind==2):
                        probability_one_occurrence_zero_order = probability_one_occurrence_zero_order * aminoacid_percentage[str(current_letter)]
                        probability_one_occurrence_first_order = probability_one_occurrence_first_order * transition_dictionary_first_order.get(str(one_previous_letter)+str(current_letter))
                        if transition_dictionary_second_order.get(str(two_previous_letters) + str(one_previous_letter) + str(current_letter)) is None:
                            probability_one_occurrence_second_order = probability_one_occurrence_second_order * 1
                            probability_one_occurrence_third_order = probability_one_occurrence_third_order * 1
                        else:
                            probability_one_occurrence_second_order = probability_one_occurrence_second_order * transition_dictionary_second_order.get(str(two_previous_letters)+str(one_previous_letter)+str(current_letter))
                            probability_one_occurrence_third_order = probability_one_occurrence_third_order * transition_dictionary_second_order.get(str(two_previous_letters)+str(one_previous_letter)+str(current_letter))
                        three_previous_letters = two_previous_letters
                        two_previous_letters = one_previous_letter
                        one_previous_letter = current_letter
                    if (ind>=3):
                        probability_one_occurrence_zero_order = probability_one_occurrence_zero_order * aminoacid_percentage[str(current_letter)]
                        probability_one_occurrence_first_order = probability_one_occurrence_first_order * transition_dictionary_first_order.get(str(one_previous_letter)+str(current_letter))                    
                        if transition_dictionary_second_order.get(str(two_previous_letters) + str(one_previous_letter) + str(current_letter)) is None:
                            probability_one_occurrence_second_order = probability_one_occurrence_second_order * 1
                        else:
                            probability_one_occurrence_second_order = probability_one_occurrence_second_order * transition_dictionary_second_order.get(str(two_previous_letters)+str(one_previous_letter)+str(current_letter))
                        if transition_dictionary_third_order.get(str(three_previous_letters) + str(two_previous_letters) + str(one_previous_letter) + str(current_letter)) is None:
                            probability_one_occurrence_third_order = probability_one_occurrence_third_order * 1
                        else:
                            probability_one_occurrence_third_order = probability_one_occurrence_third_order * transition_dictionary_third_order.get(str(three_previous_letters) + str(two_previous_letters) + str(one_previous_letter) + str(current_letter)) 
                        three_previous_letters = two_previous_letters
                        two_previous_letters = one_previous_letter
                        one_previous_letter = current_letter

                expected_number_occurrences_zero_order = probability_one_occurrence_zero_order * (total_sequence_length - motif_length + 1)
                expected_number_occurrences_first_order = probability_one_occurrence_first_order * (total_sequence_length - motif_length + 1)
                expected_number_occurrences_second_order = probability_one_occurrence_second_order * (total_sequence_length - motif_length + 1)
                expected_number_occurrences_third_order = probability_one_occurrence_third_order * (total_sequence_length - motif_length + 1)

                probability_zero_occurrence_zero_order = math.exp(-expected_number_occurrences_zero_order)
                probability_zero_occurrence_first_order = math.exp(-expected_number_occurrences_first_order)
                probability_zero_occurrence_second_order = math.exp(-expected_number_occurrences_second_order)
                probability_zero_occurrence_third_order = math.exp(-expected_number_occurrences_third_order)
                

                if (level == 'PROT'):
                    corrected_probability_zero_occurrence_zero_order = (probability_zero_occurrence_zero_order * pow(20, len(str(lines[index]))))
                    corrected_probability_zero_occurrence_first_order = (probability_zero_occurrence_first_order * pow(20, len(str(lines[index]))))
                    corrected_probability_zero_occurrence_second_order = (probability_zero_occurrence_second_order * pow(20, len(str(lines[index]))))
                    corrected_probability_zero_occurrence_third_order = (probability_zero_occurrence_third_order * pow(20, len(str(lines[index]))))
                elif (level == 'DNA'):
                    corrected_probability_zero_occurrence_zero_order = (probability_zero_occurrence_zero_order * pow(4, len(str(lines[index]))))
                    corrected_probability_zero_occurrence_first_order = (probability_zero_occurrence_first_order * pow(4, len(str(lines[index]))))
                    corrected_probability_zero_occurrence_second_order = (probability_zero_occurrence_second_order * pow(4, len(str(lines[index]))))
                    corrected_probability_zero_occurrence_third_order = (probability_zero_occurrence_third_order * pow(4, len(str(lines[index]))))
                
            
                if ((corrected_probability_zero_occurrence_zero_order < threshold) and (corrected_probability_zero_occurrence_first_order < threshold)  and (corrected_probability_zero_occurrence_second_order < threshold) and (corrected_probability_zero_occurrence_third_order < threshold)):
                    nullomer_list.append(str(lines[index]))
                    exp_num_occur_zero_order.append(expected_number_occurrences_zero_order)
                    exp_num_occur_first_order.append(expected_number_occurrences_first_order)
                    exp_num_occur_second_order.append(expected_number_occurrences_second_order)
                    exp_num_occur_third_order.append(expected_number_occurrences_third_order)
                
                    prob_corr_zero_occur_list_zero_order.append(corrected_probability_zero_occurrence_zero_order)
                    prob_corr_zero_occur_list_first_order.append(corrected_probability_zero_occurrence_first_order)
                    prob_corr_zero_occur_list_second_order.append(corrected_probability_zero_occurrence_second_order)
                    prob_corr_zero_occur_list_third_order.append(corrected_probability_zero_occurrence_third_order)
                
                
        if not (nullomer_list):    
            print("\nNo significant results found")
        else:
            print("\n** Results **\n")
            for itm1, itm2, itm3, itm4, itm5 in zip(nullomer_list, prob_corr_zero_occur_list_zero_order, prob_corr_zero_occur_list_first_order, prob_corr_zero_occur_list_second_order, prob_corr_zero_occur_list_third_order):
                max_prob = max(itm2, itm3, itm4, itm5)
                print(str(itm1) + ' - ' + str(max_prob))

#####################
        
    
    elif (correction_method == 'ADJ-BONF'): ##adjusted bonferroni correction
        print("- The selected correction method is: " + str(correction_method) + "")
        print("- The p-value threshold is: " + str(threshold) + "\n")
    
        remaining_kmers = {}
        aa_combinations = []
        for cur_len in range(min_len,max_len + 1):
            print(str(cur_len) + '-mers are now evaluated')
            aa_combinations = create_comb(cur_len)
            exclusions.update({cur_len:0})
            promising_kmers_dict = {}
            for indexx, cur_comb in enumerate(aa_combinations):
                if (indexx % 1000000) == 0 and indexx != 0:
                    print(str(indexx) + ' rows have been parsed')

                pre_probability_one_occurrence_zero_order = 1
                pre_probability_one_occurrence_first_order = 1
                pre_probability_one_occurrence_second_order = 1
                pre_probability_one_occurrence_third_order = 1

            
                for idx, cur_let in enumerate(cur_comb):
                    if (idx == 0):
                        pre_probability_one_occurrence_zero_order = pre_probability_one_occurrence_zero_order * aminoacid_percentage[str(cur_let)]
                        pre_probability_one_occurrence_first_order = pre_probability_one_occurrence_first_order * aminoacid_percentage[str(cur_let)]
                        pre_probability_one_occurrence_second_order = pre_probability_one_occurrence_second_order * aminoacid_percentage[str(cur_let)]
                        pre_probability_one_occurrence_third_order = pre_probability_one_occurrence_third_order * aminoacid_percentage[str(cur_let)]
                        one_prev_let = cur_let
                    if (idx == 1):
                        pre_probability_one_occurrence_zero_order = pre_probability_one_occurrence_zero_order * aminoacid_percentage[str(cur_let)]
                        pre_probability_one_occurrence_first_order = pre_probability_one_occurrence_first_order * transition_dictionary_first_order.get(str(one_prev_let)+str(cur_let))
                        pre_probability_one_occurrence_second_order = pre_probability_one_occurrence_second_order * transition_dictionary_first_order.get(str(one_prev_let)+str(cur_let))
                        pre_probability_one_occurrence_third_order = pre_probability_one_occurrence_third_order * transition_dictionary_first_order.get(str(one_prev_let)+str(cur_let))
                        two_prev_let = one_prev_let
                        one_prev_let = cur_let
                    if (idx == 2):
                        pre_probability_one_occurrence_zero_order = pre_probability_one_occurrence_zero_order * aminoacid_percentage[str(cur_let)]
                        pre_probability_one_occurrence_first_order = pre_probability_one_occurrence_first_order * transition_dictionary_first_order.get(str(one_prev_let)+str(cur_let))
                        if transition_dictionary_second_order.get(str(two_prev_let) + str(one_prev_let) + str(cur_let)) is None:
                            pre_probability_one_occurrence_second_order = pre_probability_one_occurrence_second_order * 1
                            pre_probability_one_occurrence_third_order = pre_probability_one_occurrence_third_order * 1
                        else:
                            pre_probability_one_occurrence_second_order = pre_probability_one_occurrence_second_order * transition_dictionary_second_order.get(str(two_prev_let)+str(one_prev_let)+str(cur_let))
                            pre_probability_one_occurrence_third_order = pre_probability_one_occurrence_third_order * transition_dictionary_second_order.get(str(two_prev_let)+str(one_prev_let)+str(cur_let))
                        three_prev_let = two_prev_let
                        two_prev_let = one_prev_let
                        one_prev_let = cur_let
                    if (idx >= 3):
                        pre_probability_one_occurrence_zero_order = pre_probability_one_occurrence_zero_order * aminoacid_percentage[str(cur_let)]
                        pre_probability_one_occurrence_first_order = pre_probability_one_occurrence_first_order * transition_dictionary_first_order.get(str(one_prev_let)+str(cur_let))                    
                        if transition_dictionary_second_order.get(str(two_prev_let) + str(one_prev_let) + str(cur_let)) is None:
                            pre_probability_one_occurrence_second_order = pre_probability_one_occurrence_second_order * 1
                        else:
                            pre_probability_one_occurrence_second_order = pre_probability_one_occurrence_second_order * transition_dictionary_second_order.get(str(two_prev_let)+str(one_prev_let)+str(cur_let))
                        if transition_dictionary_third_order.get(str(three_prev_let) + str(two_prev_let) + str(one_prev_let) + str(cur_let)) is None:
                            pre_probability_one_occurrence_third_order = pre_probability_one_occurrence_third_order * 1
                        else:
                            pre_probability_one_occurrence_third_order = pre_probability_one_occurrence_third_order * transition_dictionary_third_order.get(str(three_prev_let) + str(two_prev_let) + str(one_prev_let) + str(cur_let)) 
                        three_prev_let = two_prev_let
                        two_prev_let = one_prev_let
                        one_prev_let = cur_let

                exp_numb_occurrences_zero_order = pre_probability_one_occurrence_zero_order * (total_sequence_length - len(cur_comb) + 1)
                exp_numb_occurrences_first_order = pre_probability_one_occurrence_first_order * (total_sequence_length - len(cur_comb) + 1)
                exp_numb_occurrences_second_order = pre_probability_one_occurrence_second_order * (total_sequence_length - len(cur_comb) + 1)
                exp_numb_occurrences_third_order = pre_probability_one_occurrence_third_order * (total_sequence_length - len(cur_comb) + 1)

                prob_zero_occurrence_zero_order = math.exp(-exp_numb_occurrences_zero_order)
                prob_zero_occurrence_first_order = math.exp(-exp_numb_occurrences_first_order)
                prob_zero_occurrence_second_order = math.exp(-exp_numb_occurrences_second_order)
                prob_zero_occurrence_third_order = math.exp(-exp_numb_occurrences_third_order)


                if (len(cur_comb) == 2):
                    max_prob = max(prob_zero_occurrence_zero_order)
                elif (len(cur_comb) == 3):
                    max_prob = max(prob_zero_occurrence_zero_order, prob_zero_occurrence_first_order)
                elif (len(cur_comb) == 4):
                    max_prob = max(prob_zero_occurrence_zero_order, prob_zero_occurrence_first_order, prob_zero_occurrence_second_order)
                elif (len(cur_comb) >= 5):
                    max_prob = max(prob_zero_occurrence_zero_order, prob_zero_occurrence_first_order, prob_zero_occurrence_second_order, prob_zero_occurrence_third_order)


                promising_kmers_dict[cur_comb] = max_prob


            for key_nullo, value_prob in sorted(promising_kmers_dict.items(), key = itemgetter(1), reverse = True):
            
                if (level == 'PROT'):
                    corr_prob = value_prob * ( pow(20, cur_len) - exclusions.get(cur_len) )
                elif (level == 'DNA'):
                    corr_prob = value_prob * ( pow(4, cur_len) - exclusions.get(cur_len) )
            
                if corr_prob > threshold:
                    exclusions.update({cur_len:(exclusions.get(cur_len) + 1)})
                else:   
                    remaining_kmers[key_nullo] = corr_prob


        keys = remaining_kmers.keys() & lines
        results = {k:remaining_kmers[k] for k in keys}
        print("\n** Results **\n")
        if len(results.keys()) == 0:
            print('\nNo significant results found')
        else:
            for key, val in results.items():
                print(str(key) + " - " + str(val))


end_time = process_time() - start_time
print('\nTotal execution time: ' + str(end_time) + " seconds")
     
print("\n** The script terminated successfully **")
