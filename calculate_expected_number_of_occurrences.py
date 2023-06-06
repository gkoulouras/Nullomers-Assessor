import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import defaultdict

if __name__ == "__main__":
    fasta_file = sys.argv[1]
    nullomers_file = sys.argv[2]
else:
    print("\n**An error occured**\n")
    raise SystemExit() 


def nth_order_probs(sequences, n):
    num = defaultdict(int)
    for seq in sequences:
        for i in range(len(seq) - n):
            j = i + n + 1
            key = seq[i:j]
            num[key] += 1
    denom = defaultdict(int)
    for key, value in (num.items()):
        denom[key[0:n]] += value
    return {key: value/denom[key[0:n]] for key, value in num.items()}


full_seq_list = []
for rec in SeqIO.parse(fasta_file, "fasta"):
   x = ProteinAnalysis(str(rec.seq))
   full_seq_list.append(str(x.sequence))


aminoacid_percentage = nth_order_probs(full_seq_list, 0)
aminoacid_counter = ProteinAnalysis(''.join(full_seq_list)).count_amino_acids()
total_sequence_length = sum(aminoacid_counter.values())
print("- The frequency of residues has been calculated successfully")

transition_dictionary_first_order = nth_order_probs(full_seq_list, 1)
print("- The 1st-order transition probabilities have been calculated successfully")
    
transition_dictionary_second_order = nth_order_probs(full_seq_list, 2)
print("- The 2nd-order transition probabilities have been calculated successfully")

transition_dictionary_third_order = nth_order_probs(full_seq_list, 3)
print("- The 3rd-order transition probabilities have been calculated successfully")

##empty full_seq_list to free up memory
full_seq_list.clear()


with open(nullomers_file, encoding='utf8') as f:

    lines = f.read().splitlines()
    lines = [ x for x in lines if x and ("N" not in x) ] ##exclude blank lines and lines that include 'N's

    for index in range(len(lines)):
            
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

            expected_number_occurrences_zero_order = probability_one_occurrence_zero_order * (total_sequence_length * 2  - motif_length + 1)
            expected_number_occurrences_first_order = probability_one_occurrence_first_order * (total_sequence_length * 2 - motif_length + 1)
            expected_number_occurrences_second_order = probability_one_occurrence_second_order * (total_sequence_length * 2 - motif_length + 1)
            expected_number_occurrences_third_order = probability_one_occurrence_third_order * (total_sequence_length * 2 - motif_length + 1)

            print(str(lines[index]),'\t',str(int(round(expected_number_occurrences_zero_order))),'\t',str(int(round(expected_number_occurrences_first_order))),'\t',str(int(round(expected_number_occurrences_second_order))),'\t',str(int(round(expected_number_occurrences_third_order))))
