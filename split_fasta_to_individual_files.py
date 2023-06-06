#!/usr/bin/env python
# a script to extract fasta records from a fasta file to multiple separate fasta files based on a particular ID (time point) in a particular field, for a given delimiter
# to run, navigate to file location with command prompt and enter: python split_fasta_by_collections.py infile.fasta
from Bio import SeqIO
import os
import sys

print('started')
records = SeqIO.parse(sys.argv[1], "fasta")
collected = {}
for record in records:
    descr = record.description.split("_")[0].strip() # "_" sets the delimeter, "1" sets the field where counting starts at 0 for the first field
    try:
        collected[descr].append(record)
    except KeyError:
        collected[descr] = [record ,]

file_name = "%s.fasta" 
file_path = os.getcwd() #sets the output file path to your current working directory

for (descr, records) in collected.items():  
    with open(os.path.join(file_path,  file_name % descr), 'w') as f:
        SeqIO.write(records, f, 'fasta')
