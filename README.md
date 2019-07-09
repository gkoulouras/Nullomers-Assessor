# Nullomers Assessor
<b>Nullomers Assessor:</b> a computational method for the statistical evaluation of Nullomers (a.k.a. Minimal Absent Words)

<b>Description:</b>

Nullomers Assessor is a computational method for the evaluation and classification of a list of Nullomers as 'significantly absent' or not, depending on the distribution of residues in the specified genome/proteome. The underlying method is based on 4 different transition probabilities (orders) of the Markov chain stohastic model. The method is able to assess either nucleotide sequences or amino acid oligomers. Three different statictical correction method have been implemented and are provided with this version.

<b>Usage:</b>

Simply download and execute the <b>nullomers_assessor.py</b> script by giving the following 5 arguments.

```
--absolute-path-of-fasta-file       <str>   <mandatory>     The file should be one line fasta file
--absolute-path-of-nullomers-file   <str>   <mandatory>     The file should include a list of nullomers
--threshold                         <dbl>   <mandatory>     A float value which indicates the p-value of 
                                                            statistical correction
--level                             <str>   <mandatory>     'DNA' for nucleotide sequences, 
                                                            'PROT' for protein sequences
--statistical-correction            <str>   <mandatory>     'BONF' for a standard Bonferroni correction, 
                                                            'ADJ-BONF' for an adjusted Bonferroni correction,
                                                            'FDR' for a False Discovery Rate correction method 
```

<b>Example:</b>

$ python3 nullomers_assessor.py input/P.troglodytes_genome_oneline.fasta output/pantro_14mers.out 0.01 DNA BONF

<b>Results:</b>

Preliminary analysis has been conducted in multiple genomes and proteomes of various organisms. The results of this study are provided via <b>Nullomers Database</b> a constantly enriched web-accessed database at http://83.212.98.39/nullomersdb/default

<b>Publication:</b>

Manuscript submitted for publication. As soon as that article is accepted for publication we will update this field.
