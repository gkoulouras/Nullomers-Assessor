# Nullomers Assessor
<b>Nullomers Assessor:</b> A computational method for the statistical evaluation of nullomers (a.k.a. minimal absent words)

<b>Description:</b>

Nullomers Assessor is a probabilistic methodology for the evaluation and classification of biological nullomers as 'significant absent' or 'expected absent' based on the frequency of adjacent residues (nucleotides/amino acids) in the entire genome/proteome of a species. A 'significant absent' oligomer is a motif statistically expected to exist, but in reality it is entirely absent. The underlying method estimates the frequency of residues and subsequently calculates 3 transition probabilities, one fr each of the first three Markovian orders. The method is designed to assess either nucleotide or amino acid oligomers. Three different statistical correction methods have been implemented and provided built-in with the current version of the script.

<b>Preparatory stage:</b>

The <b>nullomers_assessor.py</b> script (which can be downloaded and run directly) requires 2 input files. The first one should be a fasta file containing the entire genome/proteome of a species. The second file should be a list of nullomers. Tools as such <b>[MAW](https://github.com/solonas13/maw)</b> or <b>[em-MAW](https://github.com/solonas13/maw/tree/master/em-maw)</b> can be used for the calculation of nullomers. The two above-mentioned tools though, require a two-line fasta file in order to calculate globally missing sequences from a given genome/proteome. The specific format can be easily achieved using the <b>fasta_formatter.py</b> script which transforms a typical fasta file to a two-line fasta file. Sample files are provided in a separate directory. Find out more at <b>[https://www.nullomers.org/](https://www.nullomers.org/)</b>.

<b>Usage:</b>

Simply download and execute the <b>nullomers_assessor.py</b> script by giving the following 6 arguments (separated by a blank space).

```
--absolute-path-of-fasta-file       <string>   <mandatory>  A typical fasta file (either DNA or protein sequences)
--absolute-path-of-nullomers-file   <string>   <mandatory>  A list of nullomers (without header)
--threshold                         <double>   <mandatory>  A float number between [0-1] indicating the threshold of 
                                                            statistical correction. Nullomers with q-values greater
                                                            than the user-specified cut-off will be discarded
--type-of-sequences                 <string>   <mandatory>  'DNA' for nucleotide sequences, 
                                                            'PROT' for protein sequences
--statistical-correction            <string>   <mandatory>  'BONF' for Bonferroni correction, 
                                                            'TARONE' for Tarone correction,
                                                            'FDR' for Benjamini-Hochberg procedure
--print-log                         <boolean>  <mandatory>  'TRUE' or 'FALSE'                                                         
```

<b>Example of usage:</b>

Unix OS:
> python3 nullomers_assessor.py input/P.troglodytes_genome_oneline.fasta output/pantro_14mers.out 0.01 DNA BONF TRUE

Windows OS:
> py nullomers_assessor.py C:\input\P.troglodytes_genome_oneline.fasta C:\output\pantro_14mers.out 0.01 DNA BONF TRUE

<b>Results:</b>

Nullomers Assessor has been applied to multiple genomes and proteomes of various organisms including <i>Homo Sapiens</i>. The results of the study are freely accessible through <b>[Nullomers Database](https://www.nullomers.org)</b>, a continuously enriched web-accessible repository of significant Minimal Absent Words.

<b>Publication:</b>

<i>Manuscript submitted for publication. We will update this section as soon as the paper is published.</i>

<b>Contact:</b>

For questions, suggestions, bug-reports or feedback, please get in touch by email:
<ul><li>gkoulouras {at symbol} gmail {dot symbol} com</li></ul>

<b>License:</b>

This project is licensed under the Apache 2.0 license, quoted below.

Copyright (c) 2020 Grigorios Koulouras

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
