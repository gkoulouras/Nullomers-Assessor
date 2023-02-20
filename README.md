# Nullomers Assessor
<b>Nullomers Assessor:</b> A computational method for statistical evaluation of nullomers and Minimal Absent Words (MAWs)

<b>Description:</b>

Nullomers Assessor is a probabilistic methodology for the evaluation of biological MAWs based on Markovian models with multiple-testing correction. A 'significant' MAW is an absent motif which is statistically expected to occur. Significant absent motifs are considered to be under negative selection. The algorithm estimates the frequency of residues and subsequently calculates 3 transition probabilities, one for each of the first three Markovian orders. The method can analyse either nucleotide or amino acid sequences. Three different statistical correction methods have been implemented and provided built-in with the current version of the script.

<b>Preparatory steps:</b>

The <b>nullomers_assessor.py</b> script requires 2 input files. The first one should be a fasta file containing the entire genome/proteome of a species. The second file should be a list of MAWs. Tools such as <b>[MAW](https://github.com/solonas13/maw)</b> or <b>[em-MAW](https://github.com/solonas13/maw/tree/master/em-maw)</b> can be used for the identification of MAWs. The two above-mentioned tools though, require a fasta file with one header and one sequence in order to calculate globally missing motifs from a given genome/proteome. The specific format can be easily achieved by using the <b>fasta_formatter.py</b> script which transforms a typical fasta file to a two-line fasta file. Sample files are provided in a separate directory. For more information visit: <b>[https://www.nullomers.org/](https://www.nullomers.org/)</b>.

<b>Usage:</b>

Simply download and execute the <b>nullomers_assessor.py</b> script by giving the following 6 arguments (separated by a blank space).

```
--absolute-path-of-fasta-file       <string>   <mandatory>  A typical fasta file (either DNA or protein sequences)
--absolute-path-of-MAWs-file        <string>   <mandatory>  A list of MAWs (without header)
--threshold                         <double>   <mandatory>  A float number between [0-1] indicating the threshold of 
                                                            statistical correction. MAWs with corrected p-values
                                                            greater than the user-specified cut-off will be discarded
--type-of-sequences                 <string>   <mandatory>  'DNA' for nucleotide sequences, 
                                                            'PROT' for protein sequences
--statistical-correction            <string>   <mandatory>  'BONF' for Bonferroni correction, 
                                                            'TARONE' for Tarone correction,
                                                            'FDR' for Benjamini-Hochberg procedure
--print-log                         <boolean>  <mandatory>  'TRUE' or 'FALSE'                                                         
```

<b>Example of usage:</b>

Unix OS:
> python3 nullomers_assessor.py input/P.troglodytes_genome.fasta output/pantro_nullomers.out 0.01 DNA BONF TRUE

Windows OS:
> py nullomers_assessor.py C:\input\P.troglodytes_genome.fasta C:\output\pantro_nullomers.out 0.01 DNA BONF TRUE

<b>Results:</b>

Nullomers Assessor has been applied to numerous genomes and proteomes of various organisms including <i>Homo Sapiens</i>. The results of the study are provided via <b>[Nullomers Database](https://www.nullomers.org)</b>, a web-based repository which hosts significant MAWs from hundreds of species and thoudands of individual virus sequences.

<b>Publication:</b>

For more information about the method and the statistical correction procedures, please consult the publication:<br>
Koulouras G, Frith MC. <b><i>Significant non-existence of sequences in genomes and proteomes.</i></b> <i>Nucleic Acids Res.</i> 2021;49(6):3139–3155. [[link]](https://academic.oup.com/nar/article/49/6/3139/6166853)

<b>Contact:</b>

For questions, suggestions, bug-reports, or feedback, please email me:
<ul><li>gkoulouras {at symbol} gmail {dot symbol} com</li></ul>

<b>License:</b>

This project is licensed under the Apache 2.0 license, quoted below.

Copyright (c) 2020 Grigorios Koulouras

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
