# Nullomers Assessor
<b>Nullomers Assessor:</b> A computational method for the statistical evaluation of Nullomers (a.k.a. Minimal Absent Words)

<b>Description:</b>

Nullomers Assessor is a computational method for the evaluation and classification of a list of Nullomers as 'significantly absent' or not, depending on the distribution of residues in a specified genome/proteome. The underlying method is based on 4 different transition probabilities (orders) of the Markov chain stohastic model. The method is able to assess either nucleotide sequences or amino acid oligomers. Three different statictical correction methods have been implemented and are provided with the current version of the <b>NullomersAssessor.py</b> script.

<b>Preparatory steps:</b>

The <b>NullomersAssessor.py</b> script (which can be directly downloaded from this repository) requires 2 input files. The first one should be a fasta file containing the entire genome/proteome of a species. The second file should be a list of nullomers. Tools such as <b>[MAW](https://github.com/solonas13/maw)</b> or <b>[em-MAW](https://github.com/solonas13/maw/tree/master/em-maw)</b> can be used for the calculation of globally missing sequences in an organism and output a list of nullomers. The two previously mentioned tools though, require a specific one-line with one-header format. This format can be easily achieved using the <b>OneHeaderFastaFormatter.py</b> script which transforms a typical fasta file to an one-line with one header file. Please find out more information on how to run the preparatory steps in the detailed [user guide](https://www.nullomers.org/Documentation_NullomersAssessor) of <b>Nullomers Database</b>.

<b>Usage:</b>

Simply download and execute the <b>NullomersAssessor.py</b> script by giving the following 5 arguments (separated by a blank space).

```
--absolute-path-of-fasta-file       <str>   <mandatory>     The file should be an one-header/one-line fasta file
--absolute-path-of-nullomers-file   <str>   <mandatory>     The file should include a list of nullomers
--threshold                         <dbl>   <mandatory>     A float number which indicates the threshold of 
                                                            statistical correction. Nullomers with q-values 
                                                            greater than the specified value are discarded
--sequences                         <str>   <mandatory>     'DNA' for nucleotide sequences, 
                                                            'PROT' for protein sequences
--statistical-correction            <str>   <mandatory>     'BONF' for standard Bonferroni correction, 
                                                            'ADJ-BONF' for adjusted Bonferroni correction,
                                                            'FDR' for False Discovery Rate correction method 
```

<b>Example of usage:</b>

Unix OS:
> python3 NullomersAssesor.py input/P.troglodytes_genome_oneline.fasta output/pantro_14mers.out 0.01 DNA BONF

Windows OS:
> py NullomersAssesor.py C:\input\P.troglodytes_genome_oneline.fasta C:\output\pantro_14mers.out 0.01 DNA BONF

<b>Results:</b>

Preliminary analysis has been conducted in multiple genomes of various organisms and the <i>Human</i> proteome. The results of this study are freely provided through <b>[Nullomers Database](https://www.nullomers.org)</b>, a continuously enriched web-accessed resource of significant Minimal Absent Words.

<b>Publication:</b>

Manuscript submitted for publication. This field will be updated as soon as the article is accepted for publication.

<b>Contact:</b>

For questions, suggestions, bug-reports or feedback, please get in touch by email:
<ul><li>gkoulouras {at} gmail {dot} com</li></ul>

<b>License:</b>

This project is licensed under the Apache 2.0 license, quoted below.

Copyright (c) 2019 Grigorios Koulouras

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
