# Nullomers Assessor
<b>Nullomers Assessor:</b> A computational method for the statistical evaluation of Nullomers (a.k.a. Minimal Absent Words)

<b>Description:</b>

Nullomers Assessor is a probabilistic method for the statistical evaluation and classification of biological nullomers into 'significantly absent' or 'insignificant', based on the frequency and distribution of residues in the entire genome/proteome of a species. A 'significant absent' oligomer is a sequence statistically foreseeable to exist, but, in reality, it is entirely absent. The underlying method estimates the frequency of residues and subsequently generates 3 different transition matrices of probabilities, based on the first, second and third order of the Markovian chains stohastic model. The method is able to assess either nucleotide or amino acid oligomers. Three different statictical correction methods have been implemented and provided built-in with the current version of <b>NullomersAssessor.py</b> script.

<b>Preparatory steps:</b>

The <b>NullomersAssessor.py</b> script (which can directly be downloaded and run) requires 2 input files. The first one should be a fasta file containing the entire genome/proteome of a species. The second file should be a list of nullomers. Tools such as <b>[MAW](https://github.com/solonas13/maw)</b> or <b>[em-MAW](https://github.com/solonas13/maw/tree/master/em-maw)</b> can be used for the calculation of globally missing sequences in an organism and export lists of nullomers. The two above-mentioned tools though, require a specific one-header fasta format in order to calculate globally missing sequences from a given genome/proteome. The specific format can easily be achieved using the <b>FastaFormatter.py</b> script which transforms a typical fasta file to a file of one-line sequence with a single header. Find out more about the underlying algorithm in the detailed [documentation](https://www.nullomers.org/Documentation_NullomersAssessor) of <b>Nullomers Database</b>.

<b>Usage:</b>

Simply download and execute the <b>NullomersAssessor.py</b> script by giving the following 5 arguments (separated by a blank space).

```
--absolute-path-of-fasta-file       <str>   <mandatory>     A typical fasta file (either DNA or protein sequences)
--absolute-path-of-nullomers-file   <str>   <mandatory>     A list of nullomers (without header)
--threshold                         <dbl>   <mandatory>     A float number between [0-1] indicating the threshold of 
                                                            statistical correction. Nullomers with corrected 
                                                            p-values greater than the specified cut-off are discarded
--sequences                         <str>   <mandatory>     'DNA' for nucleotide sequences, 
                                                            'PROT' for protein sequences
--statistical-correction            <str>   <mandatory>     'BONF' for standard Bonferroni correction, 
                                                            'ADJ-BONF' for adjusted Bonferroni correction,
                                                            'FDR' for standard False Discovery Rate correction method 
```

<b>Example of usage:</b>

Unix OS:
> python3 NullomersAssesor.py input/P.troglodytes_genome_oneline.fasta output/pantro_14mers.out 0.01 DNA BONF

Windows OS:
> py NullomersAssesor.py C:\input\P.troglodytes_genome_oneline.fasta C:\output\pantro_14mers.out 0.01 DNA BONF

<b>Results:</b>

Preliminary analysis has been conducted in multiple genomes and proteomes of various organisms including <i>Human</i>. The results of this study are freely available through <b>[Nullomers Database](https://www.nullomers.org)</b>, a continuously enriched web-accessible repository of significant Minimal Absent Words.

<b>Publication:</b>

<i>Manuscript submitted for publication. This field will be updated as soon as the article is accepted for publication.</i>

<b>Contact:</b>

For questions, suggestions, bug-reports or feedback, please get in touch by email:
<ul><li>gkoulouras {at symbol} gmail {dot symbol} com</li></ul>

<b>License:</b>

This project is licensed under the Apache 2.0 license, quoted below.

Copyright (c) 2019 Grigorios Koulouras

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of the License at [http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.
