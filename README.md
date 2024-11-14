# AA_mutationAnalysis_IAV

This directory contains (commandline) scripts to combine, clean and edit vcf files, outputted by the v-pipe (https://cbg-ethz.github.io/V-pipe/) lofreq rule.

#### detect_AAMutations.R:
Comandline accesible R script. Takes as input the v-pipe directory of interests and the ww_locations.tsv (https://github.com/cbg-ethz/cowwid/blob/master/ww_locations.tsv). The ouput of the script is a directory called MutationFrequencies, which contains the main output file *AA_mutations.tsv*.

#### mut_detect.sh:
Example of how to run detect_AAMutations.R on the commandline.

#### IAV_analysis.yml:
Environment yaml file for conda environment set-up.

#### test_data:
Small test data set including the output.


### *AA_mutations.tsv* format:

The main outputfile is basically a concatenation of all lofreq generated vcf files (vcf lists nucleotide changes), but includes additional information, for example how the nucleotide change effected the respective amino acid. Synonymous aa mutations are excluded, since they don't have a biological effect.

The file contains the following columns:
location_code, POS, #CHROM, ID, REF, ALT, QUAL, FILTER, DP, AF, SB, sample_name, date, gene_codonIDX, codonIDX, org_codon, new_codon, new_AA, org_AA, mut_lable,location

org_codon:     original codon, before nucleotide change
new_codon:     new codon after nucleotide change (mutated)
org_AA:        translated original codon in AA
new_AA         translated new codon in AA (mutated)  

gene_codonIDX: position of the codon (on AA sequence) where the nucleotide change is detected.
mut_lable:     lable of AA change after community standard (original_AA position_AA new_AA) 



