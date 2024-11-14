# AA_mutationAnalysis_IAV

This directory contains (commandline) scripts to combine, clean and edit vcf files, outputted by the v-pipe (https://cbg-ethz.github.io/V-pipe/) lofreq rule.

### detect_AAMutations.R:
Comandline accesible R script. Takes as input the v-pipe directory of interests and the ww_locations.tsv (https://github.com/cbg-ethz/cowwid/blob/master/ww_locations.tsv). The ouput of the script is a directory called MutationFrequencies, which contains the main output file AA_mutations.tsv.

### mut_detect.sh:
Example of how to run detect_AAMutations.R on the commandline.

### IAV_analysis.yml:
Environment yaml file for conda environment set-up.

### test_data:
Small test data set including the output.
