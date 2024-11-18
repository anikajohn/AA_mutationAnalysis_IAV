#!/bin/bash


#conda activate influenza_analysis_R

vpipe_dir="/cluster/work/bewi/members/anjohn/projects/wastewater/Influenza/AA_mutationAnalysis_IAV/test_data/work-IA_H3/v-pipe/"
location_dic="ww_locations.tsv"

./detect_AAMutations.R -d $vpipe_dir -l $location_dic
