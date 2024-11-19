#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(data.table))
suppressPackageStartupMessages(require(rconfig))
suppressPackageStartupMessages(require(do))

#functions to translate nucleotide muations to AA
suppressPackageStartupMessages(source("mutationTranslation.R"))
#functions to read in vcf files
suppressPackageStartupMessages(source("readingFunctions.R"))


option_list = list(
  make_option(c("-d", "--vpipe_dir"), action="store", default=NA, type='character',
              help="Path to v-pipe working directory"),
  make_option(c("-l", "--locationFile"), action="store", default=NA, type='character',
              help="Path to location translation file")
)

opt = parse_args(OptionParser(option_list=option_list))


location_translation <- fread(opt$locationFile)
dir_euler <- opt$vpipe_dir


configs <- read_ini(paste0(dir_euler,'vpipe.config'))
fasta_file <- paste0(dir_euler,configs$input$reference)
sample_file <- paste0(dir_euler,configs$input$samples_file)
results <-  paste0(dir_euler,configs$output$datadir) #output dir of v-pipe run


segment <- gsub(".fasta","",file.name(fasta_file))

##### 1. Reading of vcf files #####

samples <- read.table(sample_file)

samples_dir <- paste(samples$V1,samples$V2,sep = "/")

var_files <- paste(results,samples_dir,"variants/SNVs/snvs.vcf",sep = "/")

var_list <- var_files %>% 
  map(function(x) read_and_mark_vcf(x))

#Usually not all samples will have muation calls aka vcf files 
var_list <- var_list[!sapply(var_list, is.null)]


###### 2. Translating Nucleotide mutations to AA mutations ####

#For translating nucleotide mutations to AA mutations, the ORF start and end
#points as well as their respective fasta files are taken from nexstrain
#https://github.com/nextstrain/seasonal-flu/tree/master/config

if (segment == "H1"){
  
  #only HA1/2 because have continous indexing
  start = 72
  end = 1718
  
}else if (segment == "H3"){
  
  start = 66
  end = 1418
  
}else if (segment == "N1"){
  
  start = 9	
  end = 1418
  
}else if (segment == "N2"){
  
  start = 4	
  end = 1413
}

MutAA <- var_list %>% 
  map(function(x) identify_AA(fasta_file,start,end,x)) %>% 
  bind_rows() %>% 
  setDT()

#don't care about synonymous mutations
dt <- MutAA[new_AA != org_AA]

dt[, mut_lable := paste0(org_AA,gene_codonIDX,new_AA)]
dt[, mut_lable_nuc := paste0(REF,POS,ALT)]

dt = merge(dt,location_translation,
           by.x=c('location_code'),
           by.y=c('code'))


#formatting for Dashboard to jason format
#coverage below 100 reads is not trustworthy
dt[DP >= 100, AF := AF]
dt[DP < 100, AF := NA_real_]
colnames(dt)[colnames(dt) == "sample_name"] <- "submissionId"

#placeholder lineage information
dt[, lineage_name := "something"]
dt[, lineage_abundance := NA_real_]

#formatting for Dashboard to jason format
dt_dash_aa = dt[,c("submissionId","date","location","mut_lable","AF","lineage_name","lineage_abundance")]
dt_dash_aa[, aminoAcidMutationFrequencies := sapply(1:.N, function(i) setNames(list(AF[i]), mut_lable[i]))]
dt_dash_aa[, lineageFrequencyEstimates := sapply(1:.N, function(i) setNames(list(lineage_abundance[i]), lineage_name[i]))]

dt_dash_aa = dt_dash_aa[, .(aminoAcidMutationFrequencies = paste0("{", paste(sprintf('"%s": %.6f', mut_lable, AF), 
                                                                             collapse = ", "), "}"),
                            lineageFrequencyEstimates = paste0("{", paste(sprintf('"%s": %.6f',lineage_name,lineage_abundance), 
                                                                          collapse = ", "), "}")),
                        by = .(submissionId, date, location)]


dt_dash_nuc = dt[,c("submissionId","date","location","mut_lable_nuc","AF","lineage_name","lineage_abundance")]
dt_dash_nuc[, nucleotideMutationFrequencies := sapply(1:.N, function(i) setNames(list(AF[i]), mut_lable_nuc[i]))]
dt_dash_nuc[, lineageFrequencyEstimates := sapply(1:.N, function(i) setNames(list(lineage_abundance[i]), lineage_name[i]))]

dt_dash_nuc = dt_dash_nuc[, .(nucleotideMutationFrequencies = paste0("{", paste(sprintf('"%s": %.6f', mut_lable_nuc, AF), 
                                                                                collapse = ", "), "}"),
                              lineageFrequencyEstimates = paste0("{", paste(sprintf('"%s": %.6f', lineage_name,lineage_abundance), 
                                                                            collapse = ", "), "}")),
                          by = .(submissionId, date, location)]


dt_out = merge(dt_dash_aa,dt_dash_nuc) 
dt_out[, reference := segment]
dt_out[, primerProtocol := "EAWAG_11_24"]

#need to do it manually here, because before nummeric value is expected
dt_out[, `:=`(
  lineageFrequencyEstimates = str_replace_all(lineageFrequencyEstimates, "NA","null"),
  aminoAcidMutationFrequencies = str_replace_all(aminoAcidMutationFrequencies, "NA", "null"),
  nucleotideMutationFrequencies = str_replace_all(nucleotideMutationFrequencies, "NA", "null")
)]

###Outputting###


dir_out <-  paste0(dir_euler,"MutationFrequencies/")

if (!dir.exists(dir_out)){
  dir.create(dir_out)
}

fwrite(dt,paste0(dir_out,"Mutations.tsv"))
fwrite(dt_out,paste0(dir_out,"Mutations_Dashboard.tsv"),
       sep = "\t", quote = FALSE, na = "null")

