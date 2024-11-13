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

dt = merge(dt,location_translation,
           by.x=c('location_code'),
           by.y=c('code'))

###Outputting###


dir_out <-  paste0(dir_euler,"MutationFrequencies/")

if (!dir.exists(dir_out)){
  dir.create(dir_out)
}

fwrite(dt,paste0(dir_out,"AA_mutations.tsv"))


