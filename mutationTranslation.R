library(tidyverse)
library(Biostrings)


replace_nuc <- function(DNA_String,index,replacement){
  
  substr(DNA_String, index, index) <- replacement
  
  return(DNA_String)
  
}

AAtranslate_vcf_df <- function(fasta_path, start, end, vcf_df){
  
  genetic_code <- c(
    ATA = 'I', ATC = 'I', ATT = 'I', ATG = 'M',
    ACA = 'T', ACC = 'T', ACG = 'T', ACT = 'T',
    AAC = 'N', AAT = 'N', AAA = 'K', AAG = 'K',
    AGC = 'S', AGT = 'S', AGA = 'R', AGG = 'R',
    CTA = 'L', CTC = 'L', CTG = 'L', CTT = 'L',
    CCA = 'P', CCC = 'P', CCG = 'P', CCT = 'P',
    CAC = 'H', CAT = 'H', CAA = 'Q', CAG = 'Q',
    CGA = 'R', CGC = 'R', CGG = 'R', CGT = 'R',
    GTA = 'V', GTC = 'V', GTG = 'V', GTT = 'V',
    GCA = 'A', GCC = 'A', GCG = 'A', GCT = 'A',
    GAC = 'D', GAT = 'D', GAA = 'E', GAG = 'E',
    GGA = 'G', GGC = 'G', GGG = 'G', GGT = 'G',
    TCA = 'S', TCC = 'S', TCG = 'S', TCT = 'S',
    TTC = 'F', TTT = 'F', TTA = 'L', TTG = 'L',
    TAC = 'Y', TAT = 'Y', TAA = '*', TAG = '*',
    TGC = 'C', TGT = 'C', TGA = '*', TGG = 'W',
    X = 'X' # X = no amino acid
  )
  
  ref_DNA <- readDNAStringSet(fasta_path)
  ref_DNA <- ref_DNA[[1]]
 
  dna <-  ref_DNA[start:end] %>% 
    as.character() 
  codons <- regmatches(dna, gregexpr(".{3}", dna))[[1]] 
  
  Nuc_readingFrame <- seq(start,end,1)
  len_AAreadingFrame <- length(Nuc_readingFrame)/3
  NucToCodon <- rep(1:len_AAreadingFrame, each = 3)
  codon_index <- rep(1:3,length.out = length(Nuc_readingFrame))
  gene_trans <- data.frame(gene_nucIDX = Nuc_readingFrame, 
                           gene_codonIDX = NucToCodon,
                           codonIDX = codon_index)
  #standard vcf format
  vcf <- vcf_df %>% 
    merge(.,gene_trans,
          by.x=c('POS'),
          by.y=c('gene_nucIDX')) %>% 
    mutate(org_codon = codons[gene_codonIDX]) %>% 
    mutate(new_codon = replace_nuc(org_codon,codonIDX,ALT),
           new_AA = as.character(genetic_code[new_codon]),
           org_AA = as.character(genetic_code[org_codon]))
  
  return(vcf)
  
}




identify_AA <- function(fasta_file,
                        gene_start,
                        gene_end,
                        vcf_dt){
  dt = vcf_dt
  dt = dt[nchar(REF) == 1 & nchar(ALT) == 1]

  dt <- AAtranslate_vcf_df(fasta_file, gene_start, gene_end, dt) %>% setDT()

  # only point muyations
  dt = dt[nchar(REF) == 1 & nchar(ALT) == 1]

  return(dt)
  
}

