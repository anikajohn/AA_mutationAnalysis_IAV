library(tidyverse)
# library(Biostrings)
library(seqinr)


#' reading lofreq vcf file
#'
#' @param file_path: path to vcf file generated from lofreq run 
#'
#' @return cleaned and formatted dt, if vcf is empty NULL is returned
#' @export
#'
#' @examples
read_vcf <- function(file_path){
  
  dt = fread(file_path, 
             skip = grep("^#", readLines(file_path), value = TRUE)[length(grep("^#", readLines(file_path)))])
  
  #header has 15 lines, if 15 or less, vcf is empty
  if(nrow(dt) <= 15){
    
    return(NULL)
    
  }else{
  
  #different output format with parallelized lofreq

  #dt[, c('DP', 'AF', 'SB', 'DP4', 'V13', 'V14') := tstrsplit(INFO, ';', fixed = TRUE)]
  dt[, c('DP', 'AF', 'SB') := tstrsplit(INFO, ';', fixed = TRUE)[1:3]]
  #dt[, c('DP4', 'V13', 'V14','INFO') := NULL]
  dt[, c('INFO') := NULL]
  dt[, `:=`(
    AF = str_remove_all(AF, "AF="),
    DP = str_remove_all(DP, "DP="),
    SB = str_remove_all(SB, "SB=")
  )]
  
  # Convert columns to numeric
  dt[, c('POS', 'QUAL', 'DP', 'AF', 'SB') := lapply(.SD, as.numeric), 
     .SDcols = c('POS', 'QUAL', 'DP', 'AF', 'SB')]
  
  return(dt)
  }
}



#' reading lofreq vcf file and formatting in v-pipe context
#'
#' @param file_name: path to vcf file generated from lofreq run 
#' @param regex_sample: regex with which sample name can be extracted (v-pipe set-up)
#'
#' @return formatted vcf dt, witch additional columns location_code and date
#' @export
#'
#' @examples
#' 
read_and_mark_vcf <- function(file_path, regex_sample = "[A-Z]\\d{1}_\\d{2}_\\d{4}_\\d{2}_\\d{2}"){
  
  dt = read_vcf(file_path)
  
  if(is.null(dt)){
    
    return(NULL)
    
  }else{
    
  
    
    dt[, sample_name := str_extract(file_path,regex_sample)]
    #only care about point mutations
    dt = dt[nchar(REF) == 1 & nchar(ALT) == 1]
    dt[, c('cell','location_code', 'year','month','day') := tstrsplit(sample_name, "_", fixed=TRUE)]
    dt[, c('date') := paste(year,month,day,sep = "-")]
    dt[, c('year','month','day') := NULL]
  
  return(dt)
    
  }
}




