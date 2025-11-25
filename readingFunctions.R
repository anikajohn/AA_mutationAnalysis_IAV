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
#read_vcf <- function(file_path){
#  
#  dt = fread(file_path, 
#             skip = grep("^#", readLines(file_path), value = TRUE)[length(grep("^#", readLines(file_path)))])
#  # Check if the file exists
#  if (!file.exists(file_path)) {
#    warning(paste("File does not exist.",file_path))
#    return(NULL)
#  }
#  
#  # Check if the file is empty
#  if (file.size(file_path) == 0) {
#    warning(paste("File is empty.",file_path))
#    return(NULL)
#  }
#  # Try reading the VCF file
#  tryCatch({
#    dt = fread(file_path,
#               skip = grep("^#", readLines(file_path), value = TRUE)[length(grep("^#", readLines(file_path)))])
#
#    # Header has 15 lines, if 15 or less, VCF is empty
#    if (nrow(dt) <= 17) {
#      return(NULL)
#    } else {
#      dt[, c('DP', 'AF', 'SB', 'DP4', 'V13', 'V14') := tstrsplit(INFO, ';', fixed = TRUE)]
#      dt[, c('DP4', 'V13', 'V14','INFO') := NULL]
#      dt[, `:=`(
#        AF = str_remove_all(AF, "AF="),
#        DP = str_remove_all(DP, "DP="),
#        SB = str_remove_all(SB, "SB=")
#      )]
#      dt[, c('POS', 'QUAL', 'DP', 'AF', 'SB') := lapply(.SD, as.numeric),
#         .SDcols = c('POS', 'QUAL', 'DP', 'AF', 'SB')]
#
#      return(dt)
#    }
#  }, error = function(e) {
#    warning("Error reading VCF file: ", e$message)
#    return(NULL)
#  })
#}


read_vcf <- function(file_path){
   # Check if the file exists
  if (!file.exists(file_path)) {
    warning(paste("File does not exist.",file_path))
    return(NULL)
  }

  # Check if the file is empty
  if (file.size(file_path) == 0) {
    warning(paste("File is empty.",file_path))
    return(NULL)
  }
  # Try reading the VCF file
  tryCatch({

  #if vcf has no lines w/o # == Empty
  l = length(grep("^[^#]", readLines(file_path)))

  if(l==0){

    return(NULL)

    }else{

    dt = fread(file_path,skip="#CHROM")

    #different output format with parallelized lofre
    dt[, c('DP', 'AF', 'SB') := tstrsplit(INFO, ';', fixed = TRUE)[1:3]]
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
  },
    warning = function(w) {
    cat("Warning in file:", file_path, "\n")
    cat(conditionMessage(w), "\n\n")
    },
    error = function(e) {
    warning("Error reading VCF file: ", e$message)
    cat("Error in file:", file_path, "\n")
    cat(conditionMessage(e), "\n\n")
    return(NULL)
  })
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
read_and_mark_vcf <- function(file_path, regex_sample = "\\d{2}_\\d{4}_\\d{2}_\\d{2}"){
  
    # Check if the file exists
  if (!file.exists(file_path)) {
    warning(paste("File does not exist.",file_path))
    return(NULL)
  }

  # Check if the file is empty
  if (file.size(file_path) == 0) {
    warning(paste("File is empty.",file_path))
    return(NULL)
  }

  dt = read_vcf(file_path)
  
  if(is.null(dt)){
    
    return(NULL)
    
  }else{
    
  
    
    dt[, sample_name := str_extract(file_path,regex_sample)]
    #only care about point mutations
    print(paste("processing file:",file_path))
    dt = dt[nchar(REF) == 1 & nchar(ALT) == 1]
    dt[, c('location_code', 'year','month','day') := tstrsplit(sample_name, "_", fixed=TRUE)]
    dt[, c('date') := paste(year,month,day,sep = "-")] 
    dt[, c('year','month','day') := NULL]
  
  return(dt)
    
  }
}




