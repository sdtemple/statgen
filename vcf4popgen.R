'
File: vcf4popgen
Author: Seth Temple
Description: R functions to process and analyze VCF files for population genetics
'

library(vcfR, quietly = TRUE)
vcf_preproc <- function(file){
  '
  Preprocess a VCF file into allele matrices
  f -  character string with .vcf or .vcf.gz file name
  Return two allele matrices and a allele count/freq matrix
  '
  vcf <- read.vcfR(file, verbose = FALSE)
  
  # number of individuals in the sample population
  nsample <- ncol(vcf@gt) - 1
  
  # number of variants analyzed in the study
  nvariant <- nrow(vcf@gt)
  
  # binary matrices for allele values
  allele1 <- matrix(NA, nrow = nvariant, ncol = nsample,
                    dimnames = list(1:nvariant, 1:nsample))
  allele2 <- matrix(NA, nrow = nvariant, ncol = nsample,
                    dimnames = list(1:nvariant, 1:nsample))
  
  # populate the allele matrices
  for(i in 1:nvariant){
    for(j in 1:nsample){
      gt <- vcf@gt[i,j+1]
      allele1[i,j] <- as.integer(substr(gt,1,1))
      allele2[i,j] <- as.integer(substr(gt,3,3))
    }
  }
  
  # add a count column at the end of the allele matrices
  allele1 <- cbind(allele1, apply(allele1, 1, sum))
  allele2 <- cbind(allele2, apply(allele2, 1, sum))
  
  # determine variants that have at least 1 null allele
  del <- c()
  sumcol <- ncol(allele1)
  for(i in 1:nvariant){
    if(is.na(allele1[i,sumcol])){
      del <- c(del, i)
    } else if(is.na(allele2[i,sumcol])){
      del <- c(del, i)
    }
  }
  
  # remove variant rows that have null alleles
  allele1 <- allele1[-del,]
  allele2 <- allele2[-del,]
  
  # restate the number of variants in the allele matrices
  nvariant <- nrow(allele1)
  
  # pws stands for genome-wide statistic
  pws <- matrix(NA, nrow = nvariant, ncol = 2,
                dimnames = list(1:nvariant, c('ct','p^')))
  for(i in 1:nvariant){
    ct <- allele1[i,105] + allele2[i,105]
    pws[i,1] <- ct
    pws[i,2] <- ct / (2 * nsample)
  }
  
  # remove the counts column in the allele matrices
  allele1 <- allele1[,-sumcol]
  allele2 <- allele2[,-sumcol]
  
  dosage <- dosage_matrices(allele1, allele2)
  
  return(list(a1 = allele1, a2 = allele2, d0 = dosage$d0, d1 = dosage$d1, p = pws))
}

dosage_matrices <- function(a1, a2){
  '
  Create dosage matrices from allele matrices
  a1 - the first allele matrix
  a2 - the second allele matrix
  Return dosage matrices
  '
  d0 <- ((a1 - 1) + (a2 - 1)) * -1
  d1 <- (a1 + a2)
  return(list(d0 = d0, d1 = d1))
}
