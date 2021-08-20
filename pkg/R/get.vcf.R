'#
  Authors
Torsten Pook, torsten.pook@uni-goettingen.de

Copyright (C) 2017 -- 2020  Torsten Pook

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
'#

#' Generate vcf-file
#'
#' Generate a vcf-file for selected groups and chromosome
#' @param population Population list
#' @param path Location to save vcf-file
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param chromosomen Beschraenkung des Genotypen auf bestimmte Chromosomen (default: 1)
#' @param non.genotyped.as.missing Set to TRUE to replaced non-genotyped entries with "./."
#' @param use.id Set to TRUE to use MoBPS ids instead of Sex_Nr_Gen based names
#' @examples
#' data(ex_pop)
#' \donttest{get.vcf(path=tempdir(), ex_pop, gen=2)}
#' @return VCF-file for in gen/database/cohorts selected individuals
#' @export

get.vcf <- function(population, path=NULL, database=NULL, gen=NULL, cohorts=NULL, chromosomen="all",
                    non.genotyped.as.missing=FALSE, use.id = FALSE){

  haplo <- get.haplo(population, database=database, gen=gen, cohorts=cohorts, chromosomen=chromosomen, export.alleles=FALSE)
  # haplo <- get.haplo(population, gen=1)
  if(length(path)==0){
    path <- "population.vcf"
  } else{
    path <- paste0(path,".vcf")
  }
  chr.nr <- sum(population$info$snp)
  start <- 1
  for(index in 1:length(population$info$snp)){
    if(population$info$snp[index]>0){
      chr.nr[start:(start+population$info$snp[index]-1)] <- index
      start <- start + population$info$snp[index]
    }
  }
  bp <- population$info$bp
  snpname <- population$info$snp.name
  ref <- population$info$snp.base[1,]
  alt <- population$info$snp.base[2,]

  vcfgeno <- matrix(paste0(haplo[,(1:(ncol(haplo)/2))*2], "|", haplo[,(1:(ncol(haplo)/2))*2-1]), ncol=ncol(haplo)/2)

  if(non.genotyped.as.missing){
    is_genotyped <- get.genotyped.snp(population, gen=gen, database = database, cohorts=cohorts)

    if(sum(!is_genotyped)>0){
      vcfgeno[!is_genotyped] <- "./."
    }

  }


  ref[ref==0] <- "A"
  ref[ref==1] <- "C"
  alt[alt==0] <- "A"
  alt[alt==1] <- "C"
  options(scipen=999)
  vcfgenofull <- cbind(chr.nr, as.numeric(bp), snpname, ref, alt, ".", "PASS", ".", "GT", vcfgeno)
  vcfgenofull <- rbind(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", get.pedigree(population, database=database, gen=gen, cohorts = cohorts, id=use.id)[,1]),vcfgenofull)

  headerfile <- rbind(
    "##fileformat=VCFv4.2",
    gsub("-", "", paste0("##filedate=",  Sys.Date())),
    paste0("##source='MoBPS_", utils::packageVersion("MoBPS"),"'"),
    '##FILTER=<ID=PASS,Description="all filters passed">'
  )

  nchr <- unique(chr.nr)
  contigs <- numeric(length(nchr))
  for(index in 1:length(nchr)){
    contigs[index] <- paste0("##contig=<ID=", nchr[index],",length=", max(as.numeric(bp)[chr.nr==nchr[index]]),">")
  }

  headerfile <- rbind(headerfile, t(t(contigs)),     "##FORMAT=<ID=GT,Number=1,Type=String,Description='Genotype'>")



  utils::write.table(headerfile, file=path, quote=FALSE, col.names = FALSE, row.names = FALSE)
  utils::write.table(vcfgenofull, file=path, quote=FALSE, col.names = FALSE, row.names = FALSE, append = TRUE, sep="\t")
}


