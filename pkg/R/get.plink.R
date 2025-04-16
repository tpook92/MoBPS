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
#' @param chromosome Limit the genotype output to a selected chromosome (default: "all")
#' @param non.genotyped.as.missing Set to TRUE to replaced non-genotyped entries with "./."
#' @param use.id Set to TRUE to use MoBPS ids instead of Sex_Nr_Gen based names
#' @param fam.id Set TRUE to set the fam-ID to the individual ID
#' @param type Set 1 to only write paternal haplotype, set 2 to only write maternal haplotype, set 0 for both (default)
#' @param bve.pedigree.error Set to FALSE to ignore/correct for any pedigree errors
#' @examples
#' data(ex_pop)
#' data(ex_pop)
#' \donttest{
#' file_path <- tempdir()
#' get.vcf(path=file_path, ex_pop, gen=2)
#' file.remove(paste0(file_path, ".vcf"))
#' }
#' @return VCF-file for in gen/database/cohorts selected individuals
#' @export

get.plink <- function(population, path=NULL, database=NULL, gen=NULL, cohorts=NULL, chromosome="all",
                      non.genotyped.as.missing=FALSE, fam.id = FALSE, type = 0, use.id = TRUE,
                      bve.pedigree.error = TRUE){

  database = get.database(population, gen = gen, database = database, cohorts = cohorts)
  nindi = get.nindi(population, database = database)
  if(type ==0){
    geno <- get.geno(population, database=database, chromosome=chromosome, export.alleles=FALSE, use.id = use.id)

  } else if(type ==1){
    geno <- get.haplo(population, database=database, chromosome=chromosome, export.alleles=FALSE, use.id = use.id)[,1:nindi*2-1]

    colnames(geno) = substr(colnames(geno), start = 1, stop = nchar(colnames(geno))-5)
  } else if(type ==2){
    geno <- get.haplo(population, database=database, chromosome=chromosome, export.alleles=FALSE, use.id = use.id)[,1:nindi*2]
    colnames(geno) = substr(colnames(geno), start = 1, stop = nchar(colnames(geno))-5)

  }

  # haplo <- get.haplo(population, gen=1)
  if(length(path)==0){
    path <- "population"
  }
  chr.nr <- numeric(sum(population$info$snp))
  start <- 1
  for(index in 1:length(population$info$snp)){
    if(population$info$snp[index]>0){

      if(length(population$info$chromosome.name)>=index){
        chr.nr[start:(start+population$info$snp[index]-1)] <- population$info$chromosome.name[index]
      } else{
        chr.nr[start:(start+population$info$snp[index]-1)] <- index
      }

      start <- start + population$info$snp[index]
    }
  }
  bp <- population$info$bp
  snpname <- population$info$snp.name
  ref <- population$info$snp.base[1,]
  alt <- population$info$snp.base[2,]

  if(non.genotyped.as.missing){
    is_genotyped <- get.genotyped.snp(population, database = database)
    geno[!is_genotyped] = NA
  }


  ref[ref==0] <- "A"
  ref[ref==1] <- "C"
  alt[alt==0] <- "A"
  alt[alt==1] <- "C"

  options("scipen"=999)
  ped_tmp = get.pedigree(population, database = database, id = use.id, include.error = bve.pedigree.error)
  sex_tmp = get.pedigree(population, database = database, raw = TRUE, include.error = bve.pedigree.error)[,2]

  if(fam.id){
    fam_tmp = ped_tmp[,1]
  } else{
    fam_tmp = rep(1, nrow(ped_tmp))
  }

  fam = data.frame(fam = fam_tmp, id = ped_tmp[,1], pat = ped_tmp[,2],
                   mat = ped_tmp[,3], sex = sex_tmp, pheno = rep(0, nrow(ped_tmp)))

  map_tmp = cbind(get.map(population), ref , alt)


  bim = data.frame(chr = map_tmp[,1], id = map_tmp[,2], posg = map_tmp[,3],
                   pos = map_tmp[,4], ref = rep("A", nrow(map_tmp)), alt = rep("C", nrow(map_tmp)))

  genio::write_plink(file=path, X = geno, fam = fam, bim = bim)

}









