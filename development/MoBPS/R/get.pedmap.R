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

#' Generate plink-file (pedmap)
#'
#' Generate a ped and map file (PLINK format) for selected groups and chromosome
#' @param population Population list
#' @param path Location to save pedmap-file
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param non.genotyped.as.missing Set to TRUE to replaced non-genotyped entries with "./."
#' @param use.id Set to TRUE to use MoBPS ids instead of Sex_Nr_Gen based names
#' @examples
#' data(ex_pop)
#' \donttest{
#' file_path <- tempdir()
#' get.pedmap(path=file_path, ex_pop, gen=2)
#' file.remove(paste0(file_path, ".ped"))
#' file.remove(paste0(file_path, ".map"))
#' }
#' @return Ped and map-file for in gen/database/cohorts selected individuals
#' @export

get.pedmap <- function(population, path=NULL, database=NULL, gen=NULL, cohorts=NULL, non.genotyped.as.missing=FALSE,
                       use.id=FALSE){


  haplo <- get.haplo(population, database=database, gen=gen, cohorts=cohorts, export.alleles=FALSE)
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
  mapfile <- cbind(chr.nr, population$info$snp.name, unlist(population$info$position) ,population$info$bp)
  mapfile[is.na(mapfile)] <- 0
  haplo1 <- t(haplo[,(1:(ncol(haplo)/2))*2-1])
  haplo2 <- t(haplo[,(1:(ncol(haplo)/2))*2])
  if(non.genotyped.as.missing){

    is_genotyped <- t(get.genotyped.snp(population, gen=gen, database = database, cohorts=cohorts))
    if(sum(!is_genotyped)>0){
      haplo1[!is_genotyped] <- "N"
      haplo2[!is_genotyped] <- "N"
    }

  }
  ped <- cbind(haplo1, haplo2)
  ped <- ped[,c(0,ncol(haplo1))+ rep(1:ncol(haplo1), each=2)]
  ped[ped==0] <- "A"
  ped[ped==1] <- "C"

  family.base <- get.database(population, gen=gen, database=database, cohorts=cohorts)
  n.animals <- sum(family.base[,4]-family.base[,3]+1)
  family <- sex.s <- numeric(n.animals)
  tillnow <- 1
  for(index in 1:nrow(family.base)){
    if(diff(family.base[index,3:4])>=0){
      family[tillnow:(tillnow+diff(family.base[index,3:4]))] <- index
      sex.s[tillnow:(tillnow+diff(family.base[index,3:4]))] <- family.base[index,2]
      tillnow <- tillnow + diff(family.base[index,3:4]) +1
    }
  }
  pedi <- get.pedigree(population, database = family.base, id=use.id)
  pedfile <- cbind(family, pedi[,1],pedi[,2],pedi[,3],sex.s,0,ped)


  pedname <- paste0(path,".ped")
  utils::write.table(file=pedname, pedfile, col.names = FALSE, row.names = FALSE, quote=FALSE)
  mapname <- paste0(path,".map")
  utils::write.table(file=mapname, mapfile, col.names = FALSE, row.names = FALSE, quote=FALSE)
}

