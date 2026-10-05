'#
  Authors
Torsten Pook, torsten.pook@wur.nl

Copyright (C) 2017 -- 2025  Torsten Pook

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

#' Map generation
#'
#' Function to derive the genomic map for a given population list
#' @param population Population list
#' @param use.snp.nr Set to TRUE to display SNP number and not SNP name
#' @param morgan.position.per.chromosome Set to FALSE to Morgan position continuously over the genome
#' @param array Use only markers available on the array
#' @param mixblup.last Use only SNPs included in the last MiXBLUP evaluation (that were not filtered through MAF etc.)
#' @examples
#' data(ex_pop)
#' map <- get.map(ex_pop)
#' @return Genomic map of the population list
#' @export

get.map <- function(population, use.snp.nr=FALSE, morgan.position.per.chromosome = TRUE,
                    array = NULL,
                    mixblup.last = FALSE){
  chr.nr <- snp.nr <- cm.positions <-numeric(sum(population$info$snp))
  start <- 1
  for(index in 1:length(population$info$snp)){
    if(population$info$snp[index]>0){
      chr.nr[start:(start+population$info$snp[index]-1)] <- index
      snp.nr[start:(start+population$info$snp[index]-1)] <- 1:population$info$snp[index]
      cm.positions[start:(start+population$info$snp[index]-1)] <- population$info$position[[index]]
      start <- start + population$info$snp[index]
    }
  }

  if(morgan.position.per.chromosome==FALSE){
    cm.positions <- population$info$snp.position
  }

  if(use.snp.nr){
    mapfile <- cbind(chr.nr, snp.nr, cm.positions , as.numeric(population$info$bp))
  } else{
    mapfile <- cbind(chr.nr, population$info$snp.name, cm.positions , as.numeric(population$info$bp))
  }

  mapfile[is.na(mapfile)] <- 0

  colnames(mapfile) <- c("Chromosome", "SNP ID", "Morgan position", "BP position")

  if(length(array)>0){
    temp1 <- which(population$info$array.name == array)
    if(length(temp1)==1){
      array <- temp1
    }
    if(!is.numeric(array)){
      stop("Input for array can not be assigned! Check your input!")
    }
  }
  relevant.snps = 1:nrow(mapfile)

  if(mixblup.last & length(population$info$mixblup_lastfreq) > 0){
    relevant.snps = intersect(relevant.snps, population$info$mixblup_lastfreq[,1])
  }
  if(length(array)>0){
    relevant.snps = intersect(which(population$info$array.markers[[array]]), relevant.snps)
  }

  return(mapfile)
}

