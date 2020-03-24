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

#' Map generation
#'
#' Function to derive the genomic map for a given population list
#' @param population Population list
#' @param use.snp.nr Set to TRUE to display SNP number and not SNP name
#' @examples
#' data(ex_pop)
#' map <- get.map(ex_pop)
#' @return Genomic map of the population list
#' @export

get.map <- function(population, use.snp.nr=FALSE){
  chr.nr <- snp.nr <- numeric(sum(population$info$snp))
  start <- 1
  for(index in 1:length(population$info$snp)){
    if(population$info$snp[index]>0){
      chr.nr[start:(start+population$info$snp[index]-1)] <- index
      snp.nr[start:(start+population$info$snp[index]-1)] <- 1:population$info$snp[index]
      start <- start + population$info$snp[index]
    }
  }
  if(use.snp.nr){
    mapfile <- cbind(chr.nr, snp.nr, 0 , as.numeric(population$info$bp))
  } else{
    mapfile <- cbind(chr.nr, population$info$snp.name, 0 , as.numeric(population$info$bp))
  }

  mapfile[is.na(mapfile)] <- 0

  return(mapfile)
}

