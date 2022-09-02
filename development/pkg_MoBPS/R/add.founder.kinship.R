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

#' Add a relationship matrix for founder individuals
#'
#' Function to relationship matrix for founder individuals that is used for any calculation of the pedigree
#' @param population population list
#' @param founder.kinship Default is to use vanRaden relationship. Alternative is to enter a pedigree-matrix (order of individuals is first male then female)
#' @param gen Generation for which to enter the pedigree-matrix
#' @examples
#' data(ex_pop)
#' population <- add.founder.kinship(ex_pop)
#' @return Population list
#' @export

add.founder.kinship <- function(population, founder.kinship = "vanRaden", gen=1){

  if(length(founder.kinship)==1 && founder.kinship=="vanRaden"){
    Zt <- get.geno(population, gen=gen)
    if(population$info$miraculix){
      founder.kinship <- miraculix::relationshipMatrix(miraculix::genomicmatrix(Zt), centered = TRUE, normalized = TRUE)
    } else{
      p_i <- rowMeans(Zt) / 2
      Ztm <- Zt - p_i * 2
      founder.kinship <- crossprod(Ztm) / (2 * sum(p_i * (1-p_i)))
    }

  }

  n.animals <- sum(population$info$size[gen,])

  if(ncol(founder.kinship) != n.animals){
    stop("Dimension of founder.kinship and individual number do not match!")
  }

  population$info$founder.kinship <- sort(unique(c(population$info$founder.kinship, gen)))
  if(length(population$info$kinship)==0){
    population$info$kinship <- list()
  }

  population$info$kinship[[gen]] = founder.kinship


  return(population)
}




