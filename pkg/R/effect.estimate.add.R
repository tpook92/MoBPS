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

#' Estimation of marker effects
#'
#' Function to estimate marker effects
#' @param geno genotype dataset (marker x individuals)
#' @param pheno phenotype dataset (each phenotype in a row)
#' @param map genomic map
#' @param scaling Set FALSE to not perform variance scaling
#' @examples
#' data(ex_pop)
#' pheno <- get.pheno(ex_pop, gen=1:5)
#' geno <- get.geno(ex_pop, gen=1:5)
#' map <- get.map(ex_pop, use.snp.nr=TRUE)
#' \donttest{real.bv.add <- effect.estimate.add(geno, pheno, map)}
#' @return Empirical kinship matrix (IBD-based since Founders)
#' @export
#'

effect.estimate.add <- function(geno, pheno, map = NULL, scaling=TRUE){

  if(length(map)==0){
    map <- cbind(1:nrow(geno), 1)
  }
  if(!is.matrix(pheno) || ncol(geno)!=nrow(pheno)){
    pheno <- matrix(pheno, nrow=ncol(geno))
  }

  n <- nrow(pheno)

  if (requireNamespace("miraculix", quietly = TRUE)) {
    A <- miraculix::relationshipMatrix(miraculix::genomicmatrix(geno), centered=TRUE, normalized=TRUE)
  } else{
    p_i <- rowSums(geno)/ncol(geno)/2
    Ztm <- geno - p_i * 2
    A <- crossprod(Ztm)/ (2 * sum(p_i*(1-p_i)))
  }


  real.bv.add <- list()
  for(bven in 1:ncol(pheno)){

    if(requireNamespace("EMMREML", quietly = TRUE)){
      fm <- EMMREML::emmreml(
        pheno[,bven],
        matrix(1,nrow=n),
        diag(n),
        A)
    } else{
      stop("Usage of EMMREML without being installed!")
    }

    u_hat <- alpha_to_beta(drop(fm$uhat),A,t(geno))


    if(scaling){
      y_hat_test <- t(u_hat) %*% geno
      scaling_factor <- sqrt(fm$Vu) / sd(y_hat_test)
      u_hat <- u_hat * scaling_factor
    }

    real.bv.add[[bven]] <- cbind(map[,2:1], -u_hat, 0, u_hat)
    real.bv.add[[bven]] <- real.bv.add[[bven]][real.bv.add[[bven]][,3]!=0,]

  }



  return(real.bv.add)

}
