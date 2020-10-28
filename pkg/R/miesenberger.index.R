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

#' Miesenberger Index
#'
#' Function to selection index weights according to Miesenberger 1999
#' @param V1 Inverted phenotypic covarianz matrix
#' @param V Phenotypic covarianz matrix
#' @param G Genomic covarianz matrix
#' @param RG Genomic correlation matrix
#' @param r reliability for the breeding value estimation
#' @param w relative weighting of each trait (per genetic SD)
#' @param zw Estimated breeding value
#' @return weights of the selection index
#' @export

miesenberger.index <- function(V, G, V1=NULL, RG=NULL, r, w, zw=NULL){
  if(length(V1)==0){
    V1 <- chol2inv(chol(V))
  }
  if(length(RG)==0){
    RG <- sqrt(diag(1/diag(G))) %*% G %*% sqrt(diag(1/diag(G)))
  }

  d <- nrow(G)
  A <- RG * matrix(sqrt(diag(G)), nrow=d, ncol=d, byrow=TRUE) *
    matrix(sqrt(diag(V)), nrow=d, ncol=d, byrow=FALSE) *
    matrix(r, nrow=d, ncol=d, byrow=FALSE)
  bM <- as.numeric(V1 %*% A %*% w)
  if(length(zw)==0){
    return(bM)
  } else{
    IM <- sum( bM * zw )
    return(IM)
  }
}



