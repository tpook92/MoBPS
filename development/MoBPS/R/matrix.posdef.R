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

#' Projection positive definit
#'
#' Function to project a matrix in the space of positive definit matrices
#' @param A
#' @param verbose
#' @param matrix.name This is just for internal prints
#' @param epsilon This factor is added to the diagonal to avoid numerical issues with semi-definit matrices
#' @examples
#' @return Positive definit matrix
#' @export

matrix.posdef <- function(A, verbose = TRUE, matrix.name = "Matrix", epsilon = 0.0000001){

  test = eigen(A)

  if(min(test$values) < 0){
    test$values[test$values<0] <- 0

    if(verbose){
      warning(paste0(matrix.name, " is not positive definit."))
      cat(paste0(matrix.name, " is not positive definit.\n"))
    }
    if(verbose) cat("Generate projection on the set of positive definit matrices:")


    M <- diag(test$values)

    S <- test$vectors

    newA <- S %*% M %*% solve(S)

    diag(newA) <- diag(newA) + epsilon # Avoid numerical issues with inversion
    newA <- newA * matrix(1/sqrt(diag(newA)), nrow=nrow(newA), ncol=nrow(newA), byrow=TRUE) * matrix(1/sqrt(diag(newA)), nrow=nrow(newA), ncol=nrow(newA), byrow=FALSE)
    if(verbose) cat("new suggested genetic correlation matrix:\n")
    if(verbose) print(round(newA, digits=3))

    return(newA)
  } else{
    return(A)
  }



}
