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


#' Single Step GBLUP
#'
#' Function to perform single step GBLUP according to Legarra 2014
#' @param A11 pedigree relationship matrix of non-genotyped individuals
#' @param A12 pedigree relationship matrix between non-genotyped and genotyped individuals
#' @param A22 pedigree relationship matrix of genotyped individuals
#' @param G genomic relationship matrix of genotyped individuals
#' @return Single step relationship matrix

ssGBLUP <- function(A11, A12, A22, G){



  A22 <- add.diag(A22, 0.001) # numeric stability
  A22inv <- chol2inv(chol(A22))
  A21 <- t(A12)
  A12A22inv <- A12 %*% A22inv
  H12 <- A12A22inv %*% G
  H11 <- A11 -  A12A22inv %*% A21 + H12 %*% t(A12A22inv)
  H22 <- G

  H <- rbind(cbind(H11,H12), cbind(t(H12), H22))

  return(H)
}
