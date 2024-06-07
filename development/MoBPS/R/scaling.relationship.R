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

#' scaling.relationship
#'
#' scaling.relationship
#' @param A Population list
#' @param Z AA
#' @param p AA
#' @return scaled genomic relationship matrix
#' @export

scaling.relationship = function(A, Z, p){
  p1 = 4 * sum(p * p )
  if(sum(class(Z) == "genomicmatrix")>0){
    p2 = 2 * miraculix::vectorGeno(p, Z)
  } else{
    p2 = 2 * colSums(Z * p)
  }
  A = A + p1 -p2 - matrix(p2, nrow= nrow(A ), ncol = ncol(A), byrow = TRUE)
  A = A / (2 * sum(p*(1-p)))

  if(sum(p*(1-p))==0){
    A[is.na(A)] = 0
  }

  return(A)
}
