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

#' Row-wise Median
#'
#' Function to calculate row-wise median values
#' @param A Matrix
#' @return row-wise median values of a matrix
#' @examples
#' A <- matrix(c(1,2,3,4), ncol=2)
#' x <- rowMedian(A)
#' @return Matrix with modified diagonal entries
#' @export


rowMedian = function(A) {

  x = apply(X = A, FUN = stats::median, MARGIN = 1)
  return(x)

}
