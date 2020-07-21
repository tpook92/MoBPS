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

#' Add a genotyping array
#'
#' Function to add a genotyping array for the population
#' @param elements vector with entries to put on the diagonal of a matrix
#' @examples
#' diag(5)
#' @return Diagonal matrix
#' @export

diag.mobps <- function(elements){

  diag.matrix <- diag(elements, nrow=length(elements))
  return(diag.matrix)
}




