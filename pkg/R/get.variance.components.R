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

#' Derive age point
#'
#' Function to devide age point for each individual (Same as time.point unless copy.individual is used for aging)
#' @param population Population list
#' @examples
#' data(ex_pop)
#' get.age.point(ex_pop, gen=2)
#' @return Time point selected gen/database/cohorts-individuals are born
#' @export

get.variance.components <- function(population){

  tmp = population$info$bve.variance.components

  rownames(tmp) = c("genetic variance", "residual variance", "heritability")

  return(tmp)

}
