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


#' Estimate effective population size
#'
#' Internal function to estimate the effective population size
#' @param ld ld between markers
#' @param dist distance between markers in Morgan
#' @param n Population size
#' @return Estimated effective population size
#'

effective.size <- function(ld, dist, n){

  c <- 0
  for(index in 2*(1:max(5,dist*2))-1){
    c <- c + stats::dpois(index, lambda = dist)
  }

  return( ((1-c)^2 + c^2) / (( ld - 1/n) * (2 * c*(2-c))  ))
}
