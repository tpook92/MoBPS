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

#' Derive ID on an individual
#'
#' Function to derive the internal ID given to each individual
#' @param population Population list
#' @param plot Set to FALSE to not display overview of estimated SNP effects (default: TRUE)
#' @param extend Set to TRUE to export u_hat estimates from all breeding values instead of just the last (default: FALSE)
#' @param trait.plot Select trait for which to generate the visualization (default: 1)
#' @return matrix with estimated marker effects
#' @export

get.uhat <- function(population, extend = FALSE, plot = TRUE, trait.plot = 1){

  if(plot){

    oldpar <- graphics::par(no.readonly=TRUE)
    on.exit(graphics::par(oldpar))
    graphics::par(mfrow= c(2,1))

    u_hat = (population$info$u_hat[[length(population$info$u_hat)]][,trait.plot])

    graphics::plot(u_hat)
    graphics::lines(stats::ksmooth(1:length(u_hat), abs(u_hat), bandwidth = max(10,round(length(u_hat)/100))), col="red", lwd=2)
    graphics::abline(h=0, col="red", lty=2)
    graphics::hist(u_hat, nclass = min(100, length(u_hat)/10))
  }

  if(extend){
    return(population$info$u_hat)
  } else{
    return(population$info$u_hat[[length(population$info$u_hat)]])
  }


}



