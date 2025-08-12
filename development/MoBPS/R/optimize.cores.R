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

#' Optimize generation cores
#'
#' Function to optimize the number of cores used in the generation of new individuals
#' @param population Population list
#' @param test.size Number of individuals to generate for each core number (default: 2500)
#' @param max.cores Maximum number of cores to test (default: 10)
#' @param verbose Set to FALSE to not display any prints
#' @param plot Set to FALSE to not generate a plot of computing times per core
#' @examples
#' population = optimize.cores(max.cores=1, test.size=500)
#' @return Population-list with one or more additional new traits
#' @export

optimize.cores <- function(population = NULL, test.size = 2500, max.cores = 10, verbose= TRUE, plot=TRUE){

  if(length(population)==0){
    warning("Generate new population for the test.\nUse of the genome/trait architecture as in the real simulations is recommended to get realistic results!")
    population = creating.diploid(nsnp=10000, nindi=50, n.additive = c(100,100), chr.nr = 5, verbose = FALSE)
  }

  t = numeric(max.cores)

  for(index in 1:max.cores){
    if(verbose){cat(paste0("Test cores: ", index,"\n"))}
    t[index] = system.time({breeding.diploid(population, breeding.size = test.size, generation.cores = index, verbose=FALSE)})[3]
  }


  if(verbose) {cat(paste0("Suggested number of cores is ", which.min(t)))}
  if(plot) { plot(t, ylab="computing time in seconds", xlab="number of cores", ylim = c(0, max(t)))}
  population$info$generation.cores = which.min(t)

  return(population)
}




