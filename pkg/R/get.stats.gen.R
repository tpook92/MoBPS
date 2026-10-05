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

#' Development of genetic/breeding value
#'
#' Function to plot genetic/breeding values for multiple generation/cohorts
#' @param population population list
#' @param gen Which generation to include in the analysis (default: "all")
#' @param selection.index.weights Selection index for Genetic gain calculation
#' @param selection.index.weights.acc Weights in the overall accuracy calculation
#' @param plot Set to FALSE to not generate a visualization of outcomes
#' @param mar Margins in plot
#' @param mfrow mfrow in plot
#' @param skip.inbreeding Set to TRUE to skip inbreeding calcuation (will be much faster)
#' @examples
#' data(ex_pop)
#' get.stats.gen(ex_pop)
#' @return Estimated of avg. kinship/inbreeding based on IBD/HBD
#' @export

get.stats.gen <- function(population, gen = "all",
                              selection.index.weights = NULL,
                              selection.index.weights.acc = NULL,
                              plot = TRUE,
                              mar = c(4.1,4.1,4.1,2.1),
                              mfrow = c(1,3),
                          skip.inbreeding = FALSE){

  if(length(gen)==1 && gen[1] == "all"){
    gen = 1:get.ngen(population)
  }

  if(length(selection.index.weights)==0){
    selection.index.weights = rep(1, population$info$bv.nr)
  }

  if(length(selection.index.weights.acc)==0){
    selection.index.weights.acc = selection.index.weights
  }

  bv_index = acc_index = numeric(length(gen))
  bv = matrix(NA, nrow = length(gen), ncol = population$info$bv.nr)
  bve = matrix(NA, nrow = length(gen), ncol = population$info$bv.nr)
  pheno = matrix(NA, nrow = length(gen), ncol = population$info$bv.nr)



  inbreeding = numeric(length(gen))
  acc = matrix(NA, nrow = length(gen), ncol = population$info$bv.nr)
  for(index in 1:length(gen)){
    bv[index,] = rowMeans(get.bv(population, gen = gen[index]))
    bve[index,] = rowMeans(get.bve(population, gen = gen[index]))
    pheno[index,] = rowMeans(get.pheno(population, gen = gen[index]))
    if(!skip.inbreeding){
      inbreeding[index] = mean(inbreeding.emp(population, gen = gen[index]))
    }
    acc[index,] = analyze.bv(population, gen = gen[index])[[1]][1,]

    bv_index[index] = sum( bv[index,] * selection.index.weights)
    acc_index[index] = sum( acc[index,] * selection.index.weights / sum(selection.index.weights))

  }

  names(bv_index) = names(inbreeding) = gen

  rownames(bv) = rownames(bve) = rownames(pheno) = rownames(acc) = gen
  colnames(bv) = colnames(bve) = colnames(pheno) = colnames(acc) = population$info$trait.name


  outputs = list(
    gen,
    bv_index,
    inbreeding,
    acc_index,
    bv,
    bve,
    pheno,
    acc
  )


  if(plot){

    oldpar <- graphics::par(no.readonly=TRUE)
    on.exit(graphics::par(oldpar))

    graphics::par(mar=mar)

    graphics::par(mfrow = mfrow)

    graphics::plot( gen, outputs[[2]], ylab = "genetic values", type = "l",
          xlab = "", main = "Genetic gain")
    graphics::points( gen,outputs[[2]])



    graphics::plot( gen,outputs[[3]], ylab = "inbreeding", type = "l",  xlab = "",
          main = "inbreeding")
    graphics::points( gen,outputs[[3]])


    tmp =  outputs[[4]]
    tmp[is.na(tmp)] = 0
    graphics::plot(gen, tmp, ylab = "prediction accuracy", type = "l",
         xlab = "", main = "prediction accuracy")
    graphics::points(gen, tmp)


  }


  return(outputs)
}




