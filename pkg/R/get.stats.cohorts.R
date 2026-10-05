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
#' @param cohorts Which cohorts to include in the analysis (default: "all")
#' @param cohort.filter Only include cohorts that include this character string in their name
#' @param selection.index.weights Selection index for Genetic gain calculation
#' @param selection.index.weights.acc Weights in the overall accuracy calculation
#' @param plot Set to FALSE to not generate a visualization of outcomes
#' @param mar Margins in plot
#' @param mfrow mfrow in plot
#' @param skip.inbreeding Set to TRUE to skip inbreeding calcuation (will be much faster)
#' @examples
#' data(ex_pop)
#' get.stats.cohorts(ex_pop, cohort.filter = "_M")
#' @return Estimated of avg. kinship/inbreeding based on IBD/HBD
#' @export

get.stats.cohorts <- function(population, cohorts = "all", cohort.filter = NULL,
                              selection.index.weights = NULL,
                              selection.index.weights.acc = NULL,
                              plot = TRUE,
                              mar = c(8.1,4.1,4.1,2.1),
                              mfrow = c(1,3),
                              skip.inbreeding = FALSE){

  if(length(cohorts)==1 && cohorts[1] == "all"){
    cohorts = get.cohorts(population)
  }
  if(length(cohort.filter) > 0){

    included = rep(FALSE, length(cohorts))
    for(filter in cohort.filter){
      included[grepl(filter, cohorts)] = TRUE
    }
    cohorts = cohorts[included]
  }

  if(length(selection.index.weights)==0){
    selection.index.weights = rep(1, population$info$bv.nr)
  }

  if(length(selection.index.weights.acc)==0){
    selection.index.weights.acc = selection.index.weights
  }

  bv_index = acc_index = numeric(length(cohorts))
  bv = matrix(NA, nrow = length(cohorts), ncol = population$info$bv.nr)
  bve = matrix(NA, nrow = length(cohorts), ncol = population$info$bv.nr)
  pheno = matrix(NA, nrow = length(cohorts), ncol = population$info$bv.nr)



  inbreeding = numeric(length(cohorts))
  acc = matrix(NA, nrow = length(cohorts), ncol = population$info$bv.nr)
  for(index in 1:length(cohorts)){
    bv[index,] = rowMeans(get.bv(population, cohorts = cohorts[index]))
    bve[index,] = rowMeans(get.bve(population, cohorts = cohorts[index]))
    pheno[index,] = rowMeans(get.pheno(population, cohorts = cohorts[index]))
    if(!skip.inbreeding){
      inbreeding[index] = mean(inbreeding.emp(population, cohorts = cohorts[index]))
    }
    acc[index,] = analyze.bv(population, cohorts = cohorts[index])[[1]][1,]

    bv_index[index] = sum( bv[index,] * selection.index.weights)
    acc_index[index] = sum( acc[index,] * selection.index.weights / sum(selection.index.weights))

  }

  names(bv_index) = names(inbreeding) = cohorts

  rownames(bv) = rownames(bve) = rownames(pheno) = rownames(acc) = cohorts
  colnames(bv) = colnames(bve) = colnames(pheno) = colnames(acc) = population$info$trait.name


  outputs = list(
    cohorts,
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

    graphics::plot( outputs[[2]], ylab = "genetic values", type = "l",
          xaxt = "n", xlab = "", main = "Genetic gain")
    graphics::points( outputs[[2]])
    graphics::axis(1,
         at = seq_along(outputs[[2]]),
         labels = FALSE)

    graphics::text(x = seq_along(outputs[[2]]),
         y = graphics::par("usr")[3] - 0.05 * diff(graphics::par("usr")[3:4]),
         labels = names(outputs[[2]]),
         srt = 45, adj = 1, xpd = TRUE)


    graphics::plot( outputs[[3]], ylab = "inbreeding", type = "l", xaxt = "n", xlab = "",
          main = "inbreeding")
    graphics::points( outputs[[3]])

    graphics::axis(1,
         at = seq_along(outputs[[2]]),
         labels = FALSE)

    graphics::text(x = seq_along(outputs[[2]]),
         y = graphics::par("usr")[3] - 0.05 * diff(graphics::par("usr")[3:4]),
         labels = names(outputs[[2]]),
         srt = 45, adj = 1, xpd = TRUE)

    tmp =  outputs[[4]]
    tmp[is.na(tmp)] = 0
    graphics::plot(tmp, ylab = "prediction accuracy", type = "l",
          xaxt = "n", xlab = "", main = "prediction accuracy")
    graphics::points( tmp)

    graphics::axis(1,
         at = seq_along(outputs[[2]]),
         labels = FALSE)

    graphics::text(x = seq_along(outputs[[2]]),
         y = graphics::par("usr")[3] - 0.05 * diff(graphics::par("usr")[3:4]),
         labels = names(outputs[[2]]),
         srt = 45, adj = 1, xpd = TRUE)


  }


  return(outputs)
}




