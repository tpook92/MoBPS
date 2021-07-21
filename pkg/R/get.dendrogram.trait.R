'#
  Authors
Torsten Pook, torsten.pook@uni-goettingen.de
Azadeh Hassanpour, azadeh.hassanpour@uni-goettingen.de

Copyright (C) 2017 -- 2021  Torsten Pook

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

#' Dendrogram
#'
#' Function calculate a dendogram for the traits
#' @param population Population list
#' @param path provide a path if the dendrogram would be saved as a png-file
#' @param database Groups of individuals to consider
#' @param gen Quick-insert for database (vector of all generations to consider)
#' @param cohorts Quick-insert for database (vector of names of cohorts to consider)
#' @param traits Traits to include in the dendrogram (default: all traits)
#' @param type Which traits values to consider (default: "pheno", alt: "bv", "bve")
#' population <- creating.diploid(nsnp=1000, nindi=100, n.additive = c(100,100,100),
#'                shuffle.cor = matrix(c(1,0.8,0.2,0.8,1,0.2,0.2,0.2,1), ncol=3), shuffle.traits = 1:3)
#' population <- breeding.diploid(population, phenotyping = "all", heritability = 0.5)
#' get.dendrogram.trait(population, gen=1, type="pheno")
#' @return Dendrogram plot for traits
#' @export






get.dendrogram.trait <- function(population, path=NULL, database=NULL, gen=NULL, cohorts=NULL, traits = NULL, type="pheno"){

  if(length(traits)==0){
    traits <- 1:population$info$bv.nr
  }

  database <- get.database(population, gen, database, cohorts)

  if(type=="pheno"){
    y <- t(get.pheno(population, database = database))
  } else if(type=="bv"){
    y <- t(get.bv(population, database = database))
  } else{
    y <- t(get.bve(population, database = database))
  }

  if (length(traits) <=1){
    stop("The distance between the two phenotypic clusters needs at least two traits!")
  } else{
    y_scaled <- scale(y) ##in case we have different traits with different units
    y_dend <- stats::as.dendrogram(stats::hclust(stats::as.dist(1-stats::cor(y_scaled, method="spearman")), method="complete"))

    if(length(path)==0){
      graphics::plot(y_dend, xlab = "Dendrogram based on phenotypic distances")
    } else{
      if (requireNamespace("grDevices", quietly = TRUE)) {
        grDevices::png(file=paste0(path, ".png"), width=2000, height= 1200, res=300)
        graphics::plot(y_dend, xlab = "Dendrogram based on phenotypic distances")
        grDevices::dev.off()
      } else{
        stop("Use of grDevices without being installed!")
      }
    }
  }

  y_dend
}
