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
#' Function calculate a dendogram
#' @param population Population list
#' @param path provide a path if the dendrogram would be saved as a png-file
#' @param database Groups of individuals to consider
#' @param gen Quick-insert for database (vector of all generations to consider)
#' @param cohorts Quick-insert for database (vector of names of cohorts to consider)
#' @param method Method used to calculate genetic distances (default: "Nei", alt: "Rogers", "Prevosti", "Modified Rogers"
#' @param individual.names Names of the individuals in the database ((default are MoBPS internal names based on position))
#' data(ex_pop)
#' get.dendrogram(ex_pop, gen=2)
#' @return Dendrogram plot for genotypes
#' @export


get.dendrogram <- function(population, path=NULL, database=NULL, gen=NULL, cohorts=NULL, method=NULL, individual.names = NULL){

  if (requireNamespace("NAM", quietly = TRUE)){
    database <- get.database(population, gen, database, cohorts)

    GD <- t(get.geno(population, database = database))

    if(length(individual.names)>0){
      rownames(GD) <- individual.names
    }

    if (is.null(method) || method =="Nei"){
      gdist = NAM::Gdist(GD, method = 1)
      geno_dend <- stats::as.dendrogram(stats::hclust(gdist,method='ward.D'))
      x_lab <- "Dendrogram based on Nei's genetic distances"
    }
    else if (method == "Rogers"){
      gdist = NAM::Gdist(GD, method = 4)
      geno_dend <- stats::as.dendrogram(stats::hclust(gdist, method='ward.D'))
      x_lab <- "Dendrogram based on Rogers' genetic distances"
    }
    else if (method == "Prevosti"){
      gdist = NAM::Gdist(GD, method = 5)
      geno_dend <- stats::as.dendrogram(stats::hclust(gdist, method='ward.D'))
      x_lab <- "Dendrogram based on Prevosti's genetic distances"
    }
    else if (method == "Modified Rogers"){
      gdist = NAM::Gdist(GD, method = 6)
      geno_dend <- stats::as.dendrogram(stats::hclust(gdist, method='ward.D'))
      x_lab <- "Dendrogram based on Modified Rogers' genetic distances"
    } else{
      stop("Invalid method! Chose between: 'Nei', 'Rogers', 'Prevosti', 'Modified Rogers'")
    }

    if(length(path)==0){
      graphics::plot(geno_dend, xlab = x_lab)
    } else{
      if (requireNamespace("grDevices", quietly = TRUE)) {
        grDevices::png(file=paste0(path, ".png"), width=2000, height= 1200, res=300)
        graphics::plot(geno_dend, xlab = x_lab)
        grDevices::dev.off()
      } else{
        stop("Use of grDevices without being installed!")
      }
    }

  }

  geno_dend
}


