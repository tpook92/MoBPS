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

#' Phylogenetic Tree
#'
#' Function calculate a phylogenetic tree
#' @param population Population list
#' @param path provide a path if the dendrogram would be saved as a png-file
#' @param database Groups of individuals to consider
#' @param gen Quick-insert for database (vector of all generations to consider)
#' @param cohorts Quick-insert for database (vector of names of cohorts to consider)
#' @param method Method used to calculate genetic distances (default: "Nei", alt: "Rogers", "Prevosti", "Modified Rogers"
#' @param individual.names Names of the individuals in the database ((default are MoBPS internal names based on position))
#' @param circular Set to TRUE to generate a fan/circular layout tree
#' get.phylogenetic.tree(ex_pop, gen=1, circular=TRUE)
#' @return Dendrogram plot for traits
#' @export


get.phylogenetic.tree <- function(population, path=NULL, database=NULL,
                           gen=NULL, cohorts=NULL, method=NULL, individual.names= NULL,
                           circular = FALSE){

  if (requireNamespace("NAM", quietly = TRUE)){
    database <- get.database(population, gen, database, cohorts)

    GD <- t(get.geno(population, database = database))

    if(length(individual.names)>0){
      rownames(GD) <- individual.names
    }

    if (requireNamespace("phylogram", quietly = TRUE)){
      if (is.null(method)){
        gdist = NAM::Gdist(GD, method = 1)
      }
      else if (method == "Rogers"){
        gdist = NAM::Gdist(GD, method = 4)
      }
      else if (method == "Prevosti"){
        gdist = NAM::Gdist(GD, method = 5)
      }
      else if (method == "Modified Rogers"){
        gdist = NAM::Gdist(GD, method = 6)
      }

      fit <- stats::hclust(gdist, method='ward.D')

      if(length(path)==0){
        if (!circular){
          graphics::plot(phylogram::as.phylo(fit),cex = 0.5,show.tip.label = T)

        } else {
          graphics::plot(phylogram::as.phylo(fit),type='fan',label.offset=0.1,no.margin=TRUE,cex=0.5,show.tip.label = T)
        }
      } else{
        if (requireNamespace("grDevices", quietly = TRUE)) {
          grDevices::png(file=paste0(path, ".png"), width=2000, height= 1200, res=300)
          if (!circular){
            graphics::plot(phylogram::as.phylo(fit),cex = 0.5,show.tip.label = T)

          } else {
            graphics::plot(phylogram::as.phylo(fit),type='fan',label.offset=0.1,no.margin=TRUE,cex=0.5,show.tip.label = T)
          }
          grDevices::dev.off()
        } else{
          stop("Use of grDevices without being installed!")
        }
      }


    } else{
      stop("Generation of phylogenetic trees requires the phylogram R-package!")
    }
  }

}



