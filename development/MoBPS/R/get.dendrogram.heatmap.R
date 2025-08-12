'#
  Authors
Torsten Pook, torsten.pook@wur.nl
Azadeh Hassanpour, azadeh.hassanpour@uni-goettingen.de

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

#' Dendrogram Heatmap
#'
#' Function calculate a dendogram heat
#' @param population Population list
#' @param path provide a path if the dendrogram would be saved as a png-file
#' @param database Groups of individuals to consider
#' @param gen Quick-insert for database (vector of all generations to consider)
#' @param cohorts Quick-insert for database (vector of names of cohorts to consider)
#' @param method Method used to calculate genetic distances (default: "Nei", alt: "Rogers", "Prevosti", "Modified Rogers"
#' @param individual.names Names of the individuals in the database ((default are MoBPS internal names based on position))
#' @param traits Traits to include in the dendrogram (default: all traits)
#' @param type Which traits values to consider (default: "pheno", alt: "bv", "bve")
#' @examples
#' population <- creating.diploid(nsnp=1000, nindi=40, n.additive = c(100,100,100),
#'                trait.cor = matrix(c(1,0.8,0.2,0.8,1,0.2,0.2,0.2,1), ncol=3), shuffle.traits = 1:3)
#' population <- breeding.diploid(population, phenotyping = "all", heritability = 0.5)
#' get.dendrogram.heatmap(population, gen=1, type="pheno")
#' @return Dendrogram plot of genotypes vs phenotypes
#' @export


get.dendrogram.heatmap <- function(population, path=NULL, database=NULL, gen=NULL, cohorts=NULL,
                                   method = NULL, individual.names = NULL, traits = NULL, type="pheno"){

  if (requireNamespace("NAM", quietly = TRUE)){
    database <- get.database(population, gen, database, cohorts)

    if(length(traits)==0){
      traits <- 1:population$info$bv.nr
    }


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
    }

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


    if (requireNamespace("gplots", quietly = TRUE)){
      mycol <- gplots::colorpanel(100, "orange", "blue", "white")
      if(length(path)==0){
        gplots::heatmap.2(x=y_scaled, Rowv=geno_dend, Colv=y_dend, col=mycol,
                                          scale="none", density.info="none",
                                          trace="none",margins=c(15,2), cexRow = 0.3, cexCol = 1.6,
                                          keysize = 1)
      } else{
        if (requireNamespace("grDevices", quietly = TRUE)) {
          grDevices::png(file=paste0(path, ".png"), width=2000, height= 1200, res=300)
          gplots::heatmap.2(x=y_scaled, Rowv=geno_dend, Colv=y_dend, col=mycol,
                            scale="none", density.info="none",
                            trace="none",margins=c(15,2), cexRow = 0.3, cexCol = 1.6,
                            keysize = 1)
          grDevices::dev.off()
        } else{
          warning("Use of grDevices without being installed!")
        }
      }



    } else{
      warning("Use of gplots without being installed!")
    }
  }


}



