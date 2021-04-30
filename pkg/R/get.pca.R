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

#' Principle components analysis
#'
#' Function to perform a principle component analysis
#' @param path Location were to save the PCA-plot
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param coloring Coloring by "group", "sex", "plain"
#' @param components Default: c(1,2) for the first two principle components
#' @param plot Set to FALSE to not generate a plot
#' @examples
#' data(ex_pop)
#' get.pca(ex_pop, gen=2)
#' @return Genotype data for in gen/database/cohorts selected individuals
#' @export

get.pca <- function(population, path=NULL,  database=NULL, gen=NULL, cohorts=NULL, coloring="group",
                    components = c(1,2), plot = TRUE){

  database <- get.database(population, gen, database, cohorts)

  geno <- get.geno(population, database = database)

  if(population$info$miraculix){
    A <- miraculix::relationshipMatrix(miraculix::genomicmatrix(geno),centered=TRUE, normalized=TRUE )
  } else{
    p_i <- rowSums(geno)/ncol(geno)/2
    Ztm <- geno - p_i * 2
    A <- crossprod(Ztm)/ (2 * sum(p_i*(1-p_i)))
  }

  b <- eigen(A)
  if(length(coloring)==1){
    if(coloring=="sex"){
      col <- rep("blue", ncol(A))
      col <- numeric(ncol(A))
      start <- 1
      for(index in 1:nrow(database)){
        add <- database[index,4] - database[index,3] +1
        if(database[index,2]==2){
          col[start:(start+add)] <- "red"
        } else{
          col[start:(start+add)] <- "blue"
        }

        start <- start + add
      }
    } else if(coloring=="group"){
      col <- numeric(ncol(A))
      start <- 1
      for(index in 1:nrow(database)){
        add <- database[index,4] - database[index,3] +1
        col[start:(start+add-1)] <- index
        start <- start + add
      }
    } else if(coloring=="class"){
      col <- numeric(ncol(A))
      start <- 1
      for(index in 1:nrow(database)){
        add <- database[index,4] - database[index,3] +1
        col[start:(start+add-1)] <- get.class(population, database = database[index,,drop=FALSE])
        start <- start + add
      }
    } else if(coloring=="genclass"){
      col <- numeric(ncol(A))
      start <- 1
      for(index in 1:nrow(database)){
        actives <- population$info$cohorts[which(population$info$cohorts[,2]==database[index,1]),]
        suppressWarnings(storage.mode(actives) <-  "integer")
        add <- database[index,4] - database[index,3] +1
        for(index2 in 1:nrow(actives)){
          if(database[index,2]==1){
            overwrite <- intersect(actives[index2,6]:(actives[index2,6] + actives[index2,3]), database[index,3]:database[index,4]) - database[index,3] +1
            if(length(overwrite)>0){
              col[start:(start+add-1)][overwrite] <- actives[index2,5]
            }
          }
          if(database[index,2]==2){
            overwrite <- intersect(actives[index2,7]:(actives[index2,7] + actives[index2,4]), database[index,3]:database[index,4]) - database[index,3] +1
            if(length(overwrite)>0){
              col[start:(start+add-1)][overwrite] <- actives[index2,5]
            }
          }
        }
        start <- start + add
      }

    } else{
      col <- rep(1, ncol(A))
    }
    col[col<0] <- 0
  } else{
    col <- rep(coloring, length.out = ncol(A))
  }


  if(plot){
    if(length(path)==0){
      graphics::plot(b$vectors[,components], col=col,
                     xlab=paste0("PC",components[1]," (",round(b$values[components[1]]/sum(b$values)*100, digits=2), "%)"),
                     ylab=paste0("PC", components[2]," (",round(b$values[components[2]]/sum(b$values)*100, digits=2), "%)"))
    } else{
      if (requireNamespace("grDevices", quietly = TRUE)) {
        grDevices::png(file=paste0(path, ".png"), width=2000, height= 1200, res=300)
        graphics::plot(b$vectors[,components], col=col,
                       xlab=paste0("PC",components[1]," (",round(b$values[components[1]]/sum(b$values)*100, digits=2), "%)"),
                       ylab=paste0("PC", components[2]," (",round(b$values[components[2]]/sum(b$values)*100, digits=2), "%)"))
        grDevices::dev.off()
      } else{
        stop("Use of grDevices without being installed!")
      }
    }

  }



  b$vectors
}
