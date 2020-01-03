'#
  Authors
Torsten Pook, torsten.pook@uni-goettingen.de

Copyright (C) 2017 -- 2018  Torsten Pook

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

#' Derive class
#'
#' Function to devide the class for each individual
#' @param path Location were to save the PCA-plot
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param coloring Coloring by "group", "sex", "plain"
#' @param components Default: c(1,2) for the first two principle components
#' @export

get.pca <- function(path=NULL, population, database=NULL, gen=NULL, cohorts=NULL, coloring="group",
                    components = c(1,2)){

  database <- get.database(population, gen, database, cohorts)

  geno <- get.geno(population, database = database)

  if(population$info$miraculix){
    A <- miraculix::relationshipMatrix(miraculix::genomicmatrix(geno),centered=TRUE, normalized=TRUE )
  } else{
    A <- 1
  }

  b <- eigen(A)
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
      col[start:(start+add)] <- index
      start <- start + add
    }
  } else{
    col <- rep(1, ncol(A))
  }


  if(length(path)==0){
    plot(b$vectors[,components], col=col,
         xlab=paste0("PC",components[1]," (",round(b$values[components[1]]/sum(b$values)*100, digits=2), "%)"),
         ylab=paste0("PC", components[2]," (",round(b$values[components[2]]/sum(b$values)*100, digits=2), "%)"))
  } else{
    png(file=paste0(path, ".png"), width=2000, height= 1200, res=300)
    par(mar=c(2.5,2.5,1.5,0.1))
    plot(b$vectors[,components], col=col,
         xlab=paste0("PC",components[1]," (",round(b$values[components[1]]/sum(b$values)*100, digits=2), "%)"),
         ylab=paste0("PC", components[2]," (",round(b$values[components[2]]/sum(b$values)*100, digits=2), "%)"))
    dev.off()
  }


  return(b$vectors)
}
