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

#' Draw pedigree
#'
#' Draw a pedigree for selected individuals
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param depth.pedigree Depth of the pedigree in generations
#' @param storage.save The closer this is to 1 the more strict older animals will be filtered out of the pedigree (default: 1.1, min: 1)
#' @param use.id Set to TRUE to extract individual IDs
#' @param cex Size of individual labels
#' @param path NULL or a character value means whether the pedigree graph will be saved in a pdf file. The graph in the pdf file is a legible vector drawing, and labels isn't overlapped especially when the number of individuals is big and width of the individual label is long in one generation. It is recommended that saving a pedigree graph in the pdf file. The default value is NULL (this is taken from visPedigree documentation).
#' @param showgraph A logical value indicating whether a plot will be shown in the defaulted graphic device, such as the Plots panel of Rstudio. It is useful for quick viewing of the pedigree graph without opening the pdf file. However, the graph on the defaulted graphic device may be not legible, such as overlapped labels, aliasing lines due to the restricted width and height. It's a good choice to set showgraph = FALSE when the pedigree is large. The default value is TRUE (this is taken from visPedigree documentation).
#' @param outline A logical value indicating whether shapes without label will be shown. A graph of the pedigree without individuals' label is shown when setting outline = TRUE. It is very useful for viewing the outline of the pedigree and finding the immigrant individuals in each generation when the width of a pedigree graph is longer than the maximum width (200 inches) of the pdf file. The defaulted value is FALSE (this is taken from visPedigree documentation).
#' @param compact A logical value indicating whether IDs of full-sib individuals in one generation will be deleted and replaced with the number of full-sib individuals. For example, if there are 100 full-sib individuals in one generation, they will be deleted from the pedigree and be replaced with one individual label of "100" when compact = TRUE. The default value is FALSE (this is taken from visPedigree documentation).
#' @examples
#' population = creating.diploid(nsnp=100, nindi=10)
#' population = breeding.diploid(population, breeding.size=10)
#' population = breeding.diploid(population,  selection.m.database=cbind(c(1,2),1,1,5),
#' breeding.size=10)
#' population = breeding.diploid(population, selection.m.database=cbind(c(2,3),1,1,5),
#'  breeding.size=10)
#' get.pedigree.visual(population, gen=4)
#' @return Pedigree visualization
#' @export
#'


get.pedigree.visual = function(population, database=NULL, gen=NULL, cohorts=NULL, depth.pedigree = 3,
                    storage.save = 1.1, use.id = TRUE, cex = NULL, path=NULL,
                    showgraph = TRUE,
                    outline = FALSE, compact = FALSE){

  if (requireNamespace("visPedigree", quietly = TRUE)) {
    database <- get.database(population, gen, database, cohorts)

    if(depth.pedigree==Inf){
      pedigree.database <- get.database(population, gen=1:max(database[,1]))
    } else{
      new.pedigree.database <- pedigree.database <- database
      remaining.depth <- depth.pedigree
      while(remaining.depth>0){
        parents <- get.pedigree(population, database = new.pedigree.database, raw=TRUE)
        m_parents <- rbind(parents[parents[,5]==1,4:6], parents[parents[,8]==1,7:9])
        f_parents <- rbind(parents[parents[,5]==2,4:6], parents[parents[,8]==2,7:9])
        if(nrow(m_parents)>0){
          m_gen <- unique(m_parents[,1])
          m_data <- cbind(m_gen, 1, 0,0)
          nincluded <- numeric(length(m_gen))
          for(index in 1:length(m_gen)){
            m_data[index,3] <- min(m_parents[m_parents[,1]==m_gen[index],3])
            m_data[index,4] <- max(m_parents[m_parents[,1]==m_gen[index],3])
            nincluded[index] <- length(unique(m_parents[m_parents[,1]==m_gen[index],3]))
          }

          for(index in length(m_gen):1){
            if(nincluded[index] < (m_data[index,4]-m_data[index,3]+1)/storage.save){
              m_data <- m_data[-index,]
              activ_p <- unique(m_parents[m_parents[,1]==m_gen[index],3])
              m_data <- rbind(m_data, cbind(m_gen[index], 1, activ_p, activ_p))
            }
          }

        } else{
          m_data <- NULL
        }
        if(nrow(f_parents)>0){
          f_gen <- unique(f_parents[,1])
          f_data <- cbind(f_gen, 2, 0,0)
          nincluded <- numeric(length(f_gen))
          for(index in 1:length(f_gen)){
            f_data[index,3] <- min(f_parents[f_parents[,1]==f_gen[index],3])
            f_data[index,4] <- max(f_parents[f_parents[,1]==f_gen[index],3])
            nincluded[index] <- length(unique(f_parents[f_parents[,1]==f_gen[index],3]))
          }

          for(index in length(f_gen):1){
            if(nincluded[index] < (f_data[index,4]-f_data[index,3]+1)/storage.save){
              f_data <- f_data[-index,]
              activ_p <- unique(f_parents[f_parents[,1]==f_gen[index],3])
              f_data <- rbind(f_data, cbind(f_gen[index], 2, activ_p, activ_p))
            }
          }

        } else{
          f_data <- NULL
        }

        new.pedigree.database <- get.database(population, database=rbind(m_data,f_data))
        new.pedigree.database <- unique(new.pedigree.database)
        remaining.depth <- remaining.depth - 1
        pedigree.database <- rbind(new.pedigree.database, pedigree.database)
      }

      pedigree.database <- get.database(population, database = pedigree.database)
    }

    ped= get.pedigree(population, database = pedigree.database, use.id = use.id)

    pedt = suppressWarnings(visPedigree::tidyped(ped))
    visPedigree::visped(pedt, cex=cex, file=path, showgraph = showgraph,
                        outline = outline, compact = compact)
  } else{
    warning("R-package visPedigree needed to create visualization! \nPackage is not on CRAN and can be installed from GitHub: install_github('luansheng/visPedigree')")
  }


}
