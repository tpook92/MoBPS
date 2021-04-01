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

#' Admixture Plot
#'
#' Function to generate admixture plots
#' @param population Population list
#' @param database Groups of individuals to consider
#' @param gen Quick-insert for database (vector of all generations to consider)
#' @param cohorts Quick-insert for database (vector of names of cohorts to consider)
#' @param geno Manually provided genotype dataset to use instead of gen/database/cohorts
#' @param d dimensions to consider in admixture plot (default: automatically estimate a reasonable number)
#' @param verbose Set to FALSE to not display any prints
#' @param plot Set to FALSE to not generate an admixture plot
#' @param sort Set to TRUE to sort individuals according to contributes from the first dimension
#' @param sort.cutoff Skip individuals with contributions under this threshold (and use next dimension instead)
#' data(ex_pop)
#' get.admixture(ex_pop, gen=4:6)
#' @return Matrix with admixture proportion
#' @export


get.admixture <- function(population, geno=NULL, gen=NULL, database=NULL, cohorts= NULL, d=NULL, verbose=TRUE, plot=TRUE, sort=FALSE, sort.cutoff=0.01){

  if (requireNamespace("alstructure", quietly = TRUE)) {

    if(length(geno)==0){
      geno <- get.geno(population, gen=gen, database=database, cohorts=cohorts)
    }


    if(is.null(d)){
      d <- alstructure::estimate_d(geno) #unknown number of cluster
      if(verbose){
        cat(paste0(cat("Number of clusters are d =", d)))
      }
      if(d==1){
        stop("Estimated number of clusters is 1! Please choose dimension manually and/or check if you are including all intended individuals!")
      }
      fit <- alstructure::alstructure(geno, d = d) #estimating the ancestry coefficients
      ordered_factors_est <- alstructure::order_pops(fit$P, fit$Q_hat, method = "ave_admixture")
      Q_est <- ordered_factors_est$Q_ordered
    } else if (d >= 2){ #known number of cluster
      fit <- alstructure::alstructure(geno, d = d)
      ordered_factors_est <- alstructure::order_pops(fit$P, fit$Q_hat, method = "ave_admixture")
      Q_est <- ordered_factors_est$Q_ordered
    } else{
      stop("Incorrect input for dimensionality (d) of the Admixture plot!")
    }


  } else{
    stop("The R-package alstructure is required to generate admixture plots. \n
         Please install via github before use!\n
         install_github('storeylab/alstructure')")
  }


  if(plot){

    if(requireNamespace("RColorBrewer")){
      qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
      col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
      used_color <- sample(col_vector, nrow(Q_est))
    } else{
      used_color <- 1:nrow(Q_est)
      if(nrow(Q_est)>8){
        warning("At most 8 colors supported in admixture plot when not using RColorBrewer")
      }
    }

    if(sort){
      order <- NULL
      for(index in 1:d){
        if(length(order)< ncol(Q_est)){
          if(length(order)==0 ){
            remaining <- 1:ncol(Q_est)
          } else{
            remaining <- (1:ncol(Q_est))[-order]
          }

          sorting <- sort(Q_est[index, remaining], index.return=TRUE, decreasing=TRUE)
          order <- c(order, remaining[sorting$ix[sorting$x>sort.cutoff]])

        }

      }
    } else{
      order <- 1:ncol(Q_est)
    }

    a <- graphics::barplot(Q_est[,order], col=used_color, border=NA)

    graphics::axis(1, at=a, labels = colnames(geno)[order], las=2)
  }

  colnames(Q_est) <- colnames(geno)

  return(Q_est)
}
