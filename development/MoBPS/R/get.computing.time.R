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

#' Derive ID on an individual
#'
#' Function to derive the internal ID given to each individual
#' @param population Population list
#' @param verbose Set to FALSE to not display any prints
#' @param extend Set to TRUE to return computing times with detailled overview on generation and BVE (default: FALSE)
#' @param per.call Set to TRUE to return computing times per call of breeding.diploid / creating.diploid() (default: FALSE)
#' @examples
#' data(ex_pop)
#' get.computing.time(ex_pop)
#' @return Computing times overview
#' @export

get.computing.time <- function(population, verbose = TRUE, extend = FALSE, per.call = FALSE){

  if(verbose && length(population$info$comp.times.general)>0){
    comp.times.general = round(colSums(population$info$comp.times.general), digits =1)
    comp.times.creating = round(colSums(population$info$comp.times.creating), digits =1)

    cat(paste0("Total time spent for generation: ", comp.times.general[7] + comp.times.creating[6], " seconds.\n\n"))
    cat("Time spent per step: \n")
    if((comp.times.creating[6])>0) cat(paste0((comp.times.creating[6]), " seconds for creation of founder population.\n"))
    if((comp.times.general[1])>0) cat(paste0((comp.times.general[1]), "seconds for initialization.\n"))
    if(comp.times.general[2]>0) cat(paste0( comp.times.general[2], " seconds for calculation of true genomic values.\n"))
    if(comp.times.general[3]>0) cat(paste0( comp.times.general[3], " seconds for phenotyping.\n"))
    if(comp.times.general[4]>0) cat(paste0( comp.times.general[4]," second for breeding value estimation.\n"))
    if(comp.times.general[5]>0) cat(paste0( comp.times.general[5], " seconds for selection.\n"))
    if(comp.times.general[6]>0) cat(paste0( comp.times.general[6], " seconds for generation of new individuals.\n"))
  }

  if(verbose && length(population$info$comp.times.bve)>0 && population$info$bve){
    comp.times.bve = round(colSums(population$info$comp.times.bve), digits =1)
    cat(paste0("Total time spent for BVE: ", comp.times.bve[10], " seconds.\n\n"))
    if(sum(comp.times.bve[1:3]>0)) cat("Time spent per step: \n")
    if((comp.times.bve[1])>0) cat(paste0((comp.times.bve[1]), " seconds for calculating Z, y and preparation.\n"))
    if((comp.times.bve[2])>0) cat(paste0((comp.times.bve[2]), " seconds for calculating A / G.\n"))
    if((comp.times.bve[3])>0) cat(paste0((comp.times.bve[3]), " seconds for solving mixed model.\n"))
  }

  if(verbose && length(population$info$comp.times.generation)>0){
    comp.times.generation = round(colSums(population$info$comp.times.generation), digits =1)
    cat(paste0("Total time spent for generation of new individuals: ", comp.times.bve[6], " seconds.\n\n"))
    if(sum(comp.times.generation[1:5]>0)) cat("Time spent per step: \n")
    if((comp.times.generation[1])>0) cat(paste0((comp.times.generation[1]), " seconds for preparation.\n"))
    if((comp.times.generation[2])>0) cat(paste0((comp.times.generation[2]), " seconds for generation (meiosis).\n"))
    if((comp.times.generation[3])>0) cat(paste0((comp.times.generation[3]), " seconds for generation (genomic value calculation).\n"))
    if(sum(comp.times.generation[4:5])>0) cat(paste0(sum(comp.times.generation[4:5]), " seconds for generation (multi-core).\n"))
  }

  if(per.call){
    comp.times.general = round((population$info$comp.times.general), digits =1)
    comp.times.bve = round((population$info$comp.times.bve), digits =1)
    comp.times.generation = round((population$info$comp.times.generation), digits =1)

  }
  if(extend){
    return(list(comp.times.general, comp.times.bve, comp.times.generation))
  } else{
    return(comp.times.general)
  }

}



