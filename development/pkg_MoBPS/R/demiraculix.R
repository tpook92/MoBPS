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

#' Remove miraculix-coding for genotypes
#'
#' Internal function to decode all genotypes to non-miraculix objects
#' @param population Population list
#' @return Population list
#' @examples
#' # This is only relevant with the package miraculix is installed and used
#' population <- creating.diploid(nsnp=100, nindi=50)
#' \donttest{population <- demiraculix(population)}
#' @export

demiraculix <- function(population){

  if(!population$info$miraculix){
    warning("Miraculix is not active. Why do you use this function?")
    return(population)
  } else{
    population$info$miraculix <- FALSE
  }

  for(gen in 1:length(population$breeding)){
    for(sex in 1:2){
      if(length(population$breeding[[gen]][[sex]])){
        for(nr in 1:length(population$breeding[[gen]][[sex]])){

          if(length(population$breeding[[gen]][[sex]][[nr]][[9]])>0){
            genotype <- as.matrix(population$breeding[[gen]][[sex]][[nr]][[9]])
            population$breeding[[gen]][[sex]][[nr]][[9]] <- genotype[,1]
            population$breeding[[gen]][[sex]][[nr]][[10]] <- genotype[,2]
          }
        }
      }

    }
  }

  return(population)

}
