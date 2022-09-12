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

#' Create a phenotypic transformation
#'
#' Function to perform create a transformation of phenotypes
#' @param population Population list
#' @param phenotypic.transform.function Phenotypic transformation to apply
#' @param trait Trait for which a transformation is to be applied
#' data(ex_pop)
#' trafo <- function(x){ return(x^2)}
#' ex_pop <- creating.phenotypic.transform(ex_pop, phenotypic.transform.function=trafo)
#' @return Population-list with a new phenotypic transformation function
#' @export

creating.phenotypic.transform <- function(population, phenotypic.transform.function=NULL, trait=1){

  if(length(phenotypic.transform.function)>0){
    population$info$phenotypic.transform[trait] <- TRUE
    population$info$phenotypic.transform.function[[trait]] <- phenotypic.transform.function
  }
  return(population)

}
