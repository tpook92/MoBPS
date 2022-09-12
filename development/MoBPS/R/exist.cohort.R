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

#' Function to extract if a cohort exists
#'
#' Function to extract if a cohort exists
#' @param population Population list
#' @param cohort Cohort to check if it is contained in the population list
#' @examples
#' data(ex_pop)
#' exist.cohort(ex_pop, cohort = "StrangeName_42")
#' exist.cohort(ex_pop, cohort = "Cohort_1_M")
#' @return TRUE/FALSE
#' @export
#'
exist.cohort <- function(population, cohort){
  cohorts  <- get.cohorts(population)
  if(sum(cohorts==cohort)==1){
    return(TRUE)
  } else if(sum(cohorts==cohort)==0){
    return(FALSE)
  } else{
    stop("Cohort issue!")
  }
}
