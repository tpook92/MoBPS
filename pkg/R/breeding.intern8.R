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

#' Internal function to simulate one meiosis
#'
#' Internal function to simulate one meiosis
#' @param info.parent position of the parent in the dataset
#' @param parent list of information regarding the parent
#' @param population Population list
#' @param mutation.rate Mutation rate in each marker (default: 10^-5)
#' @param remutation.rate Remutation rate in each marker (default: 10^-5)
#' @param recombination.rate Average number of recombination per 1 length unit (default: 1M)
#' @param recom.f.indicator Use step function for recombination map (transform snp.positions if possible instead)
#' @param duplication.rate Share of recombination points with a duplication (default: 0 - DEACTIVATED)
#' @param duplication.length Average length of a duplication (Exponentially distributed)
#' @param duplication.recombination Average number of recombinations per 1 length uit of duplication (default: 1)
#' @param gene.editing If TRUE perform gene editing on newly generated individual
#' @param gen.architecture Used underlying genetic architecture (genome length in M)
#' @param nr.edits Number of edits to perform per individual
#' @param decodeOriginsU Used function for the decoding of genetic origins [[5]]/[[6]]
#' @param delete.same.origin If TRUE delete recombination points when genetic origin of adjacent segments is the same
#' @param recombination.function Function used to calculate position of recombination events (default: MoBPS::recombination.function.haldane())
#' @examples
#' data(ex_pop)
#' child_gamete <- breeding.intern8(info.parent = c(1,1,1), parent = ex_pop$breeding[[1]][[1]][[1]],
#'                                 population = ex_pop)
#' @return Inherited parent gamete
#' @export
#'
breeding.intern8 <- function(info.parent, parent,  population , mutation.rate = 10^-5, remutation.rate = 10^-5, recombination.rate=1,
                            recom.f.indicator=NULL, duplication.rate=0, duplication.length=0.01,
                            duplication.recombination=1, delete.same.origin=FALSE,
                            gene.editing=FALSE, nr.edits= 0,
                            gen.architecture=0,
                            decodeOriginsU=MoBPS::decodeOriginsR,
                            recombination.function=MoBPS::recombination.function.haldane){

  return(list(c(0,1), NULL, 1, info.parent, NULL, NULL, 0.5, NULL, NULL))
}
