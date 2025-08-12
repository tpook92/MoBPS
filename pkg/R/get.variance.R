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

#' Derive variances components (add/dom)
#'
#' Function to derive underlying variance components (add/dom)
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @return Table with realized narrow/broad-sense heritability, sigma_g,a,d
#' @examples
#' data(ex_pop)
#' get.variance(ex_pop, gen = 2)
#' @export
#'
get.variance = function(population, gen = NULL, database = NULL, cohorts = NULL){

  database = get.database(population, gen, database, cohorts)

  bv = get.bv(population, database = database)
  pheno = get.pheno(population, database = database)
  geno = get.geno(population, database = database)

  data_table = NULL

  for(trait in 1:nrow(bv)){

    sigma_g = stats::var(bv[trait,])
    sigma_e = stats::var(pheno[trait,] - bv[trait,])

    # calculation of narrow sense heritability:

    effects = get.qtl.effects(population)[[1]][[trait]]
    a = (effects[,5] - effects[,3]) /2
    d = effects[,4] - a - effects[,3]

    add_effect = a %*% geno[effects[,6],]
    dom_effect = d %*% (geno[effects[,6],] == 1)

    sigma_a = stats::var(as.numeric(add_effect))
    sigma_d = stats::var(as.numeric(dom_effect))

    h2_broad = sigma_g / (sigma_g + sigma_e)
    h2_narrow = sigma_a / (sigma_a + sigma_e)


    data_table = rbind(data_table, c(h2_narrow, h2_broad, sigma_g, sigma_e, sigma_a, sigma_d))
  }
  colnames(data_table) = c("h2_narrow", "h2_broad", "sigma_g2", "sigma_e2", "sigma_a2", "sigma_d2")

  rownames(data_table) = rownames(bv)

  return(data_table)
}
