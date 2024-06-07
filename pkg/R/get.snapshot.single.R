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

#' Derive snapshot of selected individuals
#'
#' Function to devide snapshot of genotyping/phenotyping state of selected individuals
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param phenotype.data Set to TRUE to include information of number of phenotypes generated
#' @param gain.data Set to TRUE to add information on changes in genetic level between cohorts (default: FALSE)
#' @param digits Number of digits provided for the gain.data output (default: 3)
#' @param time.diff Set to a target time interval to receive information between transitions to other cohorts (default: NA)
#' @param use.all.copy Set to TRUE to extract phenotyping
#' @param min.time Earliest time point relevant for output results (default: -Inf)
#' @param max.time Latest time point relevant for output results (default: Inf)
#' @param verbose Set to FALSE to not display any prints (default: TRUE)
#' @param time.points Use this parameter to manual provide a vector of time points on results should be generated for (default: NULL)
#' @param include.culled Set to TRUE to also include culled individuals in the statistics provided
#' @examples
#' data(ex_pop)
#' get.snapshot.single(ex_pop, cohorts = "Cohort_2_M")
#' @return Snapshot Matrix
#' @export

#save(file = "C:/Users/pook001/OneDrive - Wageningen University & Research/temp.RData", list = c("population"))
#cohorts = "E-line_piglets_5_M"
# get.snapshot(population, cohorts = "E-line_piglets_5_M")

get.snapshot.single = function(population, database=NULL, gen=NULL, cohorts=NULL, phenotype.data = FALSE, gain.data = FALSE,
                        digits = 3, time.diff = 1, min.time = -Inf, max.time = Inf, use.all.copy = TRUE, verbose = TRUE,
                        time.points = NULL, include.culled = FALSE){


  database <- get.database(population, gen, database, cohorts)

  ids = get.id(population, database = database)

  n_pheno = get.npheno(population, database = database, use.all.copy = use.all.copy)
  bv_temp = get.bv(population, database = database)

  ptime = get.pheno.time(population, database = database, use.all.copy = use.all.copy)
  gtime = get.geno.time(population, database = database, use.all.copy = use.all.copy)
  ctime = get.culling.time(population, database = database)


  # extract a matrix with information of all cohorts
  cohorts_list = get.cohorts(population, extended = TRUE)
  cohorts_ids = cohorts_list[,10:11]
  storage.mode(cohorts_ids) = "numeric"

  # only these cohorts are even candidate for including individuals from the database
  potential_cohorts = which((max(ids) > cohorts_ids[,1]) & (min(ids) < cohorts_ids[,2]))

  overlap1 = overlap2 = numeric(length(potential_cohorts))

  for(index in 1:length(potential_cohorts)){

    ids_potential  = get.id(population, cohorts = cohorts_list[potential_cohorts[index],1])
    to_analyse = which( ids %in% ids_potential)

    overlap1[index] = length(to_analyse)
    overlap2[index] = length(to_analyse) / length(ids_potential)

  }

  candidates = which(overlap1 == max(overlap1))

  candidates = candidates[which.max(overlap2[candidates])[1]]

  if(verbose){
    if(overlap2[candidates]==1){
      cat(paste0("Generate analysis starting from cohort: ", names(potential_cohorts)[candidates], "\n"))

    } else{
      cat(paste0("Generate analysis starting from subgroup of cohort: ", names(potential_cohorts)[candidates], "\n"))

    }
  }

  min_time = min(max.time, max(min.time, as.numeric(cohorts_list[potential_cohorts[candidates],8])))
  max_time = min(max.time, max(c(ptime, gtime, ctime, min_time, population$info$max.time.point), na.rm = TRUE))


  if(length(time.points)==0){
    time.points = seq(as.numeric(min_time), as.numeric(max_time), by = time.diff)
  }



  id_potential = get.id(population, cohorts = names(potential_cohorts)[candidates])

  to_analyse = which(ids %in% id_potential)
  to_analyse2 = which(id_potential %in% ids)

  if(gain.data){
    bv_base = rowMeans(get.bv(population, database = database)[,to_analyse2,drop = FALSE])
  }

  ids = ids[to_analyse]
  ptime = ptime[to_analyse]
  gtime = gtime[to_analyse]
  ctime = ctime[to_analyse]

  n_pheno = n_pheno[,to_analyse,drop = FALSE]
  bv_temp = bv_temp[,to_analyse, drop = FALSE]

  results = matrix(0, nrow = length(time.points), ncol = 5 + phenotype.data * population$info$bv.nr
                   + gain.data * population$info$bv.nr + 3 * include.culled)




  if(length(time.points)>0){
    for(index in 1:length(time.points)){


      results[index,1] = time.points[index]
      results[index,2] = paste0(names(potential_cohorts)[candidates])

      still_alive = ctime > time.points[index] | is.na(ctime)


      results[index,3] = sum(still_alive)
      results[index,4] = sum(gtime[still_alive] <= time.points[index], na.rm = TRUE)
      results[index,5] = sum(ptime[still_alive] <= time.points[index], na.rm = TRUE)

      if(phenotype.data){
        results[index, 1:population$info$bv.nr + 5 + 3 * include.culled] = rowSums(t(t(n_pheno[,still_alive,drop = FALSE]) * (ptime[still_alive] <= time.points[index])), na.rm = TRUE)
      }

      if(gain.data){
        results[index,1:population$info$bv.nr + 5 + 3 * include.culled + population$info$bv.nr * (phenotype.data)] = round(rowMeans(bv_temp[,still_alive,drop = FALSE]) - bv_base, digits = digits)

      }

      if(include.culled){
        results[index,6] = sum(!still_alive)
        results[index,7] = sum(gtime[!still_alive] <= time.points[index], na.rm = TRUE)
        results[index,8] = sum(ptime[!still_alive] <= time.points[index], na.rm = TRUE)
      }


    }
  }


  results = results[results[,3]!=0,,drop = FALSE]


  if(include.culled){
    names = c("time point", "cohort", "n_alive", "n_geno_alive", "n_pheno_alive", "n_culled", "n_geno_culled", "n_pheno_culled")
  } else{
    names = c("time point", "cohort", "n_indi", "n_geno", "n_pheno")
  }
  if(phenotype.data){
    names = c(names, paste0("Pheno_" , get.trait.name(population)))
  }
  if(gain.data){
    names = c(names, paste0("BV_diff_" , get.trait.name(population)))
  }
  colnames(results) = names

  # remove any cohort from the table with no individuals in it.


  return(results)
}
