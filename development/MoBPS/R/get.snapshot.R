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
#' @examples
#' data(ex_pop)
#' get.snapshot(ex_pop, gen = 2)
#' @return Snapshot-matrix
#' @export

#save(file = "C:/Users/pook001/OneDrive - Wageningen University & Research/temp.RData", list = c("population"))
#cohorts = "E-line_piglets_5_M"
# get.snapshot(population, cohorts = "E-line_piglets_5_M")
# get.snapshot(ex_pop, gen = 2)
get.snapshot = function(population, database=NULL, gen=NULL, cohorts=NULL, phenotype.data = FALSE, gain.data = FALSE,
                        digits = 3, use.all.copy = TRUE, time.diff = NA){


  database <- get.database(population, gen, database, cohorts)


  ids = get.id(population, database = database)

  ptime = get.pheno.time(population, database = database, use.all.copy = use.all.copy)
  gtime = get.geno.time(population, database = database, use.all.copy = use.all.copy)
  ctime = get.culling.time(population, database = database)

  if(gain.data){
    bv_base = rowMeans(get.bv(population, database = database))
  }

  # extract a matrix with information of all cohorts
  cohorts_list = get.cohorts(population, extended = TRUE)

  cohorts_list = cohorts_list[sort(as.numeric(cohorts_list[,8]), index.return = TRUE)$ix,,drop = FALSE]
  cohorts_ids = cohorts_list[,10:11]
  storage.mode(cohorts_ids) = "numeric"

  # only these cohorts are even candidate for including individuals from the database
  potential_cohorts = which((max(ids) > cohorts_ids[,1]) & (min(ids) < cohorts_ids[,2]))

  results = matrix(0, nrow = length(potential_cohorts), ncol = 5 + phenotype.data * population$info$bv.nr
                   + gain.data * population$info$bv.nr)

  for(index in 1:length(potential_cohorts)){

    ids_potential  = get.id(population, cohorts = cohorts_list[potential_cohorts[index],1])


    results[index,1] = cohorts_list[potential_cohorts[index],8]
    results[index,2] = cohorts_list[potential_cohorts[index],1]

    to_analyse = which(ids %in% ids_potential)
    to_analyse2 = which( ids_potential %in% ids)
    results[index,3] = length(to_analyse)

    results[index,4] = sum(gtime[to_analyse] <= as.numeric(results[index,1] ), na.rm = TRUE)
    results[index,5] = sum(ptime[to_analyse] <= as.numeric(results[index,1]), na.rm = TRUE)


    if(phenotype.data){
      n_pheno = get.npheno(population, cohorts = cohorts_list[potential_cohorts[index],1], use.all.copy = use.all.copy)
      n_pheno[, which(ptime[to_analyse] > as.numeric(results[index,1]))] = 0
      results[index, 1:population$info$bv.nr + 5] = rowSums(n_pheno[,to_analyse2,drop = FALSE])
    }

    if(gain.data){

      bv_temp = get.bv(population, cohorts = cohorts_list[potential_cohorts[index],1])
      results[index,1:population$info$bv.nr + 5 + population$info$bv.nr * (phenotype.data)] = round(rowMeans(bv_temp[,to_analyse2,drop = FALSE]) - bv_base, digits = digits)

    }
  }

  results = results[results[,3]!=0,,drop = FALSE]

  if(!is.na(time.diff)){

    id_list = list()

    for(index in 1:nrow(results)){

      id_list[[index]] = get.id(population, cohorts = results[index,2])
    }

    next_appearance = rep(Inf, (nrow(results)))

    for(index in 1:(nrow(results)-1)){


      for(index2 in (index+1):nrow(results)){

        if(length(intersect(id_list[[index]], id_list[[index2]])) > 0){

          next_appearance[index] = as.numeric(results[index2, 1])
          break
        }

      }


    }

    results_extend = list()

    database_indi = get.database(population, database = database, per.individual = TRUE)

    for(index in 1:nrow(results)){

      ids_potential  = get.id(population, cohorts = results[index,2])
      to_analyse = which(ids %in% ids_potential)

      database_tmp = database_indi[to_analyse,,drop = FALSE]

      tmp = get.snapshot.single(population, database = database_tmp,
                                min.time = as.numeric(results[index,1]) + time.diff,
                                max.time = next_appearance[index],
                                time.diff = time.diff, verbose = FALSE,
                                  phenotype.data = phenotype.data,
                                gain.data = gain.data)

      tmp = tmp[tmp[,1] != next_appearance[index],,drop = FALSE]

      tmp[,2] = paste0(results[index,2], "_*time",tmp[,1] ,"*")
      results_extend[[index]] = t(tmp)


    }


    results = rbind(results, matrix(unlist(results_extend), ncol = ncol(results), byrow = TRUE))
    results = results[results[,3]!=0,,drop = FALSE]

  }
  # remove any cohort from the table with no individuals in it.

  results = results[sort(as.numeric(results[,1]), index.return = TRUE)$ix,,drop = FALSE]


  names = c("time point", "cohort", "n_indi", "n_geno", "n_pheno")
  if(phenotype.data){
    names = c(names, paste0("Pheno_" , get.trait.name(population)))
  }
  if(gain.data){
    names = c(names, paste0("BV_diff_" , get.trait.name(population)))
  }
  colnames(results) = names

  return(results)
}
