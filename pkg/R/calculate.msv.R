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

#' Set class
#'
#' Function to overwrite class values for individuals
#' @param population Population list
#' @param database Groups of individuals to consider for the export
#' @param gen Quick-insert for database (vector of all generations to export)
#' @param cohorts Quick-insert for database (vector of names of cohorts to export)
#' @param id Individual IDs to search/collect in the database
#' @param msv.path Path of the msv-software executable
#' @param traits Traits (integer value) for which to compute MSV scores
#' @param msv.snpeff snpeff parameter in msv software (default: 1)
#' @param msv.nrgen nrgen parameter in msv software (default: 2)
#' @param msv.selprop selprop parameter in msv software (default: 0.01)
#' @param msv.files Directory to write msv input files in
#' @param msv.type MSV criterion to use (default: "ExpBVSelGrOff", alt: "Index5", "Index6")
#' @param ebv.remove Set TRUE to subtract the EBV from the MSV output
#' @param msv.path.paternal_plink Path to a paternal plink file (instead of auto-generation)
#' @param msv.path.maternal_plink Path to a materal plink file (instead of auto-generation)
#' @param msv.path.ebv Path to a EBV file (instead of auto-generation)
#' @param msv.path.sol Path to a SOL file (instead of auto-generation)
#' @param msv.path.allfreq Path to a allele-freq file (instead of auto-generation)
#' @param msv.paternal_plink Set FALSE to not auto-generate a new input file
#' @param msv.maternal_plink Set FALSE to not auto-generate a new input file
#' @param msv.ebv Set FALSE to not auto-generate a new input file
#' @param msv.sol Set FALSE to not auto-generate a new input file
#' @param msv.allfreq Set FALSE to not auto-generate a new input file
#' @param verbose Set to FALSE to not display any prints
#' @param array Use only markers available on the array
#' @param remove.missing Remove markers not genotyped in any individual from the export
#' @param mixblup.last Use only SNPs included in the last MiXBLUP evaluation (that were not filtered through MAF etc.)
#' @examples
#' data(ex_pop)
#' "TO-BE_ADDED"
#' @return Population-List with newly entered class values
#' @export

calculate.msv <- function(population, database=NULL, gen=NULL, cohorts=NULL, id = NULL,
                       msv.path = "/lustre/shared/JV1/progs_gitlab/hpblup/ubuntu_20_04_2024_2_ifort_amd/msv",
                       traits = NULL,
                       msv.snpeff = 1,
                       msv.nrgen = 2,
                       msv.selprop = 0.01,
                       msv.files = "MSV_files",
                       msv.type = "ExpBVSelGrOff",
                       ebv.remove = FALSE,
                       msv.path.paternal_plink = NULL,
                       msv.path.maternal_plink = NULL,
                       msv.path.ebv = NULL,
                       msv.path.sol = NULL,
                       msv.path.allfreq = NULL,
                       msv.paternal_plink = TRUE,
                       msv.maternal_plink = TRUE,
                       msv.ebv = TRUE,
                       msv.sol = TRUE,
                       msv.allfreq = TRUE,
                       verbose = TRUE,
                       array = NULL,
                       remove.missing = TRUE,
                       mixblup.last = FALSE
                       ){



  if(mixblup.last && length(population$info$mixblup_array) == 1){

    if( length(array) == 0 ||   population$info$mixblup_array != array){
      array = population$info$mixblup_array
      if(verbose){
        cat(paste0("Used array automatically set to ", population$info$array.name[array], "\n"))
      }
    }

  }


  if(length(msv.path.paternal_plink)>0){
    msv.paternal_plink = FALSE
  }
  if(length(msv.path.maternal_plink)>0){
    msv.maternal_plink = FALSE
  }
  if(length(msv.path.ebv)>0){
    msv.ebv = FALSE
  }
  if(length(msv.path.sol)>0){
    msv.sol = FALSE
  }
  if(length(msv.path.allfreq)>0){
    msv.allfreq = FALSE
  }

  if(msv.paternal_plink){
    msv.path.paternal_plink = paste0(msv.files, if(msv.files=="") {} else{"/"}, "paternal_plink")
  }
  if(msv.maternal_plink){
    msv.path.maternal_plink = paste0(msv.files, if(msv.files=="") {} else{"/"}, "maternal_plink")
  }
  if(msv.ebv){
    msv.path.ebv = paste0(msv.files, if(msv.files=="") {} else{"/"}, "EBV.txt")
  }
  if(msv.sol){
    msv.path.sol = paste0(msv.files, if(msv.files=="") {} else{"/"}, "sol.dat")
  }
  if(msv.allfreq){
    msv.path.allfreq = paste0(msv.files, if(msv.files=="") {} else{"/"}, "all_freq")
  }


  if(length(traits) == 0){
    traits = 1:population$info$bv.nr
  }
  database <- get.database(population, gen, database, cohorts, id = id)



    # extract latest estimated u_hat
    alpha_hat = population$info$u_hat[[length(population$info$u_hat)]][,traits, drop = FALSE] # estimated SNP effects

    if(sum(is.na(alpha_hat))>0){
      alpha_hat[is.na(alpha_hat)] = 0
    }

    geno = get.geno(population = population, database = database) # genotypes
    p = rowMeans(geno)/2

    # frequency of estimated "good" allele
    freq_alphahat = p
    freq_alphahat[alpha_hat==0] = 0 # these can be ignored
    freq_alphahat[alpha_hat<0] = 1 -freq_alphahat[alpha_hat<0] # these can be ignored



  ebv = get.bve(population, database = database)[traits,,drop = FALSE] # estimated breeding value


  data = NULL

  if(length(msv.path)>0){

    if(sum(dir() == msv.files)==0){
      if(verbose){ cat(paste0("Create temporary folder for MSV files: ",   msv.files, "\n"))}
      dir.create(msv.files)
    }

    get.plink(population, database = database,
              path = msv.path.paternal_plink, fam.id = TRUE, type = 3, array = array, mixblup.last = mixblup.last, remove.missing = remove.missing)

    get.plink(population, database = database,
              path = msv.path.maternal_plink, fam.id = TRUE, type = 4, array = array, mixblup.last = mixblup.last, remove.missing = remove.missing)




    for(index in 1:ncol(alpha_hat)){

      tmp33 = alpha_hat[,index]

      if(mixblup.last){
        tmp33 = tmp33[population$info$mixblup_lastfreq[,1]]
      }

      utils::write.table(file = msv.path.ebv, cbind(get.id(population, database = database), ebv[index,]), row.names = FALSE, col.names = FALSE, quote = FALSE)
      utils::write.table(file = msv.path.sol, cbind(1,1,1,1:length(tmp33), tmp33), row.names = FALSE, col.names = FALSE, quote = FALSE)
      utils::write.table(file = msv.path.allfreq, cbind(1:length(tmp33), 0), row.names = FALSE, col.names = FALSE, quote = FALSE)

      msv_cmd = paste0(msv.path, " --matgeno ", msv.path.maternal_plink ,".bed --patgeno ", msv.path.paternal_plink,".bed --freq ", msv.path.allfreq," --ebv ", msv.path.ebv," --sol ", msv.path.sol, " --snpeff ", msv.snpeff," --nrgen ",msv.nrgen," --selprop ", msv.selprop )
      system(msv_cmd)

      if(msv.type == "ExpBVSelGrOff"){
        msv_data = utils::read.table("ExpBVSelGrOff_per_individual.txt")
      } else if(msv.type =="Index5"){
        msv_data = utils::read.table("Index_5_6.txt")[,1:2]
      } else if(msv.type =="Index6"){
        msv_data = utils::read.table("Index_5_6.txt")[,1:2]
      }

      if(ebv.remove){
        msv_data[,2] = msv_data[,2] - ebv[index,as.character(msv_data[,1])]
      }

      data = rbind(data, msv_data[,2])
    }

  } else{
    # potentially add tobias stuff


  }


  colnames(data) = msv_data[,1]
  for(index in 1:nrow(database)){

    ids = get.id(population, database = database[index,])

    msv_activ = data[,as.character(ids)]

    if(length(population$breeding[[database[index,1]]][[database[index,2] + 50]])==0){
      population$breeding[[database[index,1]]][[database[index,2] + 50]] = matrix(NA, nrow = population$info$bv.nr, ncol = population$info$size[database[index,1], database[index,2]])
    }

    population$breeding[[database[index,1]]][[database[index,2] + 50]][,database[index,3]: database[index,4]] = msv_activ

  }

  return(population)

}
