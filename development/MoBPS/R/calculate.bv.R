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

#' Calculate breeding values
#'
#' Internal function to calculate the breeding value of a given individual
#' @param population Population list
#' @param gen Generation of the individual of interest
#' @param sex Sex of the individual of interest
#' @param nr Number of the individual of interest
#' @param activ_bv traits to consider
#' @param decodeOriginsU Used function for the decoding of genetic origins [[5]]/[[6]]
#' @param import.position.calculation Function to calculate recombination point into adjacent/following SNP
#' @param store.effect.freq If TRUE store the allele frequency of effect markers per generation
#' @param bit.storing Set to TRUE if the MoBPS (not-miraculix! bit-storing is used)
#' @param nbits Bits available in MoBPS-bit-storing
#' @param output_compressed Set to TRUE to get a miraculix-compressed genotype/haplotype
#' @param bv.ignore.traits Vector of traits to ignore in the calculation of the genomic value (default: NULL; Only recommended for high number of traits and experienced users!)
#' @return [[1]] true genomic value [[2]] allele frequency at QTL markers
#' @examples
#' data(ex_pop)
#' calculate.bv(ex_pop, gen=1, sex=1, nr=1, activ_bv = 1)
#' @export

calculate.bv <- function(population, gen, sex, nr, activ_bv, import.position.calculation=NULL,
                         decodeOriginsU=decodeOriginsR, store.effect.freq=FALSE,
                         bit.storing=FALSE, nbits=30, output_compressed=FALSE,
                         bv.ignore.traits=NULL){



  bv_final <- numeric(length(activ_bv))
  snp.before <- c(0,population$info$cumsnp)
  cindex <- 1
  back <- 0
  first <- TRUE
  freq_list <- "notrelevant"

  if(length(bv.ignore.traits)!=length(activ_bv)){
    # Falls noetig koennten Haplotypen hier erst bestimmt werden.
    if(population$info$miraculix){
      if (requireNamespace("miraculix", quietly = TRUE)) {
        geno_self <- miraculix::computeSNPS(population, gen, sex, nr, what="geno")
      } else{
        stop("Usage of miraculix without being installed!")
      }
    } else{
      hap <- compute.snps(population, gen, sex, nr, import.position.calculation=import.position.calculation, decodeOriginsU=decodeOriginsU, bit.storing=bit.storing, nbits=nbits, output_compressed=output_compressed)
      geno_self <- colSums(hap)
    }

    if(sum(population$info$is.maternal)>0){
      index_mother <- population$breeding[[gen]][[sex]][[nr]][[8]]
      if(population$info$miraculix){
        if (requireNamespace("miraculix", quietly = TRUE)) {
          geno_mother <- miraculix::computeSNPS(population, index_mother[1], index_mother[2], index_mother[3], what="geno")
        } else{
          stop("Usage of miraculix without being installed!")
        }
      } else{
        hap <- compute.snps(population, index_mother[1], index_mother[2], index_mother[3],
                            import.position.calculation=import.position.calculation, decodeOriginsU=decodeOriginsU,
                            bit.storing=bit.storing, nbits=nbits, output_compressed=output_compressed)
        geno_mother <- colSums(hap)
      }
    }

    if(sum(population$info$is.paternal)>0){
      index_father <- population$breeding[[gen]][[sex]][[nr]][[7]]
      if(population$info$miraculix){
        if (requireNamespace("miraculix", quietly = TRUE)) {
          geno_father <- miraculix::computeSNPS(population, index_father[1], index_father[2], index_father[3], what="geno")
        } else{
          stop("Usage of miraculix without being installed!")
        }
      } else{
        hap <- compute.snps(population, index_father[1], index_father[2], index_father[3],
                            import.position.calculation=import.position.calculation, decodeOriginsU=decodeOriginsU,
                            bit.storing=bit.storing, nbits=nbits, output_compressed=output_compressed)
        geno_father <- colSums(hap)
      }
    }


    for(bven in activ_bv){


      if(length(intersect(bv.ignore.traits, bven))==1){
        back <- back +1
        cindex <- cindex+1
        next

        # NEED TO ADJUST effect.p.add.same or this!
      } else{
        back_old <- back
        back <- 0

      }
      if(population$info$is.maternal[bven]){
        geno <- geno_mother
      } else if(population$info$is.paternal[bven]){
        geno <- geno_father
      } else{
        geno <- geno_self
      }



      bv <- population$info$base.bv[bven] # Mittelwert

      real.bv.adds <- population$info$real.bv.add[[bven]]
      # Additive Effekte -Vektorielle loesung
      if(length(real.bv.adds)>0){

        recalc <- FALSE

        for(temp1 in 0:back_old){
          recalc <- recalc || (!population$info$effect.p.add.same[bven] || population$info$is.maternal[bven]  || population$info$is.paternal[bven] || (bven>1 && (population$info$is.maternal[bven-1]  || population$info$is.paternal[bven-1])))
        }

        if(first){
          recalc <- TRUE
        }
        first <- FALSE

        if(recalc){
          position <- population$info$effect.p.add[[bven]]
          neff <- nrow(real.bv.adds)
          take <- (geno[position] + 2L ) * neff + population$info$neff[[bven]]
        }


        ## ^ seems to be extremely inefficient!
        if(population$info$bve.poly.factor[bven]==1){
          bv <- bv + population$info$bve.mult.factor[bven] * sum((real.bv.adds[take]))
        } else{
          bv <- bv + population$info$bve.mult.factor[bven] * sum((real.bv.adds[take])^population$info$bve.poly.factor[bven])
        }

      }


      real.bv.mults <- population$info$real.bv.mult[[bven]]
      # Multiplikative Effekte - Vektorielle Version
      if(length(real.bv.mults)>0){
        position1 <- population$info$effect.p.mult1[[bven]]
        position2 <- population$info$effect.p.mult2[[bven]]
        bv <- bv + sum(real.bv.mults[cbind(1:nrow(real.bv.mults), 5L + geno[position1]*3L + geno[position2])])
      }

      # Wuerfel Effekte - Nicht vektorielle da Listenzugriffe (verbesserbar wenn viele solche Effekte betrachtet werden)
      real.bv.dices <- population$info$real.bv.dice[[bven]]
      if(length(real.bv.dices)>0 && length(real.bv.dices[[1]])>0){
        for(index in 1:length(real.bv.dices[[1]])){
          posis <- numeric(nrow(real.bv.dices[[1]][[index]]))
          for(index2 in 1:length(posis)){
            posis[index2] <- snp.before[real.bv.dices[[1]][[index]][index2,2]] + real.bv.dices[[1]][[index]][index2,1L]
          }
          activ_p <- 1L
          for(index2 in 1:length(posis)){
            activ_p <- activ_p + geno[posis[index2]] * (3L^(length(posis)-index2))
          }
          bv <- bv + real.bv.dices[[2]][[index]][activ_p]
        }

      }
      bv_final[cindex] <- bv
      cindex <- cindex + 1
    }

    if(store.effect.freq){
      freq_list <- cbind(geno[population$info$effect.p]==0, geno[population$info$effect.p]==1, geno[population$info$effect.p]==2)
      storage.mode(freq_list) <- "integer"
    } else{
      freq_list <- "notrelevant"
    }

  }

  return(list(bv_final, freq_list))

}
