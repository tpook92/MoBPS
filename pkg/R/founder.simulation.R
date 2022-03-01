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

#' Founder simulation
#'
#' Function to generate founder genotypes
#' @param sex.quota Share of newly added female individuals (deterministic if sex.s="fixed", alt: sex.s="random")
#' @param n.gen Number of generations to simulate (default: 100)
#' @param dataset SNP dataset, use "random", "allhetero" "all0" when generating a dataset via nsnp,nindi
#' @param nsnp number of markers to generate in a random dataset
#' @param nindi number of inidividuals to generate in a random dataset
#' @param sex.quota Share of newly added female individuals (deterministic if sex.s="fixed", alt: sex.s="random")
#' @param nfinal Number of final individuals to include (default: nindi)
#' @param sex.quota.final Share of female individuals in the final generation
#' @param big.output Set to TRUE to export map, population list and pedigree relationship
#' @param depth.pedigree Depth of the pedigree in generations (default: 7)
#' @param plot Set to FALSE to not generate LD-decay plot and allele frequency spectrum
#' @param display.progress Set FALSE to not display progress bars. Setting verbose to FALSE will automatically deactive progress bars
#' @param freq frequency of allele 1 when randomly generating a dataset
#' @param sex.s Specify which newly added individuals are male (1) or female (2)
#' @param chromosome.length Length of the newly added chromosome (default: 5)
#' @param length.before Length before the first SNP of the dataset (default: 5)
#' @param length.behind Length after the last SNP of the dataset (default: 5)
#' @param snps.equidistant Use equidistant markers (computationally faster! ; default: TRUE)
#' @param snp.position Location of each marker on the genetic map
#' @param change.order If TRUE sort markers according to given marker positions
#' @param bit.storing Set to TRUE if the MoBPS (not-miraculix! bit-storing is used)
#' @param nbits Bits available in MoBPS-bit-storing
#' @param randomSeed Set random seed of the process
#' @param miraculix If TRUE use miraculix package for data storage, computations and dataset generation
#' @param miraculix.dataset Set FALSE to deactive miraculix package for dataset generation
#' @param template.chip Import genetic map and chip from a species ("cattle", "chicken", "pig")
#' @param position.scaling Manual scaling of snp.position
#' @param vcf Path to a vcf-file used as input genotypes (correct haplotype phase is assumed!)
#' @param vcf.maxsnp Maximum number of SNPs to include in the genotype file (default: Inf)
#' @param chr.nr Vector containing the assosiated chromosome for each marker (default: all on the same)
#' @param bp Vector containing the physical position (bp) for each marker (default: 1,2,3...)
#' @param bpcm.conversion Convert physical position (bp) into a cM position (default: 0 - not done)
#' @param snp.name Vector containing the name of each marker (default ChrXSNPY - XY chosen accordingly)
#' @param hom0 Vector containing the first allelic variant in each marker (default: 0)
#' @param hom1 Vector containing the second allelic variant in each marker (default: 1)
#' @param beta.shape1 First parameter of the beta distribution for simulating allele frequencies
#' @param beta.shape2 Second parameter of the beta distribution for simulating allele frequencies
#' @param map map-file that contains up to 5 colums (Chromsome, SNP-id, M-position, Bp-position, allele freq - Everything not provides it set to NA). A map can be imported via MoBPSmaps::ensembl.map()
#' @param verbose Set to FALSE to not display any prints
#' @examples
#' population <- founder.simulation(nindi=100, nsnp=1000, n.gen=5)
#' @export
#'

founder.simulation <- function(nindi=100, sex.quota=0.5, nsnp = 0, n.gen=100, nfinal=NULL, sex.quota.final=NULL, big.output = FALSE,
                               plot = TRUE, display.progress=TRUE, depth.pedigree = 7,
                               dataset=NULL, vcf=NULL, chr.nr=NULL, bp=NULL, snp.name=NULL, hom0=NULL, hom1=NULL,
                               bpcm.conversion=0,
                               freq="beta", sex.s="fixed",
                               chromosome.length=NULL,length.before=5, length.behind=5,
                               snps.equidistant=NULL,
                               change.order=FALSE,
                               snp.position=NULL,
                               position.scaling=FALSE,
                               bit.storing=FALSE,
                               nbits=30, randomSeed=NULL,
                               miraculix=TRUE,
                               miraculix.dataset=TRUE,
                               template.chip=NULL,
                               beta.shape1=1,
                               beta.shape2=1,
                               map=NULL,
                               verbose=TRUE,
                               vcf.maxsnp=Inf){


  if(length(nfinal)==0){
    nfinal <- nindi
  }
  if(length(sex.quota.final)==0){
    sex.quota.final <- sex.quota
  }

  if(length(map)==0 && nsnp==0){
    nsnp <- 10000
  } else if(length(map)>0){
    nsnp <- nrow(map)
  }


  n_new <- 50

  population <- creating.diploid(nindi = nindi, nsnp = nsnp, sex.quota = sex.quota,
                                 dataset=dataset, vcf=vcf, chr.nr=chr.nr, bp=bp,
                                 snp.name=snp.name, hom0= hom0, hom1=hom1,
                                 bpcm.conversion = bpcm.conversion,
                                 freq=freq, sex.s=sex.s,
                                 chromosome.length = chromosome.length,
                                 length.before = length.before, length.behind = length.behind,
                                 snps.equidistant = snps.equidistant,
                                 change.order = change.order,
                                 snp.position = snp.position,
                                 bit.storing = bit.storing,
                                 nbits = nbits, randomSeed = randomSeed,
                                 miraculix = miraculix,
                                 miraculix.dataset = miraculix.dataset,
                                 template.chip = template.chip,
                                 beta.shape1 = beta.shape1,
                                 beta.shape2 = beta.shape2,
                                 map = map,
                                 verbose=verbose,
                                 vcf.maxsnp = vcf.maxsnp)



  if(n.gen>1){
    if(verbose && display.progress){
      pb <- utils::txtProgressBar(min = 0, max = n.gen, style = 3)
    }

    for(index in 1: (n.gen-1)){

      if(verbose && display.progress){
        utils::setTxtProgressBar(pb, index)
      }

      population <- breeding.diploid(population, breeding.size = nindi, breeding.sex = sex.quota, verbose = verbose, display.progress=display.progress)

      if(index==3){
        size1 <- utils::object.size(population$breeding[[1]])
        size2 <- utils::object.size(population$breeding[[2]])
        size3 <- utils::object.size(population$breeding[[3]])
        growth <- size3 - size2
        n_new <- min(n_new, as.numeric(ceiling(sqrt(size1/growth) *1.5)))
      }

      if(index%%n_new==0){
        population <-  new.base.generation(population, base.gen=index+1)
      }

    }
  }

  if(verbose && display.progress){
    utils::setTxtProgressBar(pb, n.gen)
  }

  population <- breeding.diploid(population, breeding.size = nfinal, breeding.sex = sex.quota.final, verbose = verbose, display.progress=display.progress)

  if(verbose && display.progress){
    close(pb)
  }

  if(verbose) {cat("Start estimation of LD decay / effective population size.\n")}

  if(plot){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    graphics::par(mfrow=c(1,2))
  }

  geno <- get.geno(population, gen = nrow(population$info$size))


  ldinfo <- ld.decay(population, genotype.dataset = geno, type="cm", plot = plot)
  p_i <- rowMeans(geno)/2

  if(plot){

    tryCatch(  {
      graphics::hist(p_i[p_i!=0 & p_i !=1], nclass=20, xlab="allele frequency spectrum", main="Allele frequency spectrum")
    },
    error = function(e) {})


  }

  if(verbose){ cat(paste0(sum(p_i==0 | p_i==1), " of ", length(p_i), " markers are fixated.\n"))}

  effs <- numeric(length(ldinfo[[3]]$x))
  for(index in 1:length(ldinfo[[3]]$x)){
    effs[index] <- effective.size(ldinfo[[3]]$y[index],
                                  dist = ldinfo[[3]]$x[index],
                                  n = nfinal)
  }

  if(verbose) {cat(paste0("Effective population size is estimated to be around ", ceiling(mean(effs[-(1:10)])), ".\n"))}

  haplo <- get.haplo(population, gen = nrow(population$info$size))
  map <- get.map(population)

  storage.mode(haplo) <- "integer"

  if(big.output){
    pedigree <- kinship.exp(population, gen = nrow(population$info$size), depth.pedigree = depth.pedigree, verbose = verbose, mult=2)
    return(list(haplo, map, population, pedigree))
  } else{
    return(haplo)
  }

}


