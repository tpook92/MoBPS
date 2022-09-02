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

#' Simulation of a given pedigree
#'
#' Function to simulate a given pedigree
#' @param pedigree Pedigree-file (matrix with 3 columns (Individual ID, Father ID, Mother ID), optional forth columns with earliest generations to generate an individual)
#' @param keep.ids Set to TRUE to keep the IDs from the pedigree-file instead of the default MoBPS ids
#' @param plot Set to FALSE to not generate an overview of inbreeding and number of individuals over time
#' @param dataset SNP dataset, use "random", "allhetero" "all0" when generating a dataset via nsnp,nindi
#' @param nsnp number of markers to generate in a random dataset
#' @param freq frequency of allele 1 when randomly generating a dataset
#' @param sex.s Specify which newly added individuals are male (1) or female (2)
#' @param add.chromosome If TRUE add an additional chromosome to the dataset
#' @param chromosome.length Length of the newly added chromosome (default: 5)
#' @param length.before Length before the first SNP of the dataset (default: 5)
#' @param length.behind Length after the last SNP of the dataset (default: 5)
#' @param snps.equidistant Use equidistant markers (computationally faster! ; default: TRUE)
#' @param snp.position Location of each marker on the genetic map
#' @param bv.total Number of traits (If more than traits via real.bv.X use traits with no directly underlying QTL)
#' @param polygenic.variance Genetic variance of traits with no underlying QTL
#' @param bit.storing Set to TRUE if the MoBPS (not-miraculix! bit-storing is used)
#' @param nbits Bits available in MoBPS-bit-storing
#' @param randomSeed Set random seed of the process
#' @param miraculix If TRUE use miraculix package for data storage, computations and dataset generation
#' @param miraculix.dataset Set FALSE to deactive miraculix package for dataset generation
#' @param n.additive Number of additive QTL
#' @param n.dominant Number of dominante QTL
#' @param n.qualitative Number of qualitative epistatic QTL
#' @param n.quantitative Number of quantitative epistatic QTL
#' @param var.additive.l Variance of additive QTL
#' @param var.dominant.l Variance of dominante QTL
#' @param var.qualitative.l Variance of qualitative epistatic QTL
#' @param var.quantitative.l Variance of quantitative epistatic QTL
#' @param exclude.snps Marker were no QTL are simulated on
#' @param replace.real.bv If TRUE delete the simulated traits added before
#' @param shuffle.traits Combine different traits into a joined trait
#' @param shuffle.cor Target Correlation between shuffeled traits
#' @param real.bv.add Single Marker effects
#' @param real.bv.mult Two Marker effects
#' @param real.bv.dice Multi-marker effects
#' @param bve.mult.factor Multiplicate trait value times this
#' @param bve.poly.factor Potency trait value over this
#' @param base.bv Average genetic value of a trait
#' @param add.chromosome.ends Add chromosome ends as recombination points
#' @param new.phenotype.correlation (OLD! - use new.residual.correlation) Correlation of the simulated enviromental variance
#' @param new.residual.correlation Correlation of the simulated enviromental variance
#' @param new.breeding.correlation Correlation of the simulated genetic variance (child share! heritage is not influenced!
#' @param add.architecture Add genetic architecture (marker positions)
#' @param position.scaling Manual scaling of snp.position
#' @param name.cohort Name of the newly added cohort
#' @param template.chip Import genetic map and chip from a species ("cattle", "chicken", "pig")
#' @param vcf Path to a vcf-file used as input genotypes (correct haplotype phase is assumed!)
#' @param vcf.maxsnp Maximum number of SNPs to include in the genotype file (default: Inf)
#' @param chr.nr Vector containing the assosiated chromosome for each marker (default: all on the same)
#' @param bp Vector containing the physical position (bp) for each marker (default: 1,2,3...)
#' @param bpcm.conversion Convert physical position (bp) into a cM position (default: 0 - not done)
#' @param snp.name Vector containing the name of each marker (default ChrXSNPY - XY chosen accordingly)
#' @param hom0 Vector containing the first allelic variant in each marker (default: 0)
#' @param hom1 Vector containing the second allelic variant in each marker (default: 1)
#' @param skip.rest Internal variable needed when adding multipe chromosomes jointly
#' @param beta.shape1 First parameter of the beta distribution for simulating allele frequencies
#' @param beta.shape2 Second parameter of the beta distribution for simulating allele frequencies
#' @param time.point Time point at which the new individuals are generated
#' @param creating.type Technique to generate new individuals (usage in web-based application)
#' @param trait.name Name of the trait generated
#' @param share.genotyped Share of individuals genotyped in the founders
#' @param genotyped.s Specify with newly added individuals are genotyped (1) or not (0)
#' @param map map-file that contains up to 5 colums (Chromsome, SNP-id, M-position, Bp-position, allele freq - Everything not provides it set to NA). A map can be imported via MoBPSmaps::ensembl.map()
#' @param remove.invalid.qtl Set to FALSE to deactive the automatic removal of QTLs on markers that do not exist
#' @param bv.standard Set TRUE to standardize trait mean and variance via bv.standardization() - automatically set to TRUE when mean/var.target are used
#' @param mean.target Target mean
#' @param var.target Target variance
#' @param verbose Set to FALSE to not display any prints
#' @param is.maternal Vector coding if a trait is caused by a maternal effect (Default: all FALSE)
#' @param is.paternal Vector coding if a trait is caused by a paternal effect (Default: all FALSE)
#' @param enter.bv Internal parameter
#' @examples
#' pedigree <- matrix(c(1,0,0,
#' 2,0,0,
#' 3,0,0,
#' 4,1,2,
#' 5,1,3,
#' 6,1,3,
#' 7,1,3,
#' 8,4,6,
#' 9,4,7), ncol=3, byrow=TRUE)
#' population <- pedigree.simulation(pedigree, nsnp=1000)
#' @return Population-list
#' @export



pedigree.simulation <- function(pedigree, keep.ids=FALSE, plot=TRUE,
                              dataset=NULL, vcf=NULL, chr.nr=NULL, bp=NULL, snp.name=NULL, hom0=NULL, hom1=NULL,
                                               bpcm.conversion=0,
                                               nsnp=0, freq="beta", sex.s="fixed",
                                               chromosome.length=NULL,length.before=5, length.behind=5,
                                               real.bv.add=NULL, real.bv.mult=NULL, real.bv.dice=NULL, snps.equidistant=NULL,
                                               bv.total=0, polygenic.variance=100,
                                               bve.mult.factor=NULL, bve.poly.factor=NULL,
                                               base.bv=NULL, add.chromosome.ends=TRUE,
                                               new.phenotype.correlation=NULL,
                                               new.residual.correlation = NULL,
                                               new.breeding.correlation=NULL,
                                               add.architecture=NULL, snp.position=NULL,
                                               position.scaling=FALSE,
                                               bit.storing=FALSE,
                                               nbits=30, randomSeed=NULL,
                                               miraculix=TRUE,
                                               miraculix.dataset=TRUE,
                                               n.additive=0,
                                               n.dominant=0,
                                               n.qualitative=0,
                                               n.quantitative=0,
                                               var.additive.l=NULL,
                                               var.dominant.l=NULL,
                                               var.qualitative.l=NULL,
                                               var.quantitative.l=NULL,
                                               exclude.snps=NULL,
                                               replace.real.bv=FALSE,
                                               shuffle.traits=NULL,
                                               shuffle.cor=NULL,
                                               skip.rest=FALSE,
                                               enter.bv=TRUE,
                                               name.cohort=NULL,
                                               template.chip=NULL,
                                               beta.shape1=1,
                                               beta.shape2=1,
                                               time.point=0,
                                               creating.type=0,
                                               trait.name=NULL,
                                               share.genotyped=1,
                                               genotyped.s=NULL,
                                               map=NULL,
                                               remove.invalid.qtl=TRUE,
                                               verbose=TRUE,
                                               bv.standard=FALSE,
                                               mean.target=NULL,
                                               var.target=NULL,
                                               is.maternal = NULL, is.paternal = NULL, vcf.maxsnp=Inf){

  while(ncol(pedigree)<5){
    pedigree <- cbind(pedigree, NA)
  }
  if(sum(is.na(pedigree[,4]))>0){
    pedigree[is.na(pedigree[,4]),4] <- 1
  }

  if(sum(is.na(pedigree[,5]))>0){
    pedigree[is.na(pedigree[,5]),5] <- 1
  }


  copies <- duplicated(pedigree[,1])
  if(sum(copies)>0){
    pedigree <- pedigree[!copies,]
  }
  pedigree[is.na(pedigree)] <- 0
  n <- nrow(pedigree)
  is_founder <- rep(FALSE, n)

  for(index in 1:n){
    if(sum(pedigree[index,2]==pedigree[,1])==0 || sum(pedigree[index,3]==pedigree[,1])==0){
      is_founder[index] <- TRUE
    }
  }

  avail <- which(is_founder)
  founder_sex <- pedigree[avail,5]

  if(verbose) cat(paste0(length(avail), " individuals were identified as founders.\n"))
  if(verbose) cat(paste0("of which ", sum(founder_sex==1), "/", sum(founder_sex==2), " are male / female.\n"))
  there <- rep(FALSE, n)
  gen <- rep(1, length(avail))

  temp1 <- 1

  fixed_breeding <- list()

  pedigree_position <- matrix(0, nrow=nrow(pedigree), ncol=3)
  pedigree_position[avail,] <- cbind(1,founder_sex, 1:length(avail))
  pedigree_position[avail[founder_sex==1],3] <- 1:sum(founder_sex==1)
  pedigree_position[avail[founder_sex==2],3] <- 1:sum(founder_sex==2)

  empty <- 0

  while(length(avail)<n){


    fixed_temp <- matrix(0, nrow=n, ncol=7)
    temp1 <- temp1 + 1
    there[avail] <- TRUE
    p_there <- pedigree[avail,1:3]
    p_there1 <- pedigree_position[avail,]

    poss <- rep(FALSE, n)
    index2 <- 1
    indexf <- 1
    indexm <- 1

    for(index in (1:n)[-avail]){
      a <- which(pedigree[index,2]==p_there[,1])
      b <- which(pedigree[index,3]==p_there[,1])
      if(length(a)==1 && length(b)==1 && pedigree[index,4]<=(temp1+empty)){
        poss[index] <- TRUE
        fixed_temp[index2,] <- c(p_there1[a,], p_there1[b,], pedigree[index,5]-1)
        if(pedigree[index,5]==1){
          pedigree_position[index,1:3] <- c(temp1, 1, indexm)
          indexm <- indexm + 1
        } else{
          pedigree_position[index,1:3] <- c(temp1, 2, indexf)
          indexf <- indexf +1
        }


        index2 <- index2 +1
      }
    }
    if(verbose) cat(paste0("Generation ",temp1,": ", index2-1, " individuals. (", sum(!there), " individuals remain)\n"))
    if(verbose) cat(paste0("of which ", indexm-1, "/", indexf-1, " are male / female.\n"))

    if(index2>1){
      fixed_breeding[[temp1]] <- fixed_temp[1:(index2-1),,drop=FALSE]
    } else{
     temp1 <- temp1 -1
     empty <- empty +1
    }

    avail <- c(avail, which(poss))
    gen <- c(gen, rep(temp1, sum(poss)))

  }

  population <- creating.diploid(nindi = sum(gen==1), sex.s=founder_sex,
                                 dataset=dataset, vcf=vcf, chr.nr=chr.nr, bp=bp, snp.name=snp.name,
                                 hom0=hom0, hom1=hom1, bpcm.conversion=bpcm.conversion,
                                 nsnp=nsnp, freq=freq,
                                 chromosome.length=chromosome.length,length.before=length.before,
                                 length.behind=length.behind,
                                 real.bv.add=real.bv.add, real.bv.mult=real.bv.mult,
                                 real.bv.dice=real.bv.dice, snps.equidistant=snps.equidistant,
                                 bv.total=bv.total, polygenic.variance=polygenic.variance,
                                 bve.mult.factor=bve.mult.factor, bve.poly.factor=bve.poly.factor,
                                 base.bv=base.bv, add.chromosome.ends=add.chromosome.ends,
                                 new.phenotype.correlation=new.phenotype.correlation,
                                 new.residual.correlation = new.residual.correlation,
                                 new.breeding.correlation=new.breeding.correlation,
                                 add.architecture=add.architecture, snp.position=snp.position,
                                 position.scaling=position.scaling,
                                 bit.storing=bit.storing,
                                 nbits=nbits, randomSeed=randomSeed,
                                 miraculix=miraculix,
                                 miraculix.dataset=miraculix.dataset,
                                 n.additive=n.additive,
                                 n.dominant=n.dominant,
                                 n.qualitative=n.qualitative,
                                 n.quantitative=n.quantitative,
                                 var.additive.l=var.additive.l,
                                 var.dominant.l=var.dominant.l,
                                 var.qualitative.l=var.qualitative.l,
                                 var.quantitative.l=var.quantitative.l,
                                 exclude.snps=exclude.snps,
                                 replace.real.bv=replace.real.bv,
                                 shuffle.traits=shuffle.traits,
                                 shuffle.cor=shuffle.cor,
                                 skip.rest=skip.rest,
                                 enter.bv=enter.bv,
                                 name.cohort=name.cohort,
                                 template.chip=template.chip,
                                 beta.shape1=beta.shape1,
                                 beta.shape2=beta.shape2,
                                 time.point=time.point,
                                 creating.type=creating.type,
                                 trait.name=trait.name,
                                 share.genotyped=share.genotyped,
                                 genotyped.s=genotyped.s,
                                 map=map,
                                 remove.invalid.qtl=remove.invalid.qtl,
                                 verbose=verbose,
                                 bv.standard=bv.standard,
                                 mean.target=mean.target,
                                 var.target=var.target,
                                 is.maternal = is.maternal,
                                 is.paternal = is.paternal,
                                 vcf.maxsnp=vcf.maxsnp)

  for(index in 2:length(fixed_breeding)){
    if(length(fixed_breeding[[index]])>0){
      population <- breeding.diploid(population, fixed.breeding = fixed_breeding[[index]], breeding.size = c(sum(fixed_breeding[[index]][,7]==0), sum(fixed_breeding[[index]][,7]==1)))
    }

  }


  if(plot){

    oldpar <- graphics::par(no.readonly = TRUE)
    graphics::par(mfrow=c(1,2))
    on.exit(graphics::par(oldpar))

    inbreeding <- numeric(length(fixed_breeding))
    for(index in 1:length(fixed_breeding)){
      inbreeding[index] <- kinship.emp.fast(population = population, gen=index, ibd.obs = 0, hbd.obs = population$info$size[index,1])[2]
    }
    plot(((inbreeding-0.5) *2), col="red", type="l", lwd=2, ylab="inbreeding", xlab="generation")
    plot(table(gen), xlab="generation", ylab="number of individuals")
  }

  if(keep.ids){

    for(index in 1:nrow(pedigree_position)){
      activ <- pedigree_position[index,]
      population$breeding[[activ[1]]][[activ[2]+14]][activ[3]] <- pedigree[index,1]
    }

    if(prod(is.numeric(pedigree[,1]))==1){
      population$info$next.animal <- max(pedigree[index,1])+1
    }


  }

  if(TRUE){

    return(list(population, avail, pedigree_position))
  } else{
    return(population)
  }


}
