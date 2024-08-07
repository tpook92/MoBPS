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

#' Generation of the starting population
#'
#' Generation of the starting population
###### General
#' @param population Population list
#' @param nsnp Number of markers to generate (Split equally across chromosomes (chr.nr) unless vector is used)
#' @param nindi Number of individuals to generate (you can also provide number males / females in a vector)
#' @param name.cohort Name of the newly added cohort
#' @param generation Generation to which newly individuals are added (default: 1)
#' @param founder.pool Founder pool an individual is assign to (default: 1)
#' @param one.sex.mode Activating this will ignore all sex specific parameters and handle each individual as part of the first sex (default: FALSE)
#' @param sex.s Specify which newly added individuals are male (1) or female (2)
#' @param sex.quota Share of newly added female individuals (deterministic if sex.s="fixed", alt: sex.s="random")
#' @param class Migration level of the newly added individuals (default: 0)
#' @param verbose Set to FALSE to not display any prints
###### Genome architecture
#' @param map map-file that contains up to 5 colums (chromosome, SNP-id, M-position, Bp-position, allele freq - Everything not provides it set to NA). A map can be imported via MoBPSmaps::ensembl.map()
#' @param chromosome.length Length of the newly added chromosome (default: 5)
#' @param chr.nr Number of chromosomes (SNPs are equally split) or vector containing the associated chromosome for each marker
#' @param bp Vector containing the physical position (bp) for each marker (default: 1,2,3...)
#' @param snps.equidistant Use equidistant markers (computationally faster! ; default: TRUE)
#' @param template.chip Import genetic map and chip from a species ("cattle", "chicken", "pig")
#' @param snp.position Location of each marker on the genetic map
#' @param change.order Markers are automatically sorted according to their snp.position unless this is set to FALSE (default: TRUE)
#' @param add.chromosome If TRUE add an additional chromosome to the population
#' @param bpcm.conversion Convert physical position (bp) into a cM position (default: 0 - not done)
#' @param snp.name Vector containing the name of each marker (default ChrXSNPY - XY chosen accordingly)
#' @param hom0 Vector containing the first allelic variant in each marker (default: 0)
#' @param hom1 Vector containing the second allelic variant in each marker (default: 1)
###### Genotype data
#' @param dataset SNP dataset, use "random", "allhetero" "all0" when generating a dataset via nsnp,nindi
#' @param freq frequency of allele 1 when randomly generating a dataset (default: "beta" with parameters beta.shape1, beta.shape2; Use "same" when generating additional individuals and using the same allele frequencies)
#' @param beta.shape1 First parameter of the beta distribution for simulating allele frequencies
#' @param beta.shape2 Second parameter of the beta distribution for simulating allele frequencies
#' @param share.genotyped Share of individuals genotyped in the founders
#' @param genotyped.s Specify with newly added individuals are genotyped (1) or not (0)
#' @param vcf Path to a vcf-file used as input genotypes (correct haplotype phase is assumed!)
#' @param vcf.maxsnp Maximum number of SNPs to include in the genotype file (default: Inf)
#' @param vcf.maxindi Maximum number of individuals to include in the genotype file (default: Inf)
#' @param vcf.chromosomes Vector of chromosomes to import from vcf. Use on bgziped and tabixed vcf only. (default: NULL - all chromosomes)
#' @param vcf.VA Use the VariantAnnotation package to load in a vcf file when available (default: TRUE)
###### Traits
#' @param trait.name Name of the traits generated
#' @param mean.target Target mean for each trait
#' @param var.target Target variance for each trait
#' @param qtl.position.shared Set to TRUE to put QTL effects on the same markers for different traits
#' @param trait.cor Target correlation between QTL-based traits (underlying true genomic values)
#' @param trait.cor.include Vector of traits to be included in the modelling of corrlated traits (default: all - needs to match with trait.cor)
#' @param n.additive Number of additive QTL with effect size drawn from a gaussian distribution
#' @param n.dominant Number of dominant QTL with effect size drawn from a gaussian distribution
#' @param n.equal.additive Number of additive QTL with equal effect size (effect.size)
#' @param n.equal.dominant Number of n.equal.dominant QTL with equal effect size
#' @param n.qualitative Number of qualitative epistatic QTL
#' @param n.quantitative Number of quantitative epistatic QTL
#' @param effect.distribution Set to "gamma" for gamma distribution effects with gamma.shape1, gamma.shape2 instead of gaussian (default: "gauss")
#' @param gamma.shape1 Default: 1
#' @param gamma.shape2 Default: 1
#' @param real.bv.add Single Marker effects (list for each trait with columns for: SNP Nr, Chr Nr, Effect 00, Effect 01, Effect 11, Position (optional), Founder pool genotype (optional), Founder pool origin (optional))
#' @param real.bv.mult Two Marker effects
#' @param real.bv.dice Multi-marker effects
#' @param n.traits Number of traits (If more than traits via real.bv.X use traits with no directly underlying QTL)
#' @param base.bv Intercept of underlying true genomic values (excluding all QTL effects, default: 100)
#' @param new.residual.correlation Correlation of the simulated enviromental variance
#' @param new.breeding.correlation Correlation of the simulated genetic variance (child share! heritage is not influenced!
#' @param litter.effect.covariance Covariance matrix of the litter effect (default: no effects)
#' @param pen.effect.covariance Covariance matrix of the pen effect (default: no effects)
#' @param is.maternal Vector coding if a trait is caused by a maternal effect (Default: FALSE)
#' @param is.paternal Vector coding if a trait is caused by a paternal effect (Default: FALSE)
#' @param fixed.effects Matrix containing fixed effects (p x k -matrix with p being the number of traits and k being number of fixed effects; default: not fixed effects (NULL))
#' @param trait.pool Vector providing information for which pools QTLs of which trait are activ (default: 0 - all pools)
#' @param gxe.correlation Correlation matrix between locations / environments (default: only one location, sampled from gxe.max / gxe.min)
#' @param gxe.max Maximum correlation between locations / environments when generating correlation matrix via sampling (default: 0.85)
#' @param gxe.min Minimum correlation between locations / environments when generating correlation matrix via sampling (default: 0.70)
#' @param n.locations Number of locations / environments to consider for the GxE model
#' @param gxe.combine Set to FALSE to not view the same trait from different locations / environments as the sample trait in the prediction model (default: TRUE)
#' @param location.name Same of the different locations / environments used
#' @param dominant.only.positive Set to TRUE to always asign the heterozygous variant with the higher of the two homozygous effects (e.g. hybrid breeding); default: FALSE
#' @param exclude.snps Vector contain markers on which no QTL effects are placed
#' @param var.additive.l Variance of additive QTL
#' @param var.dominant.l Variance of dominante QTL
#' @param var.qualitative.l Variance of qualitative epistatic QTL
#' @param var.quantitative.l Variance of quantitative epistatic QTL
#' @param effect.size.equal.add Effect size of the QTLs in n.equal.additive
#' @param effect.size.equal.dom Effect size of the QTLs in n.equal.dominant
#' @param polygenic.variance Genetic variance of traits with no underlying QTL
#' @param bve.mult.factor Multiplicate trait value times this
#' @param bve.poly.factor Potency trait value over this
#' @param set.zero Set to TRUE to have no effect on the 0 genotype (or 00 for QTLs with 2 underlying SNPs)
#' @param bv.standard Set TRUE to standardize trait mean and variance via bv.standardization() - automatically set to TRUE when mean/var.target are used
#' @param replace.real.bv If TRUE delete the simulated traits added before
#' @param bv.ignore.traits Vector of traits to ignore in the calculation of the genomic value (default: NULL; Only recommended for high number of traits and experienced users!)
#' @param remove.invalid.qtl Set to FALSE to deactive the automatic removal of QTLs on markers that do not exist
###### Other
#' @param randomSeed Set random seed of the process
#' @param add.architecture Add genetic architecture (marker positions)
#' @param time.point Time point at which the new individuals are generated
#' @param creating.type Technique to generate new individuals (usage in web-based application)
#' @param size.scaling Set to value to scale all input for breeding.size / selection.size (This will not work for all breeding programs / less general than json.simulation)
#' @param progress.bar Set to FALSE to not use progress bars in any application of breeding.diploid() downstream (Keep log-files lean!)
###### Data storage
#' @param miraculix If TRUE use miraculix package for data storage, computations and dataset generation
#' @param miraculix.dataset Set FALSE to deactive miraculix package for dataset generation
#' @param add.chromosome.ends Add chromosome ends as recombination points
#' @param use.recalculate.manual Set to TRUE to use recalculate.manual to calculate genomic values (all individuals and traits jointly, default: FALSE)
#' @param store.comp.times Set to FALSE to not store computing times needed to execute creating.diploid in $info$comp.times.creating
###### Internal
#' @param skip.rest Internal variable needed when adding multipe chromosomes jointly
#' @param enter.bv Internal parameter
#' @param internal Dont touch!
#' @param internal.geno Dont touch!
#' @param internal.dataset Dont touch!
# OLD
#' @param nbits Bits available in MoBPS-bit-storing
#' @param bit.storing Set to TRUE if the MoBPS (not-miraculix! bit-storing is used)
#' @param new.phenotype.correlation (OLD! - use new.residual.correlation) Correlation of the simulated enviromental variance
#' @param length.before Length before the first SNP of the dataset (default: 5)
#' @param length.behind Length after the last SNP of the dataset (default: 5)
#' @param position.scaling Manual scaling of snp.position
#' @param shuffle.cor OLD! Use trait.cor - Target Correlation between traits
#' @param shuffle.traits OLD! Use trait.cor.include - Vector of traits to be included for modelling of correlated traits (default: all - needs to match with shuffle.cor)
#' @param bv.total OLD! Use n.traits instead. Number of traits (If more than traits via real.bv.X use traits with no directly underlying QTL)
#' @examples
#' population <- creating.diploid(nsnp=1000, nindi=100)
#' @return Population-list
#' @export


creating.diploid <- function(population=NULL,
                             #### General
                             nsnp=0,
                             nindi=0,
                             name.cohort=NULL,
                             generation=1,
                             founder.pool = 1,
                             one.sex.mode = FALSE,
                             sex.s="fixed",
                             sex.quota = 0.5,
                             class=0L,
                             verbose=TRUE,
                             #### Genome architecture
                             map=NULL,
                             chr.nr=NULL,
                             chromosome.length=NULL,
                             bp=NULL,
                             snps.equidistant=NULL,
                             template.chip=NULL,
                             snp.position=NULL,

                             change.order = TRUE,
                             add.chromosome=FALSE,
                             bpcm.conversion=0,
                             snp.name=NULL,
                             hom0=NULL,
                             hom1=NULL,
                             dataset=NULL,
                             freq="beta",
                             beta.shape1=1,
                             beta.shape2=1,
                             share.genotyped=0,
                             genotyped.s=NULL,
                             vcf=NULL,
                             vcf.maxsnp=Inf,
                             vcf.maxindi=Inf,
                             vcf.chromosomes = NULL,
                             vcf.VA = TRUE,
                             ## Traits
                             trait.name=NULL,
                             mean.target=NULL,
                             var.target=NULL,
                             qtl.position.shared = FALSE,
                             trait.cor = NULL,
                             trait.cor.include = NULL,

                             n.additive=0,
                             n.equal.additive=0,
                             n.dominant=0,
                             n.equal.dominant=0,
                             n.qualitative=0,
                             n.quantitative=0,
                             effect.distribution = "gauss",
                             gamma.shape1 = 1,
                             gamma.shape2 = 1,
                             real.bv.add=NULL,
                             real.bv.mult=NULL,
                             real.bv.dice=NULL,
                             new.residual.correlation = NULL,
                             new.breeding.correlation=NULL,
                             litter.effect.covariance = NULL,
                             pen.effect.covariance = NULL,
                             is.maternal = NULL,
                             is.paternal = NULL,
                             fixed.effects = NULL,
                             trait.pool = 0,
                             gxe.correlation = NULL,
                             n.locations = NULL,
                             gxe.max = 0.85,
                             gxe.min = 0.7,
                             location.name = NULL,
                             gxe.combine = TRUE,
                             n.traits = 0,
                             base.bv=NULL,
                             dominant.only.positive = FALSE,
                             exclude.snps=NULL,
                             var.additive.l=NULL,
                             var.dominant.l=NULL,
                             var.qualitative.l=NULL,
                             var.quantitative.l=NULL,
                             effect.size.equal.add = 1,
                             effect.size.equal.dom = 1,
                             polygenic.variance=100,
                             bve.mult.factor=NULL,
                             bve.poly.factor=NULL,
                             set.zero = FALSE,
                             bv.standard=FALSE,
                             replace.real.bv=FALSE,
                             bv.ignore.traits = NULL,
                             remove.invalid.qtl=TRUE,
                             #### Other
                             randomSeed=NULL,
                             add.architecture=NULL,
                             time.point=0,
                             creating.type=0,
                             size.scaling = 1,
                             progress.bar = TRUE,
                             #### Data storage
                             miraculix=TRUE,
                             miraculix.dataset=TRUE,
                             add.chromosome.ends=TRUE,
                             use.recalculate.manual = FALSE,
                             store.comp.times = TRUE,
                             #### Internal
                             skip.rest=FALSE,
                             enter.bv=TRUE,
                             internal=FALSE,
                             internal.geno=TRUE,
                             internal.dataset = NULL,
                             #### Old
                             nbits=30,
                             bit.storing=FALSE,
                             new.phenotype.correlation=NULL,
                             length.before=5,
                             length.behind=5,
                             position.scaling=FALSE,
                             shuffle.cor=NULL,
                             shuffle.traits=NULL,
                             bv.total=0
                             ){


  if(n.traits>0){
    if(bv.total>0){
      warning("bv.total has been overwritten with value from n.traits")
    }
    bv.total = n.traits
  }
  if(length(trait.cor)>0){
    shuffle.cor = trait.cor
  }
  if(length(trait.cor.include)>0){
    shuffle.traits = trait.cor.include
  }

  name.cohort_temp = name.cohort

  # GxE Trait generation module
  {
  if(length(n.locations)>0 && n.locations > 1 && length(gxe.correlation)==0){
    gxe.correlation = matrix(stats::runif(n.locations^2, gxe.min, gxe.max), ncol=n.locations)
    for(i in 2:nrow(gxe.correlation)) {
      for(j in 1:(i-1)) {
        gxe.correlation[i,j]=gxe.correlation[j,i]
      }
    }
    diag(gxe.correlation) = 1
    gxe.correlation = matrix.posdef(A = gxe.correlation)

    if(verbose){
      cat("Generated GxE matrix")
      print(round(gxe.correlation, digits = 3))
    }

  }

    if(length(gxe.correlation)>1 && length(n.locations)==0){
      n.locations = ncol(gxe.correlation)
    }

  if(length(location.name)==0 && length(gxe.correlation)>0){
    location.name = paste0("Location ", 1:ncol(gxe.correlation))
  }


    trait_location = NULL
    trait_nr = NULL

  if(length(gxe.correlation)>0){

    if(length(population)>0 && population$info$bv.nr >0){
      stop("GxE module is only intended for the use when no traits where previously generated")
    }

    if(length(real.bv.add)>0 || length(real.bv.mult)>0 || length(real.bv.dice)>0){
      stop("GxE module is only intended for the use with predefined MoBPS trait architectures")
    }

    # Determine total number of traits

    trait_sum <- n.additive + n.dominant + n.qualitative + n.quantitative + n.equal.additive + n.equal.dominant
    n_traits <- length(trait_sum)

    n.additive <- rep(c(n.additive, rep(0, length.out=n_traits-length(n.additive))), n.locations)
    n.dominant <- rep(c(n.dominant, rep(0, length.out=n_traits-length(n.dominant))), n.locations)
    n.equal.additive <- rep(c(n.equal.additive, rep(0, length.out=n_traits-length(n.equal.additive))), n.locations)
    n.equal.dominant <- rep(c(n.equal.dominant, rep(0, length.out=n_traits-length(n.equal.dominant))), n.locations)
    n.qualitative <- rep(c(n.qualitative, rep(0, length.out=n_traits-length(n.qualitative))), n.locations)
    n.quantitative <- rep(c(n.quantitative, rep(0, length.out=n_traits-length(n.quantitative))), n.locations)

    if(length(trait.name) < n_traits){
      trait.name = c(trait.name, paste0("Trait ", (length(trait.name)+1):n_traits))
    }

    # GxE will always result in a multi-trait model

    if(length(shuffle.cor)==0){
      shuffle.cor = diag(1, n_traits)
      shuffle.traits = 1:n_traits
    }
    n.locations = ncol(gxe.correlation)
    if(length(shuffle.cor)>0){

      shuffle.cor =   gxe.correlation  %x% shuffle.cor
      if(length(shuffle.traits)>0){
        shuffle.traits = rep(shuffle.traits, n.locations) + sort(rep(1:n.locations*n_traits - n_traits, length(shuffle.traits)))
      }
    }

    if(length(trait.name) < (n_traits * n.locations)){
      trait.name = paste0(rep(trait.name, n.locations) ," x ", rep(location.name, each = n_traits))
    }


    trait_location = rep(1:n.locations, each = n_traits)
    trait_nr =  rep(1:n_traits, n.locations)

    colnames(shuffle.cor) = rownames(shuffle.cor) = trait.name
    if(verbose && n.locations > 1 && n_traits > 1){
      cat("Used genetic correlation matrix:\n")
      print(shuffle.cor)
    }

  }
  }



  times_comp <- numeric(6)
  # Basic checks & map setup

  {



      counter.start <- counter <- c(1,1)



    if(store.comp.times){
      times_comp[1] <- as.numeric(Sys.time())
    }
    if(one.sex.mode){
      sex.quota <- 0
    }

    if(length(population)>0 && length(population$info$miraculix)>0 && !population$info$miraculix){
      miraculix = FALSE
    }
    if(length(nindi)==2){
      sex.quota <- nindi[2] / sum(nindi)
      nindi <- sum(nindi)
    }

    if(size.scaling!=1 & nindi>0){
      nindi = ceiling(nindi * size.scaling)
    }
    if(length(freq)==1 && freq=="same" && length(population$info$creating.freq)>0 ){
      freq <- population$info$creating.freq
    }
    if(length(randomSeed)>0){
      set.seed(randomSeed)
    }
    if(!miraculix && miraculix.dataset){
      if(verbose) cat("miraculix.dataset only possible when miraculix is active\n")
      miraculix.dataset <- FALSE
    }

    if(length(mean.target)>0){
      bv.standard <- TRUE
    } else{
      mean.target <- NA
    }
    if(length(var.target)>0){
      bv.standard <- TRUE
    } else{
      var.target <- NA
    }
    if(sum(set.zero)>0){
      bv.standard <- TRUE
    }

    if(length(new.phenotype.correlation)>0){
      new.residual.correlation <- new.phenotype.correlation
    }

    if(length(population)>0 & length(chr.nr)>0 & add.chromosome==FALSE){
      warning("chr.nr has automatically been set to NULL.\nThis parameter should only be used when modifiying the genome (e.g. adding chromosome / initial setup).")
      chr.nr <- NULL
    }

    preserve.bve <- length(population)==0

    if(length(map)>0){
      while(ncol(map)<5){
        map <- cbind(map, NA)
      }
      chr.nr <- map[,1]
      snp.name <- map[,2]
      chr.opt <- unique(chr.nr)

      if(sum(!is.na(map[,4]))>0){
        if(length(bp)==0){
          bp <- numeric(nrow(map))
        }
        bp[!is.na(map[,4])] <- as.numeric(map[!is.na(map[,4]),4])
      }
      if(sum(is.na(map[,3])) < nrow(map) && sum(map[,3]==0)==nrow(map)){
        if(bpcm.conversion==0){
          warning("0 Morgan is no legal position. Set position to NA")
        }

        map[map[,3]==0,3] <- NA
      } else if(sum(is.na(map[,3]))<nrow(map) && sum(map[,3]==0)>1){
        stop("0 Morgan is no legal position. Please fix!")
      }
      if(sum(!is.na(map[,3]))==nrow(map)){
        snp.position <- as.numeric(map[,3])
      } else if(sum(is.na(map[,3]))==nrow(map) && length(chromosome.length)==0){
        if(bpcm.conversion==0){
          if(verbose) cat("Assume 1 Morgan per 100.000.000bp - to change use bpcm.conversion\n")
          bpcm.conversion <- 1000000
        }
        map[,3] <- as.numeric(bp) /  bpcm.conversion / 100
        if(sum(!is.na(map[,3]))==nrow(map)){
          snp.position <- as.numeric(map[,3])
        }

      }
      if(sum(!is.na(map[,3]))==nrow(map) && length(chromosome.length)==0){

        chromosome.length <- numeric(length(chr.opt))
        for(index in 1:length(chr.opt)){
          chromosome.length[index] <- max(as.numeric(map[map[,1]==chr.opt[index],3])) + min(as.numeric(map[map[,1]==chr.opt[index],3]))
        }
      }
      if(sum(!is.na(map[,5]))>0){
        freq <- as.numeric(map[,5])
      }
      if(nsnp!=0 && sum(nsnp)!=nrow(map)){
        warning("Number of SNPs not in concordance with used map!\n")
        warning(paste0("Set number of SNPs to", nrow(map), "!\n"))

      }

      nsnp = numeric(length(chr.opt))
      for(index in 1:length(chr.opt)){
        nsnp[index] = sum(chr.nr==chr.opt[index])
      }
      #nsnp <- nrow(map)

    }
    if(length(chromosome.length)==0 || (length(chromosome.length)==1 && chromosome.length==0)){
      if(length(snp.position)>0){
        chromosome.length <- max(snp.position) + min(snp.position)
        if(chromosome.length<=0){
          warning("invalid setting for snp.position - Proceed with chromosome of 5M and equidistant markers")
          chromosome.length <- 5
          snps.equidistant <- TRUE
        }
      } else{
        chromosome.length <- 5
      }
    }
    if (requireNamespace("miraculix", quietly = TRUE)) {
      codeOriginsU <- miraculix::codeOrigins
      decodeOriginsU <- miraculix::decodeOrigins
      new_miraculix = as.numeric(substr(utils::packageVersion("miraculix"), start = 1, stop = 3))>=1.3

      if(!new_miraculix & length(population)==0){
        if(verbose){cat("Consider updating RandomFieldsUtils and miraculix for maximum performance!\n")}
      }
    } else{
      codeOriginsU <- codeOriginsR
      decodeOriginsU <- decodeOriginsR
      miraculix <- FALSE
      miraculix.dataset <- FALSE
    }

    if(add.chromosome == FALSE & length(population)>0){
      if(length(dataset)<=1 & nindi > 0){
        nsnp <- sum(population$info$snp)
      }
    }

    if(length(template.chip)==1){
      if(template.chip=="cattle"){
        target_snp <- nsnp
        chromosome.length <- MoBPS::cattle_chip[,2]
        nsnp <- round(MoBPS::cattle_chip[,3] * chromosome.length)
        if(target_snp>0){
          nsnp_temp  <- nsnp * target_snp / sum(nsnp)
          nsnp <- floor(nsnp_temp)
          if(sum(nsnp)<target_snp){
            add1 <- sort(nsnp_temp-nsnp, index.return=TRUE, decreasing = TRUE)$ix[1:(target_snp-sum(nsnp))]
            nsnp[add1] <- nsnp[add1] +1
          }
        }

        chr.nr <- numeric(sum(nsnp))
        start1 <- 1
        for(index in 1:length(nsnp)){
          if(nsnp[index]!=0){
            chr.nr[start1:(start1+nsnp[index]-1)] <- index
            start1 <- start1 + nsnp[index]
          }

        }
      } else if(template.chip=="chicken"){
        target_snp <- nsnp
        chromosome.length <- MoBPS::chicken_chip[,2]/100
        nsnp <- round(MoBPS::chicken_chip[,3] * chromosome.length)
        if(target_snp>0){
          nsnp_temp  <- nsnp * target_snp / sum(nsnp)
          nsnp <- floor(nsnp_temp)
          if(sum(nsnp)<target_snp){
            add1 <- sort(nsnp_temp-nsnp, index.return=TRUE, decreasing = TRUE)$ix[1:(target_snp-sum(nsnp))]
            nsnp[add1] <- nsnp[add1] +1
          }
        }
        chr.nr <- numeric(sum(nsnp))
        start1 <- 1
        for(index in 1:length(nsnp)){
          if(nsnp[index]!=0){
            chr.nr[start1:(start1+nsnp[index]-1)] <- index
            start1 <- start1 + nsnp[index]
          }
        }
      } else if(template.chip=="pig"){
        target_snp <- nsnp
        chromosome.length <- MoBPS::pig_chip[,2]/100
        nsnp <- round(MoBPS::pig_chip[,3] * chromosome.length)
        if(target_snp>0){
          nsnp_temp  <- nsnp * target_snp / sum(nsnp)
          nsnp <- floor(nsnp_temp)
          if(sum(nsnp)<target_snp){
            add1 <- sort(nsnp_temp-nsnp, index.return=TRUE, decreasing = TRUE)$ix[1:(target_snp-sum(nsnp))]
            nsnp[add1] <- nsnp[add1] +1
          }
        }
        chr.nr <- numeric(sum(nsnp))
        start1 <- 1
        for(index in 1:length(nsnp)){
          if(nsnp[index]!=0){
            chr.nr[start1:(start1+nsnp[index]-1)] <- index
            start1 <- start1 + nsnp[index]
          }
        }
      } else if(template.chip=="sheep"){
        target_snp <- nsnp
        chromosome.length <- MoBPS::sheep_chip[,2]
        nsnp <- round(MoBPS::sheep_chip[,3] * chromosome.length)
        if(target_snp>0){
          nsnp_temp  <- nsnp * target_snp / sum(nsnp)
          nsnp <- floor(nsnp_temp)
          if(sum(nsnp)<target_snp){
            add1 <- sort(nsnp_temp-nsnp, index.return=TRUE, decreasing = TRUE)$ix[1:(target_snp-sum(nsnp))]
            nsnp[add1] <- nsnp[add1] +1
          }
        }
        chr.nr <- numeric(sum(nsnp))
        start1 <- 1
        for(index in 1:length(nsnp)){
          if(nsnp[index]!=0){
            chr.nr[start1:(start1+nsnp[index]-1)] <- index
            start1 <- start1 + nsnp[index]
          }
        }
      } else if(template.chip=="maize"){
        target_snp <- nsnp
        chromosome.length <- MoBPS::maize_chip[,2]
        nsnp <- round(MoBPS::maize_chip[,3] * chromosome.length)
        if(target_snp>0){
          nsnp_temp  <- nsnp * target_snp / sum(nsnp)
          nsnp <- floor(nsnp_temp)
          if(sum(nsnp)<target_snp){
            add1 <- sort(nsnp_temp-nsnp, index.return=TRUE, decreasing = TRUE)$ix[1:(target_snp-sum(nsnp))]
            nsnp[add1] <- nsnp[add1] +1
          }
        }
        chr.nr <- numeric(sum(nsnp))
        start1 <- 1
        for(index in 1:length(nsnp)){
          if(nsnp[index]!=0){
            chr.nr[start1:(start1+nsnp[index]-1)] <- index
            start1 <- start1 + nsnp[index]
          }
        }
      }

    }

    if(length(dataset)>0 && sum(class(dataset) %in% "data.frame")>=1){
      dataset <- as.matrix(dataset)
    }

    if(length(dataset)>0 && sum(class(dataset) %in% "matrix")>=1){



      if(is.matrix(dataset) && is.numeric(dataset[1,1])){
        if(storage.mode(dataset)!="integer"){
          storage.mode(dataset) <- "integer"
        }
      } else if(is.matrix(dataset)){
        if(length(hom0)==0){
          hom0 <- dataset[,1]
        }
        if(length(hom1)==0){
          hom1 <- dataset[,1]
        }
        dataset <- (dataset==hom1)
        storage.mode(dataset) <- "integer"
      }


      diffs <- sum(diff(snp.position)<0)
      if(length(chr.nr)==1){
        comp <- chr.nr
      } else if (length(chr.nr)>1){
        comp <- length(unique(chr.nr))
      } else{
        comp <- 1
      }

      nsnp_temp <- nsnp
      if(sum(nsnp_temp)==0){
        nsnp_temp <- nrow(dataset)
      }
      if(diffs<comp & miraculix.dataset){
        if(length(nsnp)>1 && sum(nsnp)>0){
          dataset_temp <- dataset

          csnp <- cumsum(c(1, nsnp))
          dataset <- list()
          for(tt in 1:length(nsnp)){
            dataset[[tt]] <- miraculix::haplomatrix(dataset_temp[csnp[tt]:(csnp[tt+1]-1),,drop=FALSE])
          }
        } else if(length(chr.nr)==1 && nsnp_temp > 0){
          chr.nr <- sort(rep(1:chr.nr, length.out=nsnp_temp))
          dataset_temp <- dataset
          dataset <- list()
          for(tt in 1:max(chr.nr)){
            dataset[[tt]] <- miraculix::haplomatrix(dataset_temp[chr.nr==tt,,drop=FALSE])
          }
        } else if(length(chr.nr)==nrow(dataset)){
          dataset_temp <- dataset
          dataset <- list()
          for(tt in unique(chr.nr)){
            dataset[[tt]] <- miraculix::haplomatrix(dataset_temp[chr.nr==tt,,drop=FALSE])
          }
        } else{
          if(verbose & miraculix.dataset & length(chr.nr)>0){
            warning("Automatic miraculix sorting failed! Deactivate miraculix for dataset generation")
          }
          miraculix.dataset <- FALSE
        }
      } else{
        miraculix.dataset <- FALSE
      }
    }





    if(length(dataset)>0 && sum(class(dataset) %in% "haplomatrix")>=1){
      dataset <- list(dataset)
    }


    if(length(chr.nr)==1 && is.numeric(chr.nr)){

      if(length(nsnp)==1 && nsnp < chr.nr && nsnp != 1){
        stop("Each Chromosome must contain at least 1 SNP. Please check your inputs!")
      }
      if(length(nsnp)==1  && nsnp >= chr.nr){
        chr.nr <- sort(rep(1:chr.nr, length.out=nsnp))
        nsnp <- numeric(max(chr.nr))
        for(index in 1:length(nsnp)){
          nsnp[index] <- sum(chr.nr==index)
        }
      } else if(length(nsnp)>1 && length(nsnp)==chr.nr){
        chr.nr_temp = numeric(sum(nsnp))
        chr.nr_temp[1:nsnp[1]] = 1
        for(index in 2:length(nsnp)){
          chr.nr_temp[1:nsnp[index] + sum(nsnp[1:(index-1)])] = index
        }
        chr.nr = chr.nr_temp
      } else if(class(dataset)  %in% "matrix" && nrow(dataset)>chr.nr){
        chr.nr <- sort(rep(1:chr.nr, length.out=nrow(dataset)))
      } else if(class(dataset)  %in% "haplomatrix" && attr(dataset[[1]], "information")[2]==1){
        if(verbose) cat("Are you sure you want to generate a 1 SNP chromosome via miraculix?")
      }

    } else if(length(unique(chr.nr))>length(nsnp)){
      # in case map is provided.
      chr.unique <- unique(chr.nr)
      nsnp <- numeric(length(chr.unique))
      for(chr.c in 1:length(chr.unique)){
        nsnp[chr.c] <- sum(chr.nr==chr.unique[chr.c])
      }
    }
    if(length(chr.nr)==0 & sum(nsnp)>0){
      chr.nr <- numeric(sum(nsnp))
      for(sindex in 1:length(nsnp)){
        chr.nr[1:nsnp[sindex] + sum(nsnp[0:(sindex-1)])] <- sindex
      }
    }
  }


  # Genotype + Trait generation (only needed for chromosome 1)

  {
    if(store.comp.times){
      times_comp[2] <- as.numeric(Sys.time())
    }
    write_vcf_header <- FALSE
    if(skip.rest==FALSE){
      if(length(vcf)>0){

        exitVariantAnnotation <- 1
        if(requireNamespace("VariantAnnotation", quietly = TRUE) & vcf.VA){

          ### read (tabixed) vcf by VariantAnnotation  ---------------------------

          if(grepl("\\.gz$",vcf,perl = TRUE) && (file.exists(paste0(vcf,".tbi")) || file.exists(sub("\\.gz$",".tbi",vcf,perl = TRUE)))){
            tmp.fl <- Rsamtools::TabixFile(vcf)
          }else{
            tmp.fl <- vcf
          }
          tmp.header <- VariantAnnotation::scanVcfHeader(tmp.fl)
          write_vcf_header <- TRUE
          if(is.na(VariantAnnotation::geno(tmp.header)["GT","Description"])){
            warning("VariantAnnotation::readVcf() could not read GT field description from vcf header\n",
                    "- trying to import using vcfR::read.vcfR() without chromosome subset and header caption opportunity\n",
                    "- check your vcf header for single quotes instead of double quotes around GT FORMAT description and consider changing manually or using `bcftools reheader`")
          } else{

            #### subset by chromosome -------------------------------------------
            if(!is.null(vcf.chromosomes) && inherits(tmp.fl ,"TabixFile")){
              if(any(!vcf.chromosomes %in% rownames(VariantAnnotation::meta(tmp.header)$contig))){
                stop("Trying to subset vcf by chromosome, but some chromosome names were not found in the vcf header!\n",
                     "Potentionally, you don't have full contig info in your vcf header, try e.g. updating your sequence dictionary by \n",
                     "`bcftools reheader -f /path/to/reference.fa.fai /path/to/your/vcf.gz > new.vcf.gz && bcftools index -tf new.vcf.gz` !\n")
              }
              tmp.ranges <- GenomicRanges::GRanges(seqnames = vcf.chromosomes,
                                                   ranges = IRanges::IRanges(
                                                     start = rep(1,length(vcf.chromosomes)),
                                                     end = as.integer(VariantAnnotation::meta(tmp.header)$contig[vcf.chromosomes,])
                                                   ))
              tmp.params <- VariantAnnotation::ScanVcfParam(geno = c("GT"),
                                                            which = tmp.ranges)
            }else if(!is.null(vcf.chromosomes) && !inherits(tmp.fl ,"TabixFile")){
              stop("Trying to subset vcf by chromosome, but vcf does not seem to be bgzipped and tabix indexed!")
            }else{
              tmp.params <- VariantAnnotation::ScanVcfParam(geno = "GT") # geno = "GT"
            }

            ## read vcf
            vcf_file <- VariantAnnotation::readVcf(tmp.fl, param = tmp.params) ## needed later, or can be removed??????
            #vcf_file <- VariantAnnotation::readVcf(tmp.fl) ## needed later, or can be removed??????

            dataset <- matrix(0L, nrow = dim(vcf_file)[1], ncol = dim(vcf_file)[2]*2)
            suppressWarnings({
              dataset[,c(TRUE,FALSE)] <- as.integer(substr(VariantAnnotation::geno(vcf_file)$GT, start=1,stop=1))
              dataset[,c(FALSE,TRUE)] <- as.integer(substr(VariantAnnotation::geno(vcf_file)$GT, start=3,stop=3))
            })
            colnames(dataset) <- c(paste0(VariantAnnotation::samples(VariantAnnotation::header(vcf_file)),"_1"),
                                   paste0(VariantAnnotation::samples(VariantAnnotation::header(vcf_file)),"_2"))

            chr.nr <- as.character(MatrixGenerics::rowRanges(vcf_file)@seqnames)
            bp <- as.integer(MatrixGenerics::rowRanges(vcf_file)@ranges@start)
            snp.name <- names(MatrixGenerics::rowRanges(vcf_file))

            ## remove multivariate variants
            if(length(unlist(VariantAnnotation::fixed(vcf_file)$ALT)) > dim(vcf_file)[1] ){
              tmp.nalt <- VariantAnnotation::fixed(vcf_file)$ALT
              tmp.nalt <- which(lengths(tmp.nalt) > 1)

              warning('Currently only bivariate variants are supported. Removing ',length(tmp.nalt),' multivariate variants from vcf import!')
              dataset <- dataset[tmp.nalt * -1,]
              chr.nr <- chr.nr[tmp.nalt * -1]
              bp <- bp[tmp.nalt * -1]
              snp.name <- snp.name[tmp.nalt * -1]
              hom0 <- as.character(unlist(VariantAnnotation::ref(vcf_file)[tmp.nalt * -1]))
              hom1 <- as.character(unlist(VariantAnnotation::alt(vcf_file)[tmp.nalt * -1]))


              if(length(hom0) != length(bp)){
                hom0 <- as.character(unlist(VariantAnnotation::ref(vcf_file)[tmp.nalt * -1]))
              }
              if(length(hom1) != length(bp)){
                hom1 <- as.character(unlist(VariantAnnotation::alt(vcf_file)[tmp.nalt * -1]))
              }

            } else{
              hom0 <- as.character(unlist(VariantAnnotation::ref(vcf_file)))
              hom1 <- as.character(unlist(VariantAnnotation::alt(vcf_file)))

              if(length(hom0) != length(bp)){
                hom0 <- as.character((VariantAnnotation::ref(vcf_file)))
              }
              if(length(hom1) != length(bp)){
                hom1 <- as.character((VariantAnnotation::alt(vcf_file)))
              }
            }


            ## remove not needed objects to save space (header will be needed later!)
            #suppressWarnings(rm(vcf_file,tmp.nalt,tmp.params,tmp.ranges,tmp.fl))
            exitVariantAnnotation <- 0
          }
        }
        if(exitVariantAnnotation != 0 && requireNamespace("vcfR", quietly = TRUE)){
          vcf_file <- vcfR::read.vcfR(vcf, verbose = verbose)
          vcf_data <- vcf_file@gt[,-1]
          if(vcf.maxindi<ncol(vcf_data)){
            keep <- sort(sample(1:ncol(vcf_data), vcf.maxindi))
            keep <- 1:(vcf.maxindi) # just take the first few to not be dependent on randomness for now.
            vcf_data <- vcf_data[,keep]
          }
          dataset <- matrix(0L, nrow=nrow(vcf_data), ncol=ncol(vcf_data)*2)
          dataset[,(1:ncol(vcf_data))*2-1] <- as.integer(substr(vcf_data, start=1,stop=1))
          dataset[,(1:ncol(vcf_data))*2] <- as.integer(substr(vcf_data, start=3,stop=3))

          chr.nr <- vcf_file@fix[,1]
          bp <- as.numeric(vcf_file@fix[,2])
          snp.name <- vcf_file@fix[,3]
          hom0 <- vcf_file@fix[,4]
          hom1 <- vcf_file@fix[,5]
          #rm(vcf_file,vcf_data)

        } else if(exitVariantAnnotation != 0 ){

          vcf_file <- as.matrix(utils::read.table(vcf))
          vcf_data <- vcf_file[,-(1:9)]
          if(vcf.maxindi<ncol(vcf_data)){
            keep <- sort(sample(1:ncol(vcf_data), vcf.maxindi))
            keep <- 1:(vcf.maxindi) # just take the first few to not be dependent on randomness for now.
            vcf_data <- vcf_data[,keep]
          }
          dataset <- matrix(0L, nrow=nrow(vcf_data), ncol=ncol(vcf_data)*2)
          dataset[,(1:ncol(vcf_data))*2-1] <- as.integer(substr(vcf_data, start=1,stop=1))
          dataset[,(1:ncol(vcf_data))*2] <- as.integer(substr(vcf_data, start=3,stop=3))

          chr.nr <- vcf_file[,1]
          bp <- as.numeric(vcf_file[,2])
          snp.name <- vcf_file[,3]
          hom0 <- vcf_file[,4]
          hom1 <- vcf_file[,5]
        }


        hom0[hom0=="NA"] = NA
        hom1[hom1=="NA"] = NA

        if(change.order && length(snp.position) ==0 && sum(is.na(is.numeric(bp)))==0){
          chr_list = unique(chr.nr)
          factor_add = numeric(length(bp))
          for(index in 1:length(chr_list)){
            factor_add[chr.nr==chr_list[index]] = 10^10 * index
          }
          new_order = sort(as.numeric(bp) + factor_add, index.return=TRUE)$ix

          dataset <- dataset[new_order,,drop=FALSE]
          chr.nr <- chr.nr[new_order]
          bp <- bp[new_order]
          snp.name <- snp.name[new_order]
          hom0 <- hom0[new_order]
          hom1 <- hom1[new_order]
        }

        if(bpcm.conversion!=0){

          snp.position <- as.numeric(bp) /  bpcm.conversion / 100

        }



        if(length(population)>0){

          if(sum(is.na(population$info$snp.base))>0){
            ex_all = colSums(is.na(population$info$snp.base))
            new_all = is.na(hom0) + is.na(hom1)

            check = which(ex_all > new_all)

            replace1 = check[!is.na(population$info$snp.base[1,check]) & (hom0[check] == population$info$snp.base[1,check])]
            replace1 = replace1[!is.na(replace1)]
            population$info$snp.base[2,replace1] = hom1[replace1]


            replace2 = check[!is.na(population$info$snp.base[1,check]) & (hom1[check] == population$info$snp.base[1,check])]
            replace2 = replace2[!is.na(replace2)]
            population$info$snp.base[2,replace2] = hom0[replace2]


            if((length(replace1) + length(replace2)) >0 ){
              if(verbose) cat(paste0("Alternative alleles were added to ", length(replace1) + length(replace2), " SNPs.\n"))
            }

          }
          change = which((!(hom0 == population$info$snp.base[1,]) & (hom1 == population$info$snp.base[1,])) | (is.na(hom1) & (hom0 == population$info$snp.base[2,])))
          if(length(change)>0){
            if(verbose) cat(paste0("Reference allele is different in ", length(change), " SNPs. Automaticially fixed!\n"))
            dataset[change,] = 1L - dataset[change,]
          }
          if(sum(!(hom0 == population$info$snp.base[1,]))>length(change)){

            warning(paste0("Reference allele is different in ", sum(!(hom0 == population$info$snp.base[1,])) - length(change), " SNPs. No fix possible!"))
          }
          hom0 <- population$info$snp.base[1,]
          hom1 <- population$info$snp.base[2,]
        }

        if(length(population)>0){
          # this information is only required for the first input file!

          if(change.order && length(snp.position)>0 && length(unique(chr.nr))>1 && !add.chromosome){


            max_add = max(snp.position) + 1
            uni_chr = unique(chr.nr)
            for(index in 1:length(uni_chr)){
              snp.position[chr.nr==uni_chr[index]] = snp.position[chr.nr==uni_chr[index]] + max_add * (index-1)
              # this position is never used down-stream . this is just for sorting them
            }

          }
          chr.nr <- rep(1, length(chr.nr)) # everything is generated jointly!

          # Run checks if things match!
          if(FALSE){
            bp <- NULL
            snp.name <- NULL
            hom0 <- NULL
            hom1 <- NULL
          }

        }


        if(vcf.maxsnp< nrow(dataset)){
          keep <- sort(sample(1:nrow(dataset), vcf.maxsnp))
          dataset <- dataset[keep,,drop=FALSE]
          chr.nr <- chr.nr[keep]
          bp <- bp[keep]
          snp.name <- snp.name[keep]
          hom0 <- hom0[keep]
          hom1 <- hom1[keep]
        }


        if(length(population)==0){

          chr_tmp = unique(chr.nr)
          nsnp = length(chr_tmp)
          for(index in 1:length(chr_tmp)){
            nsnp[index] = sum(chr.nr == chr_tmp[index])
          }
        }


        if(length(dataset)>0 && sum(class(dataset) %in% "matrix")>=1){

          diffs <- sum(diff(bp)<0)
          if(length(chr.nr)==1){
            comp <- chr.nr
          } else if (length(chr.nr)>1){
            comp <- length(unique(chr.nr))
          } else{
            comp <- 1
          }

          nsnp_temp <- nsnp
          if(sum(nsnp_temp)==0){
            nsnp_temp <- nrow(dataset)
          }
          if(diffs<comp & miraculix.dataset){
            dataset_temp <- dataset
            dataset <- list()
            for(tt in unique(chr.nr)){
              dataset[[tt]] <- miraculix::haplomatrix(dataset_temp[chr.nr==tt,])
            }
          } else{
            if(verbose && miraculix.dataset && length(population)==0){
              warning("Automatic miraculix sorting failed! Deactivate miraculix for dataset generation")
            }
            miraculix.dataset <- FALSE


          }
        }

      }




      if(sum(nsnp)>0 && length(freq)==1 && freq=="beta"){
        freq <- stats::rbeta(sum(nsnp), shape1=beta.shape1, shape2=beta.shape2)
      }
      if(sum(nsnp)>0 && length(freq)<sum(nsnp)){
        freq <- rep(freq, length.out=sum(nsnp))
      }
      if(length(freq)>1 && length(freq)>sum(nsnp)){
        nsnp <- length(freq)
      }
      if(length(freq)>0 && sum(is.na(freq))>0){
        replace <- which(is.na(freq))
        freq[replace] <- stats::rbeta(sum(is.na(freq)), shape1=beta.shape1, shape2=beta.shape2)
        freq <- as.numeric(freq)
      }


      if(!is.list(var.additive.l) ){
        var.additive.l <- list(var.additive.l)
      }
      if(!is.list(var.dominant.l)){
        var.dominant.l <- list(var.dominant.l)
      }
      if(!is.list(var.qualitative.l)){
        var.qualitative.l <- list(var.qualitative.l)
      }
      if(!is.list(var.quantitative.l)){
        var.quantitative.l <- list(var.quantitative.l)
      }

      trait_sum <- n.additive + n.dominant + n.qualitative + n.quantitative + n.equal.additive + n.equal.dominant
      test <- list(NULL)

      if(length(var.additive.l) < length(trait_sum)){
        var.additive.l <- c(var.additive.l, rep(test,length.out=length(trait_sum)-length(var.additive.l)))
      }
      if(length(var.dominant.l) < length(trait_sum)){
        var.dominant.l <- c(var.dominant.l, rep(test,length.out=length(trait_sum)-length(var.dominant.l)))
      }
      if(length(var.qualitative.l) < length(trait_sum)){
        var.qualitative.l <- c(var.qualitative.l, rep(test,length.out=length(trait_sum)-length(var.qualitative.l)))
      }
      if(length(var.quantitative.l) < length(trait_sum)){
        var.quantitative.l <- c(var.quantitative.l, rep(test,length.out=length(trait_sum)-length(var.quantitative.l)))
      }

      ntraits <- length(trait_sum)
      n.additive <- c(n.additive, rep(0, length.out=ntraits-length(n.additive)))
      n.dominant <- c(n.dominant, rep(0, length.out=ntraits-length(n.dominant)))
      n.equal.additive <- c(n.equal.additive, rep(0, length.out=ntraits-length(n.equal.additive)))
      n.equal.dominant <- c(n.equal.dominant, rep(0, length.out=ntraits-length(n.equal.dominant)))
      n.qualitative <- c(n.qualitative, rep(0, length.out=ntraits-length(n.qualitative)))
      n.quantitative <- c(n.quantitative, rep(0, length.out=ntraits-length(n.quantitative)))

      if(length(unlist(c(var.qualitative.l, var.quantitative.l, var.additive.l, var.dominant.l)))>0){
        ntraits <- max(length(trait_sum), length(var.additive.l),length(var.dominant.l), length(var.qualitative.l), length(var.quantitative.l) )
        n.additive <- c(n.additive, rep(0, length.out=ntraits-length(n.additive)))
        n.dominant <- c(n.dominant, rep(0, length.out=ntraits-length(n.dominant)))
        n.equal.additive <- c(n.equal.additive, rep(0, length.out=ntraits-length(n.equal.additive)))
        n.equal.dominant <- c(n.equal.dominant, rep(0, length.out=ntraits-length(n.equal.dominant)))
        n.qualitative <- c(n.qualitative, rep(0, length.out=ntraits-length(n.qualitative)))
        n.quantitative <- c(n.quantitative, rep(0, length.out=ntraits-length(n.quantitative)))
        trait_sum <- n.additive + n.dominant + n.qualitative + n.quantitative + n.equal.additive + n.equal.dominant
        if(length(var.additive.l) < length(trait_sum)){
          var.additive.l <- rep(var.additive.l, length.out=length(trait_sum))
        }
        if(length(var.dominant.l) < length(trait_sum)){
          var.dominant.l <- rep(var.dominant.l, length.out=length(trait_sum))
        }
        if(length(var.qualitative.l) < length(trait_sum)){
          var.qualitative.l <- rep(var.qualitative.l, length.out=length(trait_sum))
        }
        if(length(var.quantitative.l) < length(trait_sum)){
          var.quantitative.l <- rep(var.quantitative.l, length.out=length(trait_sum))
        }
      }

      if(is.list(real.bv.dice) && length(unlist(real.bv.dice))==0){
        real.bv.dice = NULL
      }

      if(length(real.bv.dice)>0){
        if(is.list(real.bv.dice)){
          mdepth <- 0
          for(index in 1:length(real.bv.dice)){
            if(is.data.frame(real.bv.dice[[index]][[1]]) || is.matrix(real.bv.dice[[index]][[1]])){
              mdepth <- 1
            }
          }
          if(mdepth==0){
            for(index in 1:length(real.bv.dice)){
              if(length(real.bv.dice[[index]])>0){
                for(index2 in 1:length(real.bv.dice[[index]])){
                  if(is.data.frame(real.bv.dice[[index]][[1]][[index2]]) || is.matrix(real.bv.dice[[index]][[1]][[index2]]) ){
                    mdepth <- 2
                  }
                }
              }
            }
          }
        }
        if(mdepth == 1){
          real.bv.dice <- list(real.bv.dice)
        }
        if(mdepth == 0){
          stop("Illegal input for real.bv.dice")
        }
        for(index in 1:length(real.bv.dice)){
          if(length(real.bv.dice[[index]])>0){
            for(index2 in 1:length(real.bv.dice[[index]][[1]])){
              if(length(real.bv.dice[[index]][[2]][[index2]]) != (3 ^nrow(real.bv.dice[[index]][[1]][[index2]]))){
                warning("")
                stop("Length of effects does not match with involved effect SNPs - should be 3^(effect SNPs) (0..0, 0..01, ..., 2..2)")
              }
            }
          }
        }
      }

      if(length(real.bv.add)>0){
        if(!is.list(real.bv.add)){
          real.bv.add <- list(real.bv.add)
        }

        if(length(trait.pool)< max(length(real.bv.add),length(trait_sum))){
          trait.pool_temp = rep(trait.pool, length.out = length(real.bv.add))
        } else{
          trait.pool_temp = trait.pool
        }

        for(index3 in 1:length(real.bv.add)){

          if(ncol(real.bv.add[[index3]])==5){
            real.bv.add[[index3]] = cbind(real.bv.add[[index3]], NA, trait.pool_temp[index3])
          }

          if(ncol(real.bv.add[[index3]])==6){
            real.bv.add[[index3]] = cbind(real.bv.add[[index3]], trait.pool_temp[index3])
          }

          if(ncol(real.bv.add[[index3]])==7){
            real.bv.add[[index3]] = cbind(real.bv.add[[index3]], FALSE)
          }
        }
      }



      if(length(population)>0){
        if(length(real.bv.add)==0 && replace.real.bv==FALSE){
          real.bv.add <- population$info$real.bv.add
          real.bv.add[[population$info$bv.calc+1]] <- NULL
        } else if(replace.real.bv==FALSE){
          if(!is.list(real.bv.add)){
            real.bv.add <- list(real.bv.add)
          }
          real.bv.add <- c(population$info$real.bv.add, real.bv.add)
          real.bv.add[[population$info$bv.calc+1]] <- NULL
        }
        if(length(real.bv.mult)==0 && replace.real.bv==FALSE){
          real.bv.mult <- population$info$real.bv.mult
          real.bv.mult[[population$info$bv.calc+1]] <- NULL
        } else if(replace.real.bv==FALSE){
          if(!is.list(real.bv.mult)){
            real.bv.mult <- list(real.bv.mult)
          }
          real.bv.mult <- c(population$info$real.bv.mult, real.bv.mult)
          real.bv.mult[[population$info$bv.calc+1]] <- NULL
        }
        if(length(real.bv.dice)==0 && replace.real.bv==FALSE){
          real.bv.dice <- population$info$real.bv.dice
          real.bv.dice[[population$info$bv.calc+1]] <- NULL
        } else if(replace.real.bv==FALSE){
          if(!is.list(real.bv.dice)){
            real.bv.dice <- list(real.bv.dice)
          }
          real.bv.dice <- c(population$info$real.bv.dice, real.bv.dice)
          real.bv.dice[[population$info$bv.calc+1]] <- NULL
        }

      }
      if(length(real.bv.add)>0 && !is.list(real.bv.add)){
        real.bv.add <- list(real.bv.add)
      }
      if(length(real.bv.mult)>0 && !is.list(real.bv.mult)){
        real.bv.mult <- list(real.bv.mult)
      }
      if(length(real.bv.dice)>0 && !is.list(real.bv.dice)){
        real.bv.dice <- list(real.bv.dice)
      }

      {
        # Check for missingness in real.bvs and replace with reasonable inputs
        cum_snp <- cumsum(c(population$info$snp, nsnp))
        snpdata <- c(population$info$snp, nsnp)

        if(qtl.position.shared){

          n_qtls = max(n.additive + n.dominant + n.equal.additive+ n.equal.dominant+ n.quantitative*2+ n.qualitative*2)

          so_far <- max(length(real.bv.dice), length(real.bv.add), length(real.bv.mult))
          n_qtls_sofar = numeric(so_far)
          if(length(real.bv.add)>0){
            for(index in 1:length(real.bv.add)){
              n_qtls_sofar[index] = n_qtls_sofar[index] + nrow(real.bv.add[[index]])
            }
          }
          if(length(real.bv.mult)>0){
            for(index in 1:length(real.bv.mult)){
              n_qtls_sofar[index] = n_qtls_sofar[index] + nrow(real.bv.mult[[index]]) * 2
            }
          }

          effect_marker <- (1:sum(snpdata))
          if(length(exclude.snps)>0){
            effect_marker <- effect_marker[-exclude.snps]
          }

          max_qtls = max(c(n_qtls, n_qtls_sofar))
          effect_marker = sample(effect_marker, max_qtls)

        } else{
          so_far <- max(length(real.bv.dice), length(real.bv.add), length(real.bv.mult))
          effect_marker <- (1:sum(snpdata))
          if(length(exclude.snps)>0){
            effect_marker <- effect_marker[-exclude.snps]
          }
        }

        if(length(real.bv.add)>0){
          for(index in 1:length(real.bv.add)){
            while(sum(is.na(real.bv.add[[index]][,c(1:2,6)]))>0){

              add_marker <- sample(effect_marker, nrow(real.bv.add[[index]]), replace=if(nrow(real.bv.add[[index]])>length(effect_marker)){TRUE} else{FALSE})
              add_snp <- real.bv.add[[index]][,1]
              add_chromo <- real.bv.add[[index]][,2]

              for(index2 in (1:nrow(real.bv.add[[index]]))[is.na(add_snp) | is.na(add_chromo)]){
                add_chromo[index2] <- sum(add_marker[index2] > cum_snp) + 1
                add_snp[index2] <- add_marker[index2] - c(0,cum_snp)[add_chromo[index2]]
              }


              if(sum(is.na(real.bv.add[[index]][,6]))>0){
                add_marker = add_snp + c(0,cum_snp)[add_chromo]
              }

              enter <- add_chromo==real.bv.add[[index]][,2] | is.na(real.bv.add[[index]][,2])| is.na(real.bv.add[[index]][,1]) | is.na(real.bv.add[[index]][,6])

              real.bv.add[[index]][enter,c(1:2, 6)] <- cbind(add_snp, add_chromo, add_marker)[enter,]
            }
          }
        }

        if(length(real.bv.mult)>0){
          for(index in 1:length(real.bv.mult)){
            for(columns in c(0,2)){
              while(sum(is.na(real.bv.mult[[index]][,1:2+columns]))>0){

                add_marker <- sample(effect_marker, nrow(real.bv.mult[[index]]), replace=if(nrow(real.bv.mult[[index]])>length(effect_marker)){TRUE} else{FALSE})
                add_snp <- real.bv.mult[[index]][,1]
                add_chromo <- real.bv.mult[[index]][,2]

                for(index2 in (1:nrow(real.bv.mult[[index]]))[is.na(add_snp) | is.na(add_chromo)]){
                  add_chromo[index2] <- sum(add_marker[index2] > cum_snp) + 1
                  add_snp[index2] <- add_marker[index2] - c(0,cum_snp)[add_chromo[index2]]
                }

                enter <- add_chromo==real.bv.mult[[index]][,2+columns] | is.na(real.bv.mult[[index]][,2])

                real.bv.mult[[index]][enter,1:2+columns] <- cbind(add_snp, add_chromo)[enter,]
              }
            }

          }
        }
      }

      if(qtl.position.shared){
        sofar_positions = NULL

        if(length(real.bv.add)>0){
          for(index in 1:length(real.bv.add)){
            if(length(real.bv.add[[index]])>0){
              sofar_positions = c(sofar_positions, real.bv.add[[index]][,6])
            }
          }
        }

        sofar_positions = sofar_positions[!is.na(sofar_positions)]
        sofar_positions = setdiff(sofar_positions, effect_marker)
        if(length(sofar_positions)>0 && length(sofar_positions)<= length(effect_marker)){
          effect_marker[1:length(sofar_positions)] = sofar_positions
        }
      }


      if(length(dominant.only.positive)<length(trait_sum)){
        dominant.only.positive <- rep(dominant.only.positive, length.out = length(trait_sum))
      }
      so_far <- max(length(real.bv.dice), length(real.bv.add), length(real.bv.mult))

      if(length(trait.pool)< length(trait_sum)){
        trait.pool = rep(trait.pool, length.out = length(trait_sum))
      }

      if(length(trait_sum)>0){
        for(index_trait in 1:length(trait_sum)){
          var_additive <- var.additive.l[[index_trait]]
          var_dominant <- var.dominant.l[[index_trait]]
          var_qualitative <- var.qualitative.l[[index_trait]]
          var_quantitative <- var.quantitative.l[[index_trait]]
          if(n.additive[index_trait]>0 && length(var_additive)<n.additive[index_trait]){
            if(length(var_additive)==0){
              var_additive <- 1
            }
            var_additive <- rep(var_additive, length.out=n.additive[index_trait])
            var.additive.l[[index_trait]] <- var_additive
          }
          if(n.dominant[index_trait]>0 && length(var_dominant)<n.dominant[index_trait]){
            if(length(var_dominant)==0){
              var_dominant <- 1
            }
            var_dominant <- rep(var_dominant, length.out=n.dominant[index_trait])
            var.dominant.l[[index_trait]] <- var_dominant
          }
          if(n.qualitative[index_trait]>0 && length(var_qualitative)<n.qualitative[index_trait]){
            if(length(var_qualitative)==0){
              var_qualitative <- 1
            }
            var_qualitative <- rep(var_qualitative, length.out=n.qualitative[index_trait])
            var.qualitative.l[[index_trait]] <- var_qualitative
          }
          if(n.quantitative[index_trait]>0 && length(var_quantitative)<n.quantitative[index_trait]){
            if(length(var_quantitative)==0){
              var_quantitative <- 1
            }
            var_quantitative <- rep(var_quantitative, length.out=n.quantitative[index_trait])
            var.quantitative.l[[index_trait]] <- var_quantitative

          }

          if(length(var_additive)!= n.additive[index_trait]){
            n.additive[index_trait] <- length(var_additive)
          }
          if(length(var_dominant)!= n.dominant[index_trait]){
            n.dominant[index_trait] <- length(var_dominant)
          }
          if(length(var_qualitative)!= n.qualitative[index_trait]){
            n.qualitative[index_trait] <- length(var_qualitative)
          }
          if(length(var_quantitative)!= n.quantitative[index_trait]){
            n.quantitative[index_trait] <- length(var_quantitative)
          }






          #stop()
          #This part is only needed in creating.diploid
          if(sum(nsnp)>0){
            snpdata <- c(snpdata, nsnp)
          } else if((is.matrix(dataset) && nrow(dataset)>0 ) || is.list(dataset)){
            if(length(chr.nr)>0 && length(unique(chr.nr))>1){
              rindex <- 1
              for(chr.index in unique(chr.nr)){
                if(is.list(dataset)){

                  if(new_miraculix){
                    snpdata <- c(snpdata, attr(dataset[[rindex]], "information")[3])
                  } else{
                    snpdata <- c(snpdata, attr(dataset[[rindex]], "information")[2])
                  }

                  rindex <- rindex + 1
                } else{
                  snpdata <- c(snpdata, sum(chr.nr==chr.index))
                }
              }
            } else{
              if(sum(class(dataset)  %in% "haplomatrix")>0){
                if(new_miraculix){
                  snpdata <- c(snpdata, attr(dataset[[1]], "information")[3])
                } else{
                  snpdata <- c(snpdata, attr(dataset[[1]], "information")[2])
                }

              } else{
                snpdata <- c(snpdata, nrow(dataset))
              }

            }
          }

          # Generating additive

          add_marker <- sample(effect_marker, n.additive[index_trait], replace=if(n.additive[index_trait]>length(effect_marker)){TRUE} else{FALSE})
          dom_marker <- sample(effect_marker, n.dominant[index_trait], replace=if(n.dominant[index_trait]>length(effect_marker)){TRUE} else{FALSE})
          add_marker1 <- sample(effect_marker, n.equal.additive[index_trait], replace=if(n.equal.additive[index_trait]>length(effect_marker)){TRUE} else{FALSE})
          dom_marker1 <- sample(effect_marker, n.equal.dominant[index_trait], replace=if(n.equal.dominant[index_trait]>length(effect_marker)){TRUE} else{FALSE})
          epi1_marker <- sample(effect_marker, n.quantitative[index_trait]*2, replace=if(n.quantitative[index_trait]*2>length(effect_marker)){TRUE} else{FALSE})
          epi2_marker <- sample(effect_marker, n.qualitative[index_trait]*2, replace=if(n.qualitative[index_trait]*2>length(effect_marker)){TRUE} else{FALSE})



          cum_snp <- cumsum(snpdata)
          real.bv.add.new <- NULL
          real.bv.mult.new <- NULL
          if(n.additive[index_trait]>0){
            add_snp <- add_chromo <- numeric(n.additive[index_trait])
            for(index in 1:n.additive[index_trait]){
              add_chromo[index] <- sum(add_marker[index] > cum_snp) + 1
              add_snp[index] <- add_marker[index] - c(0,cum_snp)[add_chromo[index]]
            }
            if(effect.distribution == "gauss"){
              add_effect <- stats::rnorm(n.additive[index_trait], 0, var_additive)
            } else{
              add_effect <- stats::rgamma(n.additive[index_trait], gamma.shape1, gamma.shape2) * sample( c(-1,1), n.additive[index_trait], replace = TRUE)
            }

            real.bv.add.new <- cbind(add_snp, add_chromo, add_effect,0,-add_effect, add_marker, trait.pool[index_trait], FALSE)
          }

          if(n.equal.additive[index_trait]>0){
            add_snp1 <- add_chromo1 <- numeric(n.equal.additive[index_trait])
            for(index in 1:n.equal.additive[index_trait]){
              add_chromo1[index] <- sum(add_marker1[index] > cum_snp) + 1
              add_snp1[index] <- add_marker1[index] - c(0,cum_snp)[add_chromo1[index]]
            }
            add_effect1 <- effect.size.equal.add
            real.bv.add.new <- rbind(real.bv.add.new, cbind(add_snp1, add_chromo1,  -add_effect1, 0, add_effect1, add_marker1, trait.pool[index_trait], FALSE))

          }

          if(n.dominant[index_trait]>0){
            dom_snp <- dom_chromo <- numeric(n.dominant[index_trait])
            for(index in 1:n.dominant[index_trait]){
              dom_chromo[index] <- sum(dom_marker[index] > cum_snp) + 1
              dom_snp[index] <- dom_marker[index] - c(0,cum_snp)[dom_chromo[index]]
            }

            if(effect.distribution == "gauss"){
              dom_effect <- stats::rnorm(n.dominant[index_trait], 0, var_dominant)
            } else{
              dom_effect <- stats::rgamma(n.dominant[index_trait], gamma.shape1, gamma.shape2) * sample( c(-1,1), n.dominant[index_trait], replace = TRUE)
            }

            if(dominant.only.positive[index_trait]){
              temp1 <- dom_effect
              temp1[temp1<0] <- 0
            } else{
              temp1 <- dom_effect
            }
            real.bv.add.new <- rbind(real.bv.add.new, cbind(dom_snp, dom_chromo, 0 ,temp1,dom_effect, dom_marker, trait.pool[index_trait], FALSE))

          }

          if(n.equal.dominant[index_trait]>0){
            dom_snp1 <- dom_chromo1 <- numeric(n.equal.dominant[index_trait])
            for(index in 1:n.equal.dominant[index_trait]){
              dom_chromo1[index] <- sum(dom_marker1[index] > cum_snp) + 1
              dom_snp1[index] <- dom_marker1[index] - c(0,cum_snp)[dom_chromo1[index]]
            }
            dom_effect1 <- effect.size.equal.dom
            real.bv.add.new <- rbind(real.bv.add.new, cbind(dom_snp1, dom_chromo1, 0 ,dom_effect1, dom_effect1, dom_marker1, trait.pool[index_trait], FALSE))

          }

          if(n.quantitative[index_trait]){
            epi1_snp <- epi1_chromo <- numeric(n.quantitative[index_trait]*2)
            for(index in 1:(n.quantitative[index_trait]*2)){
              epi1_chromo[index] <- sum(epi1_marker[index] > cum_snp) + 1
              epi1_snp[index] <- epi1_marker[index] - c(0,cum_snp)[epi1_chromo[index]]
            }

            effect_matrix <- matrix(0,nrow=n.quantitative[index_trait], ncol=9)
            for(index in 1:n.quantitative[index_trait]){


              if(effect.distribution == "gauss"){
                d1 <- sort(abs(stats::rnorm(3, 0, var_quantitative[index])))
                d2 <- sort(abs(stats::rnorm(3, 0, var_quantitative[index])))
              } else{
                d1 <- sort(stats::rgamma(3, gamma.shape1, gamma.shape2))
                d2 <- sort(stats::rgamma(3, gamma.shape1, gamma.shape2))
              }

              effect_matrix[index,] <- c(d1*d2[1], d1*d2[2], d1*d2[3])
            }
            real.bv.mult.new <- cbind(epi1_snp[1:n.quantitative[index_trait]], epi1_chromo[1:n.quantitative[index_trait]],
                                      epi1_snp[-(1:n.quantitative[index_trait])], epi1_chromo[-(1:n.quantitative[index_trait])],
                                      effect_matrix)
          }

          if(n.qualitative[index_trait]>0){
            epi2_snp <- epi2_chromo <- numeric(n.qualitative[index_trait]*2)
            for(index in 1:(n.qualitative[index_trait]*2)){
              epi2_chromo[index] <- sum(epi2_marker[index] > cum_snp) + 1
              epi2_snp[index] <- epi2_marker[index] - c(0,cum_snp)[epi2_chromo[index]]
            }

            effect_matrix <- matrix(0,nrow=n.qualitative[index_trait], ncol=9)
            for(index in 1:n.qualitative[index_trait]){

              if(effect.distribution == "gauss"){
                d1 <- -abs(stats::rnorm(9, 0, var_qualitative[index]))
              } else{
                d1 <- - stats::rgamma(9, gamma.shape1, gamma.shape2)
              }

              d1[c(3,7)] <- -d1[c(3,7)]
              effect_matrix[index,] <- d1
            }
            real.bv.mult.new <- rbind(real.bv.mult.new, cbind(epi2_snp[1:n.qualitative[index_trait]], epi2_chromo[1:n.qualitative[index_trait]],
                                                              epi2_snp[-(1:n.qualitative[index_trait])], epi2_chromo[-(1:n.qualitative[index_trait])],
                                                              effect_matrix))
          }

          real.bv.add[[index_trait+so_far]] <- real.bv.add.new
          real.bv.mult[[index_trait+so_far]] <- real.bv.mult.new

        }
      }







      if(length(population)>0 && length(population$info$bitstoring)>0){
        nbits <- population$info$bitstoring
        leftover <- population$info$leftover
        bit.storing <- TRUE
      }
      if(length(population)>0 && length(population$info$miraculix)>0 && population$info$miraculix){
        miraculix <- TRUE
      }



      nbv <- max(if(is.list(real.bv.add)){length(real.bv.add)} else{as.numeric(length(real.bv.add)>0)},
                 if(is.list(real.bv.mult)){length(real.bv.mult)} else{as.numeric(length(real.bv.mult)>0)},
                 if(is.list(real.bv.dice)){length(real.bv.dice)} else{as.numeric(length(real.bv.dice)>0)},
                 if(length(trait_sum)>1){length(trait_sum)} else{0})

      if(nbv >= bv.total){
        bv.total <- nbv
        bv.calc <- nbv
        bv.random <- rep(FALSE, bv.total)
        bv.random.variance <- c(rep(0, nbv))
      }
      if(bv.total > nbv){
        if(length(polygenic.variance)< (bv.total - nbv)){
          polygenic.variance <- rep(polygenic.variance, bv.total - nbv)
        }
        bv.random <- c(rep(FALSE, nbv), rep(TRUE, bv.total - nbv))

        bv.random.variance <- c(rep(0, nbv), polygenic.variance)
        bv.calc <- nbv +1 ## WARUM STEHT HIER +1 SINN ist calc nicht gerade anzahl der ZW mit berechenbaren wert? # ALLES SO RICHTIG in breeding.diploid!
      }

      if(length(chr.nr)>0){
        chr.opt <- unique(chr.nr)
        nsnp <- numeric(length(chr.opt))
        for(index in 1:length(chr.opt)){
          nsnp[index] <- sum(chr.nr==chr.opt[index])
        }
      }


      if(length(dataset)==0 || (length(dataset)==1 && !is.list(dataset) && dataset=="random")){
        if(miraculix && length(chr.nr)==0 && miraculix.dataset){
          suppressWarnings(dataset <- list(miraculix::rhaplo(freq, indiv = nindi, loci = nsnp)))
        } else if(miraculix && miraculix.dataset){
          dataset <- list()
          for(chr_index in 1:length(chr.opt)){
            suppressWarnings(dataset[[chr_index]] <- miraculix::rhaplo(freq[which(chr.nr==chr.opt[chr_index])], indiv = nindi, loci = nsnp[chr_index]))
          }
        } else{
          dataset <- matrix((c(stats::rbinom(nindi*2*sum(nsnp),1,freq))),ncol=nindi*2, nrow=sum(nsnp))
        }
      }

      if(length(dataset)==1 && !is.list(dataset) && dataset=="all0"){
        if(miraculix && length(chr.nr)==0 && miraculix.dataset){
          suppressWarnings(dataset <- list(miraculix::rhaplo(0, indiv = nindi, loci = nsnp)))
        } else if(miraculix && miraculix.dataset){
          dataset <- list()
          for(chr_index in 1:length(chr.opt)){
            suppressWarnings(dataset[[chr_index]] <- miraculix::rhaplo(0, indiv = nindi, loci = nsnp[chr_index]))
          }
        } else{
          dataset <- matrix((c(rep(0,nindi*2*sum(nsnp)))),ncol=nindi*2, nrow=sum(nsnp))
        }

      }
      if(length(dataset)==1 && !is.list(dataset) && dataset=="homorandom"){
        if(miraculix && miraculix.dataset && length(chr.nr)==0){
          suppressWarnings(dataset <- list(miraculix::rhaplo(freq, indiv = nindi, loci = nsnp, freq2="IRGENDWAS")))
        } else if(miraculix && miraculix.dataset){
          dataset <- list()
          for(chr_index in 1:length(chr.opt)){
            suppressWarnings(dataset[[chr_index]] <- miraculix::rhaplo(freq[which(chr.nr==chr.opt[chr_index])], indiv = nindi, loci = nsnp[chr_index], freq2="IRGENDWAS"))
          }
        } else{
          dataset <- matrix((c(stats::rbinom(nindi*sum(nsnp),1,freq))),ncol=nindi, nrow=sum(nsnp))
          dataset <- dataset[,rep(1:nindi, each=2)]
        }
      }
      if(length(dataset)==1 && !is.list(dataset) && dataset=="allhetero"){
        if(miraculix && miraculix.dataset && length(chr.nr)==0){
          suppressWarnings(dataset <- list(miraculix::rhaplo(0, indiv = nindi, loci = nsnp,1)))
        } else if(miraculix && miraculix.dataset){
          dataset <- list()
          for(chr_index in 1:length(chr.opt)){
            suppressWarnings(dataset[[chr_index]] <- miraculix::rhaplo(0, indiv = nindi, loci = nsnp[chr_index],1))
          }
        } else{
          dataset <- matrix((c(rep(0,nindi*2*sum(nsnp)))),ncol=nindi*2, nrow=sum(nsnp))
          dataset[,1:nindi*2] <- 1
        }
      }

      if(length(nsnp)<2){
        if(is.list(dataset)){
          nsnp <- numeric(length(dataset))
          for(index in 1:length(dataset)){

            if(new_miraculix){
              nsnp[index] <- attr(dataset[[index]], "information")[3]
            } else{
              nsnp[index] <- attr(dataset[[index]], "information")[2]
            }

          }

        } else{
          if(length(chr.nr)>0){
            nsnp <- numeric(length(chr.opt))
            for(index in 1:length(chr.opt)){
              nsnp[index] <- sum(chr.nr==chr.opt[index])
            }
          } else{
            nsnp <- nrow(dataset)
          }
        }
      }
      if(is.list(dataset)){
        if(new_miraculix){
          nindi <- attr(dataset[[1]], "information")[5]
        } else{
          nindi <- attr(dataset[[1]], "information")[3]
        }

      } else{
        nindi <- ncol(dataset)/2
      }




      if(change.order && length(snp.position)>0){

        if(length(chr.nr)>0 && length(unique(chr.nr))>0){
          dataset_temp <- dataset
          snp.position_temp <- snp.position
          bp_temp <- bp
          snp.name_temp <- snp.name
          chr.nr_temp <- chr.nr
          so_far <- 0
          for(index in unique(chr.nr)){

            consider <- which(chr.nr_temp==index)
            nsnp_temp <- length(consider)
            order <- sort(snp.position_temp[consider], index.return=TRUE)$ix

            if(miraculix && miraculix.dataset && sum(diff(order)<0)>0){
              stop("Change order has to be executed before dataset generation // is not possible for imported datasets")
            }
            if(!miraculix || !miraculix.dataset){
              snp.position[1:nsnp_temp + so_far] <- snp.position_temp[consider][order]
              chr.nr[1:nsnp_temp + so_far] <- chr.nr_temp[consider][order]
              dataset[1:nsnp_temp + so_far,] <- dataset_temp[consider,][order,]
              if(length(bp)==length(snp.position)){
                bp[1:nsnp_temp + so_far] <- bp_temp[consider][order]
              }
              if(length(snp.name)==length(snp.position)){
                snp.name[1:nsnp_temp + so_far] <- snp.name_temp[consider][order]
              }
            }



            so_far <- so_far + nsnp_temp

          }
        } else{
          order <- sort(snp.position,index.return=TRUE)$ix

          if(miraculix && miraculix.dataset && sum(diff(order)<0)>0){
            stop("Change order has to be executed before dataset generation // is not possible for imported datasets")
          }

          if(!miraculix || !miraculix.dataset){
            snp.position <- snp.position[order]
            dataset <- dataset[order,]
            if(length(bp)==length(snp.position)){
              bp <- bp[order]
            }
            if(length(snp.name)==length(snp.position)){
              snp.name <- snp.name[order]
            }
          }

        }

      }

      if(length(genotyped.s)==0){
        genotyped.s <- stats::rbinom(nindi, 1, share.genotyped)
      }
      if(sex.s[1]=="random"){

        sex.s <- stats::rbinom(nindi, 1,sex.quota) +1
        if(add.chromosome==TRUE){
          sex.s <- population$info$sex
        }
      }
      if(sex.s[1]=="fixed"){
        sex.s <- rep(1, round(nindi*(1-sex.quota), digits = 10))
        sex.s <- c(sex.s, rep(2, nindi - length(sex.s)))
      }



      if(length(hom0)==0 && (is.list(dataset) || is.numeric(dataset[1,1]))){
        hom0 <- integer(sum(nsnp))
        hom1 <- rep(1L,sum(nsnp))
      }

    } else{
      if(is.list(dataset)){

        if(new_miraculix){
          nsnp <- attr(dataset[[1]], "information")[3]
          nindi <- attr(dataset[[1]], "information")[5]
        } else{
          nsnp <- attr(dataset[[1]], "information")[2]
          nindi <- attr(dataset[[1]], "information")[3]
        }

      } else{
        nsnp <- nrow(dataset)
        nindi <- ncol(dataset)/2
      }

    }

  }



  # Filling remaining gaps
  {
    if(store.comp.times){
      times_comp[3] <- as.numeric(Sys.time())
    }
    if(length(chr.nr)==0){
      chr.nr <- numeric(sum(nsnp))
      start1 <- 1
      for(index in 1:length(nsnp)){
        if(nsnp[index]!=0){
          chr.nr[start1:(start1+nsnp[index]-1)] <- index
          start1 <- start1 + nsnp[index]
        }

      }
    }
    if(length(bp)==0){
      bp <- numeric(sum(nsnp))
      start1 <- 1
      for(index in 1:length(nsnp)){
        if(length(chromosome.length)>=index){
          bp[start1:(start1+nsnp[index]-1)] <- ceiling((1:nsnp[index]-0.5) * chromosome.length[index] * 100000000 / nsnp[index])
        } else{
          bp[start1:(start1+nsnp[index]-1)] <- ceiling((1:nsnp[index]-0.5) * chromosome.length[1] * 100000000 / nsnp[index])
        }

        start1 <- start1 + nsnp[index]
      }
    }
    if(length(snp.name)==0){
      nsnpnr <- numeric(sum(nsnp))
      start1 <- 1
      for(index in 1:length(nsnp)){
        nsnpnr[start1:(start1+nsnp[index]-1)] <- 1:nsnp[index]
        start1 <- start1 + nsnp[index]
      }
      snp.name <- paste0("Chr", chr.nr, "SNP", nsnpnr)
    }

    chr.opt <- unique(chr.nr)

    if(length(chr.opt)==1){
      if(bpcm.conversion>0 && length(snp.position)==0){
        snp.position <- as.numeric(bp) / bpcm.conversion / 100
        chromosome.length <- max(snp.position) + min(snp.position)
      } else if(bpcm.conversion>0 && length(snp.position)>0){
        if(verbose) cat("Do not use bpcm.conversion and snp.position jointly!\n")
      }
    }
    if(length(bpcm.conversion)!=length(chr.opt)){
      bpcm.conversion <- rep(bpcm.conversion, length.out=length(chr.opt))
    }






    if((length(population)==0 || (length(population)>0 & add.chromosome)) && (length(snp.position)>0 && snp.position[1]<=0 && (length(snps.equidistant)==0 || snps.equidistant!=TRUE))){
      snp.position[1] <- snp.position[2]/2
      warning(paste("Illegal position for SNP 1 - changed to",snp.position[1]))
    }
    mindex <- max(which(chr.opt[1]==chr.nr))
    if((length(population)==0 || (length(population)>0 & add.chromosome)) && (length(snp.position)>1 && snp.position[mindex]>=chromosome.length[1] && (length(snps.equidistant)==0 || snps.equidistant!=TRUE))){
      snp.position[mindex] <- mean(c(snp.position[mindex-1], chromosome.length[1]))
      warning(paste("Illegal position for last SNP - changed to",snp.position[mindex]))
    }

    position <- snp.position
    if(length(snp.position)>0 && length(snps.equidistant)==0){
      snps.equidistant <- FALSE
    } else if(length(snps.equidistant)==0){
      snps.equidistant <- TRUE
    }

    if(length(chr.opt)==1 && (snps.equidistant || position.scaling || max(position) > sum(chromosome.length) )){
      if(snps.equidistant){
        position <- 0.5:(sum(nsnp))*10
      }
      min.p <- min(position)
      max.p <- max(position)
      #  position.scal <- (position - (min.p - length.before)) / (max.p-min.p + length.before+ length.behind) * chromosome.length
      position <- (position) / (max.p-min.p + length.before+ length.behind) * chromosome.length[1]

    }

    if(length(chr.opt)!=length(chromosome.length)){
      if(length(snp.position)==length(chr.nr)){
        chromosome.length <- numeric(length(chr.opt))
        for(index in 1:length(chr.opt)){
          chromosome.length[index] <- max(snp.position[chr.nr==chr.opt[index]]) + min(snp.position[chr.nr==chr.opt[index]])
        }
      } else{
        chromosome.length <- rep(chromosome.length, length.out=length(chr.opt))
      }

    }
  }

  # Editing of the population list - Info + genotypes

  {
    if(store.comp.times){
      times_comp[4] <- as.numeric(Sys.time())
    }
    if(length(chr.opt)==1){
      if(length(population)==0){
        population <- list()

        ## KEINE EDITS ZWISCHEN [[1]] und [[18]] snps.equidistant ++ Miraculix aenderungen sonst erforderlich!
        population$info <- list()
        population$info$schlather.slot1 <- "miraculix_not_activated"
        population$info$chromosome <- 1L
        population$info$snp <- as.integer(nsnp)
        population$info$position <- list()
        population$info$position[[1]] <- position
        population$info$snp.base <- rbind(hom0,hom1, deparse.level = 0)
        population$info$snp.position <- position
        population$info$length <- chromosome.length
        population$info$length.total <- c(0,population$info$length)
        population$info$func <- FALSE
        population$info$size <- matrix(0L,nrow=1, ncol=2)
        population$info$bve <- FALSE
        population$info$bv.calculated <- FALSE
        population$info$breeding.totals <- list()
        population$info$bve.data <- list()
        population$info$bv.nr <- bv.total
        population$info$bv.random <- bv.random
        population$info$bv.random.variance <- bv.random.variance
        population$info$snps.equidistant <- snps.equidistant
        population$info$origin.gen <- 1L
        population$info$cumsnp <- nsnp
        population$info$bp <- bp
        population$info$snp.name <- snp.name
        population$info$next.animal <- 1
        population$info$next.litter <- 1
        population$info$next.pen <- 1
        population$info$phenotypic.transform <- rep(FALSE, bv.total)
        population$info$phenotypic.transform.function <- list()
        population$info$culling.stats <- list()
        population$info$version_MoBPS <- utils::sessionInfo()$otherPkgs$MoBPS
        population$info$version_miraculix <- utils::sessionInfo()$otherPkgs$miraculix
        population$info$cohort.index <- 1
        population$info$array.name = "Full_Array"
        population$info$array.markers = list(rep(TRUE,nsnp))
        population$info$array.is_subset = FALSE
        population$info$default.parameter.name = NULL
        population$info$default.parameter.value = list()
        population$info$chromosome.name <- chr.opt
        population$info$size.scaling = size.scaling
        population$info$founder_pools = founder.pool
        population$info$founder_multi = FALSE
        population$info$founder_multi_calc = FALSE
        population$info$pool_effects = FALSE
        population$info$pool_effects_calc = FALSE

        population$info$pedigree_error = FALSE

        population$info$pen.size <- cbind(1,1)
        population$info$litter.effect.active <- FALSE
        population$info$pen.effect.active <- FALSE
        population$info$max.time.point = time.point

        if(length(litter.effect.covariance)>0 && sum( abs(litter.effect.covariance))>0){
          population$info$litter.effect.active <- TRUE
        }
        if(length(pen.effect.covariance)>0 && sum( abs(pen.effect.covariance))>0){
          population$info$pen.effect.active <- TRUE
        }


        if(length(is.maternal)==0){
          population$info$is.maternal <- rep(FALSE, bv.total)
        } else{
          if(length(is.maternal)==bv.total){
            population$info$is.maternal <- is.maternal
          } else {
            population$info$is.maternal <- c(population$info$is.materal, rep(is.maternal, length.out = bv.total - length(population$info$is.materal)))
          }

        }
        if(length(is.paternal)==0){
          population$info$is.paternal <- rep(FALSE, bv.total)
        } else{
          if(length(is.paternal)==bv.total){
            population$info$is.paternal <- is.paternal
          } else {
            population$info$is.paternal <- c(population$info$is.paternal, rep(is.paternal, length.out = bv.total - length(population$info$is.paternal)))
          }
        }



        population$info$is.combi <- rep(FALSE, bv.total)
        population$info$progress.bar = progress.bar


        if(length(bve.mult.factor)==0){
          population$info$bve.mult.factor <- rep(1L, bv.total)
        } else{

          if(length(bve.mult.factor)==bv.total){
            population$info$bve.mult.factor <- bve.mult.factor
          } else {
            population$info$bve.mult.factor <- c(population$info$bve.mult.factor, rep(bve.mult.factor, length.out = bv.total - length(population$info$bve.mult.factor)))
          }
        }

        if(length(bve.poly.factor)==0){
          population$info$bve.poly.factor <- rep(1L, bv.total)
        } else{

          if(length(bve.poly.factor)==bv.total){
            population$info$bve.poly.factor <- bve.poly.factor
          } else {
            population$info$bve.poly.factor <- c(population$info$bve.poly.factor, rep(bve.poly.factor, length.out = bv.total - length(population$info$bve.poly.factor)))
          }

        }
        if(length(base.bv)==0){
          population$info$base.bv <- rep(100L, bv.total)
        } else{

          if(length(base.bv)==bv.total){
            population$info$base.bv <- base.bv
          } else {
            population$info$base.bv <- c(population$info$base.bv, rep(base.bv, length.out = bv.total - length(population$info$base.bv)))
          }
        }



      } else if(add.chromosome){
        if(length(chr.opt)>1){
          stop("You can only add a single chromosome using add.chromosome!")
        }
        population$info$chromosome <- population$info$chromosome + 1L
        population$info$snp <- c(population$info$snp, as.integer(nsnp))
        population$info$position[[length(population$info$position)+1]] <- position
        population$info$snp.position <- c(population$info$snp.position, position + max(population$info$length.total))
        population$info$length <- c(population$info$length, chromosome.length)
        population$info$length.total <- cumsum(c(0,population$info$length))
        population$info$snp.base <- cbind(population$info$snp.base , rbind(hom0,hom1, deparse.level = 0), deparse.level = 0)
        population$info$cumsnp <- c(population$info$cumsnp, sum(population$info$snp))
        population$info$bp <- c(population$info$bp, bp)
        population$info$snp.name <- c(population$info$snp.name, snp.name)

        population$info$chromosome.name <- c(population$info$chromosome.name, chr.opt)







        if(length(population$info$array.name)>1){
          stop("New chromosomes can not be added after more than one array is entered!")
        } else{
          population$info$array.markers[[1]] <- c(population$info$array.markers[[1]], rep(TRUE, nsnp))
        }

      }
      if(generation!=1){
        take <- which(population$info$origin.gen==generation)
        if(length(take)==1){
          if(population$info$miraculix){
            origin_code <- take
          } else{
            origin_code <- population$info$origin.gen[take]
          }

        } else{
          if(population$info$miraculix){
            if(length(population$info$origin.gen)<64){
              population$info$origin.gen <- c(population$info$origin.gen, as.integer(generation))
              origin_code <- length(population$info$origin.gen)
            } else{
              warning("To many origin generation!")
              warning("Delete second lowest origin.gen")
              switch <- sort(population$info$origin.gen, index.return=TRUE)[[2]]
              population$info$origin.gen[switch] <- as.integer(generation)
              origin_code <- switch
            }
          } else{
            if(length(population$info$origin.gen)<32){
              population$info$origin.gen <- c(population$info$origin.gen, as.integer(generation))
              origin_code <- length(population$info$origin.gen)
            } else{
              warning("To many origin generation!")
              warning("Delete second lowest origin.gen")
              switch <- sort(population$info$origin.gen, index.return=TRUE)[[2]]
              population$info$origin.gen[switch] <- as.integer(generation)
              origin_code <- switch
            }
          }

        }
      } else{
        origin_code <- generation
      }



      if(length(population)==1){
        population$breeding <- list()
        population$breeding[[1]] <- list()
      }
      if(length(population$breeding)==0 || length(population$breeding[[1]])==0){
        population$breeding[[1]][[1]] <- list()
      }
      if(length(population$breeding[[1]])==1 ){
        population$breeding[[1]][[2]] <- list()
      }

      if(generation!=1){
        if(length(population$breeding)==(generation-1) || length(population$breeding[[generation]])==0){
          population$breeding[[generation]] <- list()
          population$breeding[[generation]][[1]] <- list()
          population$info$size <- rbind(population$info$size,0L, deparse.level = 0)
        }
        if(length(population$breeding[[generation]])==1){
          population$breeding[[generation]][[2]] <- list()
        }
      }

      if(internal.geno && length(internal.dataset)>0){
        dataset <- internal.dataset
      }

      counter <- c(length(population$breeding[[generation]][[1]]),length(population$breeding[[generation]][[2]]))+1L # maennlich/weibliche Tiere bisher
      counter.start <- counter

      if(add.chromosome==FALSE){

        for(index in 1:length(sex.s)){
          sex <- sex.s[index]
          genotyped <- genotyped.s[index]

          population$breeding[[generation]][[sex]][[counter[sex]]] <- list()
          population$breeding[[generation]][[sex]][[counter[sex]]][[1]] <- c(0, sum(population$info$length))
          population$breeding[[generation]][[sex]][[counter[sex]]][[2]] <- c(0,sum(population$info$length))
          if(add.chromosome.ends==TRUE){
            population$breeding[[generation]][[sex]][[counter[sex]]][[1]] <- population$info$length.total
            population$breeding[[generation]][[sex]][[counter[sex]]][[2]] <- population$info$length.total
          }
          population$breeding[[generation]][[sex]][[counter[sex]]][[3]] <- NULL
          population$breeding[[generation]][[sex]][[counter[sex]]][[4]] <- NULL
          population$breeding[[generation]][[sex]][[counter[sex]]][[5]] <- codeOriginsU(matrix(c(origin_code, sex, counter[sex], 1),nrow=(length(population$breeding[[generation]][[sex]][[counter[sex]]][[1]])-1), ncol=4, byrow=TRUE))
          population$breeding[[generation]][[sex]][[counter[sex]]][[6]] <- codeOriginsU(matrix(c(origin_code, sex, counter[sex], 2),nrow=(length(population$breeding[[generation]][[sex]][[counter[sex]]][[2]])-1), ncol=4, byrow=TRUE))
          population$breeding[[generation]][[sex]][[counter[sex]]][[7]] <- c(generation, sex, counter[sex])
          population$breeding[[generation]][[sex]][[counter[sex]]][[8]] <- c(generation, sex, counter[sex])

          if(internal.geno){
            if(miraculix && miraculix.dataset){
              if(new_miraculix){
                population$breeding[[generation]][[sex]][[counter[sex]]][[9]] <- miraculix::haplomatrix(as.matrix(dataset[[1]], sel.indiv = index))
              } else{
                population$breeding[[generation]][[sex]][[counter[sex]]][[9]] <- miraculix::haplomatrix(as.matrix(dataset[[1]],indiv = index))
              }

              population$breeding[[generation]][[sex]][[counter[sex]]][[10]] <- "Placeholder_Pointer_Martin"
            } else if(miraculix){
              population$breeding[[generation]][[sex]][[counter[sex]]][[9]] <- miraculix::haplomatrix(dataset[,(index*2-c(1,0))])
              population$breeding[[generation]][[sex]][[counter[sex]]][[10]] <- "Placeholder_Pointer_Martin"
            } else if(bit.storing){
              population$breeding[[generation]][[sex]][[counter[sex]]][[9]] <- bit.storing(dataset[,(index*2-1)], nbits)
              population$breeding[[generation]][[sex]][[counter[sex]]][[10]] <- bit.storing(dataset[,(index*2)], nbits)
            } else{
              population$breeding[[generation]][[sex]][[counter[sex]]][[9]] <- dataset[,(index*2-1)]
              population$breeding[[generation]][[sex]][[counter[sex]]][[10]] <- dataset[,(index*2)]
            }
          }


          population$breeding[[generation]][[sex]][[counter[sex]]][[11]] <- NULL
          population$breeding[[generation]][[sex]][[counter[sex]]][[12]] <- NULL
          #      population$breeding[[generation]][[sex]][[counter[sex]]][[13]] <- "test"
          population$breeding[[generation]][[sex]][[counter[sex]]][[15]] <- rep(0L, population$info$bv.nr)
          population$breeding[[generation]][[sex]][[counter[sex]]][[16]] <- genotyped
          population$breeding[[generation]][[sex]][[counter[sex]]][[21]] <- cbind(generation, sex, counter[sex], deparse.level = 0)
          storage.mode(population$breeding[[generation]][[sex]][[counter[sex]]][[21]]) <- "integer"

          population$breeding[[generation]][[sex]][[counter[sex]]][[22]] <- if(genotyped>0){1} else{NULL}

          population$breeding[[generation]][[sex]][[counter[sex]]][[23]] <- NULL ## permanent environmental effects
          population$breeding[[generation]][[sex]][[counter[sex]]][[24]] <- NULL ## random environmental effects

          population$breeding[[generation]][[sex]][[counter[sex]]][[25]] <- FALSE ## has BV been calculated
          population$breeding[[generation]][[sex]][[counter[sex]]][[26]] <- NULL ## for which BV has been calculated

          population$breeding[[generation]][[sex]][[counter[sex]]][[27]] <- list() ## List containing each individual phenotypic observation

          population$breeding[[generation]][[sex]][[counter[sex]]][[28]] <- numeric(0)

          if(population$info$litter.effect.active){
            population$breeding[[generation]][[sex]][[counter[sex]]][[29]]  <- rep(0L, population$info$bv.nr) ## Litter effect
            population$breeding[[generation]][[sex]][[counter[sex]]][[31]]  <- 0L ## Litter nr
          }

          if(population$info$pen.effect.active){
            population$breeding[[generation]][[sex]][[counter[sex]]][[30]]  <- rep(0L, population$info$bv.nr) ## Pen effect
            population$breeding[[generation]][[sex]][[counter[sex]]][[32]]  <- 0L ## Pen nr
          }

          population$breeding[[generation]][[sex]][[counter[sex]]][[33]] = numeric(0)
          population$breeding[[generation]][[sex]][[counter[sex]]][[34]] = numeric(0)

          population$breeding[[generation]][[sex]][[counter[sex]]][[35]] = numeric(0) # Recriprocal translocation (RT)
          population$breeding[[generation]][[sex]][[counter[sex]]][[36]] = numeric(0) #


          population$breeding[[generation]][[sex]][[counter[sex]]][[37]] <- c(generation, sex, counter[sex]) # Pedigree with errors
          population$breeding[[generation]][[sex]][[counter[sex]]][[38]] <- c(generation, sex, counter[sex]) #


          population$breeding[[generation]][[sex]][[counter[sex]]][[39]] <- "placeholder"
          population$info$size[generation,sex] <- population$info$size[generation,sex] +1L
          counter[sex] <- counter[sex] + 1L
        }

        if(length(population$breeding[[generation]])==2){
          population$breeding[[generation]][[3]] <- matrix(0L, nrow= population$info$bv.nr, ncol=counter[1]-1) # Selektionsfunktion
          population$breeding[[generation]][[4]] <- matrix(0L, nrow= population$info$bv.nr, ncol=counter[2]-1)
          population$breeding[[generation]][[5]] <- rep(class,counter[1]-1) # Migrationslevel
          population$breeding[[generation]][[6]] <- rep(class,counter[2]-1)
          population$breeding[[generation]][[7]] <- matrix(0L, nrow= population$info$bv.nr, ncol=counter[1]-1) # realer ZW
          population$breeding[[generation]][[8]] <- matrix(0L, nrow= population$info$bv.nr, ncol=counter[2]-1)
          population$breeding[[generation]][[9]] <- matrix(NA, nrow= population$info$bv.nr, ncol=counter[1]-1) # Phenotype
          population$breeding[[generation]][[10]] <- matrix(NA, nrow= population$info$bv.nr, ncol=counter[2]-1)
          population$breeding[[generation]][[11]] <- rep(time.point,counter[1]-1) # Time point
          population$breeding[[generation]][[12]] <- rep(time.point,counter[2]-1)
          population$breeding[[generation]][[13]] <- rep(creating.type,counter[1]-1) # Creating.type
          population$breeding[[generation]][[14]] <- rep(creating.type,counter[2]-1)
          population$breeding[[generation]][[15]] <- seq(population$info$next.animal, population$info$next.animal + counter[1] -2, length.out= counter[1] -1)
          population$info$next.animal <- population$info$next.animal + counter[1] - 1
          population$breeding[[generation]][[16]] <- seq(population$info$next.animal, population$info$next.animal + counter[2] -2, length.out= counter[2] -1)
          population$info$next.animal <- population$info$next.animal + counter[2] - 1
          population$breeding[[generation]][[17]] <- rep(NA,counter[1]-1) # Time of death point
          population$breeding[[generation]][[18]] <- rep(NA,counter[2]-1)
          population$breeding[[generation]][[19]] <- matrix(0L, nrow= population$info$bv.nr, ncol=counter[1]-1) # Reliabilities
          population$breeding[[generation]][[20]] <- matrix(0L, nrow= population$info$bv.nr, ncol=counter[2]-1)
          population$breeding[[generation]][[21]] <- matrix(0L, nrow= population$info$bv.nr, ncol=counter[1]-1) # Last applied selection index
          population$breeding[[generation]][[22]] <- matrix(0L, nrow= population$info$bv.nr, ncol=counter[2]-1)
          population$breeding[[generation]][[23]] <- rep(time.point,counter[1]-1) # Age time point
          population$breeding[[generation]][[24]] <- rep(time.point,counter[2]-1)
          population$breeding[[generation]][[25]] <- rep(NA,counter[1]-1) # Death time point
          population$breeding[[generation]][[26]] <- rep(NA,counter[2]-1)
          population$breeding[[generation]][[27]] <- matrix(0L, nrow= population$info$bv.nr, ncol=counter[1]-1) # offspring phenotype
          population$breeding[[generation]][[28]] <- matrix(0L, nrow= population$info$bv.nr, ncol=counter[2]-1)
          population$breeding[[generation]][[29]] <- matrix(0L, nrow= population$info$bv.nr, ncol=counter[1]-1) # number of offspring used
          population$breeding[[generation]][[30]] <- matrix(0L, nrow= population$info$bv.nr, ncol=counter[2]-1)
          population$breeding[[generation]][[31]] <- rep(0,counter[1]-1) # Time of death point
          population$breeding[[generation]][[32]] <- rep(0,counter[2]-1)

          population$breeding[[generation]][[33]] <- rep(0,counter[1]-1) # Litter nr male
          population$breeding[[generation]][[34]] <- rep(0,counter[2]-1) # Litter nr female

          population$breeding[[generation]][[35]] <- rep(0,counter[1]-1) # Pen nr male
          population$breeding[[generation]][[36]] <- rep(0,counter[2]-1) # Pen nr female

          population$breeding[[generation]][[37]] <- rep(founder.pool, counter[1]-1)
          population$breeding[[generation]][[38]] <- rep(founder.pool, counter[2]-1)

          population$breeding[[generation]][[39]] <- rep(NA, counter[1]-1) ## time point of death
          population$breeding[[generation]][[40]] <- rep(NA, counter[2]-1)

          population$breeding[[generation]][[41]] <- rep(NA, counter[1]-1) ## culling type
          population$breeding[[generation]][[42]] <- rep(NA, counter[2]-1)

          population$breeding[[generation]][[43]] <- rep(NA, counter[1]-1) ## time point of first pheno
          population$breeding[[generation]][[44]] <- rep(NA, counter[2]-1)

          population$breeding[[generation]][[45]] <- rep(NA, counter[1]-1) ## time point of genotyping
          population$breeding[[generation]][[46]] <- rep(NA, counter[2]-1)

          if(counter[1] >1){
            population$breeding[[generation]][[47]] <- rbind(generation, 1, 1:(counter[1]-1),
                                                             generation, 1, 1:(counter[1]-1),
                                                             generation, 1, 1:(counter[1]-1),
                                                             population$breeding[[generation]][[15]],
                                                             population$breeding[[generation]][[15]],
                                                             population$breeding[[generation]][[15]]) # pedigree (no errors)

          } else{
            population$breeding[[generation]][[47]]  = matrix(0L, nrow= 12, ncol=counter[1] -1)
          }
          if(counter[2] > 1){
            population$breeding[[generation]][[48]] <- rbind(generation, 2, 1:(counter[2]-1),
                                                             generation, 2, 1:(counter[2]-1),
                                                             generation, 2, 1:(counter[2]-1),
                                                             population$breeding[[generation]][[16]],
                                                             population$breeding[[generation]][[16]],
                                                             population$breeding[[generation]][[16]])

          } else{
            population$breeding[[generation]][[48]] = matrix(0L, nrow= 12, ncol=counter[2] -1)
          }



          # calculate Real-ZW
        } else{
          population$breeding[[generation]][[3]] <- cbind(population$breeding[[generation]][[3]], matrix(0L, nrow= population$info$bv.nr, ncol=counter[1]-counter.start[1])) # Selektionsfunktion
          population$breeding[[generation]][[4]] <- cbind(population$breeding[[generation]][[4]], matrix(0L, nrow= population$info$bv.nr, ncol=counter[2]-counter.start[2]))
          population$breeding[[generation]][[5]] <- c(population$breeding[[generation]][[5]], rep(class ,counter[1]-counter.start[1])) # Migrationslevel
          population$breeding[[generation]][[6]] <- c(population$breeding[[generation]][[6]], rep(class ,counter[2]-counter.start[2]))
          population$breeding[[generation]][[7]] <- cbind(population$breeding[[generation]][[7]] , matrix(0L, nrow= population$info$bv.nr, ncol=counter[1]-counter.start[1])) # realer ZW
          population$breeding[[generation]][[8]] <- cbind(population$breeding[[generation]][[8]] , matrix(0L, nrow= population$info$bv.nr, ncol=counter[2]-counter.start[2]))
          population$breeding[[generation]][[9]] <- cbind(population$breeding[[generation]][[9]] , matrix(NA, nrow= population$info$bv.nr, ncol=counter[1]-counter.start[1])) # geschaetzer ZW
          population$breeding[[generation]][[10]] <-cbind(population$breeding[[generation]][[10]] , matrix(NA, nrow= population$info$bv.nr, ncol=counter[2]-counter.start[2]))
          population$breeding[[generation]][[11]] <- c(population$breeding[[generation]][[11]], rep(time.point ,counter[1]-counter.start[1])) # Time point
          population$breeding[[generation]][[12]] <- c(population$breeding[[generation]][[12]], rep(time.point ,counter[2]-counter.start[2]))
          population$breeding[[generation]][[13]] <- c(population$breeding[[generation]][[13]], rep(time.point ,counter[1]-counter.start[1])) # Creating type
          population$breeding[[generation]][[14]] <- c(population$breeding[[generation]][[14]], rep(time.point ,counter[2]-counter.start[2]))

          tmp111 = seq(population$info$next.animal, population$info$next.animal + counter[1] - counter.start[1]-1, length.out= counter[1] -counter.start[1])
          population$breeding[[generation]][[15]] <- c(population$breeding[[generation]][[15]] , tmp111)
          population$info$next.animal <- population$info$next.animal + counter[1] - counter.start[1]
          tmp222 = seq(population$info$next.animal, population$info$next.animal + counter[2] - counter.start[2]-1, length.out= counter[2] -counter.start[2])
          population$breeding[[generation]][[16]] <- c(population$breeding[[generation]][[16]] , tmp222)
          population$info$next.animal <- population$info$next.animal + counter[2] - counter.start[2]
          population$breeding[[generation]][[17]] <- c(population$breeding[[generation]][[17]], rep(NA ,counter[1]-counter.start[1])) # Time of death
          population$breeding[[generation]][[18]] <- c(population$breeding[[generation]][[18]], rep(NA ,counter[2]-counter.start[2]))
          population$breeding[[generation]][[19]] <- cbind(population$breeding[[generation]][[19]], matrix(0L, nrow= population$info$bv.nr, ncol=counter[1]-counter.start[1])) # Reliabilities
          population$breeding[[generation]][[20]] <- cbind(population$breeding[[generation]][[20]], matrix(0L, nrow= population$info$bv.nr, ncol=counter[2]-counter.start[2]))
          population$breeding[[generation]][[21]] <- cbind(population$breeding[[generation]][[21]], matrix(0L, nrow= population$info$bv.nr, ncol=counter[1]-counter.start[1])) # Last applied selection index
          population$breeding[[generation]][[22]] <- cbind(population$breeding[[generation]][[22]], matrix(0L, nrow= population$info$bv.nr, ncol=counter[2]-counter.start[2]))
          population$breeding[[generation]][[23]] <- c(population$breeding[[generation]][[23]], rep(time.point ,counter[1]-counter.start[1])) # Age Time point
          population$breeding[[generation]][[24]] <- c(population$breeding[[generation]][[24]], rep(time.point ,counter[2]-counter.start[2]))
          population$breeding[[generation]][[25]] <- c(population$breeding[[generation]][[25]], rep(NA ,counter[1]-counter.start[1])) # Death Time point
          population$breeding[[generation]][[26]] <- c(population$breeding[[generation]][[26]], rep(NA ,counter[2]-counter.start[2]))
          population$breeding[[generation]][[27]] <- cbind(population$breeding[[generation]][[27]] , matrix(0L, nrow= population$info$bv.nr, ncol=counter[1]-counter.start[1])) # offspring phenotype
          population$breeding[[generation]][[28]] <-cbind(population$breeding[[generation]][[28]] , matrix(0L, nrow= population$info$bv.nr, ncol=counter[2]-counter.start[2]))
          population$breeding[[generation]][[29]] <- cbind(population$breeding[[generation]][[29]] , matrix(0L, nrow= population$info$bv.nr, ncol=counter[1]-counter.start[1])) # number of offspring used
          population$breeding[[generation]][[30]] <-cbind(population$breeding[[generation]][[30]] , matrix(0L, nrow= population$info$bv.nr, ncol=counter[2]-counter.start[2]))

          population$breeding[[generation]][[31]] <- c(population$breeding[[generation]][[31]], rep(0, counter[1]-counter.start[1])) # Last applied selection index
          population$breeding[[generation]][[32]] <- c(population$breeding[[generation]][[32]], rep(0, counter[2]-counter.start[2]))

          population$breeding[[generation]][[33]] <- c(population$breeding[[generation]][[33]], rep(0 ,counter[1]-counter.start[1])) # Migrationslevel
          population$breeding[[generation]][[34]] <- c(population$breeding[[generation]][[34]], rep(0 ,counter[2]-counter.start[2]))
          population$breeding[[generation]][[35]] <- c(population$breeding[[generation]][[35]], rep(0 ,counter[1]-counter.start[1])) # Migrationslevel
          population$breeding[[generation]][[36]] <- c(population$breeding[[generation]][[36]], rep(0 ,counter[2]-counter.start[2]))

          population$breeding[[generation]][[37]] <- c(population$breeding[[generation]][[37]], rep(founder.pool, counter[1]-counter.start[1]))
          population$breeding[[generation]][[38]] <- c(population$breeding[[generation]][[38]], rep(founder.pool, counter[2]-counter.start[2]))


          population$breeding[[generation]][[39]] <- c(population$breeding[[generation]][[39]], rep(NA, counter[1]-counter.start[1])) ## time point of death
          population$breeding[[generation]][[40]] <- c(population$breeding[[generation]][[40]], rep(NA, counter[2]-counter.start[2]))

          population$breeding[[generation]][[41]] <- c(population$breeding[[generation]][[41]], rep(NA, counter[1]-counter.start[1])) ## culling type
          population$breeding[[generation]][[42]] <- c(population$breeding[[generation]][[42]], rep(NA, counter[2]-counter.start[2]))

          population$breeding[[generation]][[43]] <- c(population$breeding[[generation]][[43]], rep(NA, counter[1]-counter.start[1])) ## time point of first pheno
          population$breeding[[generation]][[44]] <- c(population$breeding[[generation]][[44]], rep(NA, counter[2]-counter.start[2]))

          population$breeding[[generation]][[45]] <- c(population$breeding[[generation]][[45]], rep(NA, counter[1]-counter.start[1])) ## time point of first pheno
          population$breeding[[generation]][[46]] <- c(population$breeding[[generation]][[46]], rep(NA, counter[2]-counter.start[2]))

          if(counter[1] > counter.start[1]){
            population$breeding[[generation]][[47]] <- cbind(population$breeding[[generation]][[47]], rbind(generation, 1, counter.start[1]:(counter[1]-1),
                                                                                                            generation, 1, counter.start[1]:(counter[1]-1),
                                                                                                            generation, 1, counter.start[1]:(counter[1]-1),
                                                                                                            tmp111,
                                                                                                            tmp111,
                                                                                                            tmp111)) # pedigree (no errors)
          }

          if(counter[2] > counter.start[2]){
            population$breeding[[generation]][[48]] <- cbind(population$breeding[[generation]][[48]], rbind(generation, 2, counter.start[2]:(counter[2]-1),
                                                                                                            generation, 2, counter.start[2]:(counter[2]-1),
                                                                                                            generation, 2, counter.start[2]:(counter[2]-1),
                                                                                                            tmp222,
                                                                                                            tmp222,
                                                                                                            tmp222))
          }


          population$info$founder_pools = unique(c(population$info$founder_pools, founder.pool))
          population$info$founder_multi = if(length(population$info$founder_pools)>1){TRUE} else{FALSE}
          population$info$founder_multi_calc = if(length(population$info$founder_pools)>1){TRUE} else{FALSE}
        }

        population$info$sex <- c(population$info$sex, sex.s)


      } else{
        counter <- c(1,1)
        for(index in 1:length(sex.s)){
          sex <- sex.s[index]
          genotyped <- genotyped.s[index]
          population$breeding[[generation]][[sex]][[counter[sex]]][[1]][2] <- sum(population$info$length)
          population$breeding[[generation]][[sex]][[counter[sex]]][[2]][2] <- sum(population$info$length)
          if(add.chromosome.ends==TRUE){
            population$breeding[[generation]][[sex]][[counter[sex]]][[1]] <- population$info$length.total
            population$breeding[[generation]][[sex]][[counter[sex]]][[2]] <- population$info$length.total
            population$breeding[[generation]][[sex]][[counter[sex]]][[5]] <- codeOriginsU(matrix(c(origin_code, sex, counter[sex], 1),nrow=(length(population$breeding[[generation]][[sex]][[counter[sex]]][[1]])-1), ncol=4, byrow=TRUE))
            population$breeding[[generation]][[sex]][[counter[sex]]][[6]] <- codeOriginsU(matrix(c(origin_code, sex, counter[sex], 2),nrow=(length(population$breeding[[generation]][[sex]][[counter[sex]]][[2]])-1), ncol=4, byrow=TRUE))

          }
          population$breeding[[generation]][[sex]][[counter[sex]]][[7]] <- c(generation, sex, counter[sex])
          population$breeding[[generation]][[sex]][[counter[sex]]][[8]] <- c(generation, sex, counter[sex])

          if(internal.geno){
            if(miraculix && miraculix.dataset){
              if(length(population$breeding[[generation]][[sex]][[counter[sex]]][[9]])==0){

                if(new_miraculix){
                  population$breeding[[generation]][[sex]][[counter[sex]]][[9]] <- miraculix::haplomatrix(as.matrix(dataset[[1]], sel.indiv = index))
                } else{
                  population$breeding[[generation]][[sex]][[counter[sex]]][[9]] <- miraculix::haplomatrix(as.matrix(dataset[[1]], indiv = index))
                }

                population$breeding[[generation]][[sex]][[counter[sex]]][[10]] <- "Placeholder_Pointer_Martin"
              } else{

                if(new_miraculix){
                  population$breeding[[generation]][[sex]][[counter[sex]]][[9]] <- miraculix::haplomatrix(rbind(as.matrix(population$breeding[[generation]][[sex]][[counter[sex]]][[9]]),as.matrix(dataset[[1]], sel.indiv = index)))
                } else{
                  population$breeding[[generation]][[sex]][[counter[sex]]][[9]] <- miraculix::haplomatrix(rbind(as.matrix(population$breeding[[generation]][[sex]][[counter[sex]]][[9]]),as.matrix(dataset[[1]], indiv = index)))
                }

                population$breeding[[generation]][[sex]][[counter[sex]]][[10]] <- "Placeholder_Pointer_Martin"
              }

            } else if(miraculix){
              population$breeding[[generation]][[sex]][[counter[sex]]][[9]] <- miraculix::haplomatrix(rbind(as.matrix(population$breeding[[generation]][[sex]][[counter[sex]]][[9]]),dataset[,(index*2-c(1,0))]))
              population$breeding[[generation]][[sex]][[counter[sex]]][[10]] <- "Placeholder_Pointer_Martin"
            } else if(bit.storing){
              if(leftover==0){
                population$breeding[[generation]][[sex]][[counter[sex]]][[9]] <- c(population$breeding[[generation]][[sex]][[counter[sex]]][[9]], bit.storing(dataset[,(index*2-1)]),nbits)
                population$breeding[[generation]][[sex]][[counter[sex]]][[10]] <- c(population$breeding[[generation]][[sex]][[counter[sex]]][[10]], bit.storing(dataset[,(index*2)]), nbits)
              } else{
                population$breeding[[generation]][[sex]][[counter[sex]]][[9]] <- c(population$breeding[[generation]][[sex]][[counter[sex]]][[9]][-length(population$breeding[[generation]][[sex]][[counter[sex]]][[9]])],
                                                                                   bit.storing(c(bit.snps(population$breeding[[generation]][[sex]][[counter[sex]]][[9]][length(population$breeding[[generation]][[sex]][[counter[sex]]][[9]])], nbits)[(nbits-leftover+1):nbits],dataset[,(index*2-1)]),nbits))
                population$breeding[[generation]][[sex]][[counter[sex]]][[10]] <- c(population$breeding[[generation]][[sex]][[counter[sex]]][[10]][-length(population$breeding[[generation]][[sex]][[counter[sex]]][[10]])],
                                                                                    bit.storing(c(bit.snps(population$breeding[[generation]][[sex]][[counter[sex]]][[10]][length(population$breeding[[generation]][[sex]][[counter[sex]]][[10]])], nbits)[(nbits-leftover+1):nbits],dataset[,(index*2)]),nbits))

              }
            } else{
              population$breeding[[generation]][[sex]][[counter[sex]]][[9]] <- c(population$breeding[[generation]][[sex]][[counter[sex]]][[9]], dataset[,(index*2-1)])
              population$breeding[[generation]][[sex]][[counter[sex]]][[10]] <- c(population$breeding[[generation]][[sex]][[counter[sex]]][[10]], dataset[,(index*2)])

            }
          }


          counter[sex] <- counter[sex] + 1
        }

      }


      if(length(name.cohort)==0 && add.chromosome==FALSE){
        name.cohort <- paste0("Cohort_", population$info$cohort.index)
        population$info$cohort.index <- population$info$cohort.index + 1
      }

      if(length(name.cohort)==1 && ((counter-counter.start)[1]>0 && (counter-counter.start)[2]>0)){
        name.cohort = paste0(name.cohort, c("_M", "_F"))
        if(verbose) cat("Both sexes in the cohort. Added _M, _F to cohort names!\n")
      }

      if(length(name.cohort)>=1 && add.chromosome==FALSE){

        if((counter-counter.start)[1]>0 && (counter-counter.start)[2]>0){
          population$info$cohorts <- rbind(population$info$cohorts, c(name.cohort[1], generation, (counter - counter.start)[1], 0, class, counter.start[1], 0,
                                                                      time.point, creating.type, NA, NA),
                                           c(name.cohort[2], generation, 0, (counter - counter.start)[2], class, 0, counter.start[2],
                                             time.point, creating.type, NA, NA))


          temp1 = get.id(population, cohorts = name.cohort[1])
          population$info$cohorts[nrow(population$info$cohorts)-1, 10:11] = c(min(temp1), max(temp1))
          temp1 = get.id(population, cohorts = name.cohort[2])
          population$info$cohorts[nrow(population$info$cohorts), 10:11] = c(min(temp1), max(temp1))

          if(verbose){
            posi <- get.database(population, cohorts = name.cohort[1])
            cat(paste0("Successfully generated cohort: ",  name.cohort[1], "\n",
                       "Database position: ", posi[1], " (gen), ", posi[2], " (sex), ", posi[3], " (first), ", posi[4], " (last).\n" ))
            posi <- get.database(population, cohorts = name.cohort[2])
            cat(paste0("Successfully generated cohort: ",  name.cohort[2], "\n",
                       "Database position: ", posi[1], " (gen), ", posi[2], " (sex), ", posi[3], " (first), ", posi[4], " (last).\n" ))
          }

        } else{
          population$info$cohorts <- rbind(population$info$cohorts, c(name.cohort, generation, counter - counter.start, class, counter.start,
                                                                      time.point, creating.type, NA, NA))

          temp1 = get.id(population, cohorts = name.cohort)
          population$info$cohorts[nrow(population$info$cohorts), 10:11] = c(min(temp1), max(temp1))

          if(verbose){
            posi <- get.database(population, cohorts = name.cohort)
            cat(paste0("Successfully generated cohort: ", name.cohort, "\n",
                       "Database position: ", posi[1], " (gen), ", posi[2], " (sex), ", posi[3], " (first), ", posi[4], " (last).\n" ))
          }

        }
        if(nrow(population$info$cohorts)<=2){
          colnames(population$info$cohorts) <- c("name","generation", "male individuals", "female individuals", "class", "position first male", "position first female",
                                                 "time point", "creating.type", "lowest ID", "highest ID")
        }
        rownames(population$info$cohorts) <- population$info$cohorts[,1]
      }


    } else{
      if(length(population)==0 || add.chromosome==TRUE){

        if(miraculix && miraculix.dataset && length(chr.opt)>1 && (snps.equidistant ||  sum(diff(position)<0) <= (length(chr.opt)-1)) ){

          total_snp <- sum(nsnp)

          dataset_full <- matrix(0L, ncol=nindi*2, nrow=total_snp)
          prior <- 0
          for(snp_index in 1:length(dataset)){


            dataset_temp <- as.matrix(dataset[[snp_index]])
            dataset_full[1:nrow(dataset_temp)+prior,] <- dataset_temp
            prior <- prior + nrow(dataset_temp)

          }

          dataset_full <- list(miraculix::haplomatrix(dataset_full))

        } else{
          dataset_full <- NULL
        }

        for(chr_index in 1:length(chr.opt)){
          index <- chr.opt[chr_index]
          activ <- which(chr.nr==index)
          chr_activ <- chr.nr[activ]
          bp_activ <- bp[activ]
          snp.name_activ <- snp.name[activ]
          hom0_activ <- hom0[activ]
          hom1_activ <- hom1[activ]
          if(is.list(dataset)){
            dataset_activ <- dataset[[chr_index]]
          } else{
            dataset_activ <- dataset[activ,,drop=FALSE]
          }

          snp.position_activ <- position[activ]
          if(length(freq)>1){
            freq_activ <- freq[activ]
          } else{
            freq_activ <- freq
          }

          if(chr_index==1){
            skip.rest_temp <- FALSE
            add.chromosome <- add.chromosome
          } else{
            skip.rest_temp <- TRUE
            add.chromosome <- TRUE
          }
          if(add.chromosome==FALSE){
            name.cohort <- name.cohort
          } else{
            name.cohort <- NULL
          }

          enter.bv_temp <- FALSE


          shuffle.traits_temp <- NULL
          shuffle.cor_temp <- NULL
          real.bv.add_temp <- NULL
          real.bv.mult_temp <- NULL
          real.bv.dice_temp <- NULL


          population <- creating.diploid(population=population, dataset=dataset_activ,
                                         nsnp=nsnp[chr_index], nindi=nindi,
                                         add.chromosome=add.chromosome, chr.nr = chr_activ,
                                         bp= bp_activ, snp.name = snp.name_activ,
                                         hom0= hom0_activ, hom1 = hom1_activ,
                                         class = class,
                                         generation = generation,
                                         add.chromosome.ends = add.chromosome.ends,
                                         miraculix = miraculix,
                                         miraculix.dataset = miraculix.dataset,
                                         snp.position = if(bpcm.conversion[chr_index]==0){snp.position_activ} else NULL,
                                         snps.equidistant= snps.equidistant,
                                         position.scaling= position.scaling,
                                         chromosome.length= chromosome.length[chr_index],
                                         length.before = length.before,
                                         length.behind = length.behind,
                                         skip.rest = skip.rest_temp,
                                         enter.bv = enter.bv_temp,
                                         sex.s = sex.s,
                                         genotyped.s = genotyped.s,
                                         name.cohort = name.cohort,
                                         real.bv.add = real.bv.add,
                                         real.bv.mult = real.bv.mult,
                                         real.bv.dice = real.bv.dice,
                                         base.bv = base.bv,
                                         bve.mult.factor = bve.mult.factor,
                                         bve.poly.factor = bve.poly.factor,
                                         freq = freq_activ,
                                         bpcm.conversion = bpcm.conversion[chr_index],
                                         remove.invalid.qtl=FALSE,
                                         shuffle.traits = shuffle.traits,
                                         shuffle.cor = shuffle.cor,
                                         verbose = verbose,
                                         internal=TRUE,
                                         progress.bar = progress.bar,
                                   internal.geno=if(chr_index == length(chr.opt) || !(miraculix && miraculix.dataset)){TRUE} else {FALSE},
                                   internal.dataset = dataset_full,
                                   size.scaling = size.scaling)
        }
      } else{

        if(sum(is.na(as.integer(chr.nr)))==0){
          tmp1 = min(diff(as.integer((chr.nr))))<0
        } else{
          tmp1 = min(diff(as.integer(as.factor(chr.nr))))<0
        }
        if(tmp1 || !miraculix.dataset){
          dataset_temp <- dataset
          till <- 0
          for(chr_index in 1:length(chr.opt)){
            index <- chr.opt[chr_index]
            activ <- which(chr.nr==index)
            chr_activ <- chr.nr[activ]
            bp_activ <- bp[activ]
            snp.name_activ <- snp.name[activ]
            hom0_activ <- hom0[activ]
            hom1_activ <- hom1[activ]
            dataset_activ <- dataset[activ,]


            snp.position_activ <- position[activ]
            if(length(activ)>0){
              dataset[1:length(activ)+till,] <- dataset_temp[activ,]
            }
            if((dataset[1+till,1]!=hom0_activ[1] && dataset[1+till,1]!=hom1_activ[1])){
              dataset[1:length(activ)+till,][dataset[1:length(activ)+till,]==hom0_activ] <- 0
              dataset[1:length(activ)+till,][dataset[1:length(activ)+till,]==hom1_activ] <- 1
            }
            till <- till + length(activ)

          }
        } else{
          if(!miraculix.dataset && (dataset[1,1]!=hom0[1] && dataset[1,1]!=hom1[1])  ){
            dataset[dataset==hom0_activ] <- 0
            dataset[dataset==hom1_activ] <- 1
          }
        }
        skip.rest <- TRUE

        population <- creating.diploid(population=population, dataset=dataset,
                                       class = class,
                                       generation = generation,
                                       add.chromosome.ends = add.chromosome.ends,
                                       miraculix = miraculix,
                                       skip.rest = skip.rest,
                                       sex.s = sex.s,
                                       genotyped.s = genotyped.s,
                                       name.cohort = name.cohort,
                                       real.bv.add = real.bv.add,
                                       real.bv.mult = real.bv.mult,
                                       real.bv.dice = real.bv.dice,
                                       base.bv = base.bv,
                                       bve.mult.factor = bve.mult.factor,
                                       bve.poly.factor = bve.poly.factor,
                                       hom0 =population$info$snp.base[1,],
                                       hom1 =population$info$snp.base[2,],
                                       remove.invalid.qtl =FALSE,
                                       shuffle.traits = shuffle.traits,
                                       shuffle.cor = shuffle.cor,
                                       verbose = verbose,
                                       internal=TRUE,
                                       progress.bar = progress.bar,
                                       size.scaling = size.scaling,
                                       founder.pool = founder.pool)
      }

    }

  }

  # Editing of the population list - Traits
  {
    if(store.comp.times){
      times_comp[5] <- as.numeric(Sys.time())
    }
    if(enter.bv ){
      if(bv.total>0 || length(real.bv.add)>0  || length(real.bv.mult) >0 || length(real.bv.dice)>0){
        population$info$bve <- TRUE
        if(is.list(real.bv.add)){
          population$info$real.bv.add <- real.bv.add
        } else{
          population$info$real.bv.add <- list(real.bv.add)
        }
        if(is.list(real.bv.mult)){
          population$info$real.bv.mult <- real.bv.mult
        } else{
          population$info$real.bv.mult <- list(real.bv.mult)
        }
        if(is.list(real.bv.dice)){
          population$info$real.bv.dice <- real.bv.dice
        } else{
          if(length(real.bv.dice)>0){
            warning("Invalid input for real.bv.dice!")
          }
          population$info$real.bv.dice <- list(real.bv.dice)
        }


        if(skip.rest==FALSE){
          population$info$bv.nr <- bv.total
          population$info$bv.calc <- bv.calc
        } else{
          bv.total = population$info$bv.nr
          bv.calc = population$info$bv.calc
          nbv = bv.calc
          population$info$bv.calculated = FALSE
        }


        population$info$real.bv.length <- c(length(population$info$real.bv.add),
                                            length(population$info$real.bv.mult),
                                            length(population$info$real.bv.dice))

        population$info$real.bv.add[[nbv+1]] <- "placeholder" # Use nbv instead of bv.calc
        population$info$real.bv.mult[[nbv+1]] <- "placeholder"
        population$info$real.bv.dice[[nbv+1]] <- "placeholder"



      } else if(preserve.bve){
        population$info$bve <- FALSE
        population$info$bv.nr <- 0
        population$info$bv.calc <- 0
        population$info$real.bv.length <- c(0,0,0)
      }

      if(bit.storing){
        population$info$bitstoring <- nbits
        population$info$leftover <-  sum(population$info$snp)%%nbits
      }
      if(miraculix || miraculix.dataset){
        population$info$miraculix <- TRUE
      } else{
        population$info$miraculix <- FALSE
      }

      if(bv.total>0 && (length(population$info$pheno.correlation)==0 || nrow(population$info$pheno.correlation)<bv.total)){
        population$info$pheno.correlation <- diag(1L, bv.total)
      }
      if(length(new.residual.correlation)>0){
        population$info$pheno.correlation <- t(chol(new.residual.correlation))
      }
      if(bv.total>0 && (length(population$info$bv.correlation)==0 || nrow(population$info$bv.correlation)<bv.total)){
        population$info$bv.correlation <- diag(1L, bv.total)
      }
      if(length(new.breeding.correlation)>0){
        population$info$bv.correlation <- new.breeding.correlation
      }

      if(bv.total>0 && (length(population$info$litter.effect.covariance)==0 || nrow(population$info$litter.effect.covariance)<bv.total)){
        population$info$litter.effect.covariance <- diag(0L, bv.total)
      }
      if(length(litter.effect.covariance)>0){
        population$info$litter.effect.covariance <- t(chol(litter.effect.covariance))
      }
      if(bv.total>0 && (length(population$info$pen.effect.covariance)==0 || nrow(population$info$pen.effect.covariance)<bv.total)){
        population$info$pen.effect.covariance <- diag(0L, bv.total)
      }
      if(length(pen.effect.covariance)>0){
        population$info$pen.effect.covariance <- t(chol(pen.effect.covariance))
      }

      if(length( population$info$litter.effect.covariance)==0 || sum( abs(population$info$litter.effect.covariance))==0){
        population$info$litter.effect.active <- FALSE
      } else{
        population$info$litter.effect.active <- TRUE
      }

      if(length( population$info$pen.effect.covariance)==0 || sum( abs(population$info$pen.effect.covariance))==0){
        population$info$pen.effect.active <- FALSE
      } else{
        population$info$pen.effect.active <- TRUE
      }



      if(bv.total>0){
        population$info$trait.name <- trait.name
        if(length(trait.name)<bv.total){
          population$info$trait.name <- c(population$info$trait.name, paste0("Trait ", (length(trait.name)+1):bv.total))
        }
      }


      if(length(shuffle.traits)==0){
        if(length(shuffle.cor)>0){

          if(ncol(shuffle.cor)==population$info$bv.calc){
            shuffle.traits <- 1:population$info$bv.calc
          } else{
            shuffle.traits <- 1:ncol(shuffle.cor)
            warning(paste0("shuffle.traits not specified! use the first ", ncol(shuffle.cor), " traits"))
          }
        }
      }
      if(length(shuffle.traits)>0 ){
        if(length(shuffle.traits)==1){
          shuffle.traits <- which(population$info$bv.random==FALSE)
        }
        # scaling of QTL effects
        if(!population$info$bv.calculated){
          if(use.recalculate.manual){

            ttt = length(population$info$bv.random.activ)==0
            if(ttt){
              population$info$bv.random.activ = shuffle.traits
            }
            population = recalculate.manual(population, cohorts = name.cohort_temp)
            if(ttt){
              population$info$bv.random.activ = NULL

            }
            population$info$bv.calculated = TRUE
          }
        }

        if(population$info$founder_multi){
          population$info$founder_multi_calc = FALSE
          for(index in 1:population$info$bv.nr){
            if(length(population$info$real.bv.add[[index]])>1 && sum(population$info$real.bv.add[[index]][,7]!=0)>0){
              population$info$founder_multi_calc = TRUE            }
          }
        }
        population <- breeding.diploid(population, verbose = FALSE)
        population$info$founder_multi_calc = population$info$founder_multi

        bvs <- get.bv(population, gen=1)
        scalings <- sqrt(diag(stats::var(t(bvs))))
        if(sum(is.na(scalings))>0){
          stop("scaling problems!")
        }
        for(bvnr in shuffle.traits){
          if(length(population$info$real.bv.add[[bvnr]])>0){
            population$info$real.bv.add[[bvnr]][,3:5] <- population$info$real.bv.add[[bvnr]][,3:5] / scalings[bvnr] * scalings[1]
          }


          if(length(population$info$real.bv.mult[[bvnr]])>0){
            population$info$real.bv.mult[[bvnr]][,5:13] <- population$info$real.bv.mult[[bvnr]][,5:13] / scalings[bvnr] * scalings[1]
          }

          if(length(population$info$real.bv.dice[[bvnr]])>0){
            population$info$real.bv.dice[[bvnr]][[2]] <- population$info$real.bv.dice[[bvnr]][[2]] / scalings[bvnr] * scalings[1]
          }


        }
        population$info$bv.calculated <- FALSE

        if(bit.storing){
          population$info$bitstoring <- nbits
          population$info$leftover <-  sum(population$info$snp)%%nbits
        }
        if(miraculix || miraculix.dataset){
          population$info$miraculix <- TRUE
        } else{
          population$info$miraculix <- FALSE
        }

        eigen_gen <- eigen(shuffle.cor)
        if(sum(eigen_gen$values<0)>0){
          if(verbose){
            warning("Genetic covariance matrix is not positive definit.")
            cat("Genetic covariance matrix is not positive definit.\n")
          }
          if(verbose) cat("Generate projection on the set of positive definit matrices:")

          test <- eigen_gen

          test$values[test$values<0] <- 0
          M <- diag(test$values)

          S <- test$vectors

          newA <- S %*% M %*% solve(S)

          diag(newA) <- diag(newA) + 0.0000001 # Avoid numerical issues with inversion
          newA <- newA * matrix(1/sqrt(diag(newA)), nrow=nrow(newA), ncol=nrow(newA), byrow=TRUE) * matrix(1/sqrt(diag(newA)), nrow=nrow(newA), ncol=nrow(newA), byrow=FALSE)
          if(verbose) cat("new suggested genetic correlation matrix:\n")
          shuffle.cor <- newA
          if(verbose) print(round(shuffle.cor, digits=3))
        }

        LT <- chol(shuffle.cor)
        if(nrow(LT)!=length(shuffle.traits)){
          stop("Dimension of shuffle correlation matrix doesnt work with traits to shuffle")
        } else{

          population$info$bv.correlation[shuffle.traits,shuffle.traits] <- t(LT) %*% LT
          if(sum(abs(population$info$bv.correlation[shuffle.traits,shuffle.traits]- shuffle.cor))>0.0001){
            warning("No covariance matrix for genetic correlation given! Values above diagonal used.")
          }

          store.add <- population$info$real.bv.add
          store.mult <- population$info$real.bv.mult
          store.dice <- population$info$real.bv.dice


          col <- 1
          for(index in shuffle.traits){
            new.add <- new.mult <- new.dice1 <- new.dice2 <- NULL
            row <- 1
            add.list = mult.list = list()
            for(index2 in shuffle.traits){
              if(length(store.add[[index2]])>0){

                if(LT[row,col]!=0){
                  temp1 <- store.add[[index2]] %*% diag(c(1,1,rep(LT[row,col],3),1,1,1))
                } else{
                  temp1 <- NULL
                }
                if(length(temp1)>0){
                  zeros <- rowSums(abs(temp1[,3:5, drop=FALSE]))
                  temp1 <- temp1[zeros>0,,drop=FALSE]
                }
                add.list[[index2]] = temp1

              }
              if(length(store.mult[[index2]])>0){

                if(LT[row,col]!=0){
                  temp1 <- store.mult[[index2]] %*% diag(c(1,1,1,1,rep(LT[row,col],9)))
                } else{
                  temp1 <- NULL
                }
                if(length(temp1)>0){
                  zeros <- rowSums(abs(temp1[,5:13, drop=FALSE]))
                  temp1 <- temp1[zeros>0,,drop=FALSE]
                }
                mult.list[[index2]] = temp1
              }
              if(length(store.dice[[index2]])>0){
                before <- length(new.dice2)
                new.dice1 <- c(new.dice1,store.dice[[index2]][[1]])
                new.dice2 <- c(new.dice2,store.dice[[index2]][[2]])
                for(index3 in (before+1):length(new.dice2)){
                  new.dice2[[index3]] <- new.dice2[[index3]] * LT[row,col]
                }
              }
              row <- row +1
            }

            new.add = do.call(rbind, add.list)
            new.mult = do.call(rbind, mult.list)

            # DONT REMOVE NULL - MORE WORK NEEDED HERE!
            if(length(new.add)==0){

            } else{
              population$info$real.bv.add[[index]] <- new.add
            }
            if(length(new.mult)==0){

            } else{
              population$info$real.bv.mult[[index]] <- new.mult
            }
            if(length(new.add)==0){

            } else{
              population$info$real.bv.dice[[index]] <- list(new.dice1,new.dice2)
            }
            col <- col +1
          }

        }

        for(index in shuffle.traits){
          population$info$real.bv.length[1] <- max(population$info$real.bv.length[1], if(length(population$info$real.bv.add[[index]])>0){index} else{0})
          population$info$real.bv.length[2] <- max(population$info$real.bv.length[2], if(length(population$info$real.bv.mult[[index]])>0){index} else{0})
          population$info$real.bv.length[3] <- max(population$info$real.bv.length[3], if(length(population$info$real.bv.dice[[index]][[1]])>0){index} else{0})
        }


      }

      if(population$info$bv.nr > 0){
        for(index in 1:population$info$bv.nr){
          if(length(population$info$real.bv.add) >= index){
            if(length(population$info$real.bv.add[[index]])>1){
              t <- population$info$real.bv.add[[index]]
              take <- sort(t[,1]+ cumsum(c(0,population$info$snp))[t[,2]], index.return=TRUE)
              t <- t[take$ix,,drop=FALSE]
              take <- sort(t[,1]+ t[,2] * 10^10 + t[,7]*10^8 + t[,8]*10^7)
              keep <- c(0,which(diff(take)!=0), length(take))
              if(length(keep) < (nrow(t)+1)){

                for(index2 in (2:(length(keep)))[diff(keep)>1]){
                  t[keep[index2],3:5] <- colSums(t[(keep[index2-1]+1):keep[index2],3:5, drop=FALSE])
                }
              }
              population$info$real.bv.add[[index]] <- t[keep,,drop=FALSE]
            }
          }

        }
      }


      if(length(add.architecture)>0){
        population$info$gen.architecture[[length(population$info$gen.architecture)+1]] <- list()
        population$info$gen.architecture[[length(population$info$gen.architecture)]]$length.total <- cumsum(c(0,add.architecture[[1]]))
        population$info$gen.architecture[[length(population$info$gen.architecture)]]$snp.position <- add.architecture[[2]]

      }


    }


    if(write_vcf_header){
      population$info$vcf_header <- list()
      population$info$vcf_header[[1]] <- tmp.header
    }



    if(length(population$info$real.bv.add)==0){
      population$info$real.bv.add <- list()
      population$info$real.bv.mult <- list()
      population$info$real.bv.dice <- list()
      population$info$real.bv.add[[1]] <- "placeholder" # Use nbv instead of bv.calc
      population$info$real.bv.mult[[1]] <- "placeholder"
      population$info$real.bv.dice[[1]] <- "placeholder"
    }
    if(remove.invalid.qtl && length(population$info$real.bv.add)>1){
      for(index in 1:(length(population$info$real.bv.add)-1)){
        removes <- which(population$info$real.bv.add[[index]][,1] > population$info$snp[population$info$real.bv.add[[index]][,2]])
        if(length(removes)>0){
          population$info$real.bv.add[[index]] <- population$info$real.bv.add[[index]][-removes,,drop=FALSE]
          warning(paste0(length(removes), " QTL-effects entered on markers that do not exist for ", population$info$trait.name[index], "."))
          warning(paste0(nrow(population$info$real.bv.add[[index]]), " QTL-effects remain."))
        }
      }
      for(index in 1:(length(population$info$real.bv.mult)-1)){
        removes <- which(population$info$real.bv.mult[[index]][,1] > population$info$snp[population$info$real.bv.mult[[index]][,2]])
        if(length(removes)>0){
          population$info$real.bv.mult[[index]] <- population$info$real.bv.mult[[index]][-removes,,drop=FALSE]
          warning(paste0(length(removes), " QTL-effects entered on markers that do not exist for ", population$info$trait.name[index], "."))
          warning(paste0(nrow(population$info$real.bv.mult[[index]]), " QTL-effects remain."))
        }
      }
      for(index in 1:(length(population$info$real.bv.mult)-1)){
        removes <- which(population$info$real.bv.mult[[index]][,3] > population$info$snp[population$info$real.bv.mult[[index]][,4]])
        if(length(removes)>0){
          population$info$real.bv.mult[[index]] <- population$info$real.bv.mult[[index]][-removes,,drop=FALSE]
          warning(paste0(length(removes), " QTL-effects entered on markers that do not exist for ", population$info$trait.name[index], "."))
          warning(paste0(nrow(population$info$real.bv.mult[[index]]), " QTL-effects remain."))
        }
      }

    }

    population$info$neff <- list()
    if(length(population$info$real.bv.add)>1){
      for(index in 1:(length(population$info$real.bv.add)-1)){
        if(length(population$info$real.bv.add[[index]])>0){
          population$info$neff[[index]] <- 1:nrow(population$info$real.bv.add[[index]])
        }
      }
    }




    if(bv.standard){
      population <- bv.standardization(population, mean.target = mean.target, var.target = var.target, set.zero = set.zero)

    }

  }


  if(sum(population$info$length)>1000){
    warning(paste0("Chromosome added a size of ", population$info$length[length(population$info$length)], " Morgan!
This will cost massiv computing time. Are you sure this is correct?
E.g. The entire human genome has a size of ~33 Morgan."))
  }

  if(store.comp.times){
    times_comp[6] <- as.numeric(Sys.time())
  }



  if(!internal){

    if(store.comp.times){
      add <- diff(times_comp)
      add = c(add, sum(add))
      names(add) <- c("Initialization", "Genotype Generation", "Preprocessing", "Population List genotypes", "Population list traits", "Total")
      population$info$comp.times.creating <- rbind( population$info$comp.times.creating, add)
    }

    if(length(population$info$creating.freq)==0){
      population$info$creating.freq <- freq
    }

    if(length(fixed.effects)==0){
      fixed.effects <- matrix(0, nrow= population$info$bv.nr, ncol=0)
    }
    population$info$fixed.effects <- fixed.effects


    if(sum(counter.start)>2 && population$info$bv.calculated){


      if(population$info$founder_multi_calc && ((length(population$info$founder_pools) + 1) > length(population$info$bypool_list[[1]]))){
        population = breeding.diploid(population)
      }

      activ_bv <- 1:population$info$bv.nr
      if(length(bv.ignore.traits)>0){
        activ_bv <- activ_bv[-bv.ignore.traits]
      }

      if(use.recalculate.manual){
        population = recalculate.manual(population, cohorts = name.cohort_temp, store.comp.times= store.comp.times)

#        if(store.comp.times){
#          population$info$comp.times.general[nrow(population$info$comp.times.general)-1, ] = population$info$comp.times.general[nrow(population$info$comp.times.general)-1, ] +
#            population$info$comp.times.general[nrow(population$info$comp.times.general), ]
#          population$info$comp.times.general = population$info$comp.times.general[-nrow(population$info$comp.times.general), ,drop = FALSE]
#        }

      } else{
        for(sex in 1:2){
          if(counter.start[sex]<counter[sex]){
            for(nr.animal in counter.start[sex]:(counter[sex]-1)){

              temp_out <- calculate.bv(population, generation, sex, nr.animal,
                                       decodeOriginsU=decodeOriginsU,
                                       output_compressed=FALSE,
                                       activ_bv = activ_bv,
                                       bv.ignore.traits=bv.ignore.traits)
              population$breeding[[generation]][[6+sex]][,nr.animal] <- temp_out[[1]]
              population$breeding[[generation]][[sex]][[nr.animal]][[25]] <- length(bv.ignore.traits)==0
            }


          }
        }
      }



    }

    if(use.recalculate.manual){
      population = recalculate.manual(population, cohorts = name.cohort_temp)
      population$info$bv.calculated = TRUE
    }

    population <- breeding.diploid(population, verbose = verbose)


    class(population) <- "population"


    if(length(trait_location)>0){
      population$info$trait.location = trait_location
      population$info$trait.nr = trait_nr
    }
    if(gxe.combine & length(trait_nr)>0){

      traits = unique(trait_nr)
      for(index in traits){
        population <- combine.traits(population, combine.traits = which(trait_nr == index))
      }
    }

  }

  population$info$max.time.point = max(population$info$max.time.point, time.point)

  if(one.sex.mode){
    population$info$one.sex.mode <- TRUE
  } else if(length(population$info$one.sex.mode)==0){
    population$info$one.sex.mode <- FALSE
  }

  return(population)
}
