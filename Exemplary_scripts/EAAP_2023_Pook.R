###################################################################
#### Evaluation of different breeding strategies (4.2 & 4.4) ######
###################################################################

## Using args to submit different strategies within the same shell script
args <- commandArgs(TRUE)
# args = c(1, 35)
scenario <- as.numeric(args[2])
seed <- as.numeric(args[1])

title = paste0("uniqueness_v3/sc", scenario, "seed", seed, "_potential.RData")
title2 = paste0("sc", scenario, "seed", seed, "_potential.RData")
aaa = sum(dir("uniqueness_v3") == title2)
if(aaa==1){
  stop()
}

library(MoBPS)

n_additive = 1000 # Number of additive QTLs
heritability = 0.3 # heritability


pop_size = c(500, 500) # number of males / females per cycle
sel_size = c(40, 100) # number of selected males / females per cycle

nsnp = 25000 # number of SNPs
chr = 10 # number of chromosomes
chr_length = 2.5 # length of each chromosome in Morgan

burnin = 10
cycles = 50

weight_unique = 0.05

display.progress = FALSE

## Generation of the baseline population

set.seed(seed)

population = creating.diploid(nsnp = nsnp, nindi = sel_size,
                              chr.nr = chr, chromosome.length = chr_length,
                              n.additive = n_additive,
                              share.genotyped = 1)

alpha_real = numeric(nsnp)
alpha_real[population$info$real.bv.add[[1]][,6]] = population$info$real.bv.add[[1]][,5] - population$info$real.bv.add[[1]][,3]

population = breeding.diploid(population, heritability = heritability)

for(index in 1:burnin){
  population = breeding.diploid(population, phenotyping.gen = get.ngen(population))
  population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                selection.criteria = "pheno",
                                share.genotyped = 1,
                                display.progress = display.progress)
}


population = bv.standardization(population, mean.target = 100, var.target = 1,
                                gen = get.ngen(population),
                                adapt.pheno = TRUE)
population = breeding.diploid(population, heritability = heritability)

{{

  {
    # Tobias helper function for Piters Index

    ixk <- function(p) {
      x<-qnorm(1-p)		#truncation point
      i<-dnorm(x)/p		#selection intensity
      k<-i*(i-x)		#variance reduction coefficient
      return(c(i,x,k))
    }

    get.gamete.var.Musa.method <- function(population,
                                           method=1,
                                           gen=NULL,
                                           database=NULL,
                                           cohorts=NULL,
                                           R.matrix=NULL,
                                           allele1.effect=NULL,
                                           verbose=FALSE){

      database <- get.database(population, gen = gen, database = database,
                               cohorts = cohorts)

      #database <- get.database(pop2, gen = 1)

      animal.id <- as.integer(get.id(population, database = database, use.id = TRUE))

      if(is.null(allele1.effect)){
        # allele substitution effects
        av.eff <- get.average.effect(population, database = database)
        meff <- data.frame(SNPName = names(av.eff), TRAIT = av.eff)
        rownames(meff) <- c()
      } else{
        meff <- allele1.effect
        row.names(meff) = meff[,1] # addition!!!
      }

      if(is.null(allele1.effect)){
        # phased genotypes
        haplo <- get.haplo(population,
                           database = get.database(population,
                                                   id = animal.id))[sort(get.qtl(population)),]
      } else{
        haplo <- get.haplo(population,
                           database = get.database(population,
                                                   id = animal.id))

        haplo <- haplo[rownames(haplo) %in% rownames(meff),]
      }
      # the calculation is the same as in chapter 2 of Musa's dissertation
      # this matrix indicates whether a locus in an animal is heterozygous
      hetero <- (haplo[,seq(from =1, to=ncol(haplo), by = 2)] + haplo[,seq(from =2, to=ncol(haplo), by = 2)]) == 1

      # if the beneficial allele(B or 1) is on the first haplotype, then return 1.
      # If not, return -1. If the position is homozygous, just return 0
      phase.indicator <- ((haplo[,seq(from =1, to=ncol(haplo), by = 2)]*2)-1)*hetero

      # this is setting up the marker effects for every individual individually
      marker.effect.indi <- (phase.indicator*meff[,2])

      if(is.null(R.matrix) & is.null(allele1.effect)){
        R.matrix <- get.R.Musa(population, R.as.list = TRUE, verbose=verbose)
      } else {
        map <- get.map(population)
        map <- map[map[,2] %in% meff[,1],]
        R.matrix <- get.R.Musa(population, R.as.list = TRUE, verbose = verbose, qtl.map = map)
      }

      if(is.list(R.matrix)){
        if(is.null(allele1.effect)){
          gmap <- as.data.frame(get.qtl.map(population)[,c(1,2,4,3)])
        } else{
          map <- get.map(population)
          map <- map[map[,2] %in% meff[,1],]
          gmap <- as.data.frame(map[,c(1,2,4,3)])
        }
        gmap[,1] <- as.integer(gmap[,1])

        last.snp.chrom <- c(1,table(gmap[,1]))

        df.msv.chrom <- matrix(NA, nrow = length(last.snp.chrom)-1, ncol = ncol(marker.effect.indi))

        for(i in 2:length(last.snp.chrom)){
          marker.effect.indi.chrom <- marker.effect.indi[sum(last.snp.chrom[1:(i-1)]):sum(last.snp.chrom[2:(i)]),,drop=FALSE]
          R.matrix.chrom <- R.matrix[[i-1]]

          #        msv.chrom <- sapply(1:ncol(marker.effect.indi.chrom),
          #                            function(indi) (marker.effect.indi.chrom[,indi]%*%
          #                                              (R.matrix.chrom%*%marker.effect.indi.chrom[,indi]))[1,1]
          #        )

          msv.chrom = diag(t(marker.effect.indi.chrom) %*% R.matrix.chrom %*% marker.effect.indi.chrom)

          df.msv.chrom[i-1,] <- msv.chrom
          if(verbose){
            prog = max(floor(i/(length(last.snp.chrom)-1)*50),1)
            cat('\r', "|", strrep("#",prog),
                strrep(" ", 50-prog), "| ", 2*prog, "%", sep="")
          }
        }
        flush.console()
        cat(" \n\r")
      } else{
        # this is multiplying the whole R matrix with the marker effect vectors
        # since the covariance between chromosoms is 0, I can also multiply R
        # per chromsome and then later add up the chromosome variances
        msv <- sapply(1:ncol(marker.effect.indi),
                      function(indi) (marker.effect.indi[,indi]%*%(R.matrix%*%marker.effect.indi[,indi]))[1,1]
        )
      }

      msv <- colSums(df.msv.chrom)
      msv <- rbind(msv)
      colnames(msv) <- animal.id
      return(msv)

      # msv <- colSums(df.msv.chrom)
      #
      # return(data.frame(ID = animal.id,
      #                   gam_var = msv))
    }

    get.R.Musa <- function(population, R.as.list=TRUE,verbose=TRUE,qtl.map=NULL){

      # genetic map
      if(is.null(qtl.map)){
        gmap <- as.data.frame(get.qtl.map(population)[,c(1,2,4,3)])
      } else{
        gmap <- as.data.frame(qtl.map[,c(1,2,4,3)])
      }
      colnames(gmap) <- c("CHR", "SNPName", "Position", "group1")
      gmap[,1] <- as.integer(gmap[,1])
      gmap[,3] <- as.integer(gmap[,3])
      gmap[,4] <- as.numeric(gmap[,4])
      #gmap[,4] <- gmap[,4]*100


      last.snp.chrom <- c(1,table(gmap[,1]))

      # chromosome wise

      prior = NULL
      for(i in 2:length(last.snp.chrom)){

        # i is chrom num +1
        gmap.chrom <- gmap[sum(last.snp.chrom[1:(i-1)]):sum(last.snp.chrom[2:(i)]),]
        # by splitting the operation per chromosome,
        # no unnecessary distance that would be 0 anyway are calculated
        myDist <- as.matrix(sapply(1:nrow(gmap.chrom), function(x) abs(gmap.chrom[x,4] - gmap.chrom[,4])))
        R.chrom <- exp(-2*myDist)/4
        if(R.as.list){
          if(i == 2){
            R <- list(R.chrom)
          } else{
            R <- c(R,list(R.chrom))
          }
        } else{
          # this process takes the longest
          if(i == 2){
            R <- Matrix::Matrix(R.chrom, sparse = T)
          } else{
            # binding matrices is faster than inserting values
            R <- Matrix::rbind2(
              Matrix::cbind2(R,Matrix::Matrix(data = 0, nrow = nrow(R), ncol = ncol(R.chrom), sparse = T)),
              Matrix::cbind2(Matrix::Matrix(data = 0, nrow = nrow(R.chrom), ncol = ncol(R), sparse = T),R.chrom))
          }
        }

        #R[sum(last.snp.chrom[1:(i-1)]):sum(last.snp.chrom[2:(i)]),
        #sum(last.snp.chrom[1:(i-1)]):sum(last.snp.chrom[2:(i)])] <- R.chrom
        if(verbose){
          print(paste0("Finished chromosome: ", names(last.snp.chrom)[i]))
        }
        # this was an earlier solution which does not work for large number of QTL:
        # myDist <- sapply(1:nrow(gmap), function(x) abs(gmap[x,4] - gmap[,4])) # this line is from Allier Usefullnes
        #
        # R <- exp(-2*myDist)/4
        # e <- c(1,table(gmap[,1]))
        #
        # for(i in 2:length(e)){
        #   R[sum(e[1:(i-1)]):sum(e[2:(i)]),-c(sum(e[1:(i-1)]):sum(e[2:(i)]))] <- 0
        # }

      }
      return(R)
      # The Matrix package makes storage space of matrices smaller
      # because R is a sparse matrix (has mostly 0)
      #pryr::object_size(R)
    }

  }

}}

for(index in 1:cycles){

  current_gen = get.ngen(population)

  population = breeding.diploid(population, phenotyping.gen = current_gen)


  bve.gen = current_gen
  if(scenario %in% 38:48){
    bve.gen = (current_gen-2):current_gen
  }
  population = breeding.diploid(population, bve = TRUE, bve.gen = bve.gen, estimate.u = TRUE)

  alpha_hat = population$info$u_hat[[length(population$info$u_hat)]][,1] # estimated SNP effects

  if(sum(is.na(alpha_hat))>0){
    alpha_hat[is.na(alpha_hat)] = 0
  }
  geno = get.geno(population = population, gen = current_gen) # genotypes

  p = rowMeans(geno)/2

  # minor allele frequency
  maf = p
  maf[maf>0.5] = 1 - maf[maf>0.5]

  # frequency of "good" allele
  freq_alpha = p
  freq_alpha[alpha_real==0] = 0 # these can be ignored
  freq_alpha[alpha_real<0] = 1 -freq_alpha[alpha_real<0] # these can be ignored

  # frequency of estimated "good" allele
  freq_alphahat = p
  freq_alphahat[alpha_hat==0] = 0 # these can be ignored
  freq_alphahat[alpha_hat<0] = 1 -freq_alphahat[alpha_hat<0] # these can be ignored

  ebv = get.bve(population, gen = current_gen)[1,] # estimated breeding value
  pheno = get.pheno(population, gen = current_gen)[1,] # phenotype
  bv = get.bv(population, gen = current_gen)[1,] # true underlying genomic value


  if( scenario == 1){

    # Selection based on genomic breeding value
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if( scenario == 2){

    # Selection based on genomic breeding value
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bv",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if( scenario == 3){

    # Selection based on genomic breeding value
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "pheno",
                                  share.genotyped = 1,
                                  display.progress = display.progress)


  } else if(scenario == 4){

    activ = alpha_hat!= 0

    weighting = sqrt(1 / freq_alphahat)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 5){

    activ = alpha_hat!= 0

    weighting = sqrt(1 / freq_alphahat)
    weighting[weighting==Inf] = 0

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 6){

    activ = alpha_hat!= 0

    weighting = (1 / freq_alphahat)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 7){

    activ = alpha_hat!= 0

    weighting = (1 / freq_alphahat)^(1/3)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 8){

    activ = alpha_hat!= 0

    weighting = (1 - freq_alphahat)
    weighting[weighting==Inf] = 0
    weighting[weighting<0.1] = 0.1

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 9){

    activ = alpha_hat!= 0

    weighting = (1 - freq_alphahat)^(2)
    weighting[weighting==Inf] = 0
    weighting[weighting<0.1] = 0.1

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if( scenario == 10){


    cand = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                            selection.criteria = "bve", export.selected = TRUE )

    A = kinship.exp(population, gen = get.ngen(population))

    rel_top = colMeans(A[cand[[1]][,3],]) + colMeans(A[cand[[2]][,3] + pop_size[1],])

    unique = -rel_top #
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)
  } else if( scenario == 11){


    cand = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                            selection.criteria = "bve", export.selected = TRUE )

    A = kinship.exp(population, gen = get.ngen(population))
    diag(A) = (diag(A)-0.5)*2

    rel_top = colMeans(A[cand[[1]][,3],]) + colMeans(A[cand[[2]][,3] + pop_size[1],])

    unique = -rel_top #
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)
  } else if( scenario == 12){


    cand = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                            selection.criteria = "bve", export.selected = TRUE )

    A = kinship.exp(population, gen = get.ngen(population))
    unique = -diag(A)
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)
  } else if( scenario == 13){


    cand = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                            selection.criteria = "bve", export.selected = TRUE )

    unique = -inbreeding.emp(population, gen = get.ngen(population))
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)
  } else if( scenario == 14){

    weight_unique = 0.2

    cand = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                            selection.criteria = "bve", export.selected = TRUE )

    A = kinship.exp(population, gen = get.ngen(population))

    rel_top = colMeans(A[cand[[1]][,3],]) + colMeans(A[cand[[2]][,3] + pop_size[1],])

    unique = -rel_top #
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)
  } else if( scenario == 15){

    weight_unique = 0.2

    cand = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                            selection.criteria = "bve", export.selected = TRUE )

    A = kinship.exp(population, gen = get.ngen(population))
    diag(A) = (diag(A)-0.5)*2

    rel_top = colMeans(A[cand[[1]][,3],]) + colMeans(A[cand[[2]][,3] + pop_size[1],])

    unique = -rel_top #
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)
  } else if( scenario == 16){

    weight_unique = 0.2

    cand = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                            selection.criteria = "bve", export.selected = TRUE )

    A = kinship.exp(population, gen = get.ngen(population))
    unique = -diag(A)
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)
  } else if( scenario == 17){

    weight_unique = 0.2

    cand = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                            selection.criteria = "bve", export.selected = TRUE )

    unique = -inbreeding.emp(population, gen = get.ngen(population))
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)
  } else if( scenario == 18){
    # all rare alleles with a positive effect get additional weighting
    activ = (freq_alphahat < 0.2) & (freq_alphahat > 0)
    unique = as.numeric(alpha_hat[activ] %*% geno[activ,])
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if( scenario == 19){
    alpha_rare = numeric(nsnp)
    alpha_rare[p < 0.2 & p!=0] = 1
    alpha_rare[p > 0.8 & p!=1] = (-1)
    activ = alpha_rare != 0
    unique = as.numeric(alpha_rare[activ] %*% geno[activ,])
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)
  }else if( scenario == 20){
    weight_unique = 0.2
    # all rare alleles with a positive effect get additional weighting
    activ = (freq_alphahat < 0.2) & (freq_alphahat > 0)
    unique = as.numeric(alpha_hat[activ] %*% geno[activ,])
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if( scenario == 21){
    weight_unique = 0.2
    alpha_rare = numeric(nsnp)
    alpha_rare[p < 0.2 & p!=0] = 1
    alpha_rare[p > 0.8 & p!=1] = (-1)
    activ = alpha_rare != 0
    unique = as.numeric(alpha_rare[activ] %*% geno[activ,])
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)
  } else if(scenario == 22){

    activ = alpha_real!= 0

    weighting = sqrt(1 / freq_alpha)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_real[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 23){

    activ = alpha_real!= 0

    weighting = sqrt(1 / freq_alpha)
    weighting[weighting==Inf] = 0

    unique = as.numeric((alpha_real[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 24){

    activ = alpha_real!= 0

    weighting = (1 / freq_alpha)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_real[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 25){

    activ = alpha_real!= 0

    weighting = (1 / freq_alpha)^(1/3)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_real[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 26){

    activ = alpha_real!= 0

    weighting = (1 - freq_alpha)
    weighting[weighting==Inf] = 0
    weighting[weighting<0.1] = 0.1

    unique = as.numeric((alpha_real[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 27){

    activ = alpha_real!= 0

    weighting = (1 - freq_alpha)^(2)
    weighting[weighting==Inf] = 0
    weighting[weighting<0.1] = 0.1

    unique = as.numeric((alpha_real[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if( scenario == 28){
    # all rare alleles with a positive effect get additional weighting
    unique = rnorm(nrow(geno))
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if( scenario == 29){
    weight_unique = 0.2
    # all rare alleles with a positive effect get additional weighting
    unique = rnorm(nrow(geno))
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario ==30){

    db.males <- get.database(population, gen = get.ngen(population))[1,]
    db.females <- get.database(population, gen = get.ngen(population))[2,]

    x.males <- ixk(sel_size[1]/(pop_size[1]))[2]
    x.females <- ixk(sel_size[2]/(pop_size[2]))[2]
    x.non.sex.specific <- (x.males + x.females)/2

    map <- get.map(population)
    allele1.effect <- data.frame(SNPName = map[,2], TRAIT = alpha_hat)

    gamVar.males <- get.gamete.var.Musa.method(population = population,
                                               database = db.males,
                                               allele1.effect = allele1.effect)

    gamVar.females <- get.gamete.var.Musa.method(population = population,
                                                 database = db.females,
                                                 allele1.effect = allele1.effect)

    index5.males <- sqrt(2) * x.non.sex.specific * sqrt(gamVar.males)
    index5.females <-  sqrt(2) * x.non.sex.specific * sqrt(gamVar.females)

    new_ebv =  ebv + c(index5.males, index5.females)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if( scenario == 31){
    # all rare alleles with a positive effect get additional weighting
    thres = quantile(abs(alpha_hat), probs = 0.8)
    activ = (freq_alphahat < 0.2) & (freq_alphahat > 0) & (abs(alpha_hat) > thres)
    unique = as.numeric(alpha_hat[activ] %*% geno[activ,])
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if( scenario == 32){
    weight_unique = 0.2
    # all rare alleles with a positive effect get additional weighting
    thres = quantile(abs(alpha_hat), probs = 0.8)
    activ = (freq_alphahat < 0.2) & (freq_alphahat > 0) & (abs(alpha_hat) > thres)
    unique = as.numeric(alpha_hat[activ] %*% geno[activ,])
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 33){

    activ = alpha_hat!= 0

    weighting = (1 / freq_alphahat)^(1/3)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    weight_unique = 0.2

    cand = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                            selection.criteria = "bve", export.selected = TRUE )

    A = kinship.exp(population, gen = get.ngen(population))

    rel_top = colMeans(A[cand[[1]][,3],]) + colMeans(A[cand[[2]][,3] + pop_size[1],])

    unique = -rel_top #
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (new_ebv - mean(new_ebv))/sd(new_ebv)


    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 34){

    activ = alpha_hat!= 0

    # not activ!
    weighting = rep(1, length(alpha_hat))
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])

    s1 = mean(unique)
    s2 = sd(unique)

    new_ebv = (unique - mean(unique))/sd(unique)

    weight_unique = 0.2

    cand = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                            selection.criteria = "bve", export.selected = TRUE )

    A = kinship.exp(population, gen = get.ngen(population))

    rel_top = colMeans(A[cand[[1]][,3],]) + colMeans(A[cand[[2]][,3] + pop_size[1],])

    unique = -rel_top #
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (new_ebv - mean(new_ebv))/sd(new_ebv)


    db.males <- get.database(population, gen = get.ngen(population))[1,]
    db.females <- get.database(population, gen = get.ngen(population))[2,]

    x.males <- ixk(sel_size[1]/(pop_size[1]))[2]
    x.females <- ixk(sel_size[2]/(pop_size[2]))[2]
    x.non.sex.specific <- (x.males + x.females)/2

    map <- get.map(population)
    allele1.effect <- data.frame(SNPName = map[,2], TRAIT = alpha_hat)

    gamVar.males <- get.gamete.var.Musa.method(population = population,
                                               database = db.males,
                                               allele1.effect = allele1.effect)

    gamVar.females <- get.gamete.var.Musa.method(population = population,
                                                 database = db.females,
                                                 allele1.effect = allele1.effect)

    index5.males <- sqrt(2) * x.non.sex.specific * sqrt(gamVar.males)
    index5.females <-  sqrt(2) * x.non.sex.specific * sqrt(gamVar.females)

    new_ebv =  new_ebv + c(index5.males, index5.females)/s2

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 35){

    activ = alpha_hat!= 0

    weighting = (1 / freq_alphahat)^(1/3)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])

    s1 = mean(unique)
    s2 = sd(unique)

    new_ebv = (unique - mean(unique))/sd(unique)

    weight_unique = 0.2

    cand = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                            selection.criteria = "bve", export.selected = TRUE )

    A = kinship.exp(population, gen = get.ngen(population))

    rel_top = colMeans(A[cand[[1]][,3],]) + colMeans(A[cand[[2]][,3] + pop_size[1],])

    unique = -rel_top #
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (new_ebv - mean(new_ebv))/sd(new_ebv)


    db.males <- get.database(population, gen = get.ngen(population))[1,]
    db.females <- get.database(population, gen = get.ngen(population))[2,]

    x.males <- ixk(sel_size[1]/(pop_size[1]))[2]
    x.females <- ixk(sel_size[2]/(pop_size[2]))[2]
    x.non.sex.specific <- (x.males + x.females)/2

    map <- get.map(population)
    allele1.effect <- data.frame(SNPName = map[,2], TRAIT = alpha_hat)

    gamVar.males <- get.gamete.var.Musa.method(population = population,
                                               database = db.males,
                                               allele1.effect = allele1.effect)

    gamVar.females <- get.gamete.var.Musa.method(population = population,
                                                 database = db.females,
                                                 allele1.effect = allele1.effect)

    index5.males <- sqrt(2) * x.non.sex.specific * sqrt(gamVar.males)
    index5.females <-  sqrt(2) * x.non.sex.specific * sqrt(gamVar.females)

    new_ebv =  new_ebv + c(index5.males, index5.females)/s2

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if( scenario == 36){

    cand = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                            selection.criteria = "bve", export.selected = TRUE )

    A = kinship.exp(population, gen = get.ngen(population))

    rel_top = colMeans(A)

    unique = -rel_top #
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)
  } else if( scenario == 37){

    weight_unique = 0.2

    cand = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                            selection.criteria = "bve", export.selected = TRUE )

    A = kinship.exp(population, gen = get.ngen(population))

    rel_top = colMeans(A)

    unique = -rel_top #
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)
  } else   if( scenario == 38){

    # Selection based on genomic breeding value
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 39){

    activ = alpha_hat!= 0

    weighting = sqrt(1 / freq_alphahat)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 40){

    activ = alpha_hat!= 0

    weighting = sqrt(1 / freq_alphahat)
    weighting[weighting==Inf] = 0

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 41){

    activ = alpha_hat!= 0

    weighting = (1 / freq_alphahat)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 42){

    activ = alpha_hat!= 0

    weighting = (1 / freq_alphahat)^(1/3)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 43){

    activ = alpha_hat!= 0

    weighting = (1 - freq_alphahat)
    weighting[weighting==Inf] = 0
    weighting[weighting<0.1] = 0.1

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 44){

    activ = alpha_hat!= 0

    weighting = (1 - freq_alphahat)^(2)
    weighting[weighting==Inf] = 0
    weighting[weighting<0.1] = 0.1

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if( scenario == 45){
    # all rare alleles with a positive effect get additional weighting
    activ = (freq_alphahat < 0.2) & (freq_alphahat > 0)
    unique = as.numeric(alpha_hat[activ] %*% geno[activ,])
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if( scenario == 46){
    weight_unique = 0.2
    # all rare alleles with a positive effect get additional weighting
    activ = (freq_alphahat < 0.2) & (freq_alphahat > 0)
    unique = as.numeric(alpha_hat[activ] %*% geno[activ,])
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if( scenario == 47){
    # all rare alleles with a positive effect get additional weighting
    thres = quantile(abs(alpha_hat), probs = 0.8)
    activ = (freq_alphahat < 0.2) & (freq_alphahat > 0) & (abs(alpha_hat) > thres)
    unique = as.numeric(alpha_hat[activ] %*% geno[activ,])
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if( scenario == 48){
    weight_unique = 0.2
    # all rare alleles with a positive effect get additional weighting
    thres = quantile(abs(alpha_hat), probs = 0.8)
    activ = (freq_alphahat < 0.2) & (freq_alphahat > 0) & (abs(alpha_hat) > thres)
    unique = as.numeric(alpha_hat[activ] %*% geno[activ,])
    new_ebv = weight_unique * (unique - mean(unique))/sd(unique) + (1-weight_unique)* (ebv - mean(ebv))/sd(ebv)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 49){

    activ = alpha_hat!= 0

    weighting = (1 / freq_alphahat)^(1/2.5)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  } else if(scenario == 50){

    activ = alpha_hat!= 0

    weighting = (1 / freq_alphahat)^(1/4)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_hat[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  }  else if(scenario == 51){

    activ = alpha_real!= 0

    weighting = (1 / freq_alpha)^(1/2.5)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_real[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  }  else if(scenario == 52){

    activ = alpha_real!= 0

    weighting = (1 / freq_alpha)^(1/4)
    weighting[weighting==Inf] = 0
    weighting[weighting>5] = 5

    unique = as.numeric((alpha_real[activ] * weighting[activ]) %*% geno[activ,])
    new_ebv = (unique - mean(unique))/sd(unique)

    population = insert.bve(population, bves = cbind(names(ebv), new_ebv))
    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.criteria = "bve",
                                  share.genotyped = 1,
                                  display.progress = display.progress)

  }
}


#### Evaluation of the population

if(TRUE){
  ngen = get.ngen(population)

  bv = numeric(ngen)
  bv_sd = numeric(ngen)
  inb = numeric(ngen)
  share_fixed = numeric(ngen)
  share_fixed_qtl = numeric(ngen)
  share_fixed_good = numeric(ngen)
  share_fixed_bad = numeric(ngen)
  for(index in 1:ngen){
    print(index)
    bv[index] = mean(get.bv(population, gen = index))
    bv_sd[index] = sd(get.bv(population, gen = index))
    inb[index] = mean(inbreeding.emp(population, gen = index))

    geno = get.geno(population, gen = index)
    p = rowMeans(geno)/2
    share_fixed[index] = mean( p ==0 | p ==1)
    share_fixed_qtl[index] = mean( (p ==0 | p ==1)[population$info$effect.p])

    share_fixed_good[index] = mean(((alpha_real > 0 & p ==1) | (alpha_real < 0 & p ==0))[population$info$effect.p])
    share_fixed_bad[index] = mean(((alpha_real < 0 & p ==1) | (alpha_real > 0 & p ==0))[population$info$effect.p])
  }

  save(file = paste0("uniqueness_v3/sc", scenario, "seed", seed, ".RData"),
       list = c("bv", "bv_sd", "inb", "share_fixed", "share_fixed_qtl", "share_fixed_good", "share_fixed_bad"))

}

if(TRUE){
  # Generate offspring from the best animals based on traditional EBVs
  for(index in 1:cycles){
    print(index)
    current_gen = get.ngen(population) - cycles

    sel_crit = "bve"
    if(scenario %in% c(2,22:27,51,52)){
      sel_crit = "bv"
    }

    bve.gen = current_gen
    if(scenario %in% 38:48){
      bve.gen = (current_gen-2):current_gen
    }
    population = breeding.diploid(population, bve = TRUE, bve.gen = bve.gen)

    population = breeding.diploid(population, selection.size = sel_size, breeding.size = pop_size,
                                  selection.m.database = cbind(current_gen,1),
                                  selection.f.database = cbind(current_gen,2),
                                  selection.criteria = sel_crit,
                                  share.genotyped = 1,
                                  display.progress = display.progress)
  }

  #### Evaluation of the population
  ngen = get.ngen(population)

  bv = numeric(ngen)
  bv_sd = numeric(ngen)
  inb = numeric(ngen)
  share_fixed = numeric(ngen)
  share_fixed_qtl = numeric(ngen)
  share_fixed_good = numeric(ngen)
  share_fixed_bad = numeric(ngen)
  for(index in 1:ngen){
    print(index)
    bv[index] = mean(get.bv(population, gen = index))
    bv_sd[index] = sd(get.bv(population, gen = index))
    inb[index] = mean(inbreeding.emp(population, gen = index))

    geno = get.geno(population, gen = index)
    p = rowMeans(geno)/2
    share_fixed[index] = mean( p ==0 | p ==1)
    share_fixed_qtl[index] = mean( (p ==0 | p ==1)[population$info$effect.p])

    share_fixed_good[index] = mean(((alpha_real > 0 & p ==1) | (alpha_real < 0 & p ==0))[population$info$effect.p])
    share_fixed_bad[index] = mean(((alpha_real < 0 & p ==1) | (alpha_real > 0 & p ==0))[population$info$effect.p])
  }


  save(file = title,
       list = c("bv", "bv_sd", "inb", "share_fixed", "share_fixed_qtl", "share_fixed_good", "share_fixed_bad"))

}
