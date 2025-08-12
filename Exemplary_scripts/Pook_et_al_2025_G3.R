##############################################################
##### Simulation of a maize rapid cycle breeding scheme ######
##############################################################

## MoBPS 1.9.06 was used for these simulations

# Initialization of variables
## potential parameterization:
## args <- c(5,  10 , 0 ,  0 ,  0 , 30,  30 ,  1 , 10,   2, TRUE, FALSE, TRUE, 0.2)

args <- commandArgs(TRUE)
nr <- as.numeric(args[1])
max_parent <- as.numeric(args[2])
targetdegree <- as.numeric(args[3])
replacer <- as.logical(args[4])
extra_mating <- as.logical(args[5]) # Additional mating before C1F2 for more recombination
nsel <- as.numeric(args[6]) # number of selected lines per cycle
cycle <- as.numeric(args[7]) # number of cycles performed
founder_pheno <- as.numeric(args[8]) # share of DHs phenotyped
dh_sel <- as.numeric(args[9]) # number of original DHs selected
sel_index <- as.numeric(args[10]) # Index Type to use
add_pheno <- as.logical(args[11]) # phenotype additional DHs
tc_sel <- as.logical(args[12]) # phenotype crosses
backcross <- as.logical(args[13]) # reintroduce variance after 5 cycles
share_backcross <- as.numeric(args[14]) # share reintroduce variance after 5 cycles

miraculix <- TRUE

{
sim_path <- getwd()


store <- NULL
if(is.na(share_backcross)){
  share_backcross <- 0
}
  if(!backcross){
    share_backcross <- 0
  }
print(c(nr, max_parent, targetdegree, replacer, extra_mating, nsel, cycle, founder_pheno, dh_sel, sel_index, add_pheno, tc_sel, backcross))

if(tc_sel){
  activ_t <- 4:6
} else{
  activ_t <- 1:3
}

nrep <- ceiling(1035 / nsel * 2)

dh_sel_cross <- dh_sel * (dh_sel-1) / 2
dh_cross <- ceiling(1035 / dh_sel_cross)

set.seed(nr)

library(MoBPS)
if(miraculix){
  library(miraculix)
  library(RandomFieldsUtils)
}

RandomFieldsUtils::RFoptions(cores = 3)

# opt_selection courtesy of Armin Hoelker (formerly TUM, now KWS)
opt_selection <- function(data,parents,usage_thresh,n_select,direction="max", const.dir="<="){
  if (!("package:lpSolve" %in% search())) stop("load lpSolve package")
  n_data <- nrow(data)
  mat_parents <- NULL

  for (p_i in 1:length(parents)) {
    P_vec <- ifelse(data[,3]==parents[p_i]|data[,4]==parents[p_i],1,0)
    mat_parents <- rbind(mat_parents,P_vec)
  }

  sum <- rep(1,n_data)


  ind.const <- diag(n_data)
  const.mat <- rbind(mat_parents,sum,ind.const)

  const.rhs <- c(rep(usage_thresh,length(parents)),n_select,rep(1,n_data))

  if(length(const.dir)==1) const.dir <- rep(const.dir,n_data+length(parents)+1)



  optimum <- lp(direction=direction,objective.in = data[,1],const.mat=const.mat,const.rhs=const.rhs,const.dir=const.dir)

  if (any(round(optimum$solution)>1)) stop("Individuals used more than once.")
  print(paste0("If status is 0, optimum was found, if status is 2 there is no solution. Status: ", optimum$status))

  if (round(optimum$status) == 2) stop("No optimum was found.")
  selected <- data[c(which(round(optimum$solution)==1)),]

  return(selected)
}

# this selection strategy is nowadays also directely implemented in MoBPS
# selection will later be performed in a wrapper of the MoBPS framework
# this was not used in the scenarios used in the paper
miesenberger.index <- function(V, G, V1=NULL, RG=NULL, r, w, zw=NULL){

  if(length(V1)==0){
    V1 <- chol2inv(chol(V))
  }
  if(length(RG)==0){
    RG <- sqrt(diag(1/diag(G))) %*% G %*% sqrt(diag(1/diag(G)))
  }

  d <- nrow(G)
  A <- RG * matrix(sqrt(diag(G)), nrow=d, ncol=d, byrow=TRUE) *
    matrix(sqrt(diag(V)), nrow=d, ncol=d, byrow=FALSE) *
    matrix(r, nrow=d, ncol=d, byrow=FALSE)
  bM <- as.numeric(V1 %*% A %*% w)
  if(length(zw)==0){
    return(bM)
  } else{
    IM <- sum( bM * zw )
    return(IM)
  }
}

geno <- NULL
nsnp <- NULL
bp <- NULL
sommer <- FALSE

nmatings <- nsel/2
nparents <- nmatings*2
seed_val <- 1

load("PE_DH_SNP.RData")
load("PE_map.RData")

storage.mode(datas) <- "integer"
geno <- datas

map[,3] <- NA
geno2 <- geno[,sort(rep(1:ncol(geno),2))]
rm(geno)
rm(datas)
cor_perse <- matrix(c(1, 0.94, 0.50,0.94,1,0.58,0.50,0.58,1), nrow=3)
cor_tc<- matrix(c(1, 0.86, 0.21,0.86,1,0.33,0.21,0.33,1), nrow=3)
res <- matrix(c(1, 0.3, 0.1,
                0.3,1,0.15,
                0.1,0.15,1), nrow=3, byrow=TRUE)

cor_perse_tc<- c(0.62, 0.68, 0.78)
res_perse_tc <- 0.2

heritability <- c(0.95,0.95,0.96,0.76,0.77,0.87)
cor <- rbind(cbind(cor_perse, sqrt(cor_perse*cor_tc)*cor_perse_tc), cbind(t(sqrt(cor_perse*cor_tc)*cor_perse_tc), cor_tc))
res <- rbind(cbind(res, res*res_perse_tc), cbind(res*res_perse_tc, res))

pheno_perse <- c(1,1,1,0,0,0)
pheno_tc <- c(0,0,0,1,1,1)
index_perse <- c(1,1,-1,0,0,0)
index_tc <- c(0,0,0,1,1,-1)

# Taken from Su et al. 2017 ## https://www.frontiersin.org/articles/10.3389/fpls.2017.00706/full
chromosome.length <- c(3.0926, 3.8280, 2.2878, 2.7934, 1.8979, 1.5386, 1.8530, 1.6792, 2.0009, 1.3951)


# 50k array is sufficent
keep <- round(1:50000*(nrow(geno2)/50000))
geno2 <- geno2[keep,]
map <- map[keep,]
geno2[geno2==2] <- 1

# Generation of the Founder population with 402 DH-lines
population <- creating.diploid(dataset = geno2,  chromosome.length = chromosome.length,
                               map = map,
                               miraculix = miraculix,
                               sex.quota = 0, name.cohort = "DH-Founder",
                               n.additive = rep(500,6),
                               n.equal.dominant = rep(500,6),
                               n.qualitative = rep(100,6),
                               n.quantitative = rep(100,6),
                               shuffle.traits = TRUE, shuffle.cor = cor,
                               new.residual.correlation = res)

population <- bv.standardization(population, mean.target = 100, var.target = 10, gen=1)

temp1 <- population$info$real.bv.add
temp2 <- population$info$real.bv.mult
temp1[[7]] <- NULL
temp2[[7]] <- NULL

# initialize a second population for introduction of genetic diversity
pop2 <- creating.diploid(dataset = geno2, map=map,
                         chromosome.length = chromosome.length, sex.quota = 0, name.cohort = "DH-Founder",
                 real.bv.add = temp1,
                 miraculix = miraculix,
                 real.bv.mult = temp2,
                 base.bv = population$info$base.bv, verbose=FALSE)

# Generation of Phenotypes
population <- breeding.diploid(population, new.bv.observation.gen = 1, n.observation = pheno_perse, share.phenotyped= founder_pheno, heritability = heritability)

# Generation of C1F1
population <- breeding.diploid(population, bve=TRUE, rrblup.bve = TRUE, selection.size = c(dh_sel,0), breeding.all.combination=TRUE,
                               multiple.bve.weights.m = c(1,1,-1,0,0,0), selection.criteria = "bve",
                               breeding.sex=0, name.cohort = "C1F1", selection.m = "function",
                               estimate.reliability=TRUE, verbose=TRUE, display.progress = FALSE)

population <- breeding.diploid(population, breeding.size = c(1035,0),
                               selection.m.cohorts = "DH-Founder",
                               name.cohort = "C0F2", verbose=TRUE, display.progress = FALSE)

population <- breeding.diploid(population, breeding.size = c(200,0),
                               selection.m.cohorts = "C1F1",
                               dh.mating = TRUE, dh.sex = 0,
                               name.cohort = "C1-DH", verbose=TRUE, display.progress = FALSE
                               )
# Additional phenotyping of DH-lines
if(add_pheno){
  population <- breeding.diploid(population, new.bv.observation.cohorts = "C1-DH" , share.phenotyped= founder_pheno, n.observation = pheno_perse)
}


if(!extra_mating){
  # Generation of C1F2
  population <- breeding.diploid(population, selection.size=c(dh_sel_cross,0), breeding.size = c(1035,0),
                                 repeat.mating = dh_cross, max.offspring = 1, selfing.mating = TRUE, selfing.sex = 0,
                                 selection.m.cohorts = "C1F1", name.cohort = "C1F2", verbose=TRUE, display.progress = FALSE)
} else{

  # Generation of C1F2 (with an additional step to generate further recombination) 
  # no results reported in the paper as it had very limited impact
  population <- breeding.diploid(population, selection.size=c(dh_sel_cross,0), breeding.size = c(dh_sel_cross,0),
                                 repeat.mating = 1, max.offspring = 2, selfing.mating = FALSE,
                                 selection.m.cohorts = "C1F1", name.cohort = "C1F1.5", verbose=TRUE, display.progress = FALSE)

  population <- breeding.diploid(population, selection.size=c(dh_sel_cross,0), breeding.size = c(1035,0),
                                 repeat.mating = dh_cross, max.offspring = 1, selfing.mating = TRUE, selfing.sex = 0,
                                 selection.m.cohorts = "C1F1.5", name.cohort = "C1F2", verbose=TRUE, display.progress = FALSE)


}

# Additional phenotyping of test crosses
if(tc_sel){
  population <- breeding.diploid(population, new.bv.observation.cohorts = "C1F2" , share.phenotyped= founder_pheno, n.observation = pheno_tc)
}


# Breeding value estimation for C1F2
population <- breeding.diploid(population, bve=TRUE, bve.cohorts = c("C1F2", "DH-Founder"),
                               bve.insert.cohorts = "C1F2", rrblup.bve = TRUE, estimate.reliability = TRUE, bve.ignore.traits = 4:6)

## only sel_index = 2 was considering in the paper

if(sel_index==1){
  # Currently used Index in Rapid Cycle
  bves <- get.bve(population, cohorts = "C1F2")[1:3,]
  bves[3,] <- -abs(bves[3,] - mean(bves[3,]))
  selection_index <- colSums(bves)
} else if(sel_index==2){
  # Negative weight on plant height final
  bves <- get.bve(population, cohorts = "C1F2")[1:3,]
  selection_index <- colSums(bves*(c(1,1,-1)))
} else {
  # Reliability based index
  bves <- get.bve(population, cohorts = "C1F2")[1:3, ]
  reli <- get.reliabilities(population, cohorts = "C1F2")[1:3, ]
  bve.cov <- var(t(bves))

  V <- bve.cov
  V1 <- MASS::ginv(V)

  bv <- get.bv(population, cohorts = "C1F2")[1:3,]
  G_cov <- var(t(bv))
  RG <- sqrt(diag(1/diag(G_cov))) %*% G_cov %*% sqrt(diag(1/diag(G_cov)))


  bve.index <- matrix(0, nrow=nrow(bves), ncol(bves))
  for(index5 in 1:ncol(bves)){
    bve.index[,index5] <- miesenberger.index(V1=V1, V= V, G = G_cov, RG = RG, r = sqrt(reli[,index5]), w = c(1,1,-1))
  }

  selection_index <- colSums(bve.index * bves)

}

###
parent_list <- get.pedigree2(population, cohorts="C1F2")
parents <- unique(c(parent_list[,2], parent_list[,3]))
library(lpSolve)

if(!extra_mating){
  selected <- opt_selection(data=cbind(selection_index,parent_list[,1:3]), parents = parents, usage_thresh = max_parent,
                            n_select = nsel)
} else{

  order  <- sort(selection_index, decreasing = TRUE, index.return=TRUE)$ix
  selected <- cbind(selection_index[order], names(selection_index)[order])[1:nsel,]
}

nserler <- nsel
while(nrow(selected)<nsel & nserler < 150){
  nserler <- nserler + 1
  selected <- opt_selection(data=cbind(selection_index,parent_list[,1:3]), parents = parents, usage_thresh = max_parent,
                            n_select = nserler)

}

picks1 <- which(duplicated(c( selected[,2], parent_list[,1]))[-(1:nrow(selected))])

picks <- sample(picks1, nsel, replace = if(nsel>nrow(selected)){TRUE}else{FALSE})

cur_gen <- get.database(population, cohorts = "C1F2")[1,1]

fixed.breeding11 <- cbind(cur_gen,1,picks)

fixed.breeding1 <- cbind(cur_gen,1,picks[1:(nsel/2)], cur_gen,1,picks[(nsel/2+1):nsel],0)

chosen1 <- which(duplicated(c(names(selection_index)[picks], names(sort(selection_index, decreasing = TRUE))))[-(1:length(picks))])

#Generation of C2F2
population <- breeding.diploid(population, breeding.size=c(1035,0),
                               fixed.breeding = fixed.breeding1, repeat.mating = nrep,
                               name.cohort="C2F2", verbose=TRUE, display.progress = FALSE)

population <- breeding.diploid(population, breeding.size=c(200,0), dh.mating = TRUE, name.cohort = paste0("C2-DH") ,
                               selection.m.cohorts = paste0("C2F2"), verbose=TRUE, display.progress = FALSE, repeat.mating=1)

if(add_pheno){
  population <- breeding.diploid(population, new.bv.observation.cohorts = paste0("C2-DH") , share.phenotyped= founder_pheno)
}

# DH generation
population <- breeding.diploid(population, breeding.size=c(3,0), fixed.breeding = cbind(cur_gen,1,picks1[1:3],cur_gen,1,picks1[1:3],0),
                               repeat.mating = 1, dh.mating = TRUE, name.cohort = paste0("C1F2", "_DH"), verbose=TRUE, display.progress = FALSE)
population <- breeding.diploid(population, breeding.size=c(60,0), selection.size = c(4,0), breeding.all.combination = TRUE,
                               selection.m.database = cbind(cur_gen,1,picks1[1:4]), repeat.mating = 10, name.cohort = paste0("C1F2", "_TopCross"),
                               verbose=TRUE, display.progress = FALSE)


fixed.breeding <- list(fixed.breeding1)

}

## recurrent rapid cycle breeding scheme
for(index in 3:(cycle+1)){

  #### Phenotyping ######  
  activ_coh <- paste0("C", index-1, "F2")
  next_coh <- paste0("C", index, "F2")

  add_phenos <- NULL

  activ_t_temp <- activ_t
  if(index<=3){
    activ_t_temp <- 1:3
  }

  if(add_pheno){
    if(index>5 & !backcross){
      add_phenos <- c(add_phenos, paste0("C",max(1,(index-8)):(index-5),"-DH"))
    } else if(index>5 & backcross){
      add_phenos <- c(add_phenos, paste0("C",max(1,(index-6)):(index-5),"-DH"), paste0("C", (index-3):(index-2), "F2_Diversity_DH"))
    }
  }
  if(tc_sel){
    if(index>3){
      add_phenos <- c(add_phenos, paste0("C",(index-3),"F2"))
    }
  }

  if(FALSE){
    if(!tc_sel & backcross){
      if(index>3){
        add_phenos <- c(add_phenos, paste0("C",(index-2),"F2", "_Diversity_DH"))
      }
    }
  }

  if((add_pheno & index > 6 ) || (tc_sel & index > 3)){

  } else{
    add_phenos <- c(add_phenos, "DH-Founder")
  }

  #### Breeding value estimation #####

  population <- breeding.diploid(population, bve=TRUE, bve.cohorts = c(activ_coh, add_phenos),
                                 bve.insert.cohorts = activ_coh, rrblup.bve = TRUE, estimate.reliability = TRUE, bve.ignore.traits = (1:6)[-activ_t_temp])

  # phenotypes are only available after selection of the lines themselves
  if(tc_sel){
    population <- breeding.diploid(population, new.bv.observation.cohorts = activ_coh , share.phenotyped= founder_pheno, n.observation = pheno_tc)
  }
  
  #### Selection #####

  if(sel_index==1){
    # Currently used Index in Rapid Cycle
    bves <- get.bve(population, cohorts = activ_coh)[activ_t_temp,]
    bves[3,] <- -abs(bves[3,] - mean(bves[3,]))
    selection_index <- colSums(bves)
  } else if(sel_index==2){
    # Negative weight on plant height final
    bves <- get.bve(population, cohorts = activ_coh)[activ_t_temp,]
    selection_index <- colSums(bves*(c(1,1,-1)))
  } else {
    # Reliability based index

    bves <- get.bve(population, cohorts = activ_coh)[activ_t_temp,]
    reli <- get.reliabilities(population, cohorts = activ_coh)[activ_t_temp,]
    bve.cov <- var(t(bves))

    V <- bve.cov
    V1 <- MASS::ginv(V)

    bv <- get.bv(population, cohorts = activ_coh)[activ_t_temp,]
    G_cov <- var(t(bv))
    RG <- sqrt(diag(1/diag(G_cov))) %*% G_cov %*% sqrt(diag(1/diag(G_cov)))

    bve.index <- matrix(0, nrow=nrow(bves), ncol(bves))
    for(index5 in 1:ncol(bves)){
      bve.index[,index5] <- miesenberger.index(V1=V1, V= V, G = G_cov, RG = RG, r = sqrt(reli[,index5]), w = c(1,1,-1))
    }
    selection_index <- colSums(bve.index * bves)
  }

  #### Selection via AlphaMate #####
  if(targetdegree>0){

    options(scipen=999)

    addpf <- add_pheno + tc_sel*10 + backcross*100 + share_backcross*1110000 + 10000000 # have a unique directory name for all simulations running in parallel
    temp_geno <- get.geno(population, cohorts=activ_coh)
    pi1 <- rowMeans(temp_geno)/2
    temp_kin <- crossprod(temp_geno- 2*pi1) / sum(pi1*(1-pi1)) / 2

    bves <- get.bve(population, cohorts = activ_coh)[activ_t_temp,]
    bves[3,] <- - bves[3,]
    selection_index <- colSums(bves)

    diag_kin <- diag(temp_kin)
    kin_mat <- temp_kin - min(temp_kin)
    kin_mat <- kin_mat/max(as.vector(kin_mat))*2
    kin_mat_tri <- kin_mat[upper.tri(kin_mat)]
    kin_mat <- t(kin_mat)
    kin_mat[upper.tri(kin_mat)] <- kin_mat_tri

    dir.create(paste0("alphamaterun/degree_", targetdegree, "_", nr,"_", index, "_", addpf))
    print(paste0("alphamaterun/degree_", targetdegree, "_", nr,"_", index, "_", addpf))
    write.table(kin_mat, sprintf("./alphamaterun/degree_%i_%i_%i_%i/Nrm.txt", targetdegree,nr,index, addpf), sep="\t", row.names=TRUE, col.names=FALSE)

    criterion <- data.frame(Genotype = names(selection_index),
                            Value = as.vector(selection_index))
    write.table(criterion, sprintf("./alphamaterun/degree_%i_%i_%i_%i/Criterion.txt", targetdegree, nr,index, addpf), sep="\t", row.names=FALSE, quote=FALSE, col.names=FALSE)


    spec_text <- sprintf("Seed , %i\nNrmMatrixFile , Nrm.txt\nSelCriterionFile , Criterion.txt\nNumberOfMatings , %i\nNumberOfParents , %i\nTargetDegree , %i\nLimitContributions , Yes\nLimitContributionsMax , 1\nStop",
                         seed_val, nmatings, nparents, targetdegree)
    write(x = spec_text, file=sprintf("./alphamaterun/degree_%i_%i_%i_%i/AlphaMateSpec.txt", targetdegree, nr,index, addpf))

    setwd(paste0("alphamaterun/degree_", targetdegree, "_", nr,"_", index, "_", addpf))



    suppressWarnings(file.remove("MatingPlanModeOptTarget1.txt"))


    system(sprintf(paste0( sim_path, "/AlphaMate_0.2.0/bin/AlphaMate_Linux_0.2.02")))

    suppressWarnings(file.remove("Nrm.txt"))

    setwd(sim_path)

    matings <- read.table(sprintf("./alphamaterun/degree_%i_%i_%i_%i/MatingPlanModeOptTarget1.txt", targetdegree, nr, index, addpf), sep="", header=TRUE)

    write.table(matings, sprintf("./alphamaterun/degree_%i_%i_%i_%i/matingplan.tabular", targetdegree, nr, index, addpf), sep="\t", quote=FALSE, row.names=FALSE)

    picks1 <- match(c(as.character(matings$Parent1), as.character(matings$Parent2)),
                   names(selection_index))
    picks <- sample(picks1, nsel)



  }


    if(backcross & share_backcross>0 & index>3){
      nsela <- round(nsel * share_backcross)
      nselb <- nsel - nsela

      picks1a <- sort(selection_index[1:round(1035*share_backcross)], index.return=TRUE, decreasing = TRUE)$ix[1:nsela]

      if(nselb>0){
        picks1b <- sort(selection_index[-(1:round(1035*share_backcross))], index.return=TRUE, decreasing = TRUE)$ix[1:nselb] + round(1035*share_backcross)
      } else{
        picks1b <- NULL
      }

      picks1 <- c(picks1a, picks1b)

      if(length(picks1b)>0){
        toppicks <- picks1b[1:4]
      } else{
        toppicks <- picks1[1:4]
      }

    } else{
      picks1 <- sort(selection_index, index.return=TRUE, decreasing = TRUE)$ix[1:nsel]

      toppicks <- picks1[1:4]
    }

  if(targetdegree==0){
    picks <- sample(picks1, nsel)
  }

  cur_gen <- get.database(population, cohorts = activ_coh)[1,1]


  fixed.breeding[[index-1]] <- cbind(cur_gen,1,picks[1:(nsel/2)], cur_gen,1,picks[(nsel/2+1):nsel],0)

  #### Generation of external material #####
  if(backcross){

    bv_material <- rowMeans(get.bv(population, cohorts=activ_coh))

    pop1 <- pop2

    new_bvs <- rowMeans(get.bv(pop1, gen=1))

    i <- 1
    while(sum(bv_material*c(1,1,-1))> sum(new_bvs*c(1,1,-1))){

      cur <- length(pop1$breeding)

      weights = bv_material - new_bvs
      weights[c(1,2,4,5)][weights[c(1,2,4,5)]<0] = weights[c(1,2,4,5)][weights[c(1,2,4,5)]<0] * 5
      weights[c(3,6)][weights[c(3,6)]>0] = weights[c(3,6)][weights[c(3,6)]>0] * 5

      # some slow selection to get a population with similar genomic value but different genotypes

      pop1 <- breeding.diploid(pop1, breeding.size = c(ncol(geno2)/2,0), selection.size = c(200,0), selection.criteria = "bv",
                               multiple.bve.scale.m = "unit", multiple.bve.weights.m = weights,
                               verbose=FALSE, display.progress = FALSE)

      new_bvs <- rowMeans(get.bv(pop1, gen=cur+1))

      i <- i+1
    }

    store <- rbind(store, cbind(new_bvs, bv_material, i-1))
    print(paste0("Additional Founder simulation took ", i-1, " cycles"))
    print(new_bvs)
    print(bv_material)

    haplo_new <- get.haplo(pop1, gen= (i-1))
    bv_new <- get.bv(pop1, gen=(i-1))

    population <- creating.diploid(population = population, dataset = haplo_new, sex.quota = 0, name.cohort = paste0(activ_coh, "_Diversity"),
                                   miraculix = miraculix)


    if(add_pheno){

      population <- breeding.diploid(population, breeding.size = c(200,0), dh.mating = TRUE, dh.sex = 0,
                                     selection.m = "random",
                                     selection.m.cohorts = paste0(activ_coh, "_Diversity"),
                                     name.cohort = paste0(activ_coh, "_Diversity_DH"),
                                     verbose=TRUE, display.progress = FALSE)
      population <- breeding.diploid(population, new.bv.observation.cohorts = paste0(activ_coh, "_Diversity_DH") , n.observation = pheno_perse)
    }


    database_enter <- get.database(population, cohorts=paste0(activ_coh, "_Diversity"))


    fixed.breeding[[index-1]] <- cbind(fixed.breeding[[index-1]][,1], rep(1, 1035), fixed.breeding[[index-1]][,3],
                                       fixed.breeding[[index-1]][,4], fixed.breeding[[index-1]][,2], fixed.breeding[[index-1]][,6])

    fixed.breeding[[index-1]] <- fixed.breeding[[index-1]][sample(1:1035,1035),] # shuffle matings
    fixed.breeding[[index-1]][1:round(1035*share_backcross),1:3] <- cbind(1,1, sample(database_enter[3]:database_enter[4],round(1035*share_backcross), replace=TRUE))
    # make sure last mating are the ones generated with a backcross
    nrep <- 1
  }

  #### Generation of crosses #####
  population <- breeding.diploid(population, breeding.size=c(1035,0),
                                 fixed.breeding = fixed.breeding[[index-1]], repeat.mating = nrep,
                                 name.cohort=next_coh, verbose=TRUE, display.progress = FALSE)

  population <- breeding.diploid(population, breeding.size=c(3,0), fixed.breeding = cbind(cur_gen,1,toppicks[1:3],cur_gen,1,toppicks[1:3],0), repeat.mating = 1, dh.mating = TRUE, name.cohort = paste0(activ_coh, "_DH"), verbose=TRUE, display.progress = FALSE)
  population <- breeding.diploid(population, breeding.size = c(60,0), selection.size = c(4,0), breeding.all.combination = TRUE, selection.m.database = cbind(cur_gen,1,toppicks[1:4]), repeat.mating = 10, name.cohort = paste0(activ_coh, "_TopCross"), verbose=TRUE, display.progress = FALSE)


  # DH generation for BVE
  population <- breeding.diploid(population, breeding.size=c(200,0), dh.mating = TRUE, name.cohort = paste0("C",index,"-DH") ,
                                 selection.m.cohorts = paste0("C",index,"F2"), verbose=TRUE, display.progress = FALSE, repeat.mating=1)

  if(add_pheno){
    population <- breeding.diploid(population, new.bv.observation.cohorts = paste0("C",index,"-DH") , share.phenotyped= founder_pheno)
  }

}



#### Analysis of simulation results #####

priorbv <- get.bv(population, cohorts = "DH-Founder")

priorbv <- rbind(priorbv, priorbv[1,]+priorbv[2,]-priorbv[3,], priorbv[4,]+priorbv[5,]-priorbv[6,])

resultbv <- list()

resultvar <- list()
gains <- list()
accs <- list()
for(index in 1:cycle){
  temp1 <- get.bv(population, cohorts = paste0("C",index,"F2"))
  temp1 <- rbind(temp1, temp1[1,]+temp1[2,]-temp1[3,], temp1[4,]+temp1[5,]-temp1[6,])
  resultbv[[length(resultbv)+1]] <- temp1
  resultvar[[length(resultvar)+1]] <- diag(var(t(resultbv[[length(resultbv)]])))

  gains[[length(gains)+1]] <- rowMeans(resultbv[[length(gains)+1]]) - rowMeans(priorbv)

  temp1 <- get.bv(population, cohorts = paste0("C",index,"F2_DH"))
  temp1 <- rbind(temp1, temp1[1,]+temp1[2,]-temp1[3,], temp1[4,]+temp1[5,]-temp1[6,])
  resultbv[[length(resultbv)+1]] <- temp1

  resultvar[[length(resultvar)+1]] <- diag(var(t(resultbv[[length(resultbv)]])))

  gains[[length(gains)+1]] <- rowMeans(resultbv[[length(gains)+1]]) - rowMeans(priorbv)

  temp1 <- get.bv(population, cohorts = paste0("C",index,"F2_TopCross"))
  temp1 <- rbind(temp1, temp1[1,]+temp1[2,]-temp1[3,], temp1[4,]+temp1[5,]-temp1[6,])
  resultbv[[length(resultbv)+1]] <- temp1

  resultvar[[length(resultvar)+1]] <- diag(var(t(resultbv[[length(resultbv)]])))

  gains[[length(gains)+1]] <- rowMeans(resultbv[[length(gains)+1]]) - rowMeans(priorbv)

  suppressWarnings(accs[[length(accs)+1]] <- analyze.bv(population, cohorts=paste0("C",index,"F2"))[[1]][1,])
}


temp1 <- get.bv(population, cohorts = paste0("DH-Founder"))
temp1 <- rbind(temp1, temp1[1,]+temp1[2,]-temp1[3,], temp1[4,]+temp1[5,]-temp1[6,])
resultbv[[length(resultbv)+1]] <- temp1
resultvar[[length(resultvar)+1]] <- diag(var(t(resultbv[[length(resultbv)]])))

gains[[length(gains)+1]] <- rowMeans(resultbv[[length(gains)+1]]) - rowMeans(priorbv)

temp1 <- get.bv(population, cohorts = paste0("C",0,"F2"))
temp1 <- rbind(temp1, temp1[1,]+temp1[2,]-temp1[3,], temp1[4,]+temp1[5,]-temp1[6,])
resultbv[[length(resultbv)+1]] <- temp1
resultvar[[length(resultvar)+1]] <- diag(var(t(resultbv[[length(resultbv)]])))

gains[[length(gains)+1]] <- rowMeans(resultbv[[length(gains)+1]]) - rowMeans(priorbv)

save(file=paste0("Rapid_cycle/result13/Max",max_parent,"OptiSel", replacer, targetdegree, "Extra", extra_mating, "Nsel", nsel, "Cycle", cycle, "Pheno", founder_pheno, "DHs", dh_sel, "_results", sel_index,"_tc", tc_sel, "_phenoDH", add_pheno,"_back", backcross, "_", share_backcross, "_", nr,  "_small.RData"), list=c("resultvar", "gains", "accs", "store"))

if(nr==1){
  save(file=paste0("Rapid_cycle/pops13/Max",max_parent,"OptiSel", replacer, targetdegree, "Extra", extra_mating, "Nsel", nsel, "Cycle", cycle, "Pheno", founder_pheno, "DHs", dh_sel, "_sim", sel_index,"_tc", tc_sel, "_phenoDH", add_pheno,  "_back", backcross,"_", share_backcross, "_", nr, ".RData"), list=c("population"))
}

hetero1 <- NULL
hetero2 <- NULL
for(index in 1:cycle){
  resultgeno <- get.geno(population, cohorts = paste0("C",index,"F2"))

  hetero1 <- c(hetero1, mean(resultgeno==1)) ## Share of heterozygous markers
  p_i <- rowMeans(resultgeno)/2
  p_i[p_i>0.5] <- 1-p_i[p_i>0.5]
  hetero2 <- c(hetero2, sum(p_i==0)/length(p_i)) # Share fixated markers
}


save(file=paste0("Rapid_cycle/result13large/Max",max_parent,"OptiSel", replacer, targetdegree, "Extra", extra_mating, "Nsel", nsel, "Cycle", cycle, "Pheno", founder_pheno, "DHs", dh_sel, "_results", sel_index,"_tc", tc_sel, "_phenoDH", add_pheno,"_back", backcross, "_", share_backcross, "_", nr,  ".RData"), list=c("resultbv", "resultvar", "gains", "accs", "hetero1", "hetero2", "store"))




