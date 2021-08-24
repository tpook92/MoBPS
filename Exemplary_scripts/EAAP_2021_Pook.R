
# The share of genotyped animals and the correlation matrices do NOT
# reflect reality, but are not openly available
# For the presented results in the EAAP talk realistic numbers were used!

share_gen <- 1

cor_matrix1 <- diag(14) # residual covariance matrix
cor_matrix2 <- diag(14) # genetic covariance matrix

# phenotypic variance / heritablity of traits
# 13/14 are the behavior / tail biting traits
tar_var <- c(rep(1,12),0.91,1)^2
heri <- c(rep(0.3,12),0.29,0.3)

# Chose settings from a shell-script
args <- commandArgs(TRUE)
# e.g. args <- c(1,1,1) for the reference scenario
run <- as.numeric(args[1])
scenario <- as.numeric(args[2])
size <- as.numeric(args[3]) # use small numbers to scale down number of individuals for fast simulation


test1 <- FALSE # activating this will simulate the BVE instead of computing (very fast but not realistic!)

if(is.na(size)){
  size <- 1
}

set.seed(run)
library(MoBPS)
library(RandomFieldsUtils)

# Scenarios 11 - 30 consider different correlations between tail-biting and behavior trait
if(scenario<=10){
  cor_matrix2[13,14] <- cor_matrix2[14,13] <- 0.5
} else{
  cor_matrix2[13,14] <- cor_matrix2[14,13] <- (scenario - 10)/20
}

# Set up phenotyping / traits
# trait 13: behavior trait
# trait 14: tail-biting


sigma_g <- tar_var * heri
bucht <- 6

# Depending on the cohort difference traits are phenotyped
# 13/14 are the behavior / tail biting traits
pheno_OFT <- c(0,1,1,1,rep(0,10))
pheno_CTS <- c(1,0,1,1,1,0,0,0,0,1,1,1,1,1)
pheno_Betrieb <- c(rep(0,5),1,1,1, rep(0,6))
pheno_Zleistung <- c(1,0, rep(1,12))

# selection index
sel_index <- c(0.4,0,0,0.1,0.2,0,0.15,0,0.1,0,0.05,0,0,0)

# share of fellow culprits in a pen
management <- 0.5
# Should a pen only contain at least half sibs
halfsib_bucht <- TRUE

# Different scenarios presented
if(scenario == 1){
  #Reference

} else if(scenario==2){
  # Basic TB strategy
  sel_index[14] <- (-0.2)
} else if(scenario==3){
  # 4 animals per pen
  bucht <- 4
  sel_index[14] <- (-0.2)
} else if(scenario==4){
  # selection on correlated behaviour trait
  sel_index[13] <- (-0.2)
} else if(scenario==5){
  # improved management
  sel_index[14] <- (-0.2)
  management <- 0.2
} else if(scenario==6){
  # improved management / video analysis
  sel_index[14] <- (-0.2)
  management <- 0
} else if(scenario==7){
  sel_index[14] <- (-0.2)
  halfsib_bucht <- FALSE
}else if(scenario==8){
  # all animals become fellow culprits (not presented at EAAP)
  sel_index[14] <- (-0.2)
  management <- 1
} else if(scenario==9){
  # all animals become fellow culprits + pens with unrelated animals  (not presented at EAAP)
  sel_index[14] <- (-0.2)
  halfsib_bucht <- FALSE
  management <- 1
} else if(scenario==10){
  # all animals become fellow culprits + 4 animals per pen  (not presented at EAAP)
  bucht <- 4
  sel_index[14] <- (-0.2)
  management <- 1
} else if(scenario<=30){
  # selection on the correlated behaviour trait with different correlation between tail-biting and behaviour
  sel_index[13] <- (-0.2)
}

# Build up some basic LD structure for the founders
dataset <- founder.simulation(nindi=100, sex.quota = 0.5, map = MoBPSmaps::map_pig2, display.progress = FALSE,verbose=FALSE)


# Creation of the initial population list / founder population
population <- creating.diploid(dataset = dataset, map = MoBPSmaps::map_pig2, n.additive = rep(1000,14),
                               shuffle.cor = cor_matrix2,
                               new.phenotype.correlation = cor_matrix1,
                               mean.target = 0,
                               var.target = sigma_g)

# Phenotypic transformation for traits with discret realizations
f1 <- function(x){  if(x < qnorm(0.01, mean=0, sd=0.91)){    y <- 6  } else if(x < qnorm(0.06, mean=0, sd=0.91)){    y <- 5  } else if(x < qnorm(0.26, mean=0, sd=0.91)){    y <- 4  } else if(x < qnorm(0.76, mean=0, sd=0.91)){    y <- 3  } else if(x < qnorm(0.96, mean=0, sd=0.91)){    y <- 2  } else {    y <- 1  };  return(y) }
f2 <- function(x){   if(x < qnorm(0.01, mean=0, sd=1)){     y <- 1   } else {     y <- 0   };   return(y) }
population <- creating.phenotypic.transform(population, phenotypic.transform.function = f1, trait = 13)
population <- creating.phenotypic.transform(population, phenotypic.transform.function = f2, trait = 14)

population <- breeding.diploid(population, breeding.size = 10500*size, name.cohort = "piglets_0",
                               share.genotyped = 0, display.progress = FALSE)

population <- breeding.diploid(population,  heritability = heri)



tail_share <- NULL
for(index in 0:9){


  # Cohorts to use in the BVE for the selection of animals
  # to use in the on-farm-test and central-test-station
  if(index>0){
    XXX <- c(paste0("CTS_", index-1), paste0("CTSS_", index-1),
             paste0("NB_", index-1), paste0("NG_", index-1), paste0("TB_", index-1),
             paste0("OFT_", index-1, "_F"), paste0("OFT_", index-1, "_M"),
             paste0("OFTS_", index-1, "_M"),
             paste0("SB_", index-1), paste0("SS_", index-1),
             paste0("STS_", index-1))
  } else{
    XXX <- NULL
  }

  # Only enter "new" breeding value estimates for piglets
  XXX1 <- XXX
  XXX2 <- c(paste0("piglets_",index,"_M"), paste0("piglets_",index,"_F"))
  XXX3 <- c(XXX1,XXX2)


  if(index!=0 && !test1){
    population <- breeding.diploid(population,
                                   bve.ignore.traits = c(2,3,6,8,10,12),
                                   bve=TRUE,
                                   singlestep.active = TRUE,
                                   bve.cohorts = XXX3,
                                   bve.insert.cohorts = XXX2)

    # Use the R-package rrBLUP for BVE as heritability is not 0.3 due to fellow culprit behavior
    population <- breeding.diploid(population,
                                   bve.ignore.traits = c(1:13),
                                   bve=TRUE,
                                   rrblup.bve = TRUE,
                                   singlestep.active = TRUE,
                                   bve.cohorts = XXX3,
                                   bve.insert.cohorts = XXX2)
  } else{
    population <- breeding.diploid(population, bve=TRUE, bve.pseudo = TRUE,
                                   bve.cohorts = XXX3, bve.pseudo.accuracy = c(rep(0.3,12),0,0))
  }

  # Central test station
  population <- breeding.diploid(population, breeding.size = c(1500,0)*size, selection.size = c(1500,0)*size,
                                 copy.individual.m=TRUE,
                                 selection.m.cohorts = paste0("piglets_",index,"_M"),
                                 name.cohort = paste0("CTS_", index),
                                 selection.criteria = "bve",
                                 multiple.bve.weights.m = sel_index,
                                 multiple.bve.scale.m = "bve_sd", display.progress = FALSE)

  # On farm test male
  population <- breeding.diploid(population, breeding.size = c(2500,0)*size, selection.size = c(4000,0)*size,
                                 copy.individual.m=TRUE,
                                 selection.m.cohorts = paste0("piglets_",index,"_M"),
                                 ignore.best = 1500*size,
                                 name.cohort = paste0("OFT_", index, "_M"),
                                 selection.criteria = "bve",
                                 multiple.bve.weights.m = sel_index,
                                 multiple.bve.scale.m = "bve_sd", display.progress = FALSE)

  # On farm test female
  population <- breeding.diploid(population, breeding.size = c(0,4000)*size, selection.size = c(0,4000)*size,
                                 copy.individual.f=TRUE,
                                 selection.f.cohorts = paste0("piglets_",index,"_F"),
                                 name.cohort = paste0("OFT_", index, "_F"),
                                 selection.criteria = "bve",
                                 multiple.bve.weights.m = sel_index,
                                 multiple.bve.scale.m = "bve_sd", display.progress = FALSE)

  # Phenotyping
  population <- breeding.diploid(population, phenotyping.cohorts = c(paste0("OFT_", index, "_F"), paste0("OFT_", index, "_M")),
                                 n.observation = pheno_OFT)
  population <- breeding.diploid(population, phenotyping.cohorts = paste0("CTS_", index),
                                 n.observation = pheno_CTS)



  # Simulation of the fellow culprit behaviour
  {
    phenos <- get.pheno(population, cohorts= paste0("CTS_", index))
    bvs <- get.bv(population, cohorts= paste0("CTS_", index))
    pedi <- get.pedigree(population, cohorts= paste0("CTS_", index))
    before <- c(cor(phenos[14,], bvs[14,]), mean(phenos[14,]))

    if(halfsib_bucht){

      sires <- unique(pedi[,2])

      for(index2 in 1:length(sires)){
        sibs <- which(pedi[,2]==sires[index2])
        if(length(sibs)>1){
          sibs <- sample(sibs)
        }

        if(floor(length(sibs)/bucht)>0){
          for(index3 in 1:floor(length(sibs)/bucht)){
            activ_pheno <- phenos[14, sibs[1:bucht+index3*bucht-bucht]]

            if(sum(is.na(activ_pheno))>0){
              print(activ_pheno)
              print(sibs[1:bucht+index3*bucht-bucht])
              print(index3)
            }
            if(sum(activ_pheno)>0){
              activ_pheno[activ_pheno==0] <- rbinom(sum(activ_pheno==0),1, management)

            }

            phenos[14, sibs[1:bucht+index3*bucht-bucht]] <- activ_pheno

          }
        }

      }

    } else{

      sibs <- sample(1:nrow(pedi))

      if(floor(length(sibs)/bucht)>0){
        for(index3 in 1:floor(length(sibs)/bucht)){
          activ_pheno <- phenos[14, sibs[1:bucht+index3*bucht-bucht]]

          if(sum(activ_pheno)>0){
            activ_pheno[activ_pheno==0] <- rbinom(sum(activ_pheno==0),1, management)

          }
          phenos[14, sibs[1:bucht+index3*bucht-bucht]] <- activ_pheno

        }
      }
    }


    after <- c(cor(phenos[14,], bvs[14,]), mean(phenos[14,]))

    tail_share <- rbind(tail_share, c(before, after ))

    print(tail_share)

    population <- insert.bve(population, bves = cbind(colnames(phenos), t(phenos)), type="pheno")
  }

  # pre selection + genotyping of on farm test and central test station boars

  population <- breeding.diploid(population, breeding.size = c(1050,0)*size, selection.size = c(1050,0)*size,
                                 copy.individual.m = TRUE, name.cohort = paste0("CTSS_", index),
                                 selection.m.cohorts = paste0("CTS_", index),
                                 added.genotyped =  share_gen,
                                 selection.criteria = "bve",
                                 multiple.bve.weights.m = sel_index,
                                 multiple.bve.scale.m = "bve_sd", display.progress = FALSE)


  population <- breeding.diploid(population, breeding.size = c(1750,0)*size, selection.size = c(1750,0)*size,
                                 copy.individual.m = TRUE, name.cohort=paste0("OFTS_", index, "_M"),
                                 selection.m.cohorts = paste0("OFT_", index, "_M"),
                                 added.genotyped = share_gen,
                                 selection.m = "random", display.progress = FALSE)


  XXX_insert <- c(paste0("CTS_", index), paste0("CTSS_", index), paste0("OFT_", index, "_F"), paste0("OFT_", index, "_M"),
                  paste0("OFTS_", index, "_M"))

  XXX <- c(XXX, XXX_insert)


  if(!test1){
    population <- breeding.diploid(population, bve=TRUE, bve.cohorts = XXX, singlestep.active = TRUE,
                                   bve.insert.cohorts = XXX_insert,
                                   bve.ignore.traits = c(2,3,6,8,10,12))

    # Use the R-package rrBLUP for BVE as heritability is not 0.3 due to fellow culprit behavior
    population <- breeding.diploid(population, bve=TRUE, bve.cohorts = XXX, singlestep.active = TRUE,
                                   rrblup.bve = TRUE,
                                   bve.insert.cohorts = XXX_insert,
                                   bve.ignore.traits = c(1:13))
  } else{
    population <- breeding.diploid(population, bve=TRUE, bve.pseudo = TRUE,
                                   bve.cohorts = XXX_insert, bve.pseudo.accuracy = c(rep(0.6,12),-0.6,-0.4))
  }


  # nucleus gits
  population <- breeding.diploid(population, breeding.size = c(0,600)*size, selection.size = c(0,600)*size,
                                 copy.individual.f = TRUE, name.cohort = paste0("NG_", index),
                                 selection.f.cohorts = paste0("OFT_", index, "_F"),
                                 added.genotyped = 1,
                                 selection.criteria = "bve",
                                 multiple.bve.weights.m = sel_index,
                                 multiple.bve.scale.m = "bve_sd", display.progress = FALSE)

  # remaining sows ((for more phenotyping))
  database_rest <- group.diff(population, cohorts=paste0("OFT_", index, "_F"), remove.cohorts=paste0("NG_", index))

  population <- breeding.diploid(population, breeding.size = c(0,3400)*size, selection.size = c(0,3400)*size,
                                 copy.individual.f = TRUE, name.cohort=paste0("SS_", index),
                                 selection.f.database = database_rest,
                                 selection.m = "random", display.progress = FALSE)

  population <- breeding.diploid(population, phenotyping.cohorts = paste0("SS_", index),
                                 n.observation = pheno_Betrieb)

  # nucleus boars
  population <- breeding.diploid(population, breeding.size = c(50,0)*size, selection.size = c(50,0)*size,
                                 copy.individual.m = TRUE, name.cohort = paste0("NB_", index),
                                 selection.m.cohorts = c(paste0("CTS_", index), paste0("OFT_", index, "_M")),
                                 added.genotyped = 1,
                                 selection.criteria = "bve",
                                 multiple.bve.weights.m = sel_index,
                                 multiple.bve.scale.m = "bve_sd", display.progress = FALSE)

  # terminal boars
  population <- breeding.diploid(population, breeding.size = c(500,0)*size, selection.size = c(550,0)*size,
                                 copy.individual.m = TRUE, name.cohort = paste0("TB_", index),
                                 ignore.best = 50 * size,
                                 selection.m.cohorts = c(paste0("CTS_", index), paste0("OFT_", index, "_M")),
                                 added.genotyped = 1,
                                 selection.criteria = "bve",
                                 multiple.bve.weights.m = sel_index,
                                 multiple.bve.scale.m = "bve_sd", display.progress = FALSE)

  # remaining boars ((for more phenotyping))
  database_rest <- group.diff(population, cohorts=c(paste0("CTS_", index)),
                              remove.cohorts=c(paste0("NB_", index), paste0("TB_", index)))

  n_new <- sum(database_rest[,4] - database_rest[,3] + 1)

  population <- breeding.diploid(population, breeding.size = c(n_new,0), selection.size = c(n_new,0),
                                 copy.individual.m = TRUE, name.cohort=paste0("STS_", index),
                                 selection.m.database = database_rest,
                                 selection.m = "random", display.progress = FALSE)

  population <- breeding.diploid(population, phenotyping.cohorts = paste0("STS_", index),
                                 n.observation = pheno_Zleistung)

  database_rest <- group.diff(population, cohorts=c(paste0("OFT_", index, "_M")),
                              remove.cohorts=c(paste0("NB_", index), paste0("TB_", index)))

  n_new <- sum(database_rest[,4] - database_rest[,3] + 1)

  population <- breeding.diploid(population, breeding.size = c(n_new,0), selection.size = c(n_new,0),
                                 copy.individual.m = TRUE, name.cohort=paste0("SB_", index),
                                 selection.m.cohorts = paste0("OFT_", index, "_M"),
                                 selection.m = "random", display.progress = FALSE)

  population <- breeding.diploid(population, phenotyping.cohorts = paste0("SB_", index),
                                 n.observation = pheno_Betrieb)

  # piglets for the next cycle

  population <- breeding.diploid(population, breeding.size = 10500*size, name.cohort = paste0("piglets_",index+1),
                                 selection.m.cohorts = paste0("NB_", index),
                                 selection.f.cohorts = paste0("NG_", index),
                                 share.genotyped = 0, display.progress = FALSE)

}

# Derive some basic results from the simulation
# Accuracy of the BVE
# Development of genomic values
# Inbreeding development

acc <- NULL

for(index in 0:10){
  acc <- cbind(acc, analyze.bv(population, cohorts=paste0("piglets_",index,c("_M", "_F")) )[[1]][1,])
}

BVS <- NULL
for(index in 0:10){
  BVS <- cbind(BVS, rowMeans(get.bv(population, cohorts=paste0("piglets_",index,c("_M", "_F")))))
}

kin <- NULL
for(index in 0:10){
  kin <- rbind(kin, kinship.emp.fast(population=population, cohorts=paste0("piglets_",index,c("_M", "_F")), ibd.obs = 1000, hbd.obs = 250))
}

save(file=paste0("EAAP_2021/results/Scenario", scenario, "Nr", run, "SIZE", size, "_results.RData"), list=c("tail_share", "BVS", "kin", "acc"))
save(file=paste0("EAAP_2021/data/Scenario", scenario, "Nr", run, "SIZE", size, ".RData"), list=c("population", "tail_share", "BVS", "kin", "acc"))


warnings()
