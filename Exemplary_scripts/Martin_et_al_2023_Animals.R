## This program includes simulations using MoBPS v.1.9.20 software.
## The software MoBPS and a MoBPS User Manual is available at https://cran.r-project.org/web/packages/MoBPS/index.html
### or at https://github.com/tpook92/MoBPS 

args <- commandArgs(TRUE)
# exemplary setting for args:
# args = c(1,1,0.1) For first simulation of baseline scenario with 10% of the individuals
nr <- as.numeric(args[1])
scenario <- as.numeric(args[2])
n <- as.numeric(args[3])



# Loading in R - Packages
library(MoBPS)
library(RandomFieldsUtils)
library(miraculix)

# In our simulation we used the commercial software MiXBLUP for a fast
# multi-trait evaluation.
# mixblup.bve = FALSE will run separate univariate models for each trait
# With full number of individuals this require about 15 times as much computing time

mixblup.bve = FALSE
mixblup.path = "/cm/shared/apps/mixblup/current/MiXBLUP.exe"

n_rep <- 20 #  number of breeding cycles
scaling <- "bve" # traits in the index are scaled by SD of the estimated breeding values

# Setting up phenotyping classes

pheno <- list()

pheno_matrix <- matrix(c(0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0,
                         1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), byrow = TRUE, ncol =22)


bve_acc <- NULL

# setting numbers of individuals (for simulation: n = 1; only relevant when scaling down individual numbers)
N_indi <- ceiling(7654*n)
N_offspring <- ceiling(12932*n)
N_sel <- ceiling(129*n)
N_self <- ceiling(1293*n)
N_selp <- ceiling(340*n)
N_selfp <- ceiling(5892*n)
N_allmale <- ceiling(469*n)
N_allfemale <- ceiling(7185*n)
N_newM <- ceiling(469*n)
N_newF <- ceiling(7185*n)
N_licensing <- ceiling(150*n)

heritability = c(0.1,  0.2, 0.25, 0.3, 0.26, 0.13, 0.22, 0.17, 0.25, 0.4, 0.4, 0.26, 0.2, 0.5, 0.33, 0.36, 0.14, 0.12, 0.19, 0.4, 0.32, 0.15)
repeatability = c(0.3, 0.2, 0.25, 0.3, 0.26, 0.13, 0.22, 0.17, 0.25, 0.4, 0.4, 0.26, 0.2, 0.5, 0.33, 0.36, 0.14, 0.12, 0.19, 0.4, 0.32, 0.15)
selectionindex = c(10, 7.5, 7.5, 20, 10, 2.5, 2.5, 3.125, 5, 0, 15, 0, 0, 0, 2.5, 2.5, 2.5, 6.25, 3.125, 0, 0, 0)

selectionindex_licensing <- selectionindex


# setting up different scenarios
keep_pheno <- FALSE
pheno_breeding <- FALSE

BVE_offspring <- TRUE
share_pheno <- 0.25
selectionindex[9] <- 5
pheno_field <- pheno_matrix[7,]
pheno_station <- pheno_matrix[3,]
pheno_ram <- pheno_matrix[1,]

# NA50
if(scenario==2){
  share_pheno <- 0.5
}

# NA100
if(scenario==3){
  share_pheno <- 1
}

# NA100+
if(scenario==4){
  selectionindex[9] <- 10
  share_pheno <- 1
}

# ST+FT
if(scenario==5){
  pheno_field <- pheno_matrix[6,]

}

# FT
if(scenario==6){
  pheno_field <- pheno_matrix[6,]

  pheno_station <- pheno_matrix[7,]

}



set.seed(nr)


# Phenotyping of combine traits is done manually in case only some of the
# traits included in the combination are phenotyped

phenotyping_combi <- function(population, cohorts=NULL){

  to_process <- which(population$info$is.combi)

  pheno <- pheno_ori <- get.pheno(population, cohorts=cohorts)

  pop1 <- breeding.diploid(population, phenotyping.cohorts = cohorts)

  pheno_temp <- rowMeans(get.pheno(pop1, cohorts=cohorts))

  for(index in 1:nrow(pheno)){
    pheno[index,is.na(pheno[index,])] <- pheno_temp[index]
  }

  for(bven in to_process){
    activ_traits = (!is.na(pheno_ori)) * population$info$combi.weights[[bven]]
    pheno[bven,] = colSums(population$info$combi.weights[[bven]] * pheno, na.rm=TRUE)
    pheno[bven, colSums(activ_traits)==0] <- NA
  }

  pheno[1:19,][is.na(pheno_ori)[1:19,]] <- NA

  new.pheno = cbind(colnames(pheno), t(pheno))

  population <- insert.pheno(population, phenos = new.pheno)

  return(population)
}


# LD build-up for realistic founder population
# insert real genomic data

dataset_out <- founder.simulation(nindi = N_indi, vcf = "Genetic_data.vcf", n.gen=5,  ## Data is available upon request
                                  big.output = TRUE,
                                  bpcm.conversion = 1000000)

map <- dataset_out[[2]]
dataset <- dataset_out[[1]]

rm(dataset_out)

set.seed(nr)

G <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 1, 0.59, 0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0.59, 1, 0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0.4, 0.8, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 1, 0.74, 0.28, 0, 0, 0.8, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0.74, 1, 0.58, 0.29, 0, 0, 0, 0, 0, 0.8, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0.28, 0.58, 1, 0.56, 0, 0, 0, 0.8, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0.29, 0.56, 1, 0, 0, 0, 0, 0.8, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0.8, 0, 0, 0, 0, 1, -0.75, 0.12, 0.03, 0.04, 0.1, 0.15, 0.08, 0.05, -0.03,
              0, 0, 0, 0, 0, 0, 0, 0, 0, -0.75, 1, -0.01, 0.02, 0.29, 0.07, -0.05, 0.02, -0.17, 0.14,
              0, 0, 0, 0, 0, 0, 0.8, 0, 0, 0.12, -0.01, 1, 0.13, 0.03, 0.22, 0.38, 0.29, -0.1, 0.06,
              0, 0, 0, 0, 0, 0, 0, 0.8, 0, 0.03, 0.02, 0.13, 1, -0.01, 0.02, -0.09, -0.09, -0.15, 0.1,
              0, 0, 0, 0, 0, 0.8, 0, 0, 0, 0.04, 0.29, 0.03, -0.01, 1, 0.05, -0.09, 0.03, -0.11, 0.29,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0.1, 0.07, 0.22, 0.02, 0.05, 1, 0.21, 0.35, -0.27, 0.29,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0.15, -0.05, 0.38, -0.09, -0.09, 0.21, 1, 0.53, -0.02, 0.09,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0.08, 0.02, 0.29, -0.09, 0.03, 0.35, 0.53, 1, -0.08, 0.16,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0.05, -0.17, -0.1, -0.15, -0.11, -0.27, -0.02, -0.08, 1, -0.42,
              0, 0, 0, 0, 0, 0, 0, 0, 0, -0.03, 0.14, 0.06, 0.1, 0.29, 0.29, 0.09, 0.16, -0.42, 1), ncol = 19)
G[1,9] <- G[9,1] <- (-0.2)

# The G matrix is not positive definite
# Projection of the G Matrix to the space of positive definite matrices
test <- eigen(G)
test$values[test$values<0] <- 0
M <- diag(test$values)
S <- test$vectors
newG <- S %*% M %*% solve(S)
diag(newG) <- diag(newG) + 0.0000001
newG <- newG * matrix(1/sqrt(diag(newG)), nrow=nrow(newG), ncol=nrow(newG), byrow=TRUE) * matrix(1/sqrt(diag(newG)), nrow=nrow(newG), ncol=nrow(newG), byrow=FALSE)

print(round(newG, digits=3))


# Initialization of founder population

population <- creating.diploid(dataset= dataset, sex.s = c(rep(1,N_newM), rep(2,N_newF)),
                               name.cohort="breedingpopulation_0",
                               n.additive = c(rep(1000,19)),
                               map = map,
                               trait.name = c("NL", "wool", "MC", "BC", "ADGf", "FLNf", "UMDf", "UFDf", "NA", "ADGs",
                                              "FCR", "UMDs", "UFDs", "FLNs", "SW", "BMA", "WC", "SFA", "PKF"), #single trait names
                               shuffle.traits = 1:19,
                               shuffle.cor = newG)

rm(dataset)

summary(population)

# trait standardization

population <- bv.standardization(population, mean.target = 100,
                                 var.target = 10)

# combining traits into trait complexes

population <- add.combi(population, trait = 20, combi.weights = c(0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0), trait.name = "ADG")

population <- add.combi(population, trait = 21, combi.weights = c(0,0,0,0,0,1,1,0,0,0,0,1,0,1,1,1,1,0,0), trait.name = "MEAT")

population <- add.combi(population, trait = 22, combi.weights = c(0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,2,1), trait.name = "FAT")


# Accuracy of the BVE report
cohorts <- get.cohorts(population)
Accs1 <- matrix(0, nrow=23, ncol=length(cohorts))


# breeding cycles

for(index in 0:n_rep){

  population <- breeding.diploid(population,
                                 phenotyping.cohorts = paste0("breedingpopulation_",index,"_M"),
                                 n.observation = pheno_ram,
                                 heritability = heritability,
                                 repeatability = repeatability)

  population <- phenotyping_combi(population, cohorts = paste0("breedingpopulation_",index,"_M"))

  population <- breeding.diploid(population,
                                 phenotyping.cohorts = paste0("breedingpopulation_",index,"_F"),
                                 n.observation = pheno_matrix[2,],
                                 heritability = heritability,
                                 repeatability = repeatability)

  population <- phenotyping_combi(population, cohorts = paste0("breedingpopulation_",index,"_F"))


  # generating offspring

  population <- breeding.diploid(population, breeding.size=N_offspring,
                                 repeat.mating = matrix(c(1, 0.56,
                                                          2, 0.4,
                                                          3, 0.04), byrow=TRUE, ncol = 2),
                                 selection.m.cohort = paste0("breedingpopulation_",index,"_M"),
                                 selection.f.cohort = paste0("breedingpopulation_",index,"_F"),
                                 name.cohort=paste0("offspring_",index),
                                 display.progress = FALSE,
                                 time.point = index)

  # culling offspring: lamb loss rate

  population <- breeding.diploid(population, culling.cohorts = c(paste0("offspring_",index,"_M"),
                                                                 paste0("offspring_",index,"_F")),
                                 culling.share1 = 0.1)


  # phenotyping offspring

  population <- breeding.diploid(population,
                                 phenotyping.cohorts = c(paste0("offspring_",index,"_M"),
                                                         paste0("offspring_",index,"_F")),
                                 phenotyping.class = 0,
                                 n.observation = pheno_matrix[1,],
                                 heritability = heritability,
                                 repeatability = repeatability)

  population <- breeding.diploid(population,
                                 phenotyping.cohorts = c(paste0("offspring_",index,"_M"),
                                                         paste0("offspring_",index,"_F")),
                                 phenotyping.class = 0,
                                 n.observation = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                 heritability = heritability,
                                 share.phenotyped = share_pheno,
                                 repeatability = repeatability)


  population <- phenotyping_combi(population, cohorts = c(paste0("offspring_",index,"_M"),
                                                          paste0("offspring_",index,"_F")))



  # generating cohort "station testing for meat & carcass traits"

  pedigree <- get.pedigree(population, cohorts = paste0("offspring_",index,"_M"))
  pedigree_raw <- get.pedigree(population, cohorts = paste0("offspring_",index,"_M"), raw=TRUE)

  # list of all rams available
  rams <- unique(pedigree[,2])
  # take 8 male lambs per ram for station testing
  stationsheep <- matrix(0, nrow=length(rams)*8, ncol=3)
  for(index1 in 1:length(rams)){
    avail = which(pedigree[,2]==rams[index1])
    take <- sample(avail, 8, replace = if(length(avail)<8){TRUE} else {FALSE})
    stationsheep[1:8 + index1*8 -8, ] <- pedigree_raw[take,1:3]
  }

  population <- breeding.diploid(population,
                                 copy.individual.m = TRUE,
                                 fixed.breeding = cbind(stationsheep, stationsheep, 0),
                                 name.cohort=paste0("stationtesting_",index),
                                 display.progress = FALSE)


  # phenotyping station testing lambs

  population <- breeding.diploid(population,
                                 phenotyping.cohorts = paste0("stationtesting_",index),
                                 n.observation = pheno_station,
                                 heritability = heritability,
                                 repeatability = repeatability)

  population <- phenotyping_combi(population, cohorts = paste0("stationtesting_",index))



  # generation of cohort field testing

  pedigree <- get.pedigree(population, cohorts = paste0("offspring_",index,"_M"))
  pedigree_raw <- get.pedigree(population, cohorts = paste0("offspring_",index,"_M"), raw=TRUE)

  # list of all available rams
  rams <- unique(pedigree[,2])
  # take 20 male lambs per ram for field testing
  fieldsheep <- matrix(0, nrow=length(rams)*20, ncol=3)
  for(index2 in 1:length(rams)){
    avail = which(pedigree[,2]==rams[index2])
    take <- sample(avail, 20, replace = if(length(avail)<20){TRUE} else {FALSE})
    fieldsheep[1:20 + index2*20-20, ] <- pedigree_raw[take,1:3]
  }

  population <- breeding.diploid(population,
                                 copy.individual.m = TRUE,
                                 fixed.breeding = cbind(fieldsheep, fieldsheep, 0),
                                 name.cohort=paste0("fieldtesting_",index),
                                 display.progress = FALSE)

  # phenotyping field testing cohort

  population <- breeding.diploid(population,
                                 phenotyping.cohorts = paste0("fieldtesting_",index),
                                 n.observation = pheno_field,
                                 heritability = heritability,
                                 repeatability = repeatability)

  population <- phenotyping_combi(population, cohorts = paste0("fieldtesting_",index))

  # culling module

  # male
  culling_prob = c(0, 0.1208, 0.4155, 0.2633, 0.1256, 0.029, 0.0386, 0.0024, 0.0024, 0.0024, 0.0, 0.0024, 1)

  # female
  culling_probf = c(0, 0.0864, 0.2038, 0.1998, 0.2135, 0.1224, 0.0924, 0.045, 0.0217, 0.008, 0.0037, 0.0023, 0.0, 0.0003, 0, 0, 0, 0, 0, 0, 0, 0.0003, 0.0003, 1)

  for(index2 in max(0,(index-12)):index){
    population <- breeding.diploid(population, culling.cohort = paste0("offspring_",index2,"_M"), culling.share1 = culling_prob[index-index2+1])
  }
  for(index2 in max(0,(index-23)):index){
    population <- breeding.diploid(population, culling.cohort = paste0("offspring_",index2,"_F"), culling.share1 = culling_probf[index-index2+1])

  }

  # culling module for cohorts still alive

  if(index < 13){
    population <- breeding.diploid(population, culling.cohort = "breedingpopulation_0_M", culling.share1 = culling_prob[index+1])
  }
  if(index < 24){
    population <- breeding.diploid(population, culling.cohort = "breedingpopulation_0_F", culling.share1 = culling_probf[index+1])
  }


  alive_male <- sum(get.class(population, cohorts = paste0("breedingpopulation_",index,"_M"))==0)
  alive_female <- sum(get.class(population, cohorts = paste0("breedingpopulation_",index,"_F"))==0)

  # generating cohort 'licensing' selected from male offspring

  if(!BVE_offspring){
    population <- breeding.diploid(population, breeding.size=c(max(N_licensing, N_allmale-alive_male),0),
                                   selection.size=c(max(N_licensing, N_allmale-alive_male),0),
                                   selection.m.cohorts = paste0("offspring_",index,"_M"),
                                   selection.criteria = "pheno",
                                   multiple.bve.weights.m = selectionindex_licensing,
                                   copy.individual.m=TRUE,
                                   display.progress = FALSE,
                                   name.cohort = paste0("licensing_",index,"_M"))



    if(!keep_pheno){
      db <- get.database(population, cohorts=paste0("licensing_",index,"_M"))
      for(index3 in db[3]:db[4]){
        population$breeding[[db[1]]][[db[2]]][[index3]][[15]][9] <- 0
        population$breeding[[db[1]]][[8+db[2]]][9,index3] <- NA
      }
    }


    # phenotyping cohort 'licensing'

    population <- breeding.diploid(population,
                                   phenotyping.cohorts = c(paste0("licensing_",index,"_M")),
                                   n.observation = pheno_matrix[1,],
                                   heritability = heritability,
                                   repeatability = repeatability)

    population <- phenotyping_combi(population, cohorts = c(paste0("licensing_",index,"_M")))
  }





  # select cohorts used in breeding value estimation

  if( index==0) {
    Cohorten_BVE <- c(paste0("offspring_",index,"_M"),
                      paste0("offspring_",index,"_F"),
                      paste0("stationtesting_",index),
                      paste0("fieldtesting_",index),
                      paste0("breedingpopulation_",index,"_F"),
                      paste0("breedingpopulation_",index,"_M"))}



  if( index==1) {
    Cohorten_BVE <- c(paste0("offspring_",index,"_F"),
                      paste0("offspring_",index-1,"_F"),
                      paste0("stationtesting_",index),
                      paste0("stationtesting_",index-1),
                      paste0("fieldtesting_",index),
                      paste0("fieldtesting_",index-1),
                      paste0("licensing_",index-1,"_M"),
                      paste0("breedingpopulation_",index,"_F"),
                      paste0("breedingpopulation_",index-1,"_F"),
                      paste0("breedingpopulation_",index,"_M"),
                      paste0("breedingpopulation_",index-1,"_M"))}

  if( index>1) {
    Cohorten_BVE <- c(paste0("offspring_",index,"_F"),
                      paste0("offspring_",index-1,"_F"),
                      paste0("offspring_",index-2,"_F"),
                      paste0("stationtesting_",index),
                      paste0("stationtesting_",index-1),
                      paste0("stationtesting_",index-2),
                      paste0("fieldtesting_",index),
                      paste0("fieldtesting_",index-1),
                      paste0("fieldtesting_",index-2),
                      paste0("licensing_",index-1,"_M"),
                      paste0("licensing_",index-2,"_M"),
                      paste0("breedingpopulation_",index,"_F"),
                      paste0("breedingpopulation_",index-1,"_F"),
                      paste0("breedingpopulation_",index-2,"_F"),
                      paste0("breedingpopulation_",index,"_M"),
                      paste0("breedingpopulation_",index-1,"_M"),
                      paste0("breedingpopulation_",index-2,"_M"))}



  BVE_eintragen <- c(paste0("offspring_",index,"_F"),
                     paste0("offspring_",index,"_M"),
                     paste0("stationtesting_",index),
                     paste0("fieldtesting_",index),
                     paste0("breedingpopulation_",index,"_F"),
                     paste0("breedingpopulation_",index,"_M"))


  # Depending on if BVE is performed animals for 'licensing' can be considered in the BVE
  if(!BVE_offspring){
    BVE_eintragen <- c(BVE_eintragen, paste0("licensing_",index,"_M"))
    Cohorten_BVE <- c(Cohorten_BVE, paste0("licensing_",index,"_M"))
  }
  # pedigree based BVE


  population <- breeding.diploid(population, bve=TRUE,
                                 mixblup.bve = mixblup.bve,
                                 mixblup.path = mixblup.path,
                                 bve.cohorts = c(Cohorten_BVE),
                                 bve.insert.cohorts = BVE_eintragen,
                                 relationship.matrix = "pedigree")



  if(BVE_offspring){
    population <- breeding.diploid(population, breeding.size=c(max(N_licensing, N_allmale-alive_male),0),
                                   selection.size=c(max(N_licensing, N_allmale-alive_male),0),
                                   selection.m.cohorts = paste0("offspring_",index,"_M"),
                                   selection.criteria = "bve",
                                   multiple.bve.weights.m = selectionindex_licensing,
                                   copy.individual.m=TRUE,
                                   display.progress = FALSE,
                                   name.cohort = paste0("licensing_",index,"_M"))



    if(!keep_pheno){
      db <- get.database(population, cohorts=paste0("licensing_",index,"_M"))
      for(index3 in db[3]:db[4]){
        population$breeding[[db[1]]][[db[2]]][[index3]][[15]][9] <- 0
        population$breeding[[db[1]]][[8+db[2]]][9,index3] <- NA
      }
    }

    population <- breeding.diploid(population,
                                   phenotyping.cohorts = c(paste0("licensing_",index,"_M")),
                                   n.observation = pheno_matrix[1,],
                                   heritability = heritability,
                                   repeatability = repeatability)

    population <- phenotyping_combi(population, cohorts = c(paste0("licensing_",index,"_M")))
  }

  # selection of sheep for new breeding population

  # male

  bve_acc <- rbind(bve_acc, analyze.bv(population, cohorts=paste0("licensing_",index,"_M"))[[1]][1,])
  population <- breeding.diploid(population, breeding.size=c(max(N_sel, N_allmale-alive_male),0),
                                 selection.size=c(max(N_sel, N_allmale-alive_male),0),
                                 selection.m.cohorts = paste0("licensing_",index,"_M"),
                                 selection.criteria = "bve",
                                 multiple.bve.scale.m = scaling,
                                 multiple.bve.weights.m = selectionindex,
                                 copy.individual.m=TRUE,
                                 display.progress = FALSE,
                                 name.cohort = paste0("offspring_",index,"_M_sel"))



  insert_pheno <- get.pheno(population, cohorts=paste0("offspring_",index,"_M_sel"), use.all.copy = TRUE)
  insert_pheno <- cbind(colnames(insert_pheno), t(insert_pheno))
  population <- insert.pheno(population, phenos = insert_pheno)

  if(!keep_pheno){
    db <- get.database(population, cohorts=paste0("offspring_",index,"_M_sel"))
    for(index3 in db[3]:db[4]){
      population$breeding[[db[1]]][[db[2]]][[index3]][[15]][9] <- 0
      population$breeding[[db[1]]][[8+db[2]]][9,index3] <- NA
    }
  }

  population <- phenotyping_combi(population, cohorts = c(paste0("offspring_",index,"_M_sel")))

  # female

  bve_acc <- rbind(bve_acc, analyze.bv(population, cohorts=paste0("offspring_",index,"_F"))[[1]][1,])

  print(bve_acc)

  population <- breeding.diploid(population, breeding.size=c(0,max(N_self, N_allfemale-alive_female)),
                                 selection.size=c(0,max(N_self, N_allfemale-alive_female)),
                                 selection.criteria = "bve",
                                 multiple.bve.scale.m = scaling,
                                 multiple.bve.weights.f = selectionindex,
                                 selection.f.cohorts = paste0("offspring_",index,"_F"),
                                 copy.individual.f=TRUE,
                                 display.progress = FALSE,
                                 name.cohort = paste0("offspring_",index,"_F_sel"))


  if(!keep_pheno){
    db <- get.database(population, cohorts=paste0("offspring_",index,"_F_sel"))
    for(index3 in db[3]:db[4]){
      population$breeding[[db[1]]][[db[2]]][[index3]][[15]][9] <- 0
      population$breeding[[db[1]]][[8+db[2]]][9,index3] <- NA
    }
  }

  # selection of sheep remaining in the breeding population


  population <- breeding.diploid(population, breeding.size=c(min(N_selp, alive_male),0),
                                 selection.size=c(min(N_selp, alive_male),0),
                                 selection.m.cohorts = paste0("breedingpopulation_",index,"_M"),
                                 selection.criteria = "bve",
                                 multiple.bve.scale.m = scaling,
                                 multiple.bve.weights.f = selectionindex,
                                 copy.individual.m=TRUE,
                                 display.progress = FALSE,
                                 name.cohort = paste0("breedingpopulation_",index,"_M_keep"))

  population <- breeding.diploid(population, breeding.size=c(0,min(N_selfp, alive_female)),
                                 selection.size=c(0,min(N_selfp, alive_female)),
                                 selection.f.cohorts = paste0("breedingpopulation_",index,"_F"),
                                 selection.criteria = "bve",
                                 multiple.bve.scale.m = scale,
                                 multiple.bve.weights.f = selectionindex,
                                 copy.individual.f=TRUE,
                                 display.progress = FALSE,
                                 name.cohort = paste0("breedingpopulation_",index,"_F_keep"))



  # generation of the new breeding population

  population <- breeding.diploid(population, breeding.size = c(N_newM,0),
                                 selection.size=c(N_newM,0),
                                 combine = TRUE,
                                 selection.m.cohorts = c(paste0("breedingpopulation_",index,"_M_keep"),
                                                         paste0("offspring_",index,"_M_sel")),
                                 display.progress = FALSE,
                                 name.cohort = paste0("breedingpopulation_",index+1,"_M"))

  population <- breeding.diploid(population, breeding.size = c(0,N_newF),
                                 selection.size = c(0,N_newF),
                                 combine = TRUE,
                                 selection.f.cohorts = c(paste0("breedingpopulation_",index,"_F_keep"),
                                                         paste0("offspring_",index,"_F_sel")),
                                 display.progress = FALSE,
                                 name.cohort = paste0("breedingpopulation_",index+1,"_F"))

  # phenotyping new breeding population

  population <- breeding.diploid(population,
                                 phenotyping.cohorts = paste0("breedingpopulation_",index+1,"_M"),
                                 n.observation = pheno_ram,
                                 heritability = heritability,
                                 repeatability = repeatability)

  population <- phenotyping_combi(population, cohorts = paste0("breedingpopulation_",index+1,"_M"))

  population <- breeding.diploid(population,
                                 phenotyping.cohorts = paste0("breedingpopulation_",index+1,"_F"),
                                 n.observation = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                 multiple.observation = TRUE,
                                 heritability = heritability,
                                 repeatability = repeatability)



  if(pheno_breeding){

    population <- breeding.diploid(population,
                                   phenotyping.cohorts = paste0("breedingpopulation_",index+1,"_F"),
                                   n.observation = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                                   multiple.observation = TRUE,
                                   share.phenotyped = 1,
                                   heritability = heritability,
                                   repeatability = repeatability)
  }


  ## This is an evaluation of the accuracy of the breeding value estimation
  # as breeding values for the breeding population will be particially overwritten
  # this reports the first BVE for a cohort
  cohorts <- get.cohorts(population)
  Accs2 <- matrix(0, nrow=23, ncol=length(cohorts) - ncol(Accs1))
  for(index2 in (ncol(Accs1)+1):length(cohorts)){
    suppressWarnings({Accs2[1:22,index2-ncol(Accs1)] = analyze.bv(population, cohorts = cohorts[index2])[[1]][1,]})

    a = breeding.diploid(population,
                         selection.m.cohorts = cohorts[index2],
                         selection.criteria = "bve",
                         class.m = c(-1,0),
                         multiple.bve.scale.m = scaling,
                         multiple.bve.weights.m = selectionindex,
                         display.progress = FALSE,
                         verbose = FALSE,
                         export.selected = TRUE)[[1]]

    Accs2[23,index2-ncol(Accs1)] = cor(a[,4], a[,5])
  }
  Accs1 = cbind(Accs1, Accs2)

  population <- phenotyping_combi(population, cohorts = paste0("breedingpopulation_",index+1,"_F"))


}

# Analysis

cohorts <- get.cohorts(population)

BVs <- matrix(0, nrow=22, ncol=length(cohorts))
Phenos <- matrix(0, nrow=22, ncol=length(cohorts))
BVEs <- matrix(0, nrow=22, ncol=length(cohorts))
Accs <- matrix(0, nrow=22, ncol=length(cohorts))
Kins <- matrix(0, nrow=2, ncol=length(cohorts))

BVs_var <- matrix(0, nrow=22, ncol=length(cohorts))
Phenos_var <- matrix(0, nrow=22, ncol=length(cohorts))
BVEs_var <- matrix(0, nrow=22, ncol=length(cohorts))

# Analysis: genomic values & inbreeding ALL cohorts

for(index in 1:length(cohorts)){
  BVs[,index] = rowMeans(get.bv(population, cohorts = cohorts[index]))
  Phenos[,index] = rowMeans(get.pheno(population, cohorts = cohorts[index]))
  BVEs[,index] = rowMeans(get.bve(population, cohorts = cohorts[index]))
  BVs_var[,index] = diag(var(t(get.bv(population, cohorts = cohorts[index]))))
  Phenos_var[,index] = diag(var(t(get.pheno(population, cohorts = cohorts[index]))))
  BVEs_var[,index] = diag(var(t(get.bve(population, cohorts = cohorts[index]))))
  suppressWarnings({Accs[,index] = analyze.bv(population, cohorts = cohorts[index])[[1]][1,]})
  Kins[,index] = kinship.emp.fast(population = population, cohorts = cohorts[index])
}


# extract cohorts

str <- strsplit(cohorts, split="_" )
coh_type <- numeric(length(cohorts))
rep <- numeric(length(cohorts))

for(index in 1:length(cohorts)){

  rep[index] <- str[[index]][[2]]
  str[[index]][[2]] <- "X"
  name <- NULL
  for(index2 in 1:length(str[[index]])){
    if(index2 < length(str[[index]])){
      name <- paste0(name, str[[index]][index2], "_")
    } else{
      name <- paste0(name, str[[index]][index2])
    }
  }
  coh_type[index] <- name

}


# saving simulation results

save(file=paste0("schafsimulation_scenario", scenario, "nr", nr, "size", n, ".RData"), list=c("BVs", "Kins", "Accs","coh_type", "rep", "BVEs", "Phenos", "BVs_var", "BVEs_var", "Phenos_var", "bve_acc"))
if(nr==1) {save(file=paste0("schafsimulation_scenario", scenario, "nr", nr, "size", n, "_pop.RData"), list=c("population"))}

# Check if simulation encountered any warnings / problems
warnings()


