 ## This program includes simulations using MoBPS v.1.11.06 software.
## The software MoBPS and a MoBPS User Manual is available at https://cran.r-project.org/web/packages/MoBPS/index.html
### or at https://github.com/tpook92/MoBPS 

args <- commandArgs(TRUE)
# exemplary setting for args:
 #args = c(1,2,0.1) #For first simulation of baseline scenario with 10% of the individuals

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




#n_rep <- 20 #  number of breeding cycles
scaling <- "bve" # traits in the index are scaled by SD of the estimated breeding values

# Setting up phenotyping classes

pheno <- list()

pheno_matrix <- matrix(c(1, 1,  
                         1, 0, 
                         0, 1, 
                         0, 0,
                         0, 0), byrow = TRUE, ncol =2)


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


heritability = c(0.1, 0.3) 

selectionindex = c(1, 1)




# setting up different scenarios
keep_pheno <- FALSE
pheno_breeding <- FALSE



selcrit <- "bve"
share_geno_fem <- 0
share_geno_lic <- 0
matrix <- "pedigree" 
ssGBLUP <- FALSE 
BVE <- TRUE 




# Ped
if(scenario==1){
  BVE1 <- TRUE
  share_geno_fem1 <- 0
  share_geno_lic1 <- 0
  ssGBLUP1 <- FALSE
  matrix1 <- "pedigree"
  selcrit1 <- "bve"
  N_preselectmale <- ceiling(6466*n) 
  N_preselectfemale <- ceiling(6466*n)
}



# GSTop25
if(scenario==2){
 BVE1 <- TRUE
 share_geno_fem1 <- 0
share_geno_lic1 <- 1
 matrix1 <- "vanRaden"
 ssGBLUP1 <- TRUE
 selcrit1 <- "bve"
 N_preselectmale <- ceiling(1617*n)
 N_preselectfemale <- ceiling(6466*n)
}


# GSTop50
if(scenario==3){
  BVE1 <- TRUE
  share_geno_fem1 <- 0
  share_geno_lic1 <- 1
 matrix1 <- "vanRaden"
  ssGBLUP1 <- TRUE
  selcrit1 <- "bve"
  N_preselectmale <- ceiling(3233*n)
  N_preselectfemale <- ceiling(6466*n)
}

# GS100
if(scenario==4){
  BVE1 <- TRUE
  share_geno_fem1 <- 0
  share_geno_lic1 <- 1
  matrix1 <- "vanRaden"
  ssGBLUP1 <- TRUE
  selcrit1 <- "bve"
  N_preselectmale <- ceiling(6466*n)
  N_preselectfemale <- ceiling(6466*n)
}


# GS100+Top25
if(scenario==5){
  BVE1 <- TRUE
  share_geno_fem1 <- 1
  share_geno_lic1 <- 1
  matrix1 <- "vanRaden"
  ssGBLUP1 <- TRUE
  selcrit1 <- "bve"
  N_preselectmale <- ceiling(6466*n)
  N_preselectfemale <- ceiling(1617*n)
}

# GS100+Top50
if(scenario==6){
  BVE1 <- TRUE
  share_geno_fem1 <- 1
  share_geno_lic1 <- 1
  matrix1 <- "vanRaden"
  ssGBLUP1 <- TRUE
  selcrit1 <- "bve"
  N_preselectmale <- ceiling(6466*n)
  N_preselectfemale <- ceiling(3233*n)
}

# GS100+100
if(scenario==7){
  BVE1 <- TRUE
  share_geno_fem1 <- 1
  share_geno_lic1 <- 1
  matrix1 <- "vanRaden"
  ssGBLUP1 <- TRUE
  selcrit1 <- "bve" 
  N_preselectmale <- ceiling(6466*n) 
  N_preselectfemale <- ceiling(6466*n)
}


set.seed(nr)


# LD build-up for realistic founder population
# insert real genomic data

dataset_out <- founder.simulation(nindi = N_indi,
                                  vcf = dataset_path,  ## Data is available upon request
                                  n.gen=5, 
                                  big.output = TRUE,
                                  bpcm.conversion = 1000000)

map <- dataset_out[[2]]

dataset <- dataset_out[[1]]

#rm(dataset_out)

set.seed(nr)

##genetic correlation matrix
G <- matrix(c(1, -0.1, 
              -0.1, 1), ncol = 2)





# Initialization of founder population

population <- creating.diploid(dataset= dataset, sex.s = c(rep(1,N_newM), rep(2,N_newF)),
                               name.cohort="breedingpopulation_0",
                               n.additive = c(rep(1000,2)),
                               map = map,
                               trait.name = c("health", "production"), 
                               shuffle.traits = 1:2,
                               shuffle.cor = G)

rm(dataset)

summary(population)




# trait standardization

population <- bv.standardization(population, mean.target = 100,
                                 var.target = 10)


# Accuracy of the BVE report
cohorts <- get.cohorts(population)
Accs1 <- matrix(0, nrow=3, ncol=length(cohorts))





# breeding cycles


for(index in 0:20){

  if(index < 10){
    selcrit_used = selcrit
    share_geno_fem_used = share_geno_fem
    share_geno_lic_used = share_geno_lic
    matrix_used = matrix
    ssGBLUP_used = ssGBLUP
    BVE_used = BVE
  } else{
    selcrit_used = selcrit1
    share_geno_fem_used = share_geno_fem1
    share_geno_lic_used = share_geno_lic1
    matrix_used = matrix1
    ssGBLUP_used = ssGBLUP1
    BVE_used = BVE1
  }

 ###initial genotyping of breeding rams in breeding cycle 10
  if( index==10){
   population <- breeding.diploid(population,
                                 genotyped.cohorts = paste0("breedingpopulation_",index,"_M"),
                                 genotyped.share = 1)
  }

  ##generating breeding ram cohort

  population <- breeding.diploid(population,
                                 phenotyping.cohorts = paste0("breedingpopulation_",index,"_M"),
                                 n.observation = pheno_matrix[1,],
                                 heritability = heritability)


  ##generating breeding ewes cohort

  population <- breeding.diploid(population,
                                 phenotyping.cohorts = paste0("breedingpopulation_",index,"_F"),
                                 n.observation = pheno_matrix[1,],
                                 heritability = heritability)





  # generating male and female lambs

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





  # select cohorts used in breeding value estimation
  if( index==0) {
    Cohorten_BVE <- c(paste0("breedingpopulation_",index,"_F"),
                      paste0("breedingpopulation_",index,"_M"))}


  if( index > 0) {
    Cohorten_BVE <- c(paste0("offspring_",index,"_F"),
                      paste0("offspring_",index,"_M"),
                      paste0("preselectmale_",max(0,(index-10)):(index-1),"_M"),
                      paste0("preselectfemale_",max(0,(index-10)):(index-1),"_F"),
                      paste0("breedingpopulation_",max(0,(index-10)):(index),"_F"),
                      paste0("breedingpopulation_",max(0,(index-10)):(index),"_M"))
  }
  
  
  
  BVE_eintragen <- c(paste0("offspring_",index,"_F"),
                     paste0("offspring_",index,"_M"),
                     paste0("breedingpopulation_",index,"_F"),
                     paste0("breedingpopulation_",index,"_M"))
  
  
  # breeding value estimation


  population <- breeding.diploid(population, bve=BVE_used,

                                 mixblup.bve = mixblup.bve,
                                 mixblup.path = mixblup.path,

                                 bve.cohorts = c(Cohorten_BVE),
                                 bve.insert.cohorts = BVE_eintragen,
                                 relationship.matrix = matrix_used,
                                 singlestep.active = ssGBLUP_used)


  ####generating cohort 'preselected_males'

  population <- breeding.diploid(population, breeding.size=c(max(N_preselectmale, N_allmale-alive_male),0),
                                 selection.size=c(max(N_preselectmale, N_allmale-alive_male),0),
                                 selection.m.cohorts = paste0("offspring_",index,"_M"),
                                 selection.criteria = selcrit_used,
                                 selection.index.scale.m = scaling,
                                 multiple.bve.weights.m = selectionindex,
                                 copy.individual.m=TRUE,
                                 display.progress = FALSE,
                                 name.cohort = paste0("preselectmale_",index,"_M")) 
  
  ##gentotyping preselected males 
  
  population <- breeding.diploid(population,
                                 phenotyping.cohorts = c(paste0("preselectmale_",index,"_M")),
                                 phenotyping.class = 0,
                                 n.observation = pheno_matrix[5,],
                                 
                                 heritability = heritability,
                                 genotyped.share = share_geno_lic_used,
                                 genotyped.cohorts = c(paste0("preselectmale_",index,"_M")))
  
  






  ##generating cohort 'preselected females'

  population <- breeding.diploid(population, breeding.size=c(0,max(N_preselectfemale, N_allfemale-alive_female)),
                                 selection.size=c(0,max(N_preselectfemale, N_allfemale-alive_female)),
                                 selection.criteria = selcrit_used,
                                 selection.index.scale.m = scaling,
                                 multiple.bve.weights.f = selectionindex,
                                 selection.f.cohorts = paste0("offspring_",index,"_F"),
                                 copy.individual.f=TRUE,
                                 display.progress = FALSE,
                                 name.cohort = paste0("preselectfemale_",index,"_F"))
  
  # genotyping top preselected females 
  
  population <- breeding.diploid(population,
                                 phenotyping.cohorts = c(paste0("preselectfemale_",index,"_F")),
                                 phenotyping.class = 0,
                                 n.observation = pheno_matrix[5,],
                                 
                                 heritability = heritability,
                                 genotyped.share = share_geno_fem_used,
                                 genotyped.cohorts = c(paste0("preselectfemale_",index,"_F")))
  
  
  ###2nd breeding value estimation WITH genetic information of selection candidates
  
  if( index==0) {
    Cohorten_BVE <- c(paste0("breedingpopulation_",index,"_F"),
                      paste0("breedingpopulation_",index,"_M"),
                      paste0("preselectmale_",index,"_M"),
                      paste0("preselectfemale_",index,"_F"))}
  
  
  if( index > 0) {
    Cohorten_BVE <- c(paste0("offspring_",index,"_F"),
                      paste0("offspring_",index,"_M"),
                      paste0("preselectmale_",max(0,(index-10)):(index-1),"_M"),
                      paste0("preselectfemale_",max(0,(index-10)):(index-1),"_F"),
                      paste0("breedingpopulation_",max(0,(index-10)):(index),"_F"),
                      paste0("breedingpopulation_",max(0,(index-10)):(index),"_M"))
  }
  
  BVE_eintragen <- c(paste0("offspring_",index,"_F"),
                     paste0("offspring_",index,"_M"),
                     paste0("breedingpopulation_",index,"_F"),
                     paste0("breedingpopulation_",index,"_M"),
                     paste0("preselectmale_",index,"_M"),
                     paste0("preselectfemale_",index,"_F"))
  
  population <- breeding.diploid(population, bve=BVE_used,
                                 
                                 mixblup.bve = mixblup.bve,
                                 mixblup.path = mixblup.path,
                                 
                                 bve.cohorts = c(Cohorten_BVE),
                                 bve.insert.cohorts = BVE_eintragen,
                                 relationship.matrix = matrix_used,
                                 singlestep.active = ssGBLUP_used)




  # selection of sheep for new breeding population



  bve_acc <- rbind(bve_acc, analyze.bv(population, cohorts=paste0("preselectmale_",index,"_M"))[[1]][1,])

  population <- breeding.diploid(population, breeding.size=c(max(N_sel, N_allmale-alive_male),0),
                                 selection.size=c(max(N_sel, N_allmale-alive_male),0),
                                 selection.m.cohorts = paste0("preselectmale_",index,"_M"),
                                 selection.criteria = selcrit_used,
                                 selection.index.scale.m = scaling,
                                 multiple.bve.weights.m = selectionindex,
                                 copy.individual.m=TRUE,
                                 display.progress = FALSE,
                                 name.cohort = paste0("offspring_",index,"_M_sel"))
  

  ##phenotyping new breeding rams
  
  population <- breeding.diploid(population,
                                 phenotyping.cohorts = c(paste0("offspring_",index,"_M_sel")),
                                 phenotyping.class = 0,
                                 n.observation = pheno_matrix[1,],
                                 heritability = heritability)
                                 







  # female: cohort preselected females to cohort new breeding ewes

   bve_acc <- rbind(bve_acc, analyze.bv(population, cohorts=paste0("preselectfemale_",index,"_F"))[[1]][1,])

  print(bve_acc)

  population <- breeding.diploid(population, breeding.size=c(0,max(N_self, N_allfemale-alive_female)),
                                 selection.size=c(0,max(N_self, N_allfemale-alive_female)),
                                 selection.criteria = selcrit_used,
                                 selection.index.scale.m = scaling,
                                 multiple.bve.weights.f = selectionindex,
                                 selection.f.cohorts = paste0("preselectfemale_",index,"_F"),
                                 copy.individual.f=TRUE,
                                 display.progress = FALSE,
                                 name.cohort = paste0("offspring_",index,"_F_sel"))
  
  ##phenotyping new breeding ewes
  
  population <- breeding.diploid(population,
                                 phenotyping.cohorts = c(paste0("offspring_",index,"_F_sel")),
                                 phenotyping.class = 0,
                                 n.observation = pheno_matrix[1,],
                                 
                                 heritability = heritability)




  # selection of sheep remaining in the breeding population


  population <- breeding.diploid(population, breeding.size=c(min(N_selp, alive_male),0),
                                 selection.size=c(min(N_selp, alive_male),0),
                                 selection.m.cohorts = paste0("breedingpopulation_",index,"_M"),
                                 selection.criteria = selcrit_used,
                                 selection.index.scale.m = scaling,
                                 
                                 multiple.bve.weights.f = selectionindex,
                                 copy.individual.m=TRUE,
                                 display.progress = FALSE,
                                 name.cohort = paste0("breedingpopulation_",index,"_M_keep"))

  population <- breeding.diploid(population, breeding.size=c(0,min(N_selfp, alive_female)),
                                 selection.size=c(0,min(N_selfp, alive_female)),
                                 selection.f.cohorts = paste0("breedingpopulation_",index,"_F"),
                                 selection.criteria = selcrit_used,
                                 selection.index.scale.m = scaling,
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
  
  ## This is an evaluation of the accuracy of the breeding value estimation
  # as breeding values for the breeding population will be particially overwritten
  # this reports the first BVE for a cohort
  cohorts <- get.cohorts(population)
  Accs2 <- matrix(0, nrow=3, ncol=length(cohorts) - ncol(Accs1))
  for(index2 in (ncol(Accs1)+1):length(cohorts)){
    suppressWarnings({Accs2[1:2,index2-ncol(Accs1)] = analyze.bv(population, cohorts = cohorts[index2])[[1]][1,]})
    
    a = breeding.diploid(population,
                         selection.m.cohorts = cohorts[index2],
                         selection.criteria = "bve",
                         class.m = c(-1,0),
                         selection.index.scale.m = scaling,
                         multiple.bve.weights.m = selectionindex,
                         display.progress = FALSE,
                         verbose = FALSE,
                         export.selected = TRUE)[[1]]
    
    Accs2[3,index2-ncol(Accs1)] = cor(a[,4], a[,5]) 
  }
  Accs1 = cbind(Accs1, Accs2)
  
  
}


# Analysis

cohorts <- get.cohorts(population)
BVs <- matrix(0, nrow=2, ncol=length(cohorts))
Phenos <- matrix(0, nrow=2, ncol=length(cohorts))
BVEs <- matrix(0, nrow=2, ncol=length(cohorts))
Accs <- matrix(0, nrow=2, ncol=length(cohorts))
Kins <- matrix(0, nrow=2, ncol=length(cohorts))

BVs_var <- matrix(0, nrow=2, ncol=length(cohorts))
Phenos_var <- matrix(0, nrow=2, ncol=length(cohorts))
BVEs_var <- matrix(0, nrow=2, ncol=length(cohorts))



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

if(is_torsten){
  setwd("/home/WUR/pook001/Sheep_results/")
}

save(file=paste0("schafsimulation_scenario", scenario, "nr", nr, "size", n, ".RData"), list=c("BVs", "Kins", "Accs","coh_type", "rep", "BVEs", "Phenos", "BVs_var", "BVEs_var", "Phenos_var", "bve_acc", "Accs1"))
if(nr==1) {save(file=paste0("schafsimulation_scenario", scenario, "nr", nr, "size", n, "_pop.RData"), list=c("population"))}

# Check if simulation encountered any warnings / problems
warnings()


