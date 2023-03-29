set.seed(42)
library(MoBPS)
library(MoBPSmaps)

# Simulation of the LD build up
dataset <- founder.simulation(nindi=100, map = MoBPSmaps::map_sheep2,
                              nfinal = 60+120, n.gen = 50)

# Use a genomic map from the MoBPSmap R-package
# Import genomic data from the dataset-matrix ((each individual has 2 haplotypes))
population <- creating.diploid(dataset = dataset[,1:100], nindi = 50, sex.quota=0, map=MoBPSmaps::map_sheep2,
                               name.cohort = "1yearRams_0")


population <- creating.diploid(dataset = dataset[,101:120], population = population, nindi = 10,
                               sex.quota=0, name.cohort = "2yearRams_0")

population <- creating.diploid(dataset = dataset[,121:220], population = population, nindi = 50,
                               share.genotyped = 1, sex.quota=1, name.cohort = "1yearEwes_0")

population <- creating.diploid(dataset = dataset[,221:300], population = population, nindi = 40,
                               share.genotyped = 1, sex.quota=1, name.cohort = "2yearEwes_0")

population <- creating.diploid(dataset = dataset[,301:360], population = population, nindi = 30,
                               share.genotyped = 1, sex.quota=1, name.cohort = "3yearEwes_0")

summary(population)

# Use of two correlated traits

population <- creating.trait(population, n.additive = c(1000,1000), mean.target = 100, var.target = c(30,20),
                             trait.cor = matrix(c(1, 0.2, 0.2,1), nrow=2))

# first column: Size of the litter
# second column: probability of the litter size to occur

litter_size <- matrix(c(2, 0.5,
                        3, 0.3,
                        4, 0.2), ncol=2, byrow=TRUE)

### We want to execute all the code below multiple times in succession
# This code should work generally to always to the last available cycles
# in case a cohort from cycle 0 is used below replace this with  paste0(COHORTNAME, "_", index)
# in case a cohort from cycle 1 is used below replace this with  paste0(COHORTNAME, "_", index + 1)

# For coding we would recommend to not directly code this with a loop but set index to 1
# Write the code for this case and afterwards add the loop

for(index in 1:20){

  # phenotypes for the second trait are only generated for selected cohorts

  population <- breeding.diploid(population, phenotyping.gen = index,
                                 heritability = c(0.3,0.2), n.observation = c(1,0))

  population <- breeding.diploid(population,
                                 phenotyping.cohorts = paste0(c("2yearRams_", "2yearEwes_", "3yearEwes_"), index-1),
                                 heritability = c(0.3,0.2), n.observation = c(0,1))


  # As some individuals are non genotyped a genomic breeding value estimation does not work
  # Either use a pedigree-based evaluation or single-step
  # Otherwise non-genotyped individuals will automatically be removed for the estimation

  population <- breeding.diploid(population, bve = TRUE, bve.gen = index,
                                 bve.cohorts = paste0(c("1yearRams_", "2yearRams_", "1yearEwes_", "2yearEwes_", "3yearEwes_"), index-1))


  # Generation of new individuals via reproduction

  population <- breeding.diploid(population, breeding.size = c(50,50),
                                 selection.m.cohorts = paste0("2yearRams_", index-1),
                                 selection.f.cohorts = c(paste0("2yearEwes_",index-1), paste0("3yearEwes_",index-1)),
                                 name.cohort = paste0("NewOffspring_",index),
                                 repeat.mating = litter_size,
                                 share.genotyped = 0
  )

  population <- breeding.diploid(population, breeding.size = c(50,0),
                                 selection.size = c(50,0),
                                 copy.individual.m = TRUE,
                                 selection.m.cohorts = paste0("NewOffspring_",index, "_M"),
                                 name.cohort = paste0("1yearRams_", index),
                                 add.gen=index+1,
                                 added.genotyped = 1
  )

  population <- breeding.diploid(population, breeding.size = c(0,50),
                                 selection.size = c(0,50),
                                 copy.individual.f = TRUE,
                                 selection.f.cohorts = paste0("NewOffspring_",index, "_F"),
                                 name.cohort = paste0("1yearEwes_", index),
                                 add.gen=index+1,
                                 added.genotyped = 0.5
  )


  # Generation of new individuals via selection

  population <- breeding.diploid(population, breeding.size = c(10,0), selection.size = c(10,0),
                                 selection.m.cohorts = paste0("1yearRams_", index-1),
                                 copy.individual.m = TRUE,
                                 selection.criteria = "bve",
                                 add.gen=index+1,
                                 name.cohort = paste0("2yearRams_", index))


  population <- breeding.diploid(population, breeding.size = c(0,40),
                                 selection.size = c(0,40),
                                 selection.f.cohorts = paste0("1yearEwes_", index-1),
                                 copy.individual.f = TRUE,
                                 selection.criteria = "random",
                                 add.gen=index+1,
                                 name.cohort = paste0("2yearEwes_", index)
  )

  population <- breeding.diploid(population, breeding.size = c(0,30),
                                 selection.size = c(0,30),
                                 selection.f.cohorts = paste0("2yearEwes_", index-1),
                                 copy.individual.f = TRUE,
                                 selection.criteria = "random",
                                 add.gen=index+1,
                                 name.cohort = paste0("3yearEwes_", index)
  )

}

summary(population)

rowMeans(get.bv(population, gen=21))
cor(t(get.bv(population, gen=21)))

rowMeans(get.bv(population, cohorts="1yearRams_0"))
rowMeans(get.bv(population, cohorts="1yearRams_20"))
rowMeans(get.bv(population, cohorts="2yearRams_20"))


cohorts_to_analyze <- paste0("1yearRams_", 1:20)



# By increasing ibd.obs and hbd.obs you can increase the accuracy of the calculation of kinship between individuals and within//inbreeding
inbreeding <- numeric(20)
for(index in 1:20){
  inbreeding[index] <- 2* (kinship.emp.fast(population, cohorts = cohorts_to_analyze[index])[2] -0.5)
}

plot(inbreeding, ylab="inbreeding rate", xlab = "generation", type="l", lwd=2,
     main="Inbreeding rates for 1yearRams")
