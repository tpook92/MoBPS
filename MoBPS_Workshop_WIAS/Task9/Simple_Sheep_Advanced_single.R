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



# phenotypes for the second trait are only generated for selected cohorts
population <- breeding.diploid(population, phenotyping.gen = 1,
                               heritability = c(0.3,0.2), n.observation = c(1,0))


population <- breeding.diploid(population, phenotyping.cohorts = c("2yearRams_0", "2yearEwes_0", "3yearEwes_0"),
                               heritability = c(0.3,0.2), n.observation = c(0,1))



# As some individuals are non genotyped a genomic breeding value estimation does not work
# Either use a pedigree-based evaluation or single-step
# Otherwise non-genotyped individuals will automatically be removed for the estimation

population <- breeding.diploid(population, bve = TRUE,
                               bve.cohorts = c("1yearRams_0", "2yearRams_0", "1yearEwes_0", "2yearEwes_0", "3yearEwes_0"))


# Generation of new individuals via reproduction

population <- breeding.diploid(population, breeding.size = c(50,50),
                               selection.m.cohorts = "2yearRams_0",
                               selection.f.cohorts = c("2yearEwes_0", "3yearEwes_0"),
                               name.cohort = "NewOffspring_1",
                               repeat.mating = litter_size,
                               share.genotyped = 0)

population <- breeding.diploid(population, breeding.size = c(50,0),
                               selection.size = c(50,0),
                               copy.individual.m = TRUE,
                               selection.m.cohorts = "NewOffspring_1_M",
                               name.cohort = "1yearRams_1",
                               add.gen=2,
                               added.genotyped = 1)

population <- breeding.diploid(population, breeding.size = c(0,50),
                               selection.size = c(0,50),
                               copy.individual.f = TRUE,
                               selection.f.cohorts = "NewOffspring_1_F",
                               name.cohort = "1yearEwes_1",
                               add.gen=2,
                               added.genotyped = 1)


# Generation of new individuals via selection

population <- breeding.diploid(population, breeding.size = c(10,0), selection.size = c(10,0),
                               selection.m.cohorts = "1yearRams_0",
                               copy.individual.m = TRUE,
                               selection.criteria = "bve",
                               add.gen=2,
                               name.cohort = "2yearRams_1")

population <- breeding.diploid(population, breeding.size = c(0,40),
                               selection.size = c(0,40),
                               selection.f.cohorts = "1yearEwes_0",
                               copy.individual.f = TRUE,
                               selection.criteria = "random",
                               add.gen=2,
                               name.cohort = "2yearEwes_1")

population <- breeding.diploid(population, breeding.size = c(0,30),
                               selection.size = c(0,30),
                               selection.f.cohorts = "2yearEwes_0",
                               copy.individual.f = TRUE,
                               selection.criteria = "random",
                               add.gen=2,
                               name.cohort = "3yearEwes_1")

