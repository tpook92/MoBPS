# setting a random seed makes simulations reproducible
set.seed(42)
library(MoBPS)

# Generation of the founder population
# To have cohort names it is required to generate each group separately


population <- creating.diploid(nsnp = 5000, nindi = 50, sex.quota=0,
                               chr.nr = 5, chromosome.length = 1, name.cohort = "1yearRams")

population <- creating.diploid(population = population, freq = "same", nindi = 10,
                               sex.quota=0, name.cohort = "2yearRams")

population <- creating.diploid(population = population, freq = "same", nindi = 50,
                               sex.quota=1, name.cohort = "1yearEwes")

population <- creating.diploid(population = population, freq = "same", nindi = 40,
                               sex.quota=1, name.cohort = "2yearEwes")

population <- creating.diploid(population = population, freq = "same", nindi = 30,
                               sex.quota=1, name.cohort = "3yearEwes")

# get an overview of the population ((in more complex simulations also for validation in results make sense))
summary(population)
get.cohorts(population, extended = TRUE)[,1:4]

# Generation of the trait
# var.target is controlling the genetic variance! not the phenotypic variance!

population <- creating.trait(population, n.additive = 1000, mean.target = 100,
                             var.target = 30, trait.name = "Meat")

# which SNPs have underlying effects (( you can exclude QTLs from the BVE later if you dont want them in your evaluations ))
get.qtl(population)
# if output format is unclear try ?get.qtl.effects or str()
get.qtl.effects(population)
# just looking at the first 5 QTLs
get.qtl.effects(population)[[1]][[1]][1:5,]


# Generation of phenotypes
# On default the first generation will be used to approximate the required residual variance
# Other cohorts can be selection with sigma.e.gen/database/cohorts
# "all" is a quick access for phenotyping all individuals
# "non_obs" for all previously not phenotyped individuals
# alternatively direct selection via phenotyping.gen/database/cohorts

population <- breeding.diploid(population, phenotyping = "all", heritability = 0.3)

var(t(get.bv(population, gen = 1)))

# if trait is generated/standardized when only 1yearRams exists QTL effect would be scaled to obtain target in that cohort
var(t(get.bv(population, cohorts = "1yearRams")))
get.pheno(population, cohorts="1yearRams")

#population = bv.standardization(population, cohorts = "1yearRams", var.target = 30,
#                                adapt.pheno = TRUE)
var(t(get.bv(population, cohorts = "1yearRams")))
get.pheno(population, cohorts="1yearRams")

# Generation of the new individuals via reproduction:

#
#population <- breeding.diploid(population, breeding.size = c(50,50),
#                               selection.m.cohorts = "2yearRams",
#                               selection.f.cohorts = c("2yearEwes", "3yearEwes"),
#                               name.cohort = c("1yearRamsNext", "1yearEwesNext"))

population <- breeding.diploid(population, breeding.size = c(50,0),
                               selection.m.cohorts = "2yearRams",
                               selection.f.cohorts = c("2yearEwes", "3yearEwes"),
                               name.cohort = "1yearRamsNext")

# The parameter add.gen is used to make sure that all individuals of one cycle are assign to
# the same generation of individuals

population <- breeding.diploid(population, breeding.size = c(0,50),
                               selection.m.cohorts = "2yearRams",
                               selection.f.cohorts = c("2yearEwes", "3yearEwes"),
                               name.cohort = "1yearEwesNext",
                               add.gen=2)

# Generation of the new individuals via selection:

population <- breeding.diploid(population, breeding.size = c(10,0), selection.size = c(10,0),
                               selection.m.cohorts = "1yearRams",
                               copy.individual.m = TRUE,
                               selection.criteria = "pheno",
                               add.gen=2,
                               name.cohort = "2yearRamsNext")

population <- breeding.diploid(population, breeding.size = c(0,40),
                               selection.size = c(0,40),
                               selection.f.cohorts = "1yearEwes",
                               copy.individual.f = TRUE,
                               selection.criteria = "random",
                               add.gen=2,
                               name.cohort = "2yearEwesNext")

population <- breeding.diploid(population, breeding.size = c(0,30),
                               selection.size = c(0,30),
                               selection.f.cohorts = "2yearEwes",
                               copy.individual.f = TRUE,
                               selection.criteria = "random",
                               add.gen=2,
                               name.cohort = "3yearEwesNext")


get.cohorts(population, extended = TRUE)[,1:4]

#################################################################################################
####### END of the Task - everthing below is just showing of some additional MoBPS functionality
#################################################################################################

# exemplary code to use different sires for the generation of the Ewes than for the Rams
# In our example all sires were used for generation of the Rams that why I am reducing the sire_mating1 group to just the first 5 (([1:5]))
sire_mating1 = unique(get.pedigree(population, cohorts = c("1yearRamsNext"))[,2])[1:5]

sires_mating2 = group.diff(population, cohorts = "2yearRams", remove.database = get.database(population, id = sire_mating1))

population <- breeding.diploid(population, breeding.size = c(0,50),
                               selection.m.database = sires_mating2,
                               selection.f.cohorts = c("2yearEwes", "3yearEwes"),
                               name.cohort = "1yearEwesNext_alt",
                               add.gen=2)

## some example on how to use database ((instead of cohorts / gen))

# If you specifically want to use the 53th male(1) in generation 2
database = cbind(2,1,53,53)

get.pheno(population, database = database)
get.pheno(population, cohorts = "2yearRamsNext")

# A single individual can have different names ("M51_2") assigned to it when there is multiple copies of that individual
# in different cohorts. the ID is unique to a single individum and call copies have the same ID.
get.id(population, cohorts = "1yearRams")
get.id(population, cohorts = "2yearRamsNext")


## generation of a trait with 2 QTLs, with non-additive effects; position is provided in colums 1 (SNP) / 2 chromosome
## column 3: effect AA, 4: effect AB, 5: effect BB
real.bv.add = cbind(NA,NA,c(0,0), c(-2,3), c(5,6))
pop1 = creating.trait(population, real.bv.add = real.bv.add)
get.qtl.effects(pop1)[[1]][[1]][,1:5]
get.qtl.effects(pop1)[[1]][[2]][,1:5]


# Exemplary generating chromosomes with difference size and number of SNPs
pop1 <- creating.diploid(nsnp = c(100,200,500), nindi = 50, sex.quota=0,
                               chr.nr = 3, chromosome.length = c(5,3,27), name.cohort = "1yearRams")
## there was a typo here that stated chr.nr = 5 ((while using nsnp / chromosome.length with just 3 entries))

map = get.map(pop1)
pop1 = creating.diploid(map = map, nindi = 10)


# The tools sees that 360 individuals were generated
# This counts copies of individuals multiple times!
# However, each individual is also assigned with a unique id
# e.g. As 2yearRamsNext are just a copy of 10 individuals from 1yearRams the ids are matching
# In a breeding value estimation, only one copy of an individual is used
get.id(population, cohorts = "1yearRams")
get.id(population, cohorts="2yearRamsNext")

# Examplary breeding value estimation
population <- breeding.diploid(population, bve = TRUE, bve.gen=1,
                               bve.all.genotyped = TRUE)

# Set remove.effect.position to remove real QTLs from the breeding value estimation
population <- breeding.diploid(population, bve = TRUE, bve.gen=1, remove.effect.position = TRUE,
                               bve.all.genotyped = TRUE)

# Add a genotyping array with only 2500 Markers
population <- add.array(population, marker.included = rep(c(0,1),2500), array.name = "Small_Array")

# The first genotyping array always includes all markers
# When performing new genotyping we can specify which genotyping array is used
# You can do this via the array name or (1 for the array with all markers,
# 2 for the first generated array, 3 for the second etc.

population <- breeding.diploid(population, breeding.size = c(50,0),
                               selection.m.cohorts = "2yearRamsNext",
                               selection.f.cohorts = "2yearEwesNext",
                               share.genotyped = 0.5,
                               genotyped.array = "Small_Array")


# Extract genotypes for the new generation:
get.geno(population, gen=3)
get.geno(population, gen=3, non.genotyped.as.missing = TRUE)


# Genotyping already existing indidivials:
population <- breeding.diploid(population, breeding.size = c(50,0),
                               genotyped.cohorts =  "2yearRamsNext",
                               share.genotyped = 0.5,
                               genotyped.array = "Small_Array",
                               time.point = 1)

# Overview on time points whenever an individual is switching cohorts:
get.snapshot(population, cohorts = "1yearRams")
# Overview on time points whenever the status of some individuals changes:
get.snapshot.single(population, cohorts = "1yearRams")

