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

summary(population)
get.cohorts(population, extended = TRUE)

# Generation of the trait
# var.target is controlling the genetic variance! not the phenotypic variance!

population <- creating.trait(population, n.additive = 1000, mean.target = 100, var.target = 30, trait.name = "Meat")

get.qtl(population)
get.qtl.effects(population)

# Generation of phenotypes
# On default the first generation will be used to approximate the required residual variance
# Other cohorts can be selection with sigma.e.gen/database/cohorts
# "all" is a quick access for phenotyping all individuals
# "non_obs" for all previously not phenotyped individuals
# alternatively direct selection via phenotyping.gen/database/cohorts

population <- breeding.diploid(population, phenotyping = "all", heritability = 0.3)

get.pheno(population, cohorts="1yearRams")

# Generation of the new individuals via reproduction:

population <- breeding.diploid(population, breeding.size = c(50,0),
                               selection.m.cohorts = "2yearRams",
                               selection.f.cohorts = c("2yearEwes", "3yearEwes"),
                               name.cohort = "1yearRamsNext"
                               )

# The parameter add.gen is used to make sure that all individuals of one cycle are assign to
# the same generation of individuals

population <- breeding.diploid(population, breeding.size = c(0,50),
                               selection.m.cohorts = "2yearRams",
                               selection.f.cohorts = c("2yearEwes", "3yearEwes"),
                               name.cohort = "1yearEwesNext",
                               add.gen=2
)

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


get.cohorts(population, extended = TRUE)


summary(population)
# The tools sees that 360 individuals were generated
# This counts copies of individuals multiple times!
# However, each individual is also assigned with a unique id
# e.g. As 2yearRamsNext are just a copy of 10 individuals from 1yearRams the ids are matching
# In a breeding value estimation, only one copy of an individual is used
get.id(population, cohorts = "1yearRams")
get.id(population, cohorts="2yearRamsNext")

# Examplary breeding value estimation
population <- breeding.diploid(population, bve = TRUE, bve.gen=1)

# Set remove.effect.position to remove real QTLs from the breeding value estimation
population <- breeding.diploid(population, bve = TRUE, bve.gen=1, remove.effect.position = TRUE)

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

