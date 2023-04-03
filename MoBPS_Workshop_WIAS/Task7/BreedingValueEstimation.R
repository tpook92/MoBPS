set.seed(42)
library(MoBPS)

# Generation of the founder population
population <- creating.diploid(nindi=100, nsnp=25000, chr.nr=5,
                               n.additive = 1000, chromosome.length = 3)
summary(population)

array2 <- array3 <- array4 <- rep(FALSE, 25000)
array2[sample(1:25000, 10000)] <- TRUE

# These arrays were actually not put in the task - more to get a better feeling on when predictions
# are actually breaking down.
array3[sample(1:25000, 2000)] <- TRUE
array4[sample(1:25000, 100)] <- TRUE

population <- add.array(population, marker.included = array2)
population <- add.array(population, marker.included = array3)
population <- add.array(population, marker.included = array4)

# If no input for selection.size is provided MoBPS will automatically use all individuals from the previous generation
# LD build-up
# you could for example use selection.m.database = cbind(index,1) &  selection.f.database = cbind(index,2) or cohort names

for(index in 1:10){

  population <- breeding.diploid(population, breeding.size = 100,
                                 share.genotyped = 1)

}


# example of how to build a database with the best 10 individuals to subsequently genotype them
# usually run a pedigree-blup or similar before:
#best_individuals = breeding.diploid(population, selection.size = c(10,0),
#                              export.selected = TRUE)
#best_individuals[[1]][,1:3]
#population = breeding.diploid(population, genotyped.database = best_individuals[[1]][,1:3])


# Generation of non-genotyped individuals
population <- breeding.diploid(population, breeding.size = 100,
                               share.genotyped = 0.2)


### examples for naming of cohorts
# Although _M / _F are added Cohort_11 kind of still exists
get.database(population, cohorts = "Cohort_11")

pop1 <- breeding.diploid(population, breeding.size = 100,
                               share.genotyped = 1,
                         name.cohort = c("FancyName1", "ABC"))

# There is different ways of making sure your simulation does what its supposed to be doing
get.pedigree.visual(population, gen=10)

# Generation of phenotypic data
# In generation 11 & 12 all individuals are phenotyped
population <- breeding.diploid(population, phenotyping.gen = 11:12, heritability = 0.3)
# In generation 10 only males are phenotyped
population <- breeding.diploid(population, phenotyping.database = cbind(10,1), heritability = 0.3)

## Part two
# In case you are struggling with Task 7a load in an exemplary population
# save(file="C:/Users/pook001/OneDrive - Wageningen University & Research/GitHub/MoBPS/MoBPS_Workshop_WIAS/Task7/population.RData", list = c("population"))
# load("C:/Users/pook001/OneDrive - Wageningen University & Research/GitHub/MoBPS/MoBPS_Workshop_WIAS/Task7/population.RData")


# Looking at the prints provides helpful information to validate that the tool is doing what you want it do be doing!
population <- breeding.diploid(population, bve=TRUE, bve.gen=11, depth.pedigree = 3)
population <- breeding.diploid(population, bve=TRUE, bve.gen=11, relationship.matrix = "pedigree")

# BVE with various different arrays
population <- breeding.diploid(population, bve=TRUE, bve.gen=11, bve.array = 2)
population <- breeding.diploid(population, bve=TRUE, bve.gen=11, bve.array = 3)
population <- breeding.diploid(population, bve=TRUE, bve.gen=11, bve.array = 4)
population <- breeding.diploid(population, bve=TRUE, bve.gen=11, remove.effect.position=TRUE)
population <- breeding.diploid(population, bve=TRUE, bve.gen=11, remove.effect.position=TRUE, bve.array = 4)

# default BVE method is single-step GBLUP
population <- breeding.diploid(population, bve=TRUE, bve.gen=12)
# to exclude non-genotyped individuals:
population <- breeding.diploid(population, bve=TRUE, bve.gen=12, singlestep.active = FALSE)
# or perform a pedigree-based evaluation even when genotypes are partly available
population <- breeding.diploid(population, bve=TRUE, bve.gen=12, relationship.matrix = "pedigree")

# Prediction accuracy for phenotyped individuals is much better than for remaining individuals
population <- breeding.diploid(population, bve=TRUE, bve.gen=10)
analyze.bv(population, database = cbind(10,1))
analyze.bv(population, database = cbind(10,2))

# Visualization of estimated breeding values and real breeding values
# Estimated breeding values for females show much less variation (as they were not phenotyped)
bvs <- get.bv(population, gen=10)
bves <- get.bve(population, gen=10)

plot(bvs[1,], bves[1,], col=c(rep("blue",50), rep("red",50)),
     xlab="underlying true genomic values", ylab = "estimated breeding values",
     main = "Generation 10")
legend("topleft", c("Males", "Females"), lty=c(1,1), col=c("blue", "red"), lwd=2)

