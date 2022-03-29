set.seed(42)
library(MoBPS)

# All individuals of a cohort have to have the same sex
# In case a group of females and males is generated this will automatically be split
# into "Founders_M" and "Founders_F"

population <- creating.diploid(nsnp=5000, nindi = 100, sex.quota = 0.9,
                               name.cohort = "Founders")

population <- creating.trait(population, n.additive = 1000, mean.target = 100, var.target = 10)

# Generation of offspring
population <- breeding.diploid(population, breeding.size = c(45,45), name.cohort = "Offspring" )

# Generation of phenotypes for female offspring
population <- breeding.diploid(population, phenotyping.cohorts = "Offspring_F", heritability = 0.3)

# Calculate average offspring performance for our sires/dams
# As the two cohorts are both part of generation 1 we can also just use gen instead of cohorts
population <- breeding.diploid(population, offpheno.parents.gen = 1, offpheno.offspring.gen = 2)

# As female only have limited number of sires those average should be taken with more caution
get.pheno.off(population, gen=1)
get.pheno.off.count(population, gen=1)

### Alternative solution:
### You can do the same by extracting the pedigree and manually computing average offspring performance

ped <- get.pedigree(population, gen=2)
pheno <- get.pheno(population, gen=2)

parents_name <- get.pedigree(population, database = cbind(1,1))[,1]

off_pheno <- numeric(length(parents_name))
for(index in 1:10){
  off_pheno[index] <- mean( pheno[which(ped[,2]== parents_name[index])], na.rm=TRUE)
}

# Perform a breeding value estimation
# Make sure that only 45 individuals in your simulation were phenotyped!
population <- breeding.diploid(population, bve=TRUE, bve.gen = c(1,2))

# Computing a correlation requires individuals have some offspring

## Either use just the males:
pheno_off <- get.pheno.off(population, cohorts="Founders_M")
bvs <- get.bv(population, cohorts="Founders_M")
bves <- get.bve(population, cohorts="Founders_M")

# The prediction accuracy is the correlation between
# estimated and real breeding values
cor(bvs[1,], bves[1,])
cor(bvs[1,], pheno_off[1,])

## or adapt code to only calculate correlation on those individuals with available offspring phenotype

bv <- get.bv(population, gen=1)
bve <- get.bve(population, gen=1)
off <- get.pheno.off(population, gen=1)

cor(bv[1,], off[1,], use = "complete.obs")
cor(bv[1,!is.na(off)], bve[1,!is.na(off)], use = "complete.obs")


# A single simulation is not sufficient to judge which of the two estimation techniques is better
# Average offspring phenotypes should be very similar to a pedigree-based breeding value estimation
# Which typically should be inferior to the prediction of genomic values of the parents
# However if the main objective is maximizing offspring performance and traits
# are non-additive this might be a good alternative!
