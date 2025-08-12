library(MoBPS)
set.seed(42)

# Import genomic data from the first pool
population <- creating.diploid(vcf="C:/Users/pook001/OneDrive - Wageningen University & Research/MoBPS_Workshop_WIAS_2025/Task6/Pool1.vcf",
                               name.cohort = "Pool1", sex.quota = 0, founder.pool = 1)


# Import genomic data from the second pool
population <- creating.diploid(population = population, vcf="C:/Users/pook001/OneDrive - Wageningen University & Research/MoBPS_Workshop_WIAS_2025/Task6/Pool2.vcf",
                               name.cohort = "Pool2", sex.quota = 1, founder.pool = 2)

# Generation of the three different traits

population <- creating.trait(population, n.additive = c(10,1000,500),
                             n.dominant = c(0,0,500), dominant.only.positive = TRUE)

# You can use mean.target or var.target in creating.diploid or use bv.standardization
# bv.standardization can not only be applied in the initial generation but also
# at later time stages (e.g. if you first want to model some generations of drift / selection / build up)

population <- bv.standardization(population, mean.target = c(100,100,100), var.target = c(5,5,5))

# Simulation of matings between lines from the two gene pools

population <- breeding.diploid(population, breeding.size = c(100,0),
                               selection.m.cohorts = "Pool1", selection.f.cohorts = "Pool2",
                               name.cohort = "Crossbred")

summary(population)

# Principle components analysis
get.pca(population, gen=1:2, coloring=c(rep(1,10), rep(2,10), rep(3,100)))

# Extraction of the genomic values of all included lines
bv_pool1 <- get.bv(population, cohorts="Pool1")
bv_pool2 <- get.bv(population, cohorts="Pool2")
bv_crossbred <- get.bv(population, cohorts="Crossbred")

# For trait 3 we can we much higher genomic values of the crossbred lines
# Trait architecture with dominant effects
par(mfrow=c(3,3))
hist(bv_pool1[1,], xlim=c(94,115), main="Trait 1", ylab="Pool 1", xlab="genomic value")
hist(bv_pool1[2,], xlim=c(94,115), main="Trait 2", xlab="genomic value")
hist(bv_pool1[3,], xlim=c(94,115), main="Trait 3", xlab="genomic value")
hist(bv_pool2[1,], xlim=c(94,115), main="", ylab="Pool 2", xlab="genomic value")
hist(bv_pool2[2,], xlim=c(94,115), main="", xlab="genomic value")
hist(bv_pool2[3,], xlim=c(94,115), main="", xlab="genomic value")
hist(bv_crossbred[1,], xlim=c(94,115), main="", ylab="Crossbred", xlab="genomic value")
hist(bv_crossbred[2,], xlim=c(94,115), main="", xlab="genomic value")
hist(bv_crossbred[3,], xlim=c(94,115), main="", xlab="genomic value")
rowMeans(bv_pool1)
rowMeans(bv_pool2)
rowMeans(bv_crossbred)

# On average each line is used 10 times for reproduction
# There is however sampling effects
ped <- get.pedigree(population, gen=2)
table(ped[,2])
table(ped[,3])

#################################################################################################
####### END of the Task - everthing below is just showing of some additional MoBPS functionality
#################################################################################################


# Phenotyping and realized genetic variances
population = breeding.diploid(population, heritability = rep(0.3,3), phenotyping.cohorts = "Crossbred")
get.variance(population, cohorts = "Crossbred")


# Make a trait to only have variation in one pool
real.bv.add = get.qtl.effects(population)[[1]][[1]]
real.bv.add[,7] = 1

pop1 = creating.trait(population, replace.traits = TRUE,
                      real.bv.add = real.bv.add)


# Achieving overdominance 
real.bv.add = get.qtl.effects(population)[[1]][[3]]
real.bv.add[,4] = real.bv.add[,4] * 1.5

pop2 = creating.trait(population, replace.traits = TRUE,
                      real.bv.add = real.bv.add)

bv_pool1 <- get.bv(pop1, cohorts="Pool1")
bv_pool2 <- get.bv(pop1, cohorts="Pool2")
bv_crossbred <- get.bv(pop1, cohorts="Crossbred")

par(mfrow=c(3,1))
hist(bv_pool1[1,], xlim=c(94,115), main="Trait 1", ylab="Pool 1", xlab="genomic value")
hist(bv_pool2[1,], xlim=c(94,115), main="", ylab="Pool 2", xlab="genomic value")
hist(bv_crossbred[1,], xlim=c(94,115), main="", ylab="Crossbred", xlab="genomic value")


bv_pool1 <- get.bv(pop2, cohorts="Pool1")
bv_pool2 <- get.bv(pop2, cohorts="Pool2")
bv_crossbred <- get.bv(pop2, cohorts="Crossbred")

par(mfrow=c(3,1))
hist(bv_pool1[1,], xlim=c(94,125), main="Trait 1", ylab="Pool 1", xlab="genomic value")
hist(bv_pool2[1,], xlim=c(94,125), main="", ylab="Pool 2", xlab="genomic value")
hist(bv_crossbred[1,], xlim=c(94,125), main="", ylab="Crossbred", xlab="genomic value")


# Use different population means in different founder pools
pop2 = set.mean.pool(population, mean = c(110, 98), trait = 1,
                    gen = 1, reference = "pool")


get.qtl.effects(pop2)[[1]][[1]]

bv_pool1 <- get.bv(pop2, cohorts="Pool1")
bv_pool2 <- get.bv(pop2, cohorts="Pool2")
bv_crossbred <- get.bv(pop2, cohorts="Crossbred")

par(mfrow=c(3,1))
hist(bv_pool1[1,], xlim=c(94,115), main="Trait 1", ylab="Pool 1", xlab="genomic value")
hist(bv_pool2[1,], xlim=c(94,115), main="", ylab="Pool 2", xlab="genomic value")
hist(bv_crossbred[1,], xlim=c(94,115), main="", ylab="Crossbred", xlab="genomic value")



# To make sure all individuals are used at most 10 times ((max.offspring = 10))
population <- breeding.diploid(population, breeding.size = c(100,0),
                               selection.m.cohorts = "Pool1", selection.f.cohorts = "Pool2",
                               name.cohort = "Crossbred2", max.offspring = 10)

ped <- get.pedigree(population, cohorts = "Crossbred2")
table(ped[,2])
table(ped[,3])

# To make sure each mating combination is done once (( breeding.all.combination = TRUE))

population <- breeding.diploid(population, breeding.size = c(100,0),
                               selection.m.cohorts = "Pool1", selection.f.cohorts = "Pool2",
                               name.cohort = "Crossbred3", breeding.all.combination = TRUE)

ped <- get.pedigree(population, cohorts = "Crossbred3")

table(ped[,2])
table(ped[,3])

# Less strict threshold
population <- breeding.diploid(population, breeding.size = c(100,0),
                               selection.m.cohorts = "Pool1", selection.f.cohorts = "Pool2",
                               name.cohort = "Crossbred4", max.offspring = 12)

ped <- get.pedigree(population, cohorts = "Crossbred4")
table(ped[,2])
table(ped[,3])

# To manually select mating pairs fixed.breeding

# Generate 3 new individuals:
# First: mate individual from generation 1, sex 1, nr 1 with individual from generation 1, sex 2, nr 1
# Second: mate individual from generation 1, sex 1, nr 1 with individual from generation 1, sex 2, nr 2
# Third: mate individual from generation 1, sex 1, nr 4 with individual from generation 1, sex 2, nr 6
fixed <- matrix(c(1,1,1,1,2,1,
                  1,1,1,1,2,2,
                  1,1,4,1,2,6), byrow=TRUE, ncol=6)

population <- breeding.diploid(population, fixed.breeding = fixed)

fixed <- matrix(c(1,11,
                  1,12,
                  4,16), byrow=TRUE, ncol=2)
population <- breeding.diploid(population, fixed.breeding.id = fixed)

ped <- get.pedigree(population, gen=c(6,7))

table(ped[,2])
table(ped[,3])

