library(MoBPS)
set.seed(42)

# Import genomic data from the first pool
population <- creating.diploid(vcf="C:/Users/pook/Desktop/MoBPS_workshop_gft/Task6/Pool1.vcf",
                               name.cohort = "Pool1", sex.quota = 0)


# Import genomic data from the second pool
population <- creating.diploid(population = population, vcf="C:/Users/pook/Desktop/MoBPS_workshop_gft/Task6/Pool2.vcf",
                               name.cohort = "Pool2", sex.quota = 1)

# Generation of the three different traits

population <- creating.trait(population, n.additive = c(10,1000,500), n.equal.dominant = c(0,0,500))

# You can use mean.target or var.target in creating.diploid or use bv.standardization
# bv.standardization can not only be applied in the initial generation but also
# at later time stages (e.g. if you first want to model some generations of drift / selection / build up)

population <- bv.standardization(population, mean.target = c(100,100,100), var.target = c(5,5,5))

# Simulation of matings between lines from the two gene pools

population <- breeding.diploid(population, breeding.size = c(100,0),
                               selection.m.cohorts = "Pool1", selection.f.cohorts = "Pool2",
                               name.cohort = "Hybrid")

# Principle components analysis
get.pca(population, gen=1:2, coloring=c(rep(1,10), rep(2,10), rep(3,100)))

# Extraction of the genomic values of all included lines
bv_pool1 <- get.bv(population, cohorts="Pool1")
bv_pool2 <- get.bv(population, cohorts="Pool2")
bv_hybrid <- get.bv(population, cohorts="Hybrid")

# For trait 3 we can we much higher genomic values of the hybrid lines
# Trait architecture with dominant effects
par(mfrow=c(3,3))
hist(bv_pool1[1,], xlim=c(94,115), main="Trait 1", ylab="Pool 1", xlab="genomic value")
hist(bv_pool1[2,], xlim=c(94,115), main="Trait 2", xlab="genomic value")
hist(bv_pool1[3,], xlim=c(94,115), main="Trait 3", xlab="genomic value")
hist(bv_pool2[1,], xlim=c(94,115), main="", ylab="Pool 2", xlab="genomic value")
hist(bv_pool2[2,], xlim=c(94,115), main="", xlab="genomic value")
hist(bv_pool2[3,], xlim=c(94,115), main="", xlab="genomic value")
hist(bv_hybrid[1,], xlim=c(94,115), main="", ylab="Hybrids", xlab="genomic value")
hist(bv_hybrid[2,], xlim=c(94,115), main="", xlab="genomic value")
hist(bv_hybrid[3,], xlim=c(94,115), main="", xlab="genomic value")
rowMeans(bv_pool1)
rowMeans(bv_pool2)
rowMeans(bv_hybrid)

# On average each line is used 10 times for reproduction
# There is however sampling effects
ped <- get.pedigree(population, gen=2)
table(ped[,2])
table(ped[,3])



#### Use GGplot as an alternative to evaluate results

library(ggplot2)
# loop over cohorts and extract true breeding values
cohorts <- get.cohorts(population)
res <- list() # initialize empty list
for(i in 1:length(cohorts)){
  tmp <- get.bv(population = population,cohorts = cohorts[i])
  # reformat results in long format in data frame
  res[[i]] <- data.frame(
    cohort = cohorts[i],
    trait = rep(rownames(tmp), ncol(tmp)),
    ind = rep(colnames(tmp), each = nrow(tmp)),
    bv = as.numeric(tmp)
  )
}
res <- do.call(rbind,res) # rbind list elements
res$cohort <- factor(res$cohort,levels = c("Pool1","Hybrid","Pool2")) # sort factor for plotting
# create plot
ggplot(res,
       aes(x = trait,
           y = bv,
           color = cohort))+
  geom_boxplot()+
  scale_color_manual(values = c('blue','orange','red'))+
  theme_minimal()



# To make sure all individuals are used at most 10 times ((max.offspring = 10))
population <- breeding.diploid(population, breeding.size = c(100,0),
                               selection.m.cohorts = "Pool1", selection.f.cohorts = "Pool2",
                               name.cohort = "Hybrid2", max.offspring = 10)

ped <- get.pedigree(population, cohorts = "Hybrid2")
table(ped[,2])
table(ped[,3])

# To make sure each mating combination is done once (( breeding.all.combination = TRUE))

population <- breeding.diploid(population, breeding.size = c(100,0),
                               selection.m.cohorts = "Pool1", selection.f.cohorts = "Pool2",
                               name.cohort = "Hybrid3", breeding.all.combination = TRUE)

ped <- get.pedigree(population, cohorts = "Hybrid3")

table(ped[,2])
table(ped[,3])


# To make sure all individuals are used at most 10 times ((max.offspring = 10))
population <- breeding.diploid(population, breeding.size = c(100,0),
                               selection.m.cohorts = "Pool1", selection.f.cohorts = "Pool2",
                               name.cohort = "Hybrid4", max.offspring = 12)

ped <- get.pedigree(population, cohorts = "Hybrid4")
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

ped <- get.pedigree(population, gen=6)

table(ped[,2])
table(ped[,3])

