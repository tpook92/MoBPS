devtools::install_github("tpook92/MoBPS", subdir="pkg")

real.bv.add <- cbind(c(120,42,17),c(1,5,22),c(-1,0,0.1),c(0,0,0.1),c(1,2,0))
colnames(real.bv.add) <- c("SNP", "chromosome", "Effect 0", "effect 1", "effect 2")


real.bv.mult <- cbind(c(144,6,5),c(1,3,17),c(145,188,1), c(1,5,10), c(5,0,0), c(0,0,0), c(0,0,0), c(0,0,0), c(0,0,0), c(0,0,0), c(0,0,0), c(0,0,0), c(0,0,0))
set.seed(1)
real.bv.mult[2:3,5:13] <- round(rnorm(18,1,1), digits=2)
colnames(real.bv.mult) <- c("First SNP", "First chromosome", "Second SNP", "Second chromosome", "effect 00", "effect 01", "effect 02", "effect 10", "effect 11", "effect 12", "effect 20", "effect 21", "effect 22")

location <- list(matrix(c(11,12,16,1,1,4), ncol=2), matrix(c(14,77,15,2,6,9), ncol=2))
effects <- list(rnorm(27,1,1), rnorm(27,1,1))
real.bv.dice<- list()
real.bv.dice$location <- location
real.bv.dice$effects <- effects

colnames(real.bv.dice$location[[1]]) <- c("SNP", "chromosome")
colnames(real.bv.dice$location[[2]]) <- c("SNP", "chromosome")



population = creating.diploid(nsnp = 10, chr.nr = 5, nindi = 50)
population = breeding.diploid(population, breeding.size = 50,
                              breeding.sex = 0.6,
                              selection.size = c(10,10),
                              name.cohort = "Offspring")

database <- cbind(c(1,5), c(2,1))
colnames(database) <- c("Generation", "sex")


library(MoBPS)
a <- NULL
b <- NULL
c <- NULL
d <- NULL
for(index in 1:10000){
  print(index)
pop <- creating.diploid(nsnp=10000, nindi=100,
                        chr.nr=5, chromosome.length=2,
                        n.additive=50, n.dominant=10,
                        name.cohort="Founder",
                        var.target = 1,
                        share.genotyped = 1)
pop <- breeding.diploid(pop, heritability=0.5,
                        phenotyping="all")
pop <- breeding.diploid(pop, bve=TRUE)
pop1 <- breeding.diploid(pop, breeding.size=100,
                        selection.size=c(20,20),
                        selection.m="function",
                        selection.m.cohorts="Founder_M",
                        selection.f.cohorts="Founder_F",
                        name.cohort="Offspring")
pop2 <- breeding.diploid(pop, breeding.size=100,
                        selection.size=c(5,20),
                        selection.m="function",
                        selection.m.cohorts="Founder_M",
                        selection.f.cohorts="Founder_F",
                        name.cohort="Offspring")

a <- c(a,mean(get.bv(pop1, gen=2) - mean(get.bv(pop, gen=1))))
b <- c(b,mean(get.bv(pop2, gen=2) - mean(get.bv(pop, gen=1))))

c1 <- kinship.emp(population=pop1, gen=2)
d1 <- kinship.emp(population=pop2, gen=2)

c <- c(c, mean(c1-diag(diag(c1))) * 100/99)
d <- c(d, mean(d1-diag(diag(d1))) * 100/99)
# Expected gain: i * h * sigma_a = sqrt(0.5) * 0.966 * 10 = 6.83
# Expected gain: i * h * sigma_a = sqrt(0.5) * 1.755 * 10 = 12.4

}

save(file="C:/Users/pook/Desktop/Fig1_data.RData", list=c("a", "b", "c","d"))

{
png("C:/Users/pook/Desktop/Fig1.png", width=2250, height= 1100, res=300)
par(mfrow=c(2,2))
par(mar=c(4.5,1.5,0.5,2))
hist(a, xlim=c(0.2,1.65), main="", axes = FALSE, xlab="", breaks = 50)
abline(v=mean(a), col="red", lwd=2)
axis(side=4)
axis(side=1, at=seq(0.4,1.6, by=.3))
title(ylab="Scenario 1", line=0, cex.lab=1.5)
hist(c, xlim=c(0.004,0.019), main="", axes = FALSE, xlab="", breaks = 25)
abline(v=mean(c), col="red", lwd=2)
axis(side=4, ylab="Scenario 1")
axis(side=1, at=seq(0.005,0.018, by=.003))
hist(b, xlim=c(0.2,1.65), axes=FALSE, main="", xlab="", breaks = 50)
abline(v=mean(b), col="red", lwd=2)
axis(side=4)
axis(side=1, at=seq(0.4,1.6, by=.3))
title(ylab="Scenario 2", line=0, cex.lab=1.5)
title(xlab="Genomic gain", cex.lab=1.5)
hist(d, xlim=c(0.004,0.019), axes=FALSE, main="", xlab="", breaks = 40)
abline(v=mean(d), col="red", lwd=2)
axis(side=4)
axis(side=1,at=seq(0.005,0.018, by=.003))
title(xlab="Gain in avg. kinship", cex.lab=1.5)
dev.off()
}


x = seq(-1,1, by = 0.001)

rej1 = 1 - (abs(x)/0.2) ^ 1
rej2 = 1 - (abs(x)/0.2) ^ 2

rej1[rej1<0] = 0
rej2[rej2<0] = 0
par(mfrow = c(1,1))
plot(x,rej1, type = "l", lwd =2, xlab = "distance to previously sampled cross-over in Morgan",
     ylab = "rejection probability")
lines(x,rej2, lwd =2, col = 2)

rej1 = 1 - (abs(x)/0.5) ^ 1
rej2 = 1 - (abs(x)/0.5) ^ 2
rej3 =
rej1[rej1<0] = 0
rej2[rej2<0] = 0
par(mfrow = c(1,1))
lines(x,rej1, type = "l", lwd =2, col = 3)
lines(x, rej2, lwd = 2, col= 4)

legend("topleft",
       lty = 1,
       lwd =3,
       c("recombination.distance.penalty = 0.2",
         "recombination.distance.penalty2 = 0.2",
         "recombination.distance.penalty = 0.5",
         "recombination.distance.penalty2 = 0.5"),
       col = 1:4)


{
  eps("C:/Users/pook/Desktop/Fig1.png", width=2250, height= 1100, res=300)
  par(mfrow=c(2,2))
  par(mar=c(4.5,1.5,0.5,2))
  hist(a, xlim=c(0.2,1.65), main="", axes = FALSE, xlab="", breaks = 50)
  abline(v=mean(a), col="red", lwd=2)
  axis(side=4)
  axis(side=1, at=seq(0.4,1.6, by=.3))
  title(ylab="Scenario 1", line=0, cex.lab=1.5)
  hist(c, xlim=c(0.004,0.019), main="", axes = FALSE, xlab="", breaks = 25)
  abline(v=mean(c), col="red", lwd=2)
  axis(side=4, ylab="Scenario 1")
  axis(side=1, at=seq(0.005,0.018, by=.003))
  hist(b, xlim=c(0.2,1.65), axes=FALSE, main="", xlab="", breaks = 50)
  abline(v=mean(b), col="red", lwd=2)
  axis(side=4)
  axis(side=1, at=seq(0.4,1.6, by=.3))
  title(ylab="Scenario 2", line=0, cex.lab=1.5)
  title(xlab="Genomic gain", cex.lab=1.5)
  hist(d, xlim=c(0.004,0.019), axes=FALSE, main="", xlab="", breaks = 40)
  abline(v=mean(d), col="red", lwd=2)
  axis(side=4)
  axis(side=1,at=seq(0.005,0.018, by=.003))
  title(xlab="Gain in avg. kinship", cex.lab=1.5)
  dev.off()
}
{
  png("C:/Users/pook/Desktop/Fig1b.png", width=2250, height= 1100, res=300)
  par(mfrow=c(2,2))
  par(mar=c(4.5,1.5,0.5,2))
  plot(density(a), xlim=c(0.2,1.65), main="", axes = FALSE, xlab="", breaks = 50)
  abline(v=mean(a), col="red", lwd=2)
  axis(side=4)
  axis(side=1, at=seq(0.4,1.6, by=.3))
  title(ylab="Scenario 1", line=0, cex.lab=1.5)
  plot(density(c), xlim=c(0.004,0.019), main="", axes = FALSE, xlab="", breaks = 25)
  abline(v=mean(c), col="red", lwd=2)
  axis(side=4, ylab="Scenario 1")
  axis(side=1, at=seq(0.005,0.018, by=.003))
  plot(density(b), xlim=c(0.2,1.65), axes=FALSE, main="", xlab="", breaks = 50)
  abline(v=mean(b), col="red", lwd=2)
  axis(side=4)
  axis(side=1, at=seq(0.4,1.6, by=.3))
  title(ylab="Scenario 2", line=0, cex.lab=1.5)
  title(xlab="Genomic gain", cex.lab=1.5)
  plot(density(d), xlim=c(0.004,0.019), axes=FALSE, main="", xlab="", breaks = 40)
  abline(v=mean(d), col="red", lwd=2)
  axis(side=4)
  axis(side=1,at=seq(0.005,0.018, by=.003))
  title(xlab="Gain in avg. kinship", cex.lab=1.5)
  dev.off()
}
summary(pop1)


a <- NULL
b <- NULL
for(index in 1:1000){
  pop <- creating.diploid(nsnp=10000, nindi=100,
                          chr.nr=5, chromosome.length=2,
                          n.additive=50, n.dominant=10,
                          name.cohort="Founder",
                          bv.standard = TRUE,
                          var.target = 100)
  pop <- breeding.diploid(pop, heritability=0.5,
                          phenotyping="all")
  pop <- breeding.diploid(pop, bve=TRUE)
  pop1 <- breeding.diploid(pop, breeding.size=100,
                           selection.criteria = "pheno",
                           selection.size=c(20,20),
                           selection.m="function",
                           selection.m.cohorts="Founder_M",
                           selection.f.cohorts="Founder_F",
                           name.cohort="Offspring")
  pop2 <- breeding.diploid(pop, breeding.size=100,
                           selection.criteria = "pheno",
                           selection.size=c(5,5),
                           selection.m="function",
                           selection.m.cohorts="Founder_M",
                           selection.f.cohorts="Founder_F",
                           name.cohort="Offspring")


  a <- c(a,mean(get.bv(pop1, gen=2) - mean(get.bv(pop, gen=1))))
  b <- c(b,mean(get.bv(pop2, gen=2) - mean(get.bv(pop, gen=1))))

}

# Expected gain: i * h * sigma_a = sqrt(0.5) * 0.966 * 10 = 6.83
# Expected gain: i * h * sigma_a = sqrt(0.5) * 1.755 * 10 = 12.4
summary(pop1)
analyze.bv(pop, gen=1)

str(pop$breeding[[1]][[1]][[1]])

ex_json <- jsonlite::read_json(path="C:/Users/pook/Downloads/Ex_json (4).json")
save(file="C:/Users/pook/Desktop/R-Stuff/MoBPS/data/ex_json.RData", list=c("ex_json"), version = 2)
genos <- get.geno(population,gen=3)
genos[1:5,1:10]

haplos <- get.haplo(population,gen=3)
haplos[1:5,1:10]

bv <- get.bv(population,gen=3)
bve <- get.bve(population,gen=3)
pheno <- get.pheno(population,gen=3)
bve[1:5]
bv[1:5]
pheno[1:5]

recombi <- get.recombi(population, gen=3)

ped <- get.pedigree(population, gen=5)
ped[1:5,]
ped <- get.pedigree(population, gen=5, raw=TRUE)
ped[1:5,]
ped <- get.pedigree2(population, gen=5)
ped[1:5,]
ped <- get.pedigree3(population, gen=5)
ped[1:5,]

population$breeding[[1]][[1]][[1]]


dataset <- matrix(rbinom(100,1,0.5), nrow=10)
rownames(dataset) <- paste0("SNP", 1:10)
colnames(dataset) <- paste0("Indi", sort(rep(1:5,2)), "Haplo", rep(1:2,5))

dataset


bves <- cbind(c("M1_1", "M2_1", "M3_1", "M4_1", "M5_1"), c(101.5, 102, 99.7, 98.2,103.8), c(104.2,98.9, 98.4, 101.2, 101.1))
colnames(bves) <- c("Individual Name", "Trait 1", "Trait 2")

bves
population <- insert.bve(population, bves)

population <- json.simulation(total = ex_json)
analyze.population(population,5,2, gen=1:8)

population <- json.simulation(total = ex_json)
analyze.bv(population, cohorts="New_bulls_4")
analyze.bv(population, cohorts="New_cows_4")

data("ex_pop")
get.pca(ex_pop, gen=1:2)

# Simulate Phenotypes for generation 4 with heritability 0.4
population <- breeding.diploid(population, heritability = 0.4,
                               sigma.e.gen = 4,
                               phenotyping.gen=4)
# Export genotypes and phenotypes for generation 4
genos <- get.geno(population, gen=4)
phenos <- get.pheno(population, gen=4)

# Here you perform your own method to assign breeding values to each individual
bve <- runif(ncol(genos)) # This is probably not the best technique for this =)

# Import breeding values estimated for generation 4
bves <- cbind(colnames(genos), bve)
population <- insert.bve(population, bves=bves)


# get.geno, get.haplo, get.bv, get.bve, get.pheno, get.recombi, get.pedigree
# get.class, get.cohorts, get.creating.type, get.database, get.individual.loc, get.time.point, get.vcf

get.cohorts(population)[1:5]
get.cohorts(population, extended=TRUE)[1:5,]
get.class(population, gen=1:2)
get.time.point(population , gen=1:3)
get.creating.type(population , gen=1:3)
get.vcf(population, path="C:/Users/pook/Desktop/Like_to_save_on_Desktop.vcf", gen=3)
get.pedmap(population, path="C:/Users/pook/Desktop/Like_to_save_on_Desktop", gen=2)

# Mating design MAGIC
{
{
  library(MoBPS)

# Generation of 20 fully-homozygous founders lines
# All plants are stored as male individuals (sex=0)
population <- creating.diploid(nindi = 20, sex.quota = 0, template.chip = "maize",
                               dataset = "homorandom")

# Simulate matings between all founders.
# Each plan is involved in exactly 19 matings.
population <- breeding.diploid(population, breeding.size = c(190,0),
                               breeding.all.combination = TRUE,
                               selection.size = c(20,0), max.offspring = 19)

# Simulate matings between plants of the last generation.
# Each plant is involved in exactly 2 matings.

population <- breeding.diploid(population, breeding.size = c(190,0),
                               selection.size = c(190,0), same.sex.activ = TRUE,
                               same.sex.sex = 0, max.offspring = 2)
population <- breeding.diploid(population, breeding.size = c(190,0),
                               selection.size = c(190,0), same.sex.activ = TRUE,
                               same.sex.sex = 0, max.offspring = 2)
population <- breeding.diploid(population, breeding.size = c(190,0),
                               selection.size = c(190,0), same.sex.activ = TRUE,
                               same.sex.sex = 0, max.offspring = 2)
}

# Same design in cohort modus
{
library(MoBPS)

# Generation of 20 fully-homozygous founders lines
# All plants are stored as male individuals (sex=0)
population <- creating.diploid(nindi = 20, sex.quota = 0, template.chip = "maize",
                               dataset = "homorandom", name.cohort = "F0")

# Simulate matings between all founders.
# Each plan is involved in exactly 19 matings.
population <- breeding.diploid(population, breeding.size = c(190,0),
                               breeding.all.combination = TRUE,
                               selection.size = c(20,0),
                               selection.m.cohort = "F0", name.cohort = "F1")

# Simulate matings between plants of the last generation.
# Each plant is involved in exactly 2 matings.

population <- breeding.diploid(population, breeding.size = c(190,0),
                               selection.size = c(190,0), same.sex.activ = TRUE,
                               same.sex.sex = 0, max.offspring = c(2,0),
                               selection.m.cohort = "F1", name.cohort = "F2")
population <- breeding.diploid(population, breeding.size = c(190,0),
                               selection.size = c(190,0), same.sex.activ = TRUE,
                               same.sex.sex = 0, max.offspring = c(2,0),
                               selection.m.cohort = "F2", name.cohort = "F3")
population <- breeding.diploid(population, breeding.size = c(190,0),
                               selection.size = c(190,0), same.sex.activ = TRUE,
                               same.sex.sex = 0, max.offspring = c(2,0),
                               selection.m.cohort = "F3", name.cohort = "F4")
}
}


# Introgression
{
set.seed(1)
library(MoBPS)

# Generate an input SNP-dataset
# 10 White-Layer (0) (20 haplotypes, 5'000 SNPs)
# 10 Wild population (1) (20 haplotypes, 5'000 SNPs)
dataset1 <- matrix(0, nrow = 5000, ncol = 20)
dataset2 <- matrix(1, nrow = 5000, ncol = 20)

# Generation of a trait
# Colums code: SNP, chromosome, effect 00, effect 01, effect 11
# Blue Eggshell QTL is positioned on SNP 2000, chromosome 1
major_qtl <- c(2000, 1, 0, 10000, 20000)
# In all other positions the white layer genome is assumed to be favorable
# All marker effects combiened are smaller than the blue eggshell QTL
rest <- cbind(1:5000, 1, 1, 0.5, 0)
trait <- rbind(major_qtl, rest)

# Generation of the base-population
# First 10 individuals are female (sex=2)
# Next 10 individuals are male (sex=1)
population <- creating.diploid(dataset = cbind(dataset1, dataset2),
                               real.bv.add = trait, name.cohort = "Founders",
                               sex.s = c(rep(2,10), rep(1,10)))

# Simulate random mating:
population <- breeding.diploid(population, breeding.size = c(100,100),
                               selection.size = c(10,10),
                               selection.m.cohorts = "Founders_M",
                               selection.f.cohorts = "Founders_F",
                               name.cohort = "F1")

# Simuation of matings with selection:
# Top 50 cocks are mated to the 10 founder hens
# Selection of the cocks based on their genomic value ("bv")
# Target: Increase share of white layer while preserving blue egg shell QTL

population <- breeding.diploid(population, breeding.size = c(100,100),
                               selection.size = c(50,10),
                               selection.m.cohorts = "F1_M",
                               selection.f.cohorts = "Founders_F",
                               name.cohort = "BC1", selection.m = "function",
                               selection.criteria = "bv")
population <- breeding.diploid(population, breeding.size = c(100,100),
                               selection.size = c(50,10),
                               selection.m.cohorts = "BC1_M",
                               selection.f.cohorts = "Founders_F",
                               name.cohort = "BC2", selection.m = "function",
                               selection.criteria = "bv")
population <- breeding.diploid(population, breeding.size = c(100,100),
                               selection.size = c(50,10),
                               selection.m.cohorts = "BC2_M",
                               selection.f.cohorts = "Founders_F",
                               name.cohort = "BC3", selection.m = "function",
                               selection.criteria = "bv")

# Mating of cocks and hens that are heterozygous in blue egg shell QTL
# 25% of resulting offspring should be homozygous in blue egg shell QTL

population <- breeding.diploid(population, breeding.size = c(100,100),
                               selection.size = c(50,50),
                               selection.m.cohorts = "BC3_M",
                               selection.f.cohorts = "BC3_F",
                               name.cohort = "IC",
                               selection.criteria = "bv")

# Check genomic share of wild race in the final generation
genoIC <- get.geno(population, cohorts = "IC_F")
plot(rowSums(genoIC)/200, xlab = "genome", ylab = "frequency of wild allele", type = "l")
abline(v = 2000, lwd = 2, col = "red")

png("C:/Users/pook/Desktop/wild_allele_freq.png", width=2250, height= 960, res=300)
par(mar=c(4.1,4.1,1.6,0.6))
plot(rowSums(genoIC)/200, xlab="genome", ylab="frequency of wild allele", type="l")
abline(v=2000, lwd=2, col="red")
dev.off()
}


library(miraculix)
a <- computeSNPS(population, rep(1,10), rep(1,10), 1:10, output_compressed=TRUE)

b <- copyGeno(a)
c <- zeroNthGeno(b, 1:10*10)
e <- decodeGeno(a)
f <- decodeGeno(c, N=1:10)
prod(decodeGeno(a)[,11]==decodeGeno(c)[,11])

dim(decodeGeno(a))
dim(decodeGeno(c))


args <- commandArgs(TRUE)

{
set.seed(1)
library(MoBPS)

# Generation of a base population:
# 1'000 Founder individuals
# 5'000 SNPs
# 100 additive single marker QTL
population <- creating.diploid(nindi = 1000, nsnp = 5000,
                               n.additive = 100, name.cohort = "Founders")

# Simulation of a random mating generation
# 100 bulls (sex=1), 1'000 cows (sex=2) are generated
population <- breeding.diploid(population, breeding.size = c(100,1000),
                               share.genotyped = 1,
                               selection.size = c(500,500),
                               selection.m.cohorts = "Founders_M",
                               selection.f.cohorts = "Founders_F",
                               name.cohort = "Random")

# Generate 200 offspring of both from the top 5 bulls / 200 cows
# Heritability of the trait is set to 0.5
# only phenotypes previously unobserved cows are generated
population <- breeding.diploid(population, breeding.size = 200,
                               selection.size = c(5,200), bve = TRUE,
                               heritability = 0.5,
                               phenotyping = "non_obs_f",
                               selection.criteria = "bve",
                               name.cohort = "Top",
                               selection.m.cohorts = "Random_M",
                               selection.f.cohorts = "Random_F")

# Generate additional cows using all cows of the previous generation
# Cows are added to the same generation as the previous simulation
population <- breeding.diploid(population, breeding.size = c(0,900),
                               selection.size = c(5,1000),
                               selection.criteria = "bve",
                               name.cohort = "Sec_F",
                               share.genotyped = 1,
                               selection.m.cohorts = "Random_M",
                               selection.f.cohorts = "Random_F",
                               add.gen = 3)

# Same cycle as before with additional genome editing
# Edits are chosen based on highest effects in rrBLUP
population <- breeding.diploid(population, breeding.size = c(100,100),
                               selection.size = c(5,200), bve = TRUE,
                               phenotyping = "non_obs_f",
                               selection.criteria = "bve",
                               name.cohort = "Top_Edit",
                               selection.m.cohorts = "Top_M",
                               selection.f.cohorts = c("Top_F","Sec_F"),
                               nr.edits = 20, estimate.u = TRUE)

population <- breeding.diploid(population, breeding.size = c(0,900),
                               selection.size = c(5,1000),
                               selection.criteria = "bve",
                               name.cohort = "Sec_Edit",
                               selection.m.cohorts = "Top_M",
                               selection.f.cohorts = c("Top_F","Sec_F"),
                               add.gen = 4)



bv.development.box(population, cohorts = c("Founders_F", "Random_F", "Sec_F",
                               "Top_F", "Sec_Edit", "Top_Edit_F"), display = "bv")

png("C:/Users/pook/Desktop/cattle_example.png", width=3250, height= 1660, res=300)
par(mar=c(7,3,1,1))
bv.development.box(population, cohorts = c("Founders_F", "Random_F", "Sec_F",
                                           "Top_F", "Sec_Edit", "Top_Edit_F"), display = "bv")
dev.off()
}




# Seletion Hardsweep
{
library(MoBPS)
set.seed(2)

# Generate a starting population with 5000 SNPs and 200 individuals
# and a single chromosome of length 2 Morgan.
population <- creating.diploid(nsnp = 5000, nindi = 200, chromosome.length = 2)

# LD build up via 100 generations of random mating
# Each generation contains 200 individuals
for(index in 1:100){
  population <- breeding.diploid(population, breeding.size = 200,
                                 selection.size = c(100,100))
}

# Derive allele frequency and check LD for the last generation:
genotype.check <- get.geno(population, gen = length(population$breeding))
p_i <- rowMeans(genotype.check)/2
ld.decay(population, genotype.dataset = genotype.check, step = 10, max = 500, plot = TRUE)

# Simulate a favorable mutation in a previously fixed marker
fixated_markers <- which(p_i==0) # Which markers are fixated
qtl_posi <- sample(fixated_markers, 1) # Selected a fixated marker in A
trait <- cbind(qtl_posi, 1, 0, 1, 2) # SNP, Chromosome, Effect AA, Effect AB, Effect BB
population <- creating.trait(population, real.bv.add = trait)

# Generate a mutation in the first male individual
population <- mutation.intro(population, 101, 1, 1, qtl_posi)

# Simulate generations with selection pressure
# Individuals with the favorable SNP are picked 5 times as often
for(index in 1:25){
  population <- breeding.diploid(population, breeding.size = 200,
                                 selection.size = c(100,100),
                                 best.selection.ratio.m = 5,
                                 best.selection.ratio.f = 5)
}

analyze.population(population, gen = 98:115, chromosome = 1, snp = qtl_posi)

  png("C:/Users/pook001/OneDrive - Wageningen University & Research/ld_decay_sweep.png", width=2250, height= 960, res=300)
  par(mar=c(4.1,4.1,0.6,0.6))
  ld.decay(population, genotype.dataset = genotype.check, step=5, max=500, plot = TRUE)
  dev.off()

  png("C:/Users/pook001/OneDrive - Wageningen University & Research/allele_freq_sweep.png", width=2250, height= 960, res=300)
  par(mar=c(4.1,4.1,0.6,0.6))
  analyze.population(population, gen=98:115, chromosome = 1, snp=qtl_posi )
  dev.off()

}

{

  X11()
     par(mfrow=c(1,5))
     plot(get.distance(population, gen1 = 100, gen2=110, per.marker = TRUE))
     abline(v=qtl_posi, col="red")
     plot(get.distance(population, gen1 = 100, gen2=110, per.marker = TRUE, type="reynold"))
     abline(v=qtl_posi, col="red")
     plot(get.distance(population, gen1 = 100, gen2=110, per.marker = TRUE, type="cavalli"))
     abline(v=qtl_posi, col="red")
     plot(get.distance(population, gen1 = 100, gen2=110, per.marker = TRUE, type="nei_distance"))
     abline(v=qtl_posi, col="red")
   plot(get.distance(population, gen1 = 100, gen2=110, per.marker = TRUE, type="nei_minimum"))
     abline(v=qtl_posi, col="red")


}


# Cock-Rotation

set.seed(1)
{

# Generate initial boxes with 5 hens (sex=2) and 1 cock (sex=1) each
population <- NULL
for(index in 1:7){
  population <- creating.diploid(population = population, nindi = 6,
                                 nsnp = 5000, sex.s = c(1, 2, 2, 2, 2, 2),
                                 name.cohort = paste0("Box_", index, "gen_0"))
}

# Simulate 25 generations of matings.
# Hens are rotated by one box per generation.
# selection.m.cohorts is the cohort used as sires
# selection.f.cohorts is the cohort used as dams
for(gen in 1:25){
    for(index in 1:7){
    population <- breeding.diploid(population, breeding.size = c(1,5),
                                   selection.size = c(1,5),
              selection.m.cohorts = paste0("Box_",
                if(index==1){7} else {index-1},"gen_", gen-1,"_M"),
              selection.f.cohorts = paste0("Box_", index,"gen_", gen-1,"_F"),
              name.cohort = paste0("Box_", index, "gen_", gen),
              add.gen=gen+1
    )
  }
}


# Generate a population of same size without cock rotation
pop1 <- creating.diploid(nindi = 42, nsnp = 5000,
                         sex.s = c(rep(1,7), rep(2,35)))

# Simulate 25 generations of random mating
for(gen in 1:25){
  pop1 <- breeding.diploid(pop1, breeding.size = c(7,35),
                           selection.size = c(7,35))
}

kin <- kinship.development(population, gen = 1:26, ibd.obs = 1000, hbd.obs = 100)
kin1 <- kinship.development(pop1, gen = 1:26, ibd.obs = 1000, hbd.obs = 100)

png(file="C:/Users/pook/Desktop/inbreed_rota.png", width=2450, height= 1060, res=300)
par(mar=c(4.1,4.1,1.6,0.6))
plot(0:25,kin[,1], xlab="generation", ylab="kinship", type="l", lwd=2, col="red", ylim=c(0, max(kin[,1], kin1[,1])))
lines(0:25,kin1[,1], col="blue", lwd=2)
legend("topleft", c("Cock-rotation", "Random mating"), lty=c(1,1), lwd=c(2,2), col=c("red", "blue"))
dev.off()

#  population$info$size
#a <- kinship.emp(population = population, gen = 26, sym = TRUE)
#b <- kinship.emp(population = pop1, gen = 26, sym = TRUE)


}
{
  ## Genotyping/Phenotyping of Subcohorts

# Generate a starting population with 5000 SNPs and 500 male individuals
# 3 Traits with 500 purely additive QTL each
# No genotypes / phenotypes are generate
population <- creating.diploid(nsnp=5000, nindi = 500,
                               n.additive = c(500,500,500),
                               share.genotyped = 0,
                               sex.quota = 0,
                               name.cohort="Founder")

# Generate a copy of those individuals that are supposed to be
# genotyped and/or phenotyped
population <- breeding.diploid(population, selection.size = c(250,0),
                               copy.individual.m = TRUE,
                               selection.m.cohorts = "Founder",
                               name.cohort = "PartlyPhenotyped")

population <- breeding.diploid(population, selection.size = c(100,0),
                               copy.individual.m = TRUE,
                               selection.m.cohorts = "PartlyPhenotyped",
                               name.cohort = "Phenotyped+Genotyped")

# Generate phenotypes for trait 1 & 2. Each trait is assumed to have a heritabilty of 0.5
population <- breeding.diploid(population, phenotyping.cohorts = "PartlyPhenotyped",
                               n.observation = c(1,1,0), heritability = c(0.5,0.5,0.5),
                               sigma.e.cohorts = "Founder")

# Generate genotypes / phenotypes for all traits
population <- breeding.diploid(population, phenotyping.cohorts = "Phenotyped+Genotyped",
                               genotyped.cohorts = "Phenotyped+Genotyped")

# MoBPS default: single-step GBLUP
# All 500 individuals are used with 250/250/100 providing phenotype information
population <- breeding.diploid(population, bve=TRUE,
                               bve.cohorts = c("Founder", "PartlyPhenotyped",
                                               "Phenotyped+Genotyped"))

# BVE with pedigree BLUP
population <- breeding.diploid(population, bve=TRUE,
                               bve.cohorts = c("Founder", "PartlyPhenotyped",
                                               "Phenotyped+Genotyped"),
                               relationship.matrix = "kinship")

# Not include non-genotyped individuals in the evaluation
population <- breeding.diploid(population, bve=TRUE, singlestep.active = FALSE,
                               bve.cohorts = c("Founder", "PartlyPhenotyped",
                                               "Phenotyped+Genotyped"))

# Act as if all individuals were genotyped (although they are not!)
population <- breeding.diploid(population, bve=TRUE, bve.all.genotyped = TRUE,
                               bve.cohorts = c("Founder", "PartlyPhenotyped",
                                               "Phenotyped+Genotyped"))

}

{
  # Trait estimation

  set.seed(1)
# Load in some data (Ideally real-world but here from an simulated MoBPS population
data(ex_pop)
pheno <- get.pheno(ex_pop, gen=1:5)
geno <- get.geno(ex_pop, gen=1:5)
haplo <- get.haplo(ex_pop, gen=1:5)
map <- get.map(ex_pop, use.snp.nr=TRUE)

# Estimate marker effects
real.bv.add <- effect.estimate.add(geno, pheno, map)

# Generate new population (or in this case exactly the same) with estimated marker effects
population <- creating.diploid(dataset = haplo, map = map, real.bv.add = real.bv.add)

# First run of breeding.diploid to calculate true genomic values
population <- breeding.diploid(population)

# Check how good the approximation works here:
cor(get.bv(population, gen=1)[1,], get.pheno(ex_pop, gen=1:5)[1,])
cor(get.bv(population, gen=1)[1,], get.bv(ex_pop, gen=1:5)[1,])
}

{
# 1. Trait is direct effect
# 2. Trait is maternal effect
population <- creating.diploid(nsnp = 1000, nindi=10,
                               n.additive = c(100,100),
                               var.target = c(100,20),
                               trait.cor = matrix(c(1,0.5, 0.5, 1), nrow=2),
                               is.maternal = c(FALSE, TRUE),
                               trait.name = c("direct_effect", "maternal_effect"))

# 3. Trait is the sum of both previous traits
population <- add.combi(population, trait = 3, combi.weights = c(1,1), trait.name = "combined_trait")

# Generate some offspring
population <- breeding.diploid(population, breeding.size = 100)

# Looking at trait 2 and mother id
bvs <- get.bv(population, gen=2)
pedi <- get.pedigree(population, gen=2)
cbind(pedi[,3], bvs[2,])
}




# Generation of an inbred founder population
# All plants are saves as sex 0 (male)
dhfounder_population <- creating.diploid(nindi = 200, nsnp = 500, dataset = "homorandom",
                                         sex.quota = 0)

# Generation of a heterozygous founder population
population <- creating.diploid(nindi = 50, nsnp = 500, sex.quota = 0)

# Generation of inbreds via DH technology
dh_population <- breeding.diploid(population, breeding.size = c(200,0), dh.mating = TRUE)

# Generation of inbreds via five generation of selfing
selfing_population <- population
for(index in 1:5){
  selfing_population <- breeding.diploid(selfing_population, breeding.size = c(200,0),
                                         selfing.mating = TRUE)
}

table(get.geno(dhfounder_population, gen=1))
table(get.geno(dh_population, gen=2))
for(gen in 1:6){
  print(table(get.geno(selfing_population, gen=gen)))
}
# Share of heterozygous markers in the selfing population reduces by ~50% per cycle




# Generation of an exemplary VCF-file to load it
# File should contain: 500 SNPs, 30 Individuals
# SNPs are placed on two chromosome
population <- creating.diploid(nsnp = 500, nindi = 30,
                               chr.nr = c(rep(1,250), rep(2, 250)),
                               bp = c(1:250*1000000, 1:250*50000))
get.vcf(population, path="import_test", gen=1)

# Convert provided base-pair position into positions in Morgan
# 1 Morgan = 100.000.000 bp is a common conversion rate
# This can differ between species (e.g. chicken: ~30.000.000 bp)
# or be dependent between telomere / centromere
population <- creating.diploid(vcf="import_test.vcf", bpcm.conversion = 1000000)
summary(population)
# Assume markers to be equidistant and generate a genome of fixed size:
population <- creating.diploid(vcf="import_test.vcf", chromosome.length = c(2,1))
summary(population)


# Generation of an exemplary Ped/map-files to load it
# File should contain: 500 SNPs, 30 Individuals
# SNPs are placed on two chromosome
population <- creating.diploid(nsnp = 500, nindi = 30,
                               chr.nr = c(rep(1,250), rep(2, 250)),
                               bp = c(1:250*1000000, 1:250*50000))
get.pedmap(population, path="import_test", gen=1)

map <- as.matrix(read.table("import_test.map"))
ped <- as.matrix(read.table("import_test.ped"))

# Convert PED-file into haplotype dataset (one haplotype per colum)
nsnp <- (ncol(ped)-6)/2
haplo1 <- ped[,1:nsnp*2+6-1]
haplo2 <- ped[,1:nsnp*2+6]
haplo <- t(rbind(haplo1, haplo2)[c(0,nrow(haplo1)) + sort(rep(1:nrow(haplo1),2)),])

population <- creating.diploid(dataset = haplo, map = map,  bpcm.conversion = 1000000)
summary(population)
# Assume markers to be equidistant and generate a genome of fixed size:
# This will overwrite snp positions from the map file
population <- creating.diploid(dataset = haplo, map = map,
                               snps.equidistant = TRUE, chromosome.length = c(2,1))
summary(population)


# Baseline

# Generation of a founder population
# with 100 individuals and 10,000 SNPs
# The genome contains 5 chromosomes with a length of 2 Morgan each
# Generation of one trait with 60 underlying QTL
# of which 50 are purely additive and 10 have a dominate effect
# The genomic variance of the trait simulated to be 1
# The generated cohort is named "Founder"

pop <- creating.diploid(nsnp=10000, nindi=100,
                        chr.nr=5, chromosome.length=2,
                        n.additive=50, n.dominant=10,
                        name.cohort="Founder",
                        share.genotyped = 1,
                        var.target = 1)

# Generate phenotypic observations for all individuals
# Residual variance is set to result in a heritablity of 0.5
pop <- breeding.diploid(pop, heritability=0.5,
                        phenotyping="all")

# Perform a breeding value estimation using all individuals
# of generation 1 (which are all).
pop <- breeding.diploid(pop, bve=TRUE,
                        bve.gen = 1)

# Generate 100 offspring
# Use the top 20 male and top 20 female based on their BVE
# All male individuals from the "Founder" cohort are used
# as potential sires.
# All female individuals from the "Founder cohort are used
# as potential dams.
# The resulting cohort in named "Offspring".
pop1 <- breeding.diploid(pop, breeding.size=100,
                         selection.size=c(20,20),
                         selection.criteria = "bve",
                         selection.m.cohorts="Founder_M",
                         selection.f.cohorts="Founder_F",
                         name.cohort="Offspring")

# Same procedure, just with a higher selection intensity
# on the male side.
pop2 <- breeding.diploid(pop, breeding.size=100,
                         selection.size=c(5,20),
                         selection.criteria = "bve",
                         selection.m.cohorts="Founder_M",
                         selection.f.cohorts="Founder_F",
                         name.cohort="Offspring")


# Initialize objects to store simulation outputs in
gain1 <- gain2 <- kinship1 <- kinship2 <- numeric(100)
# Run simulation 100 times
for(index in 1:100){
  pop <- creating.diploid(nsnp=10000, nindi=100,
                          chr.nr=5, chromosome.length=2,
                          n.additive=50, n.dominant=10,
                          share.genotyped = 1,
                          name.cohort="Founder",
                          var.target = 1)

  pop <- breeding.diploid(pop, heritability=0.5,
                          phenotyping="all")

  pop <- breeding.diploid(pop, bve=TRUE,
                          bve.gen = 1)

  pop1 <- breeding.diploid(pop, breeding.size=100,
                           selection.size=c(20,20),
                           selection.criteria = "bve",
                           selection.m.cohorts="Founder_M",
                           selection.f.cohorts="Founder_F",
                           name.cohort="Offspring")

  pop2 <- breeding.diploid(pop, breeding.size=100,
                           selection.size=c(5,20),
                           selection.m="function",
                           selection.m.cohorts="Founder_M",
                           selection.f.cohorts="Founder_F",
                           name.cohort="Offspring")

  # store the resulting genomic gains per run
  gain1[index] <- mean(get.bv(pop1, gen=2) - mean(get.bv(pop, gen=1)))
  gain2[index] <- mean(get.bv(pop2, gen=2) - mean(get.bv(pop, gen=1)))

  temp1 <- kinship.emp(population=pop1, gen=2)
  temp2 <- kinship.emp(population=pop2, gen=2)

  # Calculate the average of all off-diagonal values
  kinship1[index] <- mean(temp1-diag(diag(temp1))) * 100/99
  kinship2[index] <- mean(temp2-diag(diag(temp2))) * 100/99
}

save(file="C:/Users/pook/Desktop/Fig1_data.RData", list=c("a", "b", "c","d"))

# Expected gain: i * h * sigma_a = sqrt(0.5) * 0.966 * 10 = 6.83
# Expected gain: i * h * sigma_a = sqrt(0.5) * 1.755 * 10 = 12.4

{
  png("C:/Users/pook/Desktop/Fig1.png", width=2250, height= 1100, res=300)
  par(mfrow=c(2,2))
  par(mar=c(4.5,1.5,0.5,2))
  hist(gain1, xlim=c(0.2,1.65), main="", axes = FALSE, xlab="", breaks = 50)
  abline(v=mean(gain1), col="red", lwd=2)
  axis(side=4)
  axis(side=1, at=seq(0.4,1.6, by=.3))
  title(ylab="Scenario 1", line=0, cex.lab=1.5)
  hist(kinship1, xlim=c(0.004,0.019), main="", axes = FALSE, xlab="", breaks = 25)
  abline(v=mean(kinship1), col="red", lwd=2)
  axis(side=4, ylab="Scenario 1")
  axis(side=1, at=seq(0.005,0.018, by=.003))
  hist(gain2, xlim=c(0.2,1.65), axes=FALSE, main="", xlab="", breaks = 50)
  abline(v=mean(gain2), col="red", lwd=2)
  axis(side=4)
  axis(side=1, at=seq(0.4,1.6, by=.3))
  title(ylab="Scenario 2", line=0, cex.lab=1.5)
  title(xlab="Genomic gain", cex.lab=1.5)
  hist(kinship2, xlim=c(0.004,0.019), axes=FALSE, main="", xlab="", breaks = 40)
  abline(v=mean(kinship2), col="red", lwd=2)
  axis(side=4)
  axis(side=1,at=seq(0.005,0.018, by=.003))
  title(xlab="Gain in avg. kinship", cex.lab=1.5)
  dev.off()
}

{


}


## Json Example


population <- json.simulation(file="C:/Users/pook/Desktop/Simple_cattle.json")
bv.development(population, json = TRUE, bvrow = 1, confidence = 1, development = 1,
               display.creating.type = TRUE, display.sex = TRUE,
               display.cohort.name = TRUE)

png(file="C:/Users/pook/Desktop/bv_development_plot.png", width=2450, height= 1460, res=300)
bv.development(population, json = TRUE, bvrow=1, confidence = 1, development = 1,
               display.creating.type = TRUE, display.sex = TRUE,
               display.cohort.name = TRUE)
dev.off()


population <- json.simulation(file="C:/Users/pook/Desktop/Simple_cattle.json")
bv.development.box(population, json=TRUE, bvrow=1)

png(file="C:/Users/pook/Desktop/bv_development_box_plot.png", width=2450, height= 1460, res=300)
bv.development.box(population, json=TRUE, bvrow=1)
dev.off()

population <- json.simulation(total = ex_json)
kinship.development(population, json=TRUE, display.cohort.name = TRUE)

png(file="C:/Users/pook/Desktop/kinship_development_plot.png", width=2450, height= 1460, res=300)
kinship.development(population, json=TRUE, display.cohort.name = TRUE)
dev.off()

png(file="C:/Users/pook/Desktop/allelefreq_development_plot.png", width=2450, height= 1460, res=300)
analyze.population(population,5,2, gen=1:8)
dev.off()

##### Further tests before submitting new versions:

# Generate multiple correlated traits
population <- creating.diploid(nsnp=1000, nindi = 50, chr.nr = 5)
population <- creating.trait(population, n.additive = c(100), n.qualitative = c(0,50),
                             shuffle.traits = TRUE, shuffle.cor =  cbind(c(1,0.4), c(0.4,1)))

population <- breeding.diploid(population, heritability = 0.3, phenotyping = "non_obs_f")
print(sum(is.na(get.pheno(population, gen=1)))==50)

# Single Step
population <- creating.diploid(nsnp=5000, nindi=100, n.additive = 100, chr.nr=5)
for(index in 1:5){
  population <- breeding.diploid(population, breeding.size = 100, selection.size = c(50,50), share.genotyped = 0.5)
}
population <- breeding.diploid(population, heritability = 0.3, phenotyping = "all")
population <- breeding.diploid(population, bve=TRUE, singlestep.active = TRUE)
if(FALSE){
  plot(get.bve(population, gen=6)[1,])
  points((1:100)[get.genotyped(population, gen=6)], get.bve(population, gen=6)[1,get.genotyped(population, gen=6)], col="red")
  cor(get.bve(population, gen=6)[1,], get.bv(population, gen=6)[1,])
  cor(get.bve(population, gen=6)[1,get.genotyped(population, gen=6)], get.bv(population, gen=6)[1,get.genotyped(population, gen=6)])
  cor(get.bve(population, gen=6)[1,!get.genotyped(population, gen=6)], get.bv(population, gen=6)[1,!get.genotyped(population, gen=6)])
}
print(sum(get.genotyped(population, gen=6))<100 && sum(get.genotyped(population, gen=6)) >0)

#
get.genotyped(population, gen=6)


# BVE
library(MoBPS)
library(RandomFieldsUtils)
library(miraculix)
population <- creating.diploid(nsnp=10000,nindi = 1000, n.additive = c(100,100,100))
population <- breeding.diploid(population, heritability = 0.5, phenotyping = "all")
system.time({pop1 <- breeding.diploid(population, bve=TRUE, miraculix.chol=FALSE)})
system.time({pop1 <- breeding.diploid(population, bve=TRUE, miraculix.chol = TRUE)})
system.time({pop1 <- breeding.diploid(population, bve=TRUE, miraculix.chol = TRUE, miraculix.cores=5)})
pop1 <- breeding.diploid(population, bve=TRUE, sommer.bve = TRUE)
pop1 <- breeding.diploid(population, bve=TRUE, rrblup.bve = TRUE)
pop1 <- breeding.diploid(population, bve=TRUE, emmreml.bve = TRUE)
pop1 <- breeding.diploid(population, bve=TRUE, BGLR.bve = TRUE)
pop1 <- breeding.diploid(population, bve=TRUE, BGLR.bve = TRUE, BGLR.model = "BayesA")


pop1 <- breeding.diploid(population, bve=TRUE, miraculix.chol = TRUE, bve.ignore.traits =c(1,3))

## Generation of ex_pop

ex_pop <- creating.diploid(nsnp=100, nindi=50, miraculix=FALSE, n.additive = 10)
ex_pop <- breeding.diploid(ex_pop, heritability = 0.5)

for(index in 1:5){
  ex_pop <- breeding.diploid(ex_pop, phenotyping = "all")
  ex_pop <- breeding.diploid(ex_pop, breeding.size=50, selection.size = c(5,10),
                             selection.criteria =  "pheno")
}

save(file="C:/Users/pook001/OneDrive - Wageningen University & Research/GitHub/MoBPS/development/MoBPS/data/ex_pop.RData", list=c("ex_pop"), version=2)


library(MoBPS)


dataset <- founder.simulation(nsnp=10000, nindi=100)
population <- creating.diploid(dataset, share.genotyped = 0, n.additive = 100)

population$breeding[[1]][[1]][[1]][[22]]

population <- add.array(population, marker.included = rep(c(0,0,0,0,0,0,0,0,0,1), 1000))


population <- breeding.diploid(population, phenotyping = "all", genotyped.gen = 1,
                               genotyped.array = 2)

population$breeding[[1]][[1]][[1]][[22]]

population <- breeding.diploid(population, phenotyping = "all", genotyped.gen = 1, genotyped.share = 0.5,
                               genotyped.array = 1)

for(index in 1:10){
  print(population$breeding[[1]][[1]][[index]][[22]])

}

data <- get.genotyped.snp(population, gen=1)

get.vcf(population, path="C:/Users/pook/Desktop/test111.vcf", gen=1, non.genotyped.as.missing = TRUE)

get.pedmap(population, path="C:/Users/pook/Desktop/test111", gen=1, non.genotyped.as.missing = TRUE)

get.geno(population, gen=1, non.genotyped.as.missing = TRUE)
get.haplo(population, gen=1, non.genotyped.as.missing = TRUE)


kinship.emp(population = population, database = cbind(5,1,1,6), sym = TRUE)

kinship.emp.fast(population = population, database = cbind(5,1,1,6))

kinship.exp(population = population, database = cbind(5,1,1,6))

data(ex_pop)
get.admixture(ex_pop, gen=4:6, d=3, sort=TRUE)


data(ex_pop)
get.qtl(ex_pop)



kinship.emp(population = population, gen=5)


library(MoBPS)
set.seed(1)

# Generation of a small population with 10 individuals from two distinct pools
population = creating.diploid(nsnp=100, nindi=4, founder.pool = 1)
population = creating.diploid(population = population , nsnp=100, nindi=4,
                              founder.pool = 2)

# Generation of a trait that will only have variablity in pool 1
population = creating.trait(population, n.additive = 100, trait.pool = 1)

# Generation of a trait that will only have variablity in pool 2
population = creating.trait(population, n.additive = 100, trait.pool = 2)

old_bv = get.bv(population, gen = 1)
# Extraction of QTL effects to combine the traits into a joint trait
qtl_effects = get.qtl.effects(population)
real.bv.add = rbind(qtl_effects[[1]][[1]], qtl_effects[[1]][[2]])

# Generate a trait that combines the two pool specific traits into a joint trait
# Replace the old two traits
population = creating.trait(population, real.bv.add = real.bv.add,
                            replace.traits = TRUE)

# Simulate a couple of generations of random mating
population = breeding.diploid(population, breeding.size = 50)
population = breeding.diploid(population, breeding.size = 50)
population = breeding.diploid(population, breeding.size = 4)

new_bv = get.bv(population, gen=1)

old_bv
#M1_1      M2_1      M3_1      M4_1      F1_1      F2_1      F3_1      F4_1
#Trait 1  94.98455  97.95988 100.0000 100.000  98.01574  94.72518 100.0000 100.0000
#Trait 2 100.00000 100.00000 114.7768 108.875 100.00000 100.00000 100.3462 101.3287

new_bv
#M1_1     M2_1     M3_1     M4_1     F1_1     F2_1     F3_1     F4_1
#Trait 1 94.98455 97.95988 114.7768 108.875 98.01574 94.72518 100.3462 101.3287

# Which segment stems from which founder pool?
pools = get.pool(population, gen = 4, plot = TRUE)



library(MoBPS)
set.seed(1)
# Generate a founder population with 5 traits
population = creating.diploid(nindi=50, nsnp=1000, n.additive = rep(100,5),
                              mean.target = 100, var.target = 25)

# Generate some phenotypes and perform a BVE
population = breeding.diploid(population, heritability = 0.3, phenotyping = "all")
population = breeding.diploid(population, bve = TRUE, bve.gen=1)

# Only consider individuals that
# 1. Have a minimum BVE of 99 for trait 1
# 2. Have a minimum underlying true genomic value of 98 for trait 3
# 3. Have a minimum BVE of 199 for trait 4 + trait 5
threshold.selection.index = matrix(c(1,0,0,0,0,
                                     0,0,1,0,0,
                                     0,0,0,1,1), ncol=5, byrow=TRUE)

threshold.selection.value = c(99, 98, 199)
threshold.selection.sign = rep(">", 3)
threshold.selection.criteria = c("bve", "bv", "bve")

# Generate 100 offspring from 3 males / females that fullfil these contrains
# selection otherwise is random
population = breeding.diploid(population,
                              breeding.size = 100,
                              selection.size = c(3,3),
                              selection.criteria = "random",
                              threshold.selection.index = threshold.selection.index,
                              threshold.selection.value = threshold.selection.value,
                              threshold.selection.sign = threshold.selection.sign,
                              threshold.selection.criteria = threshold.selection.criteria)

#Start selection procedure.
#Selection male size:
#  Threshold selection 1 is fullfilled by 21 out of 25.
#Threshold selection 2 is fullfilled by 22 out of 25.
#Threshold selection 3 is fullfilled by 14 out of 25.
#After threshold selection 11 individuals remain to select 3.
#Selection female size:
#  Threshold selection 1 is fullfilled by 19 out of 25.
#Threshold selection 2 is fullfilled by 19 out of 25.
#Threshold selection 3 is fullfilled by 16 out of 25.
#After threshold selection 10 individuals remain to select 3.

# Exemplary further analysis:
# Extract the used parents and look at their estimated breeding values
parent_ids = unique(as.numeric(get.pedigree(population, gen=2, id=TRUE)[,c(2,3)]))
parent_database = get.database(population, id=parent_ids)
get.bve(population, database = parent_database)

#M8_1     M24_1     M25_1      F2_1     F10_1     F24_1
#Trait 1  99.42298 100.46841  99.93170 103.15000 104.41123 101.03186
#Trait 2 100.02388 103.03817  96.80039  99.05030  98.33225 100.87097
#Trait 3 103.56059 101.58174 102.27021 100.29314 101.98670  98.70362
#Trait 4 101.32436 101.60426 101.98240  99.23875 101.63607 102.05797
#Trait 5  98.89409  98.72613  98.91435 102.66770 101.44858  99.52546


library(MoBPS)
set.seed(1)

# Setting up the initial population
population = creating.diploid( nsnp=500, chr.nr =5, nindi = 100)


# distribution of litter size
litter.size = matrix(c(1, 0.70,
                       2, 0.15,
                       3, 0.10,
                       4, 0.05
), byrow=TRUE, ncol=2)

n_litter = 5
#population = breeding.diploid(population, repeat.mating = litter.size,
#                              breeding.size.litter = n_litter)

# Construction of the fertility trait with discrete phenotypic realisations
# This is a trait with only purely dominate effects
# (to make it highly related to inbreeding)

population = creating.trait(population, trait.name = "Fertility",
                            n.additive = 250,
                            mean.target = 100, var.target = 100)

population = breeding.diploid(population, heritability = 1/10)

# The trait on default has a gaussian distribution
# here I am adding code to make phenotypic observations be discrete natural numbers

sigma_ferti = sqrt(1000)
# This is the overall phenotypic variance of the fertility trait
prob_cum = cumsum(litter.size[,2])

litter.function = function(x){
  if(x < qnorm(prob_cum[1], mean=100, sd=sigma_ferti)){
    y = 1
  } else if(x < qnorm(prob_cum[2], mean=100, sd=sigma_ferti)){
    y = 2
  } else if(x < qnorm(prob_cum[3], mean=100, sd=sigma_ferti)){
    y = 3
  } else {
    y = 4
  }
  # IF you dont like If loops and want it efficiently do this instead
  # y = sum(  qnorm(prob_cum, mean=100, sd=sigma_ferti) < x) +1
  return(y)
}

# Traits that control the litter size
# should only have natural numbers as potential realisations

population = creating.phenotypic.transform(population,
                          phenotypic.transform.function = litter.function,
                                           trait = 1)

population = breeding.diploid(population,
                              repeat.mating = "genetic",
                              repeat.mating.trait = 1,
                              repeat.mating.max = 4)

n_litter = 50
# Simulate 100 years of random mating
for(generation in 2:100){
  population = breeding.diploid(population, breeding.size.litter = n_litter,
                                selection.criteria = "random")
}

# As individuals with better fertility have more offspring
# fertility increases over time

bv = numeric(100)
for(index in 1:100){
  bv[index] = mean(get.bv(population, gen=index))
}
plot(bv, type="l", xlab = "generation", ylab = "genomic value")

plot(2:get.ngen(population), rowSums(get.size(population))[-1]/n_litter, main="Litter size",
     xlab = "generation", ylab = "avg. litter ")

# The increase in the trait reduces over time as the
# difference in litter size between animals reduces over time

lines(ksmooth(2:get.ngen(population), rowSums(get.size(population))[-1]/n_litter,
              bandwidth = 10), col="red", lwd=2)



library(MoBPS)
set.seed(1)
# Generation of an initial breeding population
population <- creating.diploid(nsnp = 1000, nindi = 200, name="Breeding_Population_1")

for(index in 1:20){

  # Generation of new offspring
  population <- breeding.diploid(population, breeding.size = 100,
                          name=paste0("Offspring_",index),
                          selection.m.cohorts = paste0("Breeding_Population_", index, "_M"),
                          selection.f.cohorts = paste0("Breeding_Population_", index, "_F"),
                          time.point = index)

  # Calculate how many animals are alive
  alive_m <- sum(get.class(population,
                           cohorts = paste0("Breeding_Population_", index, "_M"))==0) + 50
  alive_f <- sum(get.class(population,
                           cohorts = paste0("Breeding_Population_", index, "_F"))==0) + 50

  # individuals from the offspring and breeding population are merged into a joint
  # breeding population for the next cycle

  population <- breeding.diploid(population, combine = TRUE,
                        selection.m.cohorts = c(paste0("Breeding_Population_", index, "_M"),
                                               paste0("Offspring_", index, "_M")),
                        name.cohort = paste0("Breeding_Population_", index+1, "_M"))

  population <- breeding.diploid(population, combine = TRUE,
                        selection.f.cohorts = c(paste0("Breeding_Population_", index, "_F"),
                                               paste0("Offspring_", index, "_F")),
                        name.cohort = paste0("Breeding_Population_", index+1, "_F"))

  # All individuals are culled with a probability of 0.3
  ### IF different culling probabilites per age are wanted culling could be
  ### applied on the Offspring cohorts
  ### if one copy of an individuals is culled, all other copies are culled as well
  ### culling.share1 can also be a vector to assign each individual will a different
  ### culling probability

  population <- breeding.diploid(population, culling.cohorts =
                                   paste0("Breeding_Population_", index+1, c("_M", "_F")),
                                 culling.share1 = 0.3)

}

# Extract the Age distribution in the final breeding population
age <- 21 - get.age.point(population,
                          cohorts=paste0("Breeding_Population_", 21, c("_M", "_F")) )
table(age)
# age
#1   2   3   4   5   6   7   8   9  10  11  16  17  21
#100  73  54  32  23  18  17   9   6   7   6   1   1   1
hist(age, main="Distribution of the age of individuals in the breeding population", nclass=50)


set.seed(1)

# Generation of a baseline population
# Generation of a trait with slightly different genetic architecture in two environments

population <- creating.diploid(nsnp=5000, nindi=100, n.additive = c(100,100),
                               share.genotyped = 1,
                               trait.cor = matrix(c(1,0.9,0.9,1), ncol=2))

# Linking of the two traits (environments)
population <- combine.traits(population, combine.traits = 1:2)

# Collection of phenotypic data
# 50 lines are phenotypes in either of the two environments
population <- breeding.diploid(population, phenotyping.database = cbind(1,1),
                               heritability = 0.3, n.observation = c(1,0))
population <- breeding.diploid(population, phenotyping.database = cbind(1,2),
                               heritability = 0.3, n.observation = c(0,1))


# Breeding value estimation on the entire set
population <- breeding.diploid(population, bve=TRUE, bve.gen = 1)

# Details on the accuracy of the BVE on the individuals sets
#(for comparability later)
analyze.bv(population, database=cbind(1,1))
analyze.bv(population, database=cbind(1,2))

# Breeding value estimation when only using one of the two sets
population <- breeding.diploid(population, bve=TRUE, bve.database = cbind(1,1))
population <- breeding.diploid(population, bve=TRUE, bve.database = cbind(1,2))


# Generation of a baseline population
# Generation of a trait with slightly different genetic architecture in two environments

set.seed(1)

population <- creating.diploid(nsnp=5000, nindi=100, n.additive = c(100,100),
                               share.genotyped = 1,
                               trait.cor = matrix(c(1,0.9,0.9,1), ncol=2))

# Collection of phenotypic data
# 50 lines are phenotypes in either of the two environments
population <- breeding.diploid(population, phenotyping.database = cbind(1,1),
                               heritability = 0.3, n.observation = c(1,0))
population <- breeding.diploid(population, phenotyping.database = cbind(1,2),
                               heritability = 0.3, n.observation = c(0,1))


# Breeding value estimation using a multi trait model
population <- breeding.diploid(population, bve=TRUE,
                               sommer.multi.bve = TRUE, bve.gen = 1)

# Details on the accuracy of the BVE on the individuals sets
#(for comparability later)
analyze.bv(population, database=cbind(1,1))
analyze.bv(population, database=cbind(1,2))

# Breeding value estimation when only using one of the two sets
population <- breeding.diploid(population, bve=TRUE, bve.database = cbind(1,1))
population <- breeding.diploid(population, bve=TRUE, bve.database = cbind(1,2))



##

# Generation of a basic breeding population
set.seed(1)
population = creating.diploid(nindi = 100, nsnp = 100, n.additive = c(100, 200),
                              mean.target = 100, var.target = 100)

# Simulate 10 generations of selection high double weighted on the second trait
for(index in 1:10){
  population = breeding.diploid(population, phenotyping.gen = index,
                                heritability = 0.3,
                                bve = TRUE, bve.gen = index)
  population = breeding.diploid(population, breeding.size = 100,
                                selection.size = c(10,10),
                                multiple.bve.weights.m = c(1,2))
}

# Visualization of genetic gains
bv.development(population, gen=1:11)

# Generation of new individuals from founders that is on a similar level than last generation
# Introduce three difference pools on similar genetic level
population = add.diversity(population, pool.gen = 1, target.gen = 11, add.gen = 12)

#Required 10 generations to generate introduced material.
#Avg. genomic value for traits:
#  Trait 1  Trait 2
#132.2063 145.4372
#Compared to target reference:
#  Trait 1  Trait 2
#134.6005 148.0267

population = add.diversity(population, pool.gen = 1, target.gen = 11, add.gen = 13)
population = add.diversity(population, pool.gen = 1, target.gen = 11, add.gen = 14)

# Get an overview of the population structure
get.pca(population, gen = c(11:14), coloring = "gen")




#### GxE

# Generation of a basic population
set.seed(1)
population = creating.diploid(nsnp=1000, nindi=50)

# Generation of two traits that are both evaluated in three locations
population = creating.trait(population, n.additive = c(50,30),
                            n.locations = 2,
                            gxe.max = 0.6, gxe.min = 0.4,
                            shuffle.cor = cbind(c(1,0.3), c(0.3, 1)))

# Check for correlation between different traits
round(cor(t(get.bv(population, gen=1))), digits = 2)
#Trait 1 x Location 1 Trait 2 x Location 1 Trait 1 x Location 2 Trait 2 x Location 2
#Trait 1 x Location 1                 1.00                 0.24                 0.39                 0.02
#Trait 2 x Location 1                 0.24                 1.00                 0.03                 0.45
#Trait 1 x Location 2                 0.39                 0.03                 1.00                 0.13
#Trait 2 x Location 2                 0.02                 0.45                 0.13                 1.00

# Phenotyping of all traits in location 2

locs = get.index(population, locations=2)
population = breeding.diploid(population, n.observation = locs, phenotyping.gen = 1,
                              heritability = 0.3)
get.pheno(population, gen = 1)[,1:3]
#M1_1     M2_1     M3_1
#Trait 1 x Location 1       NA       NA       NA
#Trait 2 x Location 1       NA       NA       NA
#Trait 1 x Location 2 103.6944 103.7062 100.8905
#Trait 2 x Location 2 118.9703 125.8905 105.0064

# selecting individuals based on location 2 and double weighting on trait 2
trait_index = get.index(population, location.weights = c(0,1,0), trait.weights = c(1,2))
selected_indi = breeding.diploid(population, selection.size = c(5,5),
                                 selection.criteria = "pheno",
                                 export.selected = TRUE,
                                 multiple.bve.weights.m = trait_index)


{

  # Generate a founder population
  population <- creating.diploid(nsnp=1000, nindi = 100, n.additive = 100)

  # Make sure the simulated trait has a set mean and variance
  population <- bv.standardization(population, gen=1, mean.target = 100, var.target = 10)

  # Generate offspring
  population <- breeding.diploid(population, breeding.size = 500, selection.size = c(50,50))

  ### Plain application of culling with a set probability

  # Apply culling with a probability of 0.3 to all individuals
  pop1 <- breeding.diploid(population, culling.gen=2, culling.share1 = 0.3)

  # Death animals put into class (-1)
  # This means that these individuals on default will not be used for reproduction anymore
  # These individuals can still be used in a breeding value estimation
  table(get.class(pop1, gen=2))
  #-1   0
  #170 330

  ### Apply culling with a probability depending on the genomic value

  pop2 <- breeding.diploid(population, culling.gen=2, culling.bv1 = 100, culling.share1 = 0.5, culling.bv2 = 105, culling.share2 = 0.1,
                           culling.index = 1)

  class <- get.class(pop2, gen=2)
  bv <- get.bv(pop2, gen=2)[1,]

  par(mfrow=c(2,1))
  hist(bv[class==0], xlim=c(90,110), main="Genomic value distribution of alive individuals", xlab="genomic value")
  hist(bv[class==(-1)], xlim=c(90,110), main="Genomic value distrition of death individuals", xlab="genomic value")

}

{


  # Initialize covariance matrix of pen & litter effects
  cov_litter = cbind(c(100,-50), c(-50,100))
  cov_pen = cbind(c(100,-50), c(-50, 100))

  # Generation of the baseline population
  # Founders of the population do not have pen or litter effects
  population = creating.diploid(nsnp = 100, nindi = 100, n.additive = c(100,100),
                                mean.target = 100, var.target = 100,
                                litter.effect.covariance = cov_litter,
                                pen.effect.covariance = cov_pen)


  population = breeding.diploid(population, breeding.size = 100,
                                repeat.mating = matrix(c(2, 0.5,
                                                         3, 0.5), ncol = 2, byrow = TRUE),
                                pen.size = matrix(c(5, 0.4,
                                                    6, 0.6), ncol=2, byrow = TRUE),
                                pen.by.sex = TRUE)


  population = breeding.diploid(population, heritability = 0.3, phenotyping = "all")

  # Extract which individual is in which pen / litter
  get.litter(population, gen = 2)
  get.pen(population, gen = 2)

  # size of the pen & litter effects
  get.litter.effect(population, gen = 2)
  get.pen.effect(population, gen = 2)

  # As a piece of warning!
  # Non of the R-package packages software for breeding value estimation include litter & pen effects
  # wont be include in calculated sigma_e / sigma_g in mobps.bve either!
  population = breeding.diploid(population, bve = TRUE)


}

library(drat)
?drat
drat::insertPackage("C:/Users/pook/Documents/GitHub/MoBPS/RandomFieldsUtils_1.0.6.tar.gz", "C:/Users/pook/Documents/GitHub/drat/")
drat::insertPackage("C:/Users/pook/Documents/GitHub/MoBPS/miraculix_1.0.5.tar.gz", "C:/Users/pook/Documents/GitHub/drat/")
drat::insertPackage("C:/Users/pook/Documents/GitHub/alstructure_0.1.0.tar.gz", "C:/Users/pook/Documents/GitHub/drat/")
drat::insertPackage("C:/Users/pook/Documents/GitHub/MoBPS/MoBPSmaps_0.1.12.tar.gz", "C:/Users/pook/Documents/GitHub/drat/")
devtools::check_win_release()


population <- creating.diploid(nsnp = 100, nindi = 100)
population <- breeding.diploid(population, breeding.size = 100, selection.size=c(25,0), selection.m.gen = 2)

## CRAN Notes 1
This is a re-submission with following changes:
* The description in the DESCRIPTION file has been updated and doi/ref have been added
* All uses of cat/print (except in summary()) can now be deactivated by use of verbose in the respective function
* Some of the uses of cat/print have been changed to warning() / stop()
* Use of par() has been reduced when possible but if still required use of on.exit() has been added
* requireNamespace() has been added for all uses of suggested packages
* graphics, stats, utils have been set to Imports (instead of Depends)
* ex_pop has been replaced by a larger dataset that is leading to nicer figures for some examples
* examples have been revisioned and use of donttest{} has been added for lengthly examples (instead of #)
* Minor updates to functions (mainly breeding.diploid) to improve run-time, reduce number of manditory parameters
*** An exact list of these minor changes can be found at https://github.com/tpook92/MoBPS/commits/master
* An updated version of RandomFieldsUtils has been added to the Additional_repositories list
)


population <- creating.diploid(nsnp=5000, nindi=100, n.additive = 100, chr.nr=5)
population <- breeding.diploid(population, breeding.size = 100, selection.size = c(50,50), share.genotyped = 0.5)


## CRAN Notes 1
This is a resubmission of the MoBPS package
The new version fixes the issues brought up by Brian Ripley on April 28
Additionally dependency issues to the miraculix package from our submission on May 10 were fixed
In additional various additional updates to the functionality of the package have been done
that require new version of enhancing R-packages RandomFieldsUtils / miraculix / MoBPSmaps (see drat repository / Additional_repositories)
An exact list of all changes can be found at https://github.com/tpook92/MoBPS/commits/master
however non of the content-related changes should be critical for CRAN submission.

## CRAN Notes 2
This is a resubmission of the MoBPS package.
This new version fixes the issues brought up by Brian Ripley on September 5.
This mostly means that the example in json.simulation() now runs without the MoBPSmaps package.
Additionally, temporary files generated in get.pedmap and get.vcf are now removed in the examples (issed by Uwe Ligges on October 30).
In addition, various additional updates to the functionality of the package have been done.
An exact list of all changes can be found at https://github.com/tpook92/MoBPS/commits/master
however non of the content-related changes should be critical for CRAN submission.

We decided to not wait till miraculix reappears on CRAN since the submission process of the associated packages
RandomFieldsUtils and RandomFields will most likely still take quite a while on those packages are only
enhancing MoBPS but not manditory in any way.
