

############# Example Selfing / DHs

library(MoBPS)
set.seed(42)

# Import genomic data from the first pool
population <- creating.diploid(vcf="C:/Users/pook/Desktop/MoBPS_workshop_gft/Task6/Pool1.vcf",
                               name.cohort = "Pool1", sex.quota = 0)


# Import genomic data from the second pool
population <- creating.diploid(population = population, vcf="C:/Users/pook/Desktop/MoBPS_workshop_gft/Task6/Pool2.vcf",
                               name.cohort = "Pool2", sex.quota = 1)


population <- breeding.diploid(population, breeding.size = c(100,0),
                               selection.m.cohorts = "Pool1",
                               name.cohort = "Selfed_plant_1",
                               selfing.mating = TRUE,
                               selfing.sex = 0)

population <- breeding.diploid(population, breeding.size = c(100,0),
                               selection.m.cohorts = "Pool1",
                               name.cohort = "DHs",
                               dh.mating = TRUE,
                               dh.sex = 0)

get.pedigree(population, cohorts = "Selfed_plant_1")


mean(get.geno(population, cohorts="Pool1")==1)

# Selfed plant have about half the number of heterozygous variants
mean(get.geno(population, cohorts="Selfed_plant_1")==1)

# DHs are fully homozygous
mean(get.geno(population, cohorts="DHs")==1)



############ Example:  Group Phenotypes

population <- creating.diploid(nsnp=100, nindi=10, n.additive = 100)

population <- breeding.diploid(population, heritability = 0.3, phenotyping = "all")

pheno <- get.pheno(population, gen=1)

# Individuals 1-3 are group phenotyped
pheno[1:3] <- mean(pheno[1:3])

# Enter new phenotypic values for individuals
new_phenos <- cbind(colnames(pheno), pheno[1,])
population <- insert.pheno(population, phenos = new_phenos)

get.geno(population, gen=1)





############# Example: Traits with non-gaussian phenotypes // Phenotypic transformation

set.seed(42)
population <- creating.diploid(nsnp = 1000, nindi = 100, n.additive = 100,
                               mean.target = 100, var.target = 50)

# Goal is to have a trait with binary phenotypes (YES / NO)
# 30% of all individuals should have the YES realisation
# Heritability = 0.5 - genetic variance is 50 -> phenotypic variance is 100

trafo <- function(x){
  if(x < qnorm(0.7, mean=100,sd=10)){
    y <- 0
  } else  {
    y <- 1
  }
  return(y)
}

population <- creating.phenotypic.transform(population, trait = 1,
                                            phenotypic.transform.function = trafo)


population <- breeding.diploid(population, heritability = 0.5,
                               phenotyping = "all")

# Look at the distribution of the phenotypes
# Performing a BVE with our direct model can be problematic
# as genomic SD is not in line with phenotypic SD
# You will get a warning - and should use a method that estimates variance components instead
# e.g. rrBLUP

table(get.pheno(population, gen=1))

population <- breeding.diploid(population, bve = TRUE, bve.gen=1)
population <- breeding.diploid(population, bve = TRUE, rrblup.bve = TRUE, bve.gen=1)

########### Example: Use more SNPs than in VCF
# This does not provide a good method to generate additional SNPs in linkage
# This could for example be done with imputation , haplotype blocks etc.

vcffile <- as.matrix(read.table("C:/Users/pook/Desktop/MoBPS_workshop_gft/Task6/Pool1.vcf"))

# Extract haplotypes from the dataset ((use the R-package vcfR to make this easier!))
haplo1 <- substr(vcffile[,-(1:9)], start=1, stop=1)
haplo2 <- substr(vcffile[,-(1:9)], start=3, stop=3)

# Order or haplotypes need to be changes so that column 1&2 are first individual, 3&4 are second etc.
dataset <- cbind(haplo1, haplo2)[, c(rep(c(0,10),10) + sort(rep(1:10,2)))]
storage.mode(dataset) <- "integer"

# add 5000 additional SNP - these are in no linkage to existing variants
dataset <- rbind(dataset, matrix(rbinom(5000*10*2,1, prob=0), ncol=20))

# Map file for real markers
map <- cbind(vcffile[,1], vcffile[,3], NA,  as.numeric(vcffile[,2]) )

# Map file for additional markers
# all markers are put on chromosome 1, BP positions are drawn at random

map <- rbind(map, cbind(1, paste0("AdditionalSNP", 1:5000), NA, sort(round(runif(5000, min=1, max=1e8)))))

# In case a map with bp position but no morgan position is provided the Morgan position will be calculated based on the BP position
# conversion rate can be given in bpcm.conversion
population <- creating.diploid(dataset = dataset,
                               map = map)

summary(population)

# Looking at the SNP positions
map <- get.map(population)
plot(map[,3])

# In case no map is provided SNPs will be placed to be equidistant (which is computationally faster downstream)
population <- creating.diploid(dataset = dataset,
                               map = map, snps.equidistant = TRUE)

map <- get.map(population)
plot(map[,3])
plot(rowMeans(get.geno(population, gen=1)))

############ Example Own GWAS methodology
library(MoBPS)

set.seed(1)
dataset <- founder.simulation(nsnp=1000)
population <- creating.diploid(dataset = dataset)

geno <- get.geno(population, gen=1)
hist(rowMeans(geno))

rowMeans(geno)[c(126,577,806)]
# here I am making sure that the QTLs I am selecting are variable

QTLs <- matrix(c(126,1,0,1,2,
                 577,1,0,1,2,
                 806, 1,0,1,2), byrow=TRUE, ncol=5)

population <- creating.trait(population, real.bv.add = QTLs)

population <- breeding.diploid(population, heritability = 0.5, phenotyping = "all")

pheno <- get.pheno(population, gen=1)

# Use any software for GWAS analysis

# Putting the data from MoBPS in to the syntax needed by the external software
# get.vcf/ get.pedmap can be very useful to write output files
# If you need additional output format that are commonly use please notify us =)

ge <- data.frame(marker = 1:1000, chrom = rep(1,1000), pos = 1:1000,
                 geno, check.names = FALSE)

ph <- data.frame(line = colnames(geno), y = pheno[1,])

a <- rrBLUP::GWAS(pheno = ph, geno = ge, plot = FALSE)

par(mfrow=c(1,1))
plot(a$y)

abline(v = c(126, 577, 806), col="red")

