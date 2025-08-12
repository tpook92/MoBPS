set.seed(42)
library(MoBPS)
library(MoBPSmaps)

# Down-sample maize array to 10k SNPs
map <- MoBPSmaps::map_maize1[sort(sample(1:nrow(MoBPSmaps::map_maize1), 10000)),]
# Sample target allele frequencies for the flint and dent pool
allele_freq_flint <- runif(10000,0,1)
allele_freq_dent <- runif(10000,0,1)

# Looking at the allele frequency spectrum
hist(allele_freq_dent)

# Other allele frequency spectra you could think of
hist(stats::rbeta(1000, 1,1))
hist(stats::rbeta(1000, 0.3,1))

# Generation of the founder population

# Sex can usually be neglected in plant breeding simulations
# Still, each individual / animal / line is assigned with a sex
# We would highly encourage simulate all lines with the same sex OR use the sex to differantiate between gene pools etc
# to give the simulation more structure
population <- creating.diploid(nindi=100, map = map, sex.quota=0,
                               freq = allele_freq_flint,
                               name = "Flint_material",
                               n.additive = 500, n.equal.dominant = 500)

# Make sure Dent tester plants are fully inbred
population <- creating.diploid(population = population, dataset="homorandom", nindi=1,
                               sex.quota=1, freq = allele_freq_dent, name="Dent_tester_1")
population <- creating.diploid(population = population, dataset="homorandom", nindi=1,
                               sex.quota=1, freq = allele_freq_dent, name="Dent_tester_2")
population <- creating.diploid(population = population, dataset="homorandom", nindi=1,
                               sex.quota=1, freq = allele_freq_dent, name="Dent_tester_3")

# Having a look at the population structure
get.pca(population, cohorts="Flint_material")
get.pca(population, gen=1)
get.pca(population, gen=1, components = c(3,4))

# Generate crosses from the founding material
population <- breeding.diploid(population,
                               selection.m.cohorts = "Flint_material",
                               selection.size = c(100,0),
                               breeding.size = c(1000,0),
                               name.cohort = "Flint_cross")

# Generate DH lines from the crosses - this includes the entire DH generation process including meiosis
# and doubling of one of the haplotypes
population <- breeding.diploid(population,
                               selection.m.cohorts = "Flint_cross",
                               selection.size = c(1000,0),
                               breeding.size = c(1000,0),
                               dh.mating = TRUE,
                               share.genotyped = 1,
                               max.offspring = 1,
                               name.cohort = "Cross_DHs")

# Generate lines in a yield trial
population <- breeding.diploid(population, selection.m.cohorts = "Cross_DHs",
                               selection.f.cohorts = "Dent_tester_1",
                               breeding.size = c(1000,0), name.cohort = "Yield_trial_1",
                               max.offspring = c(1,1000))

# Generate phenotypes within the yield trial
population <- breeding.diploid(population, heritability = 0.3,
                               phenotyping.cohorts = "Yield_trial_1")

get.pheno(population, cohorts = "Yield_trial_1")

# Extract the average phenotypes of the offspring for each cross
population <- breeding.diploid(population,
                               offpheno.parents.cohorts = "Cross_DHs",
                               offpheno.offspring.cohorts = "Yield_trial_1")

get.pheno.off(population, cohorts = "Cross_DHs")
get.pheno.off.count(population, cohorts = "Cross_DHs")


# Perform a breeding value estimation for DHs
# use rrblup to estimate variance components as Direct estimation in MoBPS will just take
# variance components from last phenotyping

population <- breeding.diploid(population, bve=TRUE,
                               bve.cohorts = "Cross_DHs",
                               rrblup.bve = TRUE,
                               input.phenotype = "off")

plot(get.pheno.off(population, cohorts = "Cross_DHs")[1,],
     get.bve(population, cohorts = "Cross_DHs")[1,])
plot(get.bve(population, cohorts = "Cross_DHs")[1,])
breeding.diploid(population, selection.m.cohorts = "Cross_DHs",
                 selection.criteria = "bve",
                 copy.individual.m = TRUE,
                 selection.size=c(10,0),
                 name.cohort = "Cross_DHs_sel",
                 export.selected = TRUE)

# Select the top 250 DHs for the second yield trial - these are the same lines as before!
# BTW: it is not necessary to generate this cohorts. It just can make it easier for downstream analysis
population <- breeding.diploid(population, selection.m.cohorts = "Cross_DHs",
                               selection.criteria = "bve",
                               copy.individual.m = TRUE,
                               selection.size=c(250,0),
                               name.cohort = "Cross_DHs_sel")

# Yield Trial 2
population <- breeding.diploid(population, breeding.size = c(750,0),
                               selection.f.cohorts = c("Dent_tester_1", "Dent_tester_2", "Dent_tester_3"),
                               selection.m.cohorts =  "Cross_DHs_sel",
                               breeding.all.combination = TRUE,
                               name.cohort = "Yield_trial_2")


population <- breeding.diploid(population, heritability = 0.3,
                               phenotyping.cohorts = "Yield_trial_2")

population <- breeding.diploid(population, offpheno.parents.cohorts = "Cross_DHs_sel",
                               offpheno.offspring.cohorts = "Yield_trial_2")

# Perform a breeding value estimation for DHs
# use rrblup to estimate variance components as Direct estimation in MoBPS will just take
# variance components from last phenotyping
population <- breeding.diploid(population, bve=TRUE, bve.cohorts = "Cross_DHs_sel",
                               rrblup.bve = TRUE, input.phenotype = "off")

# Select the top 5 DHs based on the second yield trial
population <- breeding.diploid(population, selection.m.cohorts = "Cross_DHs_sel",
                               copy.individual.m = TRUE,
                               selection.size=c(5,0), name.cohort = "Cross_DHs_final")


summary(population)

# This is the genomic value of the DHs itself - not for a cross to a dent line
par(mfrow=c(3,1))
hist(get.bv(population, cohorts="Cross_DHs")[1,], xlim=c(250,440), xlab="genomic value", main ="Genomic value in initial DHs")
hist(get.bv(population, cohorts="Cross_DHs_sel")[1,], xlim=c(250,440), xlab="genomic value", main ="Genomic value in selected DHs")
hist(get.bv(population, cohorts="Cross_DHs_final")[1,], xlim=c(250,440), xlab="genomic value", main ="Genomic value in final DHs")

#################################################################################################
####### END of the Task - everthing below is just showing of some additional MoBPS functionality
#################################################################################################

# In case you want to perform specific mating you can also manually tell the tool who to make with whom
# E.g. perform generation of the Yield_trial_2 without use of breeding.all.combination

# figure out position of the first parent
get.database(population, cohorts="Cross_DHs_sel")
#[1,]    5    1    1  250

# figure out position of the first parent
get.database(population, cohorts=c("Dent_tester_1", "Dent_tester_2", "Dent_tester_3"))
#[1,]    1    2    1    3

# Each row codes the mating between two individuals (column 1-3; column 4-6)
# This will mate the individual stored in generation 5, sex 1, nr 1 with individual generation 1, sex 2, nr 1
fixed.breeding <- cbind(5,1,1, 1,2,1)
population <- breeding.diploid(population, fixed.breeding = fixed.breeding, name.cohort = "Fixed_example")

get.id(population, cohorts=c("Dent_tester_1", "Dent_tester_2", "Dent_tester_3"))
get.id(population, cohorts="Cross_DHs_sel")
fixed.breeding.id = cbind(1337, 101)

# This is how you could set up your mating structure with fixed breeding

fixed.breeding <- matrix(0, nrow=750, ncol=6)
# Dent tester 1
fixed.breeding[1:250,] <- cbind(5,1,1:250,1,2,1)
# Dent tester 2
fixed.breeding[251:500,] <- cbind(5,1,1:250,1,2,2)
# Dent tester 3
fixed.breeding[501:750,] <- cbind(5,1,1:250,1,2,3)

population <- breeding.diploid(population, breeding.size = c(750,0),
                               fixed.breeding = fixed.breeding, name.cohort = "Yield_trial_2_alt")

