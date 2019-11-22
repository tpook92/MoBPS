# This will install the most recent version of the R-package
# install.packages("devtools")
# This is the offical release of MoBPS (1.4.14)
# devtools::install_github("tpook92/MoBPS", subdir="pkg")
# This is the version of MoBPS including most recent updates (1.4.22)
# devtools::install_github("tpook92/MoBPS", subdir="development_version")

install.packages("C:/Users/pook/Desktop/R-Stuff/MoBPS_1.4.22.tar.gz", repos = NULL, type = "source")

# To import json-files to R you need the package "jsonlite"
install.packages("jsonlite")

# This will load in the R-package
library(MoBPS)
population <- json.simulation(file="C:/Users/pook/Downloads/Simple_Cattle.json")

summary(population)
'#
Population size:
Total: 8860 Individuals
Of which 3110 are male and 5750 are female.
There are 23 generations
and 46 unique cohorts.

Genome Info:
There are 3 unique chromosomes.
In total there are 3000 SNPs.
The genome has a total length of 3 Morgan.
The genome has a physical size of about: 0.2998 GB

Trait Info:
There are 3 modelled traits.
Of which 3 have underlying QTL.
Trait names are:Milk Fat Protein
Highest correlation between genetics of traits is 0.4.
Highest correlation between enviromental effects is 0.3
'#

phenotypes <- get.pheno(population, cohorts = "CowsSecondYear")
# Each colum contains the phenotypes for one individual
hist(phenotypes[1,], main="Distribution Milk Yield", xlab="milk yield")

# All markers were simulated independently - no LD for the first generation
ld.decay(population, cohorts="Bull", chromosome= 1, step=5)
# After 5 generations there is some LD build-up
ld.decay(population, cohorts="NewBulls_5", chromosome=1, step=5)
# LD for the group of selected Bulls is much higher
ld.decay(population, cohorts="SelectedBulls_5", step=5)

get.pedmap(population, path="C:/Users/pook/Desktop/temp1", cohorts="NewBulls_5")
