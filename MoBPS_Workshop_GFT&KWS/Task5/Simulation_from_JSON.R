set.seed(42)
library(MoBPS)

# Simulation of the Breeding scheme entered in the Web-interface
population <- json.simulation(file="C:/Users/pook/Desktop/MoBPS_workshop_gft/Task2/Simple Sheep Advanced_baseline.json")

# The population is a LONG list
# to look at it via str(population) for similar is a mess
# This list contains information on all cohorts ever generated in your breeding scheme
# Typically you would overwrite your previous version with the new object after each breeding.diploid call

# Better use:
# - summary() for a brief overview
# - get.XXX with XXX being whatever you are interested in (get.geno, get.recombi, get.pedigree)
summary(population)

get.cohorts(population)

cohorts_to_analyze <- paste0("1yearEwesNext_", 1:20)

BVS <- numeric(20)
for(index in 1:20){
  BVS[index] <- mean(get.bv(population, cohorts = cohorts_to_analyze[index])[1,])
}

plot(BVS, ylab="genomic value", xlab = "generation", type="l", lwd=2,
     main="Genomic value for 1yearEwes")
