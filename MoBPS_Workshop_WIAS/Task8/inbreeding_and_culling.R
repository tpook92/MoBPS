
## Code for additions 1 & 2 is very similar
# use a loop to avoid having the same code multiple times
# initialize: scenario = 1 to run individual simulation line-by-line


library(MoBPS)

# scenario 1: base-line
# scenario 2: 40% of animals are culled before selection
# scenario 3: at max 2 animal per half-sib family are selected

inbreeding = NULL
genetic_gain = NULL

for(scenario in 1:3){

  set.seed(42)
  ### Random Seed is set to be the same for the three simulations
  # to start with the same founder population
  # founder genotypes & QTLs can have a major impact on subsequent results
  # Traits are standardized to more easily interprete plots ( + 1  = 1 genetic SD)
  population = creating.diploid(nsnp = 1000, nindi = 100, n.additive = 100,
                                chr.nr = 5, chromosome.length = 3,
                                mean.target = 100, var.target = 1)


  for(index in 1:10){

    # phenotyping of current population
    current.gen = get.ngen(population)
    population = breeding.diploid(population, heritability = 0.3,
                                  phenotyping.gen = current.gen)

    # culling only happens in Addition 2:
    if(scenario == 2){
      population = breeding.diploid(population, culling.share1 = 0.4,
                                    culling.gen = current.gen)
    }

    # to reduce inbreeding avoid matings between siblings
    if(scenario == 3){
      max.selection.halfsib = c(Inf,Inf)
      max.selection.fullsib = c(Inf,Inf)
      avoid.mating.halfsib = TRUE
    } else{
      max.selection.halfsib = c(Inf,Inf)
      max.selection.fullsib = c(Inf,Inf)
      avoid.mating.halfsib = FALSE
    }

    population = breeding.diploid(population, breeding.size = 100,
                                  selection.size = c(10,10),
                                  max.selection.halfsib = max.selection.halfsib,
                                  max.selection.fullsib = max.selection.fullsib,
                                  avoid.mating.halfsib = avoid.mating.halfsib,
                                  selection.m.database = cbind(current.gen,1),
                                  selection.f.database = cbind(current.gen,2),
                                  selection.criteria = "pheno")


  }

  # evaluate the outcomes of a given simulation
  inb = bv = numeric(11)
  for(index in 1:11){
    bv[index] = mean(get.bv(population, gen = index))
    inb[index] = mean(inbreeding.emp(population, gen = index))
  }

  inbreeding = rbind(inbreeding, inb)
  genetic_gain = rbind(genetic_gain, bv)


}

# Comparision of the different scenarios (this is just based on a single replicate!):
par(mfrow = c(1,2)) # to generate two plots next to each other

plot(genetic_gain[1,], lwd = 2, type = "l", xlab = "generation", ylab = "genetic gain", main = "Genetic gain",
     ylim = c(min(genetic_gain), max(genetic_gain)))
lines(genetic_gain[2,],lwd = 2,  col = "red")
lines(genetic_gain[3,], lwd = 2, col = "green")

plot(inbreeding[1,], lwd = 2, type = "l", xlab = "generation", ylab = "inbreeding level", main = "Inbreeding",
     ylim = c(min(inbreeding), max(inbreeding)))
lines(inbreeding[2,],lwd = 2,  col = "red")
lines(inbreeding[3,], lwd = 2, col = "green")

legend("topleft", c("Baseline", "40% Culling", "Avoid Half/Fullsib selection & mating"),
       col = c("black", "red", "green"), lty = 1, lwd = 3)



#################################################################################################
####### END of the Task - everthing below is just showing of some additional MoBPS functionality
#################################################################################################


# addressing a specific individiual (hint for the next task)
sort(get.bv(population, gen = 11), index.return = TRUE)
colnames(get.bv(population, gen = 11))[88]

which(max(get.bv(population, gen = 11)) == get.bv(population, gen = 11))
get.database(population, id = 1088)



breeding.diploid(population,
                 selection.criteria = "bv",
                 selection.size = c(10,0),
                 selection.m.database = cbind(11,1),
                 export.selected = TRUE)

breeding.diploid(population,
                 selection.criteria = "bv",
                 selection.size = c(10,0),
                 selection.m.database = cbind(11,1),
                 export.selected.database =  TRUE)


