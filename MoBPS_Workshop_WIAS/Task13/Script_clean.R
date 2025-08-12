###########################################################
################# WARNING!!! ##############################
###########################################################
#### This is a simulation script that contains errors! ####
###########################################################

# write a log-file
if(FALSE){
  zz <- file("C:/Users/pook001/OneDrive - Wageningen University & Research/MoBPS_Workshop_WIAS_2025/Task13/log_file.log", open="wt")
  sink(zz, append = TRUE, type = c("output"), split=TRUE)
}


library(MoBPS)
set.seed(40)

population = creating.diploid(nsnp = 25000, nindi = 100,
                              chr.nr = 10, progress.bar = FALSE)

economic_weight = c(0.30, 50) # 30cent per litre, 50 euro per F%

# generate a two trait
# first with epistasis, second purely additive

population = creating.trait(population,
                            n.additive = c(1000,1000),
                            n.quantitative = c(1000,0),
                            ## TP: This was leading to wrong trait architectures
                            # but c(1000,0) makes it easier to read
                            mean.target = c(9000,4),
                            var.target = c(400,0.5)^2,
                            trait.name = c("MilkYield",
                                           "F%"))

# Build up some LD
for(index in 1:3){
  suppressWarnings({
    population = breeding.diploid(population,
                                  breeding.size = 200,
                                  selection.size = c(20,100))
  })
  ## TP: in case you are sure something is doing what you want it
  # to be doing you can also suppress warnings being generated
  # First generation in the LD build-up has to little individuals
}

# two offspring per cow
population = breeding.diploid(population,
                              breeding.size = 200,
                              selection.size = c(20,100),
                              max.offspring = c(Inf,2),
                              ## TP: limit number of offspring per cow to 2
                              name.cohort = c("Bulls", "Cows"))

population = breeding.diploid(population,
                              phenotyping.cohorts = c("Cows"),
                              heritability = c(0.3,0.3))

## TP: no heritability was provided. Both traits were simulated with
# a residual variance of 100
# Bulls to not give milk!

population = breeding.diploid(population, bve = TRUE,
                              bve.cohorts = c("Bulls", "Cows"))

population = breeding.diploid(population, selection.size = c(10,100),
                              breeding.size = 200,
                              selection.m.cohorts = "Bulls",
                              selection.f.cohorts = "Cows",
                              ## TP: MoBPS on default takes last generation as sires/dams
                              # so this led to right outcomes but explicitely listing
                              # cohorts makes it easier to read
                              selection.index.weights.m = economic_weight,
                              selection.index.scale.m = "unit",
                              ## TP: index weights were provided per unit not relative weights
                              name.cohort = c("Bulls_next", "Cows_next"))

population = breeding.diploid(population,
                              phenotyping.cohorts = c("Bulls_next", "Cows_next"))
# TP: Bulls to not give milk!

# include both generations in the breeding value estimation
population = breeding.diploid(population, bve = TRUE,
                              bve.cohorts = c("Bulls", "Cows",
                                              "Bulls_next", "Cows_next"))

## TP: Breeding value estimation was written to a new population
# subsequent selection was on animals without estimated breeding values

population = breeding.diploid(population, selection.size = c(10,100),
                              selection.m.cohorts = "Bulls_next",
                              selection.f.cohorts = "Cows_next",
                              ## TP: Bulls were used as dams
                              breeding.size = 200,
                              selection.index.weights.m = economic_weight,
                              selection.index.scale.m = "unit",
                              name.cohort = c("Bulls_next2", "Cows_next2"))
## TP: new cohorts had the same name as the previous generation


# see how much progress we made - does not work?!
rowMeans(get.bv(population, cohorts =  "Bulls"))
rowMeans(get.bv(population, cohorts = "Bulls_next"))
rowMeans(get.bv(population, cohorts = "Bulls_next2"))

if(FALSE){
  sink(zz, append = TRUE, type = c("message"))
  warnings()
  sink(NULL)
  sink(NULL, type=c("message"))
}
