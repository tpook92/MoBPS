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
                              chromosome.length = 1,
                              chr.nr = 29, progress.bar = FALSE)

economic_weight = c(0.30, 50) # 30cent per litre, 50 euro per F%

# generate a two trait
# first with epistasis, second purely additive
# both have heritability 0.3

population = creating.trait(population,
                            n.additive = c(1000,1000),
                            n.quantitative = c(1000,0),
                            mean.target = c(9000,4),
                            var.target = c(400,0.5)^2,
                            
                            trait.cor = matrix(c(1, -0.3,
                                                 -0.3, 1), ncol = 2),
                            
                            trait.name = c("MilkYield",
                                           "F%"))

# Build up some LD
for(index in 1:3){
  suppressWarnings({
    population = breeding.diploid(population,
                                  breeding.size = 200,
                                  selection.size = c(20,100))
  })

}

# two offspring per cow
population = breeding.diploid(population,
                              breeding.size = 200,
                              selection.size = c(20,100),
                              max.offspring = c(Inf,2),
                              name.cohort = c("Bulls", "Cows"))

population = breeding.diploid(population,
                              phenotyping.cohorts = c("Cows"),
                              heritability = c(0.3, 0.3))

population = breeding.diploid(population, bve = TRUE,
                              bve.cohorts = c("Bulls", "Cows"))

population = breeding.diploid(population, selection.size = c(10,10),
                              breeding.size = 200,
                              selection.index.scale.m = "unit",
                              selection.index.weights.m = economic_weight,
                              name.cohort = c("Bulls_next", "Cows_next"))

population = breeding.diploid(population,
                              phenotyping.cohorts = c("Cows_next"))

# include both generations in the breeding value estimation
population = breeding.diploid(population, bve = TRUE,
                              bve.cohorts = c("Bulls", "Cows",
                                              "Bulls_next", "Cows_next"))

population = breeding.diploid(population, selection.size = c(10,100),
                              selection.m.cohorts = "Bulls_next",
                              selection.f.cohorts = "Cows_next",
                              max.offspring = c(Inf,2),
                              breeding.size = 200,
                              selection.index.weights.m = economic_weight,
                              name.cohort = c("Bulls_next2", "Cows_next2"))


get.database(population, cohorts = "Bulls_next")

# see how much progress we made - does not work?!
rowMeans(get.bv(population, cohorts =  "Bulls"))
rowMeans(get.bv(population, cohorts = "Bulls_next"))

if(FALSE){
  sink(zz, append = TRUE, type = c("message"))
  warnings()
  sink(NULL)
  sink(NULL, type=c("message"))
}
