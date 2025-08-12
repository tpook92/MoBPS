
## 1. Step: Generation Founder:

population = creating.diploid(nindi = c(250,250),
                 name.cohort = c("Bulls", "Cows"),
                 map = MoBPSmaps::map_cattle1)

summary(population)

population = creating.trait(population,
                            n.additive = 1000,
                            trait.name = "MilkYield",
                            mean.target = 100,
                            var.target = 1)

summary(population)

## 2. Breeding Program

# Phenotype Cows

population = breeding.diploid(population,
                              phenotyping.cohorts = "Cows",
                              heritability = 0.3)

# Estimate Breeding values

population = breeding.diploid(population,
                              genotyped.cohorts = c("Bulls", "Cows"))

population = breeding.diploid(population,
                              bve = TRUE,
                              bve.cohorts = c("Bulls", "Cows"))

get.bve(population, cohorts = "Bulls")
analyze.bv(population, cohorts = c("Cows"))
analyze.bv(population, cohorts = c("Bulls"))

# Select bulls
# Generate the next generation

population = breeding.diploid(population, 
                              selection.size = c(10,250),
                              selection.criteria = "bve",
                              selection.m.cohorts = "Bulls",
                              selection.f.cohorts = "Cows",
                              breeding.size = 500,
                              name.cohort = c("BullsNext",
                                              "CowsNext"))