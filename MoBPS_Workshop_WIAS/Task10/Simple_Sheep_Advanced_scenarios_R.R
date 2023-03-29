# Initialize variables to run different scenarios and runs of the simulation script
# Here this is done in a for loop
# In practise you would typically use a server

# Practical R tip:
# If you enter following command in your console:
# Rscript Script_you_want_to_run.R Possible Input Variable 123
# This will start an R session and run the selected script

# Typing "args <- commandArgs(TRUE)" in your Rscript will initialize
# a variable named args which contains c("Possible", "Input", "Variable", "123")
# "123" is no numeric number but a character string!

# In this example this could be:
# Rscript Simple_Sheep_Advanced_scenario.R 1 2
# To perform run 1; scenario 2

#args <- commandArgs(TRUE)
#run <- as.numeric(args[1])
#scenario <- as.numeric(args[2])

####################################
########### WARNING ################
####### This code takes a while! ###
####################################

# These are 40 simulation of the breeding program!

for(run in 1:10){
  for(scenario in 1:4){

    # Add additional parameter to the code at all stages that we want to make
    # changes to the breeding program.
    # The else statement provides the setting in our baseline.

    if(scenario==2){
      selected_rams <- 5
    } else{
      selected_rams <- 10
    }

    if(scenario==3){
      share_genotyped_female <- 0.5
    } else{
      share_genotyped_female <- 1
    }

    if(scenario==4){
      sel_index <- c(3,1)
    } else{
      sel_index <- c(1,1)
    }

    set.seed(run)

    # Simulation of the LD build up
    dataset <- founder.simulation(nindi=100, map = MoBPSmaps::map_sheep2,
                                  nfinal = 60+120, n.gen = 5)

    # Use a genomic map from the MoBPSmap R-package
    # Import genomic data from the dataset-matrix ((each individual has 2 haplotypes))
    population <- creating.diploid(dataset = dataset[,1:100], nindi = 50, sex.quota=0, map=MoBPSmaps::map_sheep2,
                                   name.cohort = "1yearRams_0")


    # Make sure female founders are only genotyped with a probability of 50%

    population <- creating.diploid(dataset = dataset[,101:(100+selected_rams*2)], population = population, nindi = selected_rams,
                                   sex.quota=0, name.cohort = "2yearRams_0")

    population <- creating.diploid(dataset = dataset[,121:220], population = population, nindi = 50,
                                   share.genotyped = share_genotyped_female, sex.quota=1, name.cohort = "1yearEwes_0")

    population <- creating.diploid(dataset = dataset[,221:300], population = population, nindi = 40,
                                   share.genotyped = share_genotyped_female, sex.quota=1, name.cohort = "2yearEwes_0")

    population <- creating.diploid(dataset = dataset[,301:360], population = population, nindi = 30,
                                   share.genotyped = share_genotyped_female, sex.quota=1, name.cohort = "3yearEwes_0")

    summary(population)

    # Use of two correlated traits

    set.seed(run)

    population <- creating.trait(population, n.additive = c(1000,1000), mean.target = 100, var.target = c(30,20),
                                 trait.cor = matrix(c(1, 0.2, 0.2,1), nrow=2))

    # first column: Size of the litter
    # second column: probability of the litter size to occur

    litter_size <- matrix(c(2, 0.5,
                            3, 0.3,
                            4, 0.2), ncol=2, byrow=TRUE)


    for(index in 1:20){

      # phenotypes for the second trait are only generated for selected cohorts

      population <- breeding.diploid(population, phenotyping.gen = index,
                                     heritability = c(0.3,0.2), n.observation = c(1,0))

      population <- breeding.diploid(population, phenotyping.cohorts = paste0(c("2yearRams_", "2yearEwes_", "3yearEwes_"), index-1),
                                     heritability = c(0.3,0.2), n.observation = c(0,1))


      # As some individuals are non genotyped a genomic breeding value estimation does not work
      # Either use a pedigree-based evaluation or single-step
      # Otherwise non-genotyped individuals will automatically be removed for the estimation

      population <- breeding.diploid(population, bve = TRUE, bve.gen = index,
                                     bve.cohorts = paste0(c("1yearRams_", "2yearRams_", "1yearEwes_", "2yearEwes_", "3yearEwes_"), index-1))


      # Generation of new individuals via reproduction

      population <- breeding.diploid(population, breeding.size = c(50,50),
                                     selection.m.cohorts = paste0("2yearRams_", index-1),
                                     selection.f.cohorts = c(paste0("2yearEwes_",index-1), paste0("3yearEwes_",index-1)),
                                     name.cohort = paste0("NewOffspring_",index),
                                     repeat.mating = litter_size,
                                     share.genotyped = 0
      )

      population <- breeding.diploid(population, breeding.size = c(50,0),
                                     selection.size = c(50,0),
                                     copy.individual.m = TRUE,
                                     selection.m.cohorts = paste0("NewOffspring_",index, "_M"),
                                     name.cohort = paste0("1yearRams_", index),
                                     add.gen=index+1,
                                     added.genotyped = 1
      )

      population <- breeding.diploid(population, breeding.size = c(0,50),
                                     selection.size = c(0,50),
                                     copy.individual.f = TRUE,
                                     selection.f.cohorts = paste0("NewOffspring_",index, "_F"),
                                     name.cohort = paste0("1yearEwes_", index),
                                     add.gen=index+1,
                                     added.genotyped = share_genotyped_female
      )


      # Generation of new individuals via selection

      population <- breeding.diploid(population, breeding.size = c(selected_rams,0), selection.size = c(selected_rams,0),
                                     selection.m.cohorts = paste0("1yearRams_", index-1),
                                     copy.individual.m = TRUE,
                                     selection.criteria = "bve",
                                     multiple.bve.weights.m = sel_index,
                                     add.gen=index+1,
                                     name.cohort = paste0("2yearRams_", index))


      population <- breeding.diploid(population, breeding.size = c(0,40),
                                     selection.size = c(0,40),
                                     selection.f.cohorts = paste0("1yearEwes_", index-1),
                                     copy.individual.f = TRUE,
                                     selection.criteria = "random",
                                     add.gen=index+1,
                                     name.cohort = paste0("2yearEwes_", index))

      population <- breeding.diploid(population, breeding.size = c(0,30),
                                     selection.size = c(0,30),
                                     selection.f.cohorts = paste0("2yearEwes_", index-1),
                                     copy.individual.f = TRUE,
                                     selection.criteria = "random",
                                     add.gen=index+1,
                                     name.cohort = paste0("3yearEwes_", index))

    }

    # Analysis of the outcomes

    cohorts <- get.cohorts(population)

    genomic_values <- accuracies <- kinships <- matrix(0, nrow=2, ncol=length(cohorts))

    for(index in 1:length(cohorts)){
      # This extracts the genomic values for animals from the cohort and calculates the mean
      genomic_values[,index] <- rowMeans(get.bv(population, cohorts=cohorts[[index]]))
      # This extracts accuracies of the breeding value estimation for selected cohorts
      accuracies[,index] <- analyze.bv(population, cohorts=cohorts[index])[[1]][1,]
      # This approximates avg. kinship between animals and within ((kinship -0.5) *2 is inbreeding)
      kinships[,index] <- kinship.emp.fast(population=population, cohorts = cohorts[index])

    }

    save(file=paste0("Sheep_simulation_scenario", scenario, "run", run,".RData"), list=c("cohorts", "genomic_values", "accuracies", "kinships"))
  }
}


