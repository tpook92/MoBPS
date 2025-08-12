
library(MoBPS)
set.seed(43)

population = creating.diploid(nsnp = 1000, nindi = c(20,80), n.additive = 100,
                              chr.nr = 5, chromosome.length = 3)

for(index in 1:10){
  population = breeding.diploid(population, breeding.size = c(20,80))
}

current_gen = get.ngen(population)
population = breeding.diploid(population, phenotyping.database = cbind(current_gen,2),
                              #, matrix(c(11,2), ncol = 2)
                              heritability = 0.3)

population = breeding.diploid(population, genotyped.database = cbind(current_gen,1))

# genotype the top 20 females from the last generation
population = breeding.diploid(population, selection.size = c(0,20),
                              selection.f.database = cbind(current_gen,2),
                              genotyped.selected = TRUE,
                              selection.criteria = "pheno")

get.database(population, gen = 3, cohorts = "Cohort_1_F",
             database = cbind(c(7,8,11),2))

# genotype the top 20 females from the last generation
population = breeding.diploid(population, selection.size = c(0,20),
                              selection.f.database = cbind(11,2),
                              genotyped.selected = TRUE,
                              selection.criteria = "pheno")

database_ordered = breeding.diploid(population, selection.size = c(20,80),
                                selection.criteria = c('bv', "pheno"),
                                export.selected = TRUE)

sire = database_ordered[[1]][3,1:3]
dams = database_ordered[[2]][c(3,7,17,42,79),1:3]

# generate a matrix with 6 colums:
# Generation Sire / Sex Sire / Number Sire / Generation Dam / Sex dam / Number dam
# Below is one way of doing this but R provides various other options
fixed.breeding = cbind(matrix(sire, nrow = 5, ncol = 3,
                              byrow = TRUE),
                       dams)

population = breeding.diploid(population, fixed.breeding = fixed.breeding)

# Extract pedigree and if an individual is genotyped
# raw pedigree has 9 columns:
# Individual (generation,sex,nr), sire (generation,sex,nr), dam (generation,sex,nr)
pedigree = get.pedigree(population, gen = current_gen, raw = TRUE)
is_genotyped = get.genotyped(population, gen = current_gen)


# number of genotyped offspring
tab1 = table(pedigree[is_genotyped,6])

# number of non-genotyped offspring
tab2 = table(pedigree[!is_genotyped,6])

# which sires have no genotyped offspring
to_genotype = which(!(names(tab2) %in% names(tab1)))

# for each sire without genotyped offspring add one offspring to the genotype_database
# that will later be the input to genotype more individuals

genotype_database = NULL
for(index in to_genotype){

  genotype_database = rbind(genotype_database,
                            pedigree[sample(which(pedigree[,6] == index),1), 1:3]
  )
}


# genotyping of additional animal
population = breeding.diploid(population, genotyped.database = genotype_database)

# How many individuals were now genotyped in total?
is_genotyped = get.genotyped(population, gen = current_gen)
sum(is_genotyped)

# number of genotyped offspring
table(pedigree[is_genotyped,6])

# number of non-genotyped offspring
table(pedigree[!is_genotyped,6])

#################################################################################################
####### END of the Task - everthing below is just showing of some additional MoBPS functionality
#################################################################################################

pedigree = get.pedigree(population, gen = current_gen)
exemplary_database = get.database(population, id = pedigree[1:5,2])
