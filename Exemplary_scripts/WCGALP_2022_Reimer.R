seed <- 1

library(MoBPS)
library(MoBPSmaps)

set.seed(seed)
map <- MoBPSmaps::map_pig1


# Sample the position of
del <- sample(1:nrow(map), 1000)

# manually set allele frequency for the deleterious variants
freq <- rnorm(1000, mean=0.01, sd=0.005)
freq[freq<0] <- 0
map[,5] <- runif(nrow(map), 0,1)
map[del,5] <- freq
population <- creating.diploid(nindi = 1000, map=map)


# Generation of the performance trait
population <- creating.trait(population, n.additive = 1000, trait.name = "6 Month Weight")

# Generation of a trait to model deleterious variants
real.bv.add <- cbind(rep(0,1000),0,0,0,-100)
for(index in 1:length(del)){
  chr <- sum(del[index] > c(0,population$info$cumsnp))
  snp <- del[index] - c(0,population$info$cumsnp)[chr]
  real.bv.add[index,1:2] <- c(snp,chr)
}

population <- creating.trait(population, real.bv.add = real.bv.add, trait.name = "Deleterious")

population <- breeding.diploid(population)

# Check which individuals carry at least one homozygous lethal
# Remove all animals that have on homozygous lethal

population <- breeding.diploid(population, culling.gen = 1, culling.bv1 = 0, culling.bv2 = 100, culling.share1 = 1, culling.share2 = 0, culling.index = c(0,1))

## 15 Generations of moderate selection on a single breeding nucleus
alive_stats <- NULL # number of animals without lethal variant
for(index in 1:15){
  population <- breeding.diploid(population, heritability = c(0.3, 0), phenotyping.gen = length(population$breeding))
  population <- breeding.diploid(population, bve=TRUE, bve.gen = length(population$breeding), relationship.matrix = "kinship")
  population <- breeding.diploid(population, breeding.size = 1000, selection.size = c(167, 167), multiple.bve.weights.m = c(1,0), display.progress = FALSE)
  population <- breeding.diploid(population, culling.gen = length(population$breeding), culling.bv1 = 0, culling.bv2 = 100, culling.share1 = 1, culling.share2 = 0, culling.index = c(0,1))
  alive_stats <- c(alive_stats, sum(get.class(population, gen=length(population$breeding))==0))
}

# Simulation of three breeding stocks of different size
alive_stats2 <- NULL # number of animals without lethal variant
total_stats2 <- NULL # total number of generated animals
for(index in 1:20){
  population <- breeding.diploid(population, heritability = c(0.3, 0), phenotyping.gen = length(population$breeding))

  population <- breeding.diploid(population, bve=TRUE, bve.gen = length(population$breeding), relationship.matrix = "kinship")

  add.gen <- length(population$breeding)+1

  # Stock 1 (Large)
  if(index==1 || (alive_stats2[index-1,1]>0 & alive_stats2[index-1,5]>0)){
    population <- breeding.diploid(population, breeding.size = c(333,667), selection.size = c(111, 222), multiple.bve.weights.m = c(1,0),
                                   class.m = if(index==1){0} else{1}, class.f = if(index==1){0} else{1}, new.class = 1,
                                   selection.m.database = cbind(add.gen-1,1),
                                   selection.f.database = cbind(add.gen-1,1),
                                   add.gen = add.gen, display.progress = FALSE)
  }



  # Stock 2 (Medium)
  if(index ==1){
    size <- c(264, 546)
  } else if(index<=5){
    size <- c(132, 268)
  } else if(index<=10){
    size <- c(66, 134)
  } else{
    size <- c(33, 67)
  }
  if(index==1 || (alive_stats2[index-1,2]>0 & alive_stats2[index-1,6]>0)){
    population <- breeding.diploid(population, breeding.size = size, selection.size = ceiling(size/3), multiple.bve.weights.m = c(1,0),
                                   class.m = if(index==1){0} else{2}, class.f = if(index==1){0} else{2}, new.class = 2,
                                   selection.m.database = cbind(add.gen-1,1),
                                   selection.f.database = cbind(add.gen-1,1),
                                   add.gen = add.gen, display.progress = FALSE)
  }

  # Stock 3 (Endangered)
  if(index ==1){
    size <- c(233, 467)
  } else if(index<=5){
    size <- c(66, 134)
  } else if(index<=10){
    size <- c(26, 54)
  } else{
    size <- c(12, 24)
  }

  if(index==1 || (alive_stats2[index-1,3]>0 & alive_stats2[index-1,7]>0)){
    population <- breeding.diploid(population, breeding.size = size, selection.size = ceiling(size/3), multiple.bve.weights.m = c(1,0),
                                   class.m = if(index==1){0} else{3}, class.f = if(index==1){0} else{3}, new.class = 3,
                                   selection.m.database = cbind(add.gen-1,1),
                                   selection.f.database = cbind(add.gen-1,1),
                                   add.gen = add.gen, display.progress = FALSE)
  }


  # Stock 4 (Critical)
  if(index ==1){
    size <- c(100, 200)
  } else if(index<=5){
    size <- c(25, 50)
  } else if(index<=10){
    size <- c(13, 27)
  } else{
    size <- c(5, 10)
  }

  if(index==1 || (alive_stats2[index-1,4]>0 & alive_stats2[index-1,8]>0)){
    population <- breeding.diploid(population, breeding.size = size, selection.size = ceiling(size/3), multiple.bve.weights.m = c(1,0),
                                   class.m = if(index==1){0} else{4}, class.f = if(index==1){0} else{4}, new.class = 4,
                                   selection.m.database = cbind(add.gen-1,1),
                                   selection.f.database = cbind(add.gen-1,1),
                                   add.gen = add.gen, display.progress = FALSE)
  }


  total_stats2 <- rbind(total_stats2, c(sum(get.class(population, database=cbind(add.gen,1))==1), sum(get.class(population, database=cbind(add.gen,1))==2),
                                        sum(get.class(population, database=cbind(add.gen,1))==3),sum(get.class(population, database=cbind(add.gen,1))==4),
                                        sum(get.class(population, database=cbind(add.gen,2))==1),sum(get.class(population, database=cbind(add.gen,2))==2),
                                        sum(get.class(population, database=cbind(add.gen,2))==3),sum(get.class(population, database=cbind(add.gen,2))==4)))
  population <- breeding.diploid(population, culling.gen = length(population$breeding), culling.bv1 = 0, culling.bv2 = 100, culling.share1 = 1, culling.share2 = 0, culling.index = c(0,1))
    alive_stats2 <- rbind(alive_stats2, c(sum(get.class(population, database=cbind(add.gen,1))==1), sum(get.class(population, database=cbind(add.gen,1))==2),
                                        sum(get.class(population, database=cbind(add.gen,1))==3),sum(get.class(population, database=cbind(add.gen,1))==4),
                                        sum(get.class(population, database=cbind(add.gen,2))==1),sum(get.class(population, database=cbind(add.gen,2))==2),
                                        sum(get.class(population, database=cbind(add.gen,2))==3),sum(get.class(population, database=cbind(add.gen,2))==4)))
}

# Take 2 males from the four different stocks + 10 generation old animals from stock 3/4 as a cryo reserve and mate to current females
add.gen <- length(population$breeding)+1

if(alive_stats2[20,1]>0 & alive_stats2[20,7]){
  population <- breeding.diploid(population, breeding.size = 36, selection.size=c(2,24),
                                 selection.m.database = cbind(add.gen-1,1),
                                 class.m = 1,
                                 selection.f.database = cbind(add.gen-1,2),
                                 class.f = 3,
                                 new.class = 1,
                                 add.gen = add.gen, display.progress = FALSE)
}

if(alive_stats2[20,2]>0 & alive_stats2[20,7]){
  population <- breeding.diploid(population, breeding.size = 36, selection.size=c(2,24),
                                 selection.m.database = cbind(add.gen-1,1),
                                 class.m = 2,
                                 selection.f.database = cbind(add.gen-1,2),
                                 class.f = 3,
                                 new.class = 2,
                                 add.gen = add.gen, display.progress = FALSE)
}

if(alive_stats2[20,3]>0 & alive_stats2[20,7]){
  population <- breeding.diploid(population, breeding.size = 36, selection.size=c(2,24),
                                 selection.m.database = cbind(add.gen-1,1),
                                 class.m = 3,
                                 selection.f.database = cbind(add.gen-1,2),
                                 class.f = 3,
                                 new.class = 3,
                                 add.gen = add.gen, display.progress = FALSE)
}

if(alive_stats2[10,3]>0 & alive_stats2[20,7]){
  population <- breeding.diploid(population, breeding.size = 36, selection.size=c(2,24),
                                 selection.m.database = cbind(add.gen-11,1),
                                 class.m = 3,
                                 selection.f.database = cbind(add.gen-1,2),
                                 class.f = 3,
                                 new.class = 4,
                                 add.gen = add.gen, display.progress = FALSE)
}


if(alive_stats2[20,1]>0 & alive_stats2[20,8]){
  population <- breeding.diploid(population, breeding.size = 36, selection.size=c(2,10),
                                 selection.m.database = cbind(add.gen-1,1),
                                 class.m = 1,
                                 selection.f.database = cbind(add.gen-1,2),
                                 class.f = 4,
                                 new.class = 5,
                                 add.gen = add.gen, display.progress = FALSE)
}

if(alive_stats2[20,3]>0 & alive_stats2[20,8]){
  population <- breeding.diploid(population, breeding.size = 36, selection.size=c(2,10),
                                 selection.m.database = cbind(add.gen-1,1),
                                 class.m = 3,
                                 selection.f.database = cbind(add.gen-1,2),
                                 class.f = 4,
                                 new.class = 6,
                                 add.gen = add.gen, display.progress = FALSE)
}

if(alive_stats2[20,4]>0 & alive_stats2[20,8]){
  population <- breeding.diploid(population, breeding.size = 36, selection.size=c(2,10),
                                 selection.m.database = cbind(add.gen-1,1),
                                 class.m = 4,
                                 selection.f.database = cbind(add.gen-1,2),
                                 class.f = 4,
                                 new.class = 7,
                                 add.gen = add.gen, display.progress = FALSE)
}

if(alive_stats2[10,4]>0 & alive_stats2[20,8]){
  population <- breeding.diploid(population, breeding.size = 36, selection.size=c(2,10),
                                 selection.m.database = cbind(add.gen-11,1),
                                 class.m = 4,
                                 selection.f.database = cbind(add.gen-1,2),
                                 class.f = 4,
                                 new.class = 8,
                                 add.gen = add.gen, display.progress = FALSE)
}

total_stats3 <- NULL
alive_stats3 <- NULL

if(nrow(population$info$size)<add.gen){
  total_stats3 <- alive_stats3 <- matrix(0, nrow=30, ncol=16)
} else{
  total_stats3 <- rbind(total_stats3, c(sum(get.class(population, database=cbind(add.gen,1))==1), sum(get.class(population, database=cbind(add.gen,1))==2),
                                        sum(get.class(population, database=cbind(add.gen,1))==3),sum(get.class(population, database=cbind(add.gen,1))==4),
                                        sum(get.class(population, database=cbind(add.gen,1))==5),sum(get.class(population, database=cbind(add.gen,1))==6),
                                        sum(get.class(population, database=cbind(add.gen,1))==7),sum(get.class(population, database=cbind(add.gen,1))==8),
                                        sum(get.class(population, database=cbind(add.gen,2))==1),sum(get.class(population, database=cbind(add.gen,2))==2),
                                        sum(get.class(population, database=cbind(add.gen,2))==3),sum(get.class(population, database=cbind(add.gen,2))==4),
                                        sum(get.class(population, database=cbind(add.gen,2))==5),sum(get.class(population, database=cbind(add.gen,2))==6),
                                        sum(get.class(population, database=cbind(add.gen,2))==7),sum(get.class(population, database=cbind(add.gen,2))==8)))
population <- breeding.diploid(population, culling.gen = length(population$breeding), culling.bv1 = 0, culling.bv2 = 100, culling.share1 = 1, culling.share2 = 0, culling.index = c(0,1))

alive_stats3 <- rbind(alive_stats3, c(sum(get.class(population, database=cbind(add.gen,1))==1), sum(get.class(population, database=cbind(add.gen,1))==2),
                                      sum(get.class(population, database=cbind(add.gen,1))==3),sum(get.class(population, database=cbind(add.gen,1))==4),
                                      sum(get.class(population, database=cbind(add.gen,1))==5),sum(get.class(population, database=cbind(add.gen,1))==6),
                                      sum(get.class(population, database=cbind(add.gen,1))==7),sum(get.class(population, database=cbind(add.gen,1))==8),
                                      sum(get.class(population, database=cbind(add.gen,2))==1),sum(get.class(population, database=cbind(add.gen,2))==2),
                                      sum(get.class(population, database=cbind(add.gen,2))==3),sum(get.class(population, database=cbind(add.gen,2))==4),
                                      sum(get.class(population, database=cbind(add.gen,2))==5),sum(get.class(population, database=cbind(add.gen,2))==6),
                                      sum(get.class(population, database=cbind(add.gen,2))==7),sum(get.class(population, database=cbind(add.gen,2))==8)))
}


# Random mating of the newly generation 4 breeding nucleus
for(index in 1:30){
  print(index)
  add.gen <- length(population$breeding)+1
  for(group in 1:8){
    if(alive_stats3[index,group]>0 & alive_stats3[index, group+8]){
      population <- breeding.diploid(population, breeding.size = 36, selection.size = c(18,18),
                                     new.class = group, class.m = group, class.f=group,
                                     selection.m.database = cbind(add.gen-1,1),
                                     selection.f.database = cbind(add.gen-1,2),
                                     add.gen= add.gen, display.progress = FALSE)
    }
  }
  print(!(nrow(population$info$size)<add.gen))
  if(!(nrow(population$info$size)<add.gen)){
    total_stats3 <- rbind(total_stats3, c(sum(get.class(population, database=cbind(add.gen,1))==1), sum(get.class(population, database=cbind(add.gen,1))==2),
                                          sum(get.class(population, database=cbind(add.gen,1))==3),sum(get.class(population, database=cbind(add.gen,1))==4),
                                          sum(get.class(population, database=cbind(add.gen,1))==5),sum(get.class(population, database=cbind(add.gen,1))==6),
                                          sum(get.class(population, database=cbind(add.gen,1))==7),sum(get.class(population, database=cbind(add.gen,1))==8),
                                          sum(get.class(population, database=cbind(add.gen,2))==1),sum(get.class(population, database=cbind(add.gen,2))==2),
                                          sum(get.class(population, database=cbind(add.gen,2))==3),sum(get.class(population, database=cbind(add.gen,2))==4),
                                          sum(get.class(population, database=cbind(add.gen,2))==5),sum(get.class(population, database=cbind(add.gen,2))==6),
                                          sum(get.class(population, database=cbind(add.gen,2))==7),sum(get.class(population, database=cbind(add.gen,2))==8)))
    population <- breeding.diploid(population, culling.gen = length(population$breeding), culling.bv1 = 0, culling.bv2 = 100, culling.share1 = 1, culling.share2 = 0, culling.index = c(0,1))
    alive_stats3 <- rbind(alive_stats3, c(sum(get.class(population, database=cbind(add.gen,1))==1), sum(get.class(population, database=cbind(add.gen,1))==2),
                                          sum(get.class(population, database=cbind(add.gen,1))==3),sum(get.class(population, database=cbind(add.gen,1))==4),
                                          sum(get.class(population, database=cbind(add.gen,1))==5),sum(get.class(population, database=cbind(add.gen,1))==6),
                                          sum(get.class(population, database=cbind(add.gen,1))==7),sum(get.class(population, database=cbind(add.gen,1))==8),
                                          sum(get.class(population, database=cbind(add.gen,2))==1),sum(get.class(population, database=cbind(add.gen,2))==2),
                                          sum(get.class(population, database=cbind(add.gen,2))==3),sum(get.class(population, database=cbind(add.gen,2))==4),
                                          sum(get.class(population, database=cbind(add.gen,2))==5),sum(get.class(population, database=cbind(add.gen,2))==6),
                                          sum(get.class(population, database=cbind(add.gen,2))==7),sum(get.class(population, database=cbind(add.gen,2))==8)))
  }
}



## Analysis: allele frequencies of deleterious variants

p_del_phase0 <- NULL
p_purged_phase0 <- NULL

for(index in 1:16){
  print(index)
  geno <- get.geno(population, gen= index)
  class <- get.class(population, gen= index)
  p_del1 <- rowMeans(geno[del,, drop=FALSE])
  p_purged1 <- rowMeans(geno[del,, drop=FALSE]>0)
  p_del_phase0 <- cbind(p_del_phase0, p_del1)
  p_purged_phase0 <- cbind(p_purged_phase0, p_purged1)
}

p_del <- NULL
p_purged_phase1 <- NULL
for(index in 1:20){
  print(index)
  geno <- get.geno(population, gen= index +16)
  class <- get.class(population, gen= index +16)
  p_del1 <- cbind(rowMeans(geno[del, class==1, drop=FALSE]), rowMeans(geno[del, class==2, drop=FALSE]), rowMeans(geno[del, class==3, drop=FALSE]), rowMeans(geno[del, class==4, drop=FALSE]))
  p_purged1 <- cbind(rowMeans(geno[del, class==1, drop=FALSE]>0), rowMeans(geno[del, class==2, drop=FALSE]>0), rowMeans(geno[del, class==3, drop=FALSE]>0), rowMeans(geno[del, class==4, drop=FALSE]>0))
  p_del <- cbind(p_del, p_del1)
  p_purged_phase1 <- cbind(p_purged_phase1, p_purged1)
}

p_del_phase1 <- p_del

p_del_phase2 <- NULL
p_purged_phase2 <- NULL


for(index in 1:31){
  print(index)
  if(length(population$breeding)>=index+36){
    geno <- get.geno(population, gen= index +36)
    class <- get.class(population, gen= index +36)
    p_del1 <- cbind(rowMeans(geno[del, class==1, drop=FALSE]), rowMeans(geno[del, class==2, drop=FALSE]), rowMeans(geno[del, class==3, drop=FALSE]), rowMeans(geno[del, class==4, drop=FALSE]),
                    rowMeans(geno[del, class==5, drop=FALSE]), rowMeans(geno[del, class==6, drop=FALSE]), rowMeans(geno[del, class==7, drop=FALSE]), rowMeans(geno[del, class==8, drop=FALSE]))

    p_purged1 <- cbind(rowMeans(geno[del, class==1, drop=FALSE]>0), rowMeans(geno[del, class==2, drop=FALSE]>0), rowMeans(geno[del, class==3, drop=FALSE]>0), rowMeans(geno[del, class==4, drop=FALSE]>0),
                       rowMeans(geno[del, class==5, drop=FALSE]>0), rowMeans(geno[del, class==6, drop=FALSE]>0), rowMeans(geno[del, class==7, drop=FALSE]>0), rowMeans(geno[del, class==8, drop=FALSE]>0))

  } else{
    p_del1 <- matrix(NA, nrow=1000, ncol=8)
    p_purged1 <- matrix(NA, nrow=1000, ncol=8)
  }

  p_del_phase2 <- cbind(p_del_phase2, p_del1)
  p_purged_phase2 <- cbind(p_purged_phase2, p_purged1)
}

# Analysis: Inbreeding and kinship
cohorts_ex <- get.cohorts(population, extended=TRUE)

cohorts_ex <- cohorts_ex[cohorts_ex[,4]>0,]

kin <- matrix(0, nrow=nrow(cohorts_ex), ncol=2)
for(index in 1:nrow(cohorts_ex)){
  kin[index,] <- kinship.emp.fast(population = population, cohorts = cohorts_ex[index,1], ibd.obs = 100, hbd.ob=50)
}

kin <- cbind(as.numeric(cohorts_ex[,2]), as.numeric(cohorts_ex[,5]), kin)
kin[,4] <- 2*(kin[,4]-0.5)
purged <-rbind(colMeans(p_del==0), colMeans(p_del<0.05 & p_del>0), colMeans(p_del<0.2 & p_del>=0.05), colMeans( p_del>=0.2))
purged_phase0 <-rbind(colMeans(p_del_phase0==0), colMeans(p_del_phase0<0.05 & p_del_phase0>0), colMeans(p_del_phase0<0.2 & p_del_phase0>=0.05), colMeans( p_del_phase0>=0.2))
purged_phase1 <-rbind(colMeans(p_del_phase1==0), colMeans(p_del_phase1<0.05 & p_del_phase1>0), colMeans(p_del_phase1<0.2 & p_del_phase1>=0.05), colMeans( p_del_phase1>=0.2))
purged_phase2 <-rbind(colMeans(p_del_phase2==0), colMeans(p_del_phase2<0.05 & p_del_phase2>0), colMeans(p_del_phase2<0.2 & p_del_phase2>=0.05), colMeans( p_del_phase2>=0.2))

save(file=paste0("wcgalp_reimer2/reimer", seed, ".RData"), list=c("alive_stats3", "alive_stats2", "alive_stats","p_del", "purged", "kin", "p_del_phase0", "p_del_phase1", "p_del_phase2",
                                                                 "purged_phase0", "purged_phase1", "purged_phase2", "total_stats2", "total_stats3",
                                                                 "p_purged_phase0", "p_purged_phase1", "p_purged_phase2"))

save(file=paste0("wcgalp_reimer2/reimer", seed, "_pop.RData"), list=c("alive_stats3", "alive_stats2", "alive_stats","p_del", "purged", "kin", "p_del_phase0", "p_del_phase1", "p_del_phase2",
                                                                    "purged_phase0", "purged_phase1", "purged_phase2", "total_stats2", "total_stats3",
                                                                    "p_purged_phase0", "p_purged_phase1", "p_purged_phase2", "population"))




