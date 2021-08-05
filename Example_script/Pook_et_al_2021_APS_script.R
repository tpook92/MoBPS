
## This program was simulated using
## MoBPS v.1.6.19

args <- commandArgs(TRUE)
nr <- as.numeric(args[1])
scenario <- as.numeric(args[2])
ogc <- as.logical(args[3])
ogc_target <- as.numeric(args[4])

set.seed(nr)

library(MoBPS)


sel <- c(50,500)
share_geno <- 0
heritability <- c(0.3,0.3,0.1)
weights <- c(1,1,1)
rel_matrix <- "vanRaden"
single.step <- TRUE
selection.criteria <- c("bve", "bve")
bve <- TRUE


if(scenario ==1){
  single.step <- FALSE
  selection.criteria <- c("pheno", "pheno")
  bve <- FALSE
}
if(scenario == 2){
  rel_matrix <- "kinship"
  single.step <- FALSE
}

if(scenario == 3){
  share_geno <- 0.1
}

if(scenario == 4){
  share_geno <- 0.25
}

if(scenario == 5){
  share_geno <- 1
}

if(ogc)
{
  ogc = TRUE
  sel = c(2500,500)
}

population <- creating.diploid(vcf= "GMP_filt_imp_vcf.vcf.gz",
                               sex.s = c(rep(1,17), rep(2,150), rep(1,24)),
                               n.additive = c(1000,1000,1000,1000),
                               vcf.maxsnp = 50000, mean.target = 100,
                               var.target = c(100,100/133,100,100))

trafo <- function(x){
  if(x < qnorm(0.02, mean=100,sd=sqrt(1000/3/133))){
    y <- 9
  } else   if(x < qnorm(0.42, mean=100,sd=sqrt(1000/3/133))){
    y <- 10
  } else   if(x < qnorm(0.62, mean=100,sd=sqrt(1000/3/133))){
    y <- 11
  } else   if(x < qnorm(0.82, mean=100,sd=sqrt(1000/3/133))){
    y <- 12
  } else   if(x < qnorm(0.92, mean=100,sd=sqrt(1000/3/133))){
    y <- 13
  } else   if(x < qnorm(0.97, mean=100,sd=sqrt(1000/3/133))){
    y <- 14
  } else   if(x < qnorm(0.99, mean=100,sd=sqrt(1000/3/133))){
    y <- 15
  } else  if(x < qnorm(0.9999, mean=100,sd=sqrt(1000/3/133))){
    y <- 16
  } else if(x < qnorm(0.99999, mean=100,sd=sqrt(1000/3/133))){
    y <- 17
  } else  {
    y <- 18
  }
}


population <- creating.phenotypic.transform(population, trait = 2, phenotypic.transform.function = trafo)

population <- add.founder.kinship(population)


# Split Dalmose / Relihausen

population <- breeding.diploid(population, repeat.mating = matrix(c(0.05, 2,
                                                                    0.10, 3,
                                                                    0.10, 4,
                                                                    0.20, 5,
                                                                    0.25, 6,
                                                                    0.15, 7,
                                                                    0.10, 8,
                                                                    0.05, 9), ncol=2, byrow=TRUE)) ## Litter size

ngen <- length(population$breeding)
population <- breeding.diploid(population, breeding.size = 5000, add.gen= ngen+1,
                               selection.m.database=cbind(1,2),
                               selection.f.database=cbind(1,2),
                               display.progress=FALSE)

for(index in 1:3){
  print(index)
  ngen <- length(population$breeding)
  population <- breeding.diploid(population, breeding.size = 5000, selection.size = c(50,500),
                                 add.gen= ngen+1,
                                 share.genotyped = 0, display.progress=FALSE)

}



for(index in 1:10){
  print(index)
  ngen <- length(population$breeding)

  # BVE
  population <- breeding.diploid(population,
                                 phenotyping.gen = ngen,
                                 heritability = heritability,
                                 sigma.e.gen = 2,
                                 n.observation = c(1,1,1,0),
                                 display.progress=FALSE)
  population <- breeding.diploid(population,
                                 bve = bve,
                                 bve.gen = c(length(population$breeding)-1, length(population$breeding)),
                                 bve.insert.gen = length(population$breeding),
                                 singlestep.active = single.step,
                                 relationship.matrix = rel_matrix,
                                 display.progress=FALSE)


  population <- breeding.diploid(population, breeding.size = 5000, selection.size = sel,
                                 add.gen= ngen+1,
                                 selection.criteria = selection.criteria,
                                 share.genotyped = share_geno,
                                 multiple.bve.weights.m = weights,
                                 display.progress=FALSE,
                                 ogc=ogc,
                                 ogc.target = "max.BV",
                                 ogc.uniform = "female",
                                 ogc.ub.sKin = ogc_target,
                                 relationship.matrix.ogc = "vanRaden"
                                 )

}



BVE1 <- matrix(0, nrow=4, ncol=11)
acc1 <- matrix(0, nrow=4, ncol=11)
kin1 <- matrix(0, nrow=2, ncol=11)
nsel <- matrix(0, nrow=2, ncol=11)
t2_list <- list()

for(index in 1:11){
  BVE1[,index] <- rowMeans(get.bv(population, gen=index+4))
  kin1[,index] <- kinship.emp.fast(population=population, database = cbind(index+4,1))
  if(scenario ==1){
    acc1[,index] <- analyze.bv(population=population, database = cbind(index+4,1))[[1]][2,]
  } else{
    acc1[,index] <- analyze.bv(population=population, database = cbind(index+4,1))[[1]][1,]
  }

  a <- get.pedigree(population = population, gen = index+4)
  nsel[,index] <- c(length(unique(a[,2])), length(unique(a[,3])))
  t2_list[[index]] <- get.pheno(population, gen=index+4)[2,]
}



save(file=paste0("Minipig/Sim", nr, "scenario", scenario, "ogc", ogc, "target", ogc_target,"_4.RData"), list=c("population"))
save(file=paste0("Minipig/Sim", nr, "scenario", scenario, "ogc", ogc, "target", ogc_target,"_4_result.RData"), list=c("BVE1", "kin1", "acc1", "nsel", "t2_list"))

