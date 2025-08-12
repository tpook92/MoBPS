# The RData files used here can be generated
# with the script from Task 11
# For each simulation, this includes information on genetic values, inbreeding levels and the accuracy of the BVE

bvs <- list(numeric(100), numeric(100))
n_run <- 100
for(scenario in 1:2){
  for(run in 1:n_run){

    load(paste0(paste0("C:/Users/pook001/OneDrive - Wageningen University & Research/MoBPS_Workshop_WIAS_2025/Task12/Sheep_simulation_scenario", scenario, "run", run,".RData")))
    bvs[[scenario]][run] <-  genomic_values[1,which(cohorts=="1yearRams_20")]

  }
}

# Scenarios 1 & 2 are statistically highly significantly different
# Differnce between the two scenarios is a higher selection intensity. Thus, this should not
# be that surprising.
par(mfrow=c(2,1))
hist(bvs[[1]], xlim=c(110,130), main="Scenario 1", xlab="genomic value")
hist(bvs[[2]], xlim=c(110,130), main="Scenario 2", xlab="genomic value")
t.test(bvs[[1]], bvs[[2]])




BV_avg <- list()
ACC_avg <- list()
n_run <- 100
for(scenario in 1:4){
  for(run in 1:n_run){

    load(paste0(paste0("C:/Users/pook001/OneDrive - Wageningen University & Research/MoBPS_Workshop_WIAS_2025/Task12/Sheep_simulation_scenario", scenario, "run", run,".RData")))

    if(run==1){
      BV_avg[[scenario]] <-   genomic_values / n_run
      ACC_avg[[scenario]] <-  accuracies /n_run
    } else{
      BV_avg[[scenario]] <-  BV_avg[[scenario]] + genomic_values/ n_run
      ACC_avg[[scenario]] <- ACC_avg[[scenario]] + accuracies/ n_run
    }

  }
}

what_to_plot <- paste0("1yearRams_", 0:20)
to_plot <- which(cohorts %in% what_to_plot)

# Visualization of the underlying true genomic values of the 1yearRams in the different cycles
par(mfrow=c(1,2))
traitname <- c("Meat", "Fertility")
for(trait in 1:2){
  plot(0:20, BV_avg[[1]][trait,to_plot], type="l", main=traitname[trait],xlab="cycle", ylab="genomic value", lwd=2, ylim=c(100,122))
  for(sc in 2:4){
    lines(0:20, BV_avg[[sc]][trait,to_plot], col=sc, type="l", xlab="cycle", ylab="genomic value", lwd=2)
  }
}


legend("topleft", c("Baseline", "Higher selection intensity", "Lower Genotyping share", "Higher weighting on Meat"), lty=1, lwd=2, col=1:4)

# Visualization of the prediction accuracies of the 1yearEwes in the different cycles
what_to_plot <- paste0("1yearEwes_", 0:20)
to_plot <- which(cohorts %in% what_to_plot)
traitname <- c("Meat", "Fertility")
par(mfrow=c(1,2))
for(trait in 1:2){
  plot(0:20, ACC_avg[[1]][trait,to_plot], type="l", main=traitname[trait], xlab="cycle", ylab="prediction accuracy", lwd=2, ylim=c(0,1))
  for(sc in 2:4){
    lines(0:20, ACC_avg[[sc]][trait,to_plot], col=sc, type="l", xlab="cycle", ylab="prediction accuracy", lwd=2)
  }
}

legend("topleft", c("Baseline", "Higher selection intensity", "Lower Genotyping share", "Higher weighting on Meat"), lty=1, lwd=2, col=1:4)


# The accuracy of the BVE is very similar across years.
# To get more statistical power / reduce noise look at the averages across all years instead of the plots

mean(ACC_avg[[1]][2,to_plot], na.rm=TRUE)
mean(ACC_avg[[2]][2,to_plot], na.rm=TRUE)
# less genotyping leads to lower accuracies
mean(ACC_avg[[3]][2,to_plot], na.rm=TRUE)
mean(ACC_avg[[4]][2,to_plot], na.rm=TRUE)

mean(ACC_avg[[1]][1,to_plot], na.rm=TRUE)
mean(ACC_avg[[2]][1,to_plot], na.rm=TRUE)
# less genotyping leads to lower accuracies
mean(ACC_avg[[3]][1,to_plot], na.rm=TRUE)
mean(ACC_avg[[4]][1,to_plot], na.rm=TRUE)

