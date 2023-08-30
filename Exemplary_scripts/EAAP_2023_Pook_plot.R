setwd("C:/Users/pook001/OneDrive - Wageningen University & Research/")
n = 500
nsc = 48

bv_total = bv_sd_total = inb_total= share_fixed_total = share_fixed_qtl_total =
  share_fixed_good_total = share_fixed_bad_total = matrix(0, nrow=nsc, ncol=111)
for(scenario in (1:nsc)){
  n_temp = 0
  for(seed in c(1:n)){


    tryCatch(  {


      load(paste0("uniqueness_v3/sc", scenario , "seed", seed, "_potential.RData"))



      bv_total[scenario,] = bv_total[scenario,] + bv
      bv_sd_total[scenario,] = bv_sd_total[scenario,] + bv_sd
      inb_total[scenario,] = inb_total[scenario,] + inb
      share_fixed_total[scenario,] = share_fixed_total[scenario,] + share_fixed
      share_fixed_qtl_total[scenario,] = share_fixed_qtl_total[scenario,] + share_fixed_qtl
      share_fixed_good_total[scenario,] = share_fixed_good_total[scenario,] + share_fixed_good
      share_fixed_bad_total[scenario,] = share_fixed_bad_total[scenario,] + share_fixed_bad

      n_temp = n_temp + 1
    },
    error = function(e) {})


  }

  bv_total[scenario,] = bv_total[scenario,]/n_temp
  bv_sd_total[scenario,] = bv_sd_total[scenario,]/n_temp
  inb_total[scenario,] = inb_total[scenario,]/n_temp
  share_fixed_total[scenario,] = share_fixed_total[scenario,]/n_temp
  share_fixed_qtl_total[scenario,] = share_fixed_qtl_total[scenario,] /n_temp
  share_fixed_good_total[scenario,] = share_fixed_good_total[scenario,]/n_temp
  share_fixed_bad_total[scenario,] = share_fixed_bad_total[scenario,]/n_temp

  print(n_temp)
}

save(file = "total.RData", list = c("bv_total",
                                                  "bv_sd_total",
                                                  "inb_total",
                                                  "share_fixed_total",
                                                  "share_fixed_qtl_total",
                                                  "share_fixed_good_total",
                                                  "share_fixed_bad_total"))

load("total.RData")

legend =     c("EBV", "BV", "Pheno", "min(5,p^(-1/2))", "p^(-1/2)", "p^-1", "p^(-1/3)","(1-p)", "(1-p)^2",
               "FAILED", "kinship to top individuals_0.05",  "inbreeding (pedigree)_0.05", "inbreeding (genetic)_0.05",
               "FAILED", "kinship to top individuals_0.2", "inbreeding (pedigree)_0.2", "inbreeding (genetic)_0.2",
               "Number of rare benefical_0.05", "Number of rare_0.05", "Number of rare benefical_0.2", "Number of rare_0.2",  "p^(-1/2)_limit", "p^(-1/2)", "p^-1", "p^(-1/3)","(1-p)",
               "(1-p)^2", "random_error_0.05", "random_error_0.2", "Index5",
               "Number of rare strongly ben_0.05", "Number of rare strongly ben_0.2",
               "p^(-1/3) & kinship to top individuals_0.2",
               "kinship to top individuals_0.2 & Index5",
               "p^(-1/3) & kinship to top individuals_0.2 & Index5",
               "kinship to population_0.05",
               "kinship to population_0.2"

)

used_color = NULL

#################
### Plots #######


png(filename = "EAAP_1.png", res = 300,
    height = 1500, width = 3000)
limits = c(-1,3.5)

par(mfrow=c(1,2))
include = c(2,23)
used_color[include] = 1:8

plot(0,0, xlim=c(0,50), ylim=limits, main = "QTL effects known", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
  lines(0:50, bv_total[index,11:61] - bv_total[2,11:61], col=used_color[index], lwd=2)

}

include = c(1,5)
used_color[include] = 1:8

plot(0,0, xlim=c(0,50), ylim=limits, main = "QTL effects estimated", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
  lines(0:50, bv_total[index,11:61] - bv_total[1,11:61], col=used_color[index], lwd=2)

}

dev.off()


png(filename = "EAAP_2.png", res = 300,
    height = 1500, width = 3000)
limits = c(-1,3.5)

par(mfrow=c(1,2))
include = c(2,23)
used_color[include] = 1:8

plot(0,0, xlim=c(0,50), ylim=limits, main = "QTL effects known", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
  lines(0:50, bv_total[index,11:61] - bv_total[2,11:61], col=used_color[index], lwd=2)
  lines(0:50, bv_total[index,c(11,62:111)] - bv_total[2,c(11,62:111)], lty = 2, col=used_color[index], lwd=2)

}

include = c(1,5)
used_color[include] = 1:8

plot(0,0, xlim=c(0,50), ylim=limits, main = "QTL effects estimated", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
  lines(0:50, bv_total[index,11:61] - bv_total[1,11:61], col=used_color[index], lwd=2)
  lines(0:50, bv_total[index,c(11,62:111)] - bv_total[1,c(11,62:111)], lty = 2, col=used_color[index], lwd=2)

}

dev.off()

png(filename = "EAAP_1b.png", res = 300,
    height = 1500, width = 3000)

par(mfrow = c(1,2))
plot(((bv_total[23,11:61] - bv_total[2,11:61])) / (bv_total[2,11:61] - 100)*100,
     xlab = "generation", ylab = "in percent", main = "Relative gains to BVE selection", ylim = c(-12,10))
abline(h = 0, lwd = 2)
plot(((bv_total[5,11:61] - bv_total[1,11:61])) / (bv_total[1,11:61] - 100)*100,
     xlab = "generation", ylab = "in percent", main = "Relative gains to BVE selection", ylim = c(-12,10))
abline(h = 0, lwd = 2)

dev.off()

png(filename = "EAAP_2b.png", res = 300,
    height = 1500, width = 3000)

par(mfrow = c(1,2))
plot(((bv_total[23,11:61] - bv_total[2,11:61])) / (bv_total[2,11:61] - 100)*100,
     xlab = "generation", ylab = "in percent", main = "Relative gains to BVE selection", ylim = c(-12,10))

points(((bv_total[23,c(11,62:111)] - bv_total[2,c(11,62:111)])) / (bv_total[2,c(11,62:111)] - 100)*100,
     xlab = "generation", ylab = "in percent", main = "Relative gains to BVE selection", col = 2, ylim = c(-12,10))
abline(h = 0, lwd = 2)
plot(((bv_total[5,11:61] - bv_total[1,11:61])) / (bv_total[1,11:61] - 100)*100,
     xlab = "generation", ylab = "in percent", main = "Relative gains to BVE selection", ylim = c(-12,10))

points(((bv_total[5,c(11,62:111)] - bv_total[1,c(11,62:111)])) / (bv_total[1,c(11,62:111)] - 100)*100,
     xlab = "generation", col = 2, ylab = "in percent", main = "Relative gains to BVE selection", ylim = c(-12,10))
abline(h = 0, lwd = 2)

dev.off()

plot(((bv_total[23,11:61] - bv_total[2,11:61])) )
plot(((bv_total[5,11:61] - bv_total[1,11:61])) )

plot(((bv_total[23,11:61 + 50] - bv_total[2,11:61 + 50])) / (bv_total[2,11:61 + 50] - 100))
plot(((bv_total[5,11:61 + 50] - bv_total[1,11:61 + 50])) / (bv_total[1,11:61 + 50] - 100))

plot(((bv_total[23,11:61 + 50] - bv_total[2,11:61 + 50])))
plot(((bv_total[5,11:61 + 50] - bv_total[1,11:61 + 50])))
## weighting ebv


png(filename = "EAAP_3.png", res = 300,
    height = 1500, width = 3000)
limits = c(-0.5,0.5)

par(mfrow=c(1,2))

include = c(1,5,4,6:9)
used_color[include] = 1:8

plot(0,0, xlim=c(0,50), ylim=limits, main = "genetic gain", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
#  lines(0:50, bv_total[index,11:61] - bv_total[1,11:61], col=used_color[index], lwd=2)
  lines(0:50, bv_total[index,c(11,62:111)] - bv_total[1,c(11,62:111)], lty = 1, col=used_color[index], lwd=2)

}
legend("topleft", legend[include],
       col=used_color[include], lty = 1, lwd=2, cex = 0.7)

plot(-10000,0, xlim=c(0,50), ylim=c(0,1), main = "inbreeding", ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(0:50, inb_total[index,11:61], col=used_color[index], lwd=2)
}

dev.off()

png(filename = "EAAP_4.png", res = 300,
    height = 1500, width = 3000)
limits = c(-3,4)

par(mfrow=c(1,2))

include = c(2,23,22,24:27)
used_color[include] = 1:8

plot(0,0, xlim=c(0,50), ylim=limits, main = "genetic gain", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
#  lines(0:50, bv_total[index,11:61] - bv_total[2,11:61], col=used_color[index], lwd=2)
  lines(0:50, bv_total[index,c(11,62:111)] - bv_total[2,c(11,62:111)], lty = 1, col=used_color[index], lwd=2)

}
legend("topleft", legend[include],
       col=used_color[include], lty = 1, lwd=2)

plot(-10000,0, xlim=c(0,50), ylim=c(0,1), main = "inbreeding", ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(0:50, inb_total[index,11:61], col=used_color[index], lwd=2)
}

dev.off()



png(filename = "EAAP_5.png", res = 300,
    height = 1500, width = 3000)
limits = c(-0.5,0.5)

par(mfrow=c(1,2))

include = c(1,7,19,18,31)
used_color[include] = c(1,5,3,7,4,6)

plot(0,0, xlim=c(0,50), ylim=limits, main = "genetic gain", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
#  lines(0:50, bv_total[index,11:61] - bv_total[1,11:61], col=used_color[index], lwd=2)
  lines(0:50, bv_total[index,c(11,62:111)] - bv_total[1,c(11,62:111)], lty = 1, col=used_color[index], lwd=2)

}
legend("topleft", legend[include],
       col=used_color[include], lty = 1, lwd=2, cex = 0.7)

plot(-10000,0, xlim=c(0,50), ylim=c(0,1), main = "inbreeding", ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(0:50, inb_total[index,11:61], col=used_color[index], lwd=2)
}

dev.off()

png(filename = "EAAP_5b.png", res = 300,
    height = 1500, width = 3000)
limits = c(-0.5,0.5)

par(mfrow=c(1,2))

include = c(1,7,21,20,32)
used_color[include] = c(1,5,3,7,4,6)

plot(0,0, xlim=c(0,50), ylim=limits, main = "genetic gain", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
  #  lines(0:50, bv_total[index,11:61] - bv_total[1,11:61], col=used_color[index], lwd=2)
  lines(0:50, bv_total[index,c(11,62:111)] - bv_total[1,c(11,62:111)], lty = 1, col=used_color[index], lwd=2)

}
legend("topleft", legend[include],
       col=used_color[include], lty = 1, lwd=2, cex = 0.7)

plot(-10000,0, xlim=c(0,50), ylim=c(0,1), main = "inbreeding", ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(0:50, inb_total[index,11:61], col=used_color[index], lwd=2)
}

dev.off()

png(filename = "EAAP_6.png", res = 300,
    height = 1500, width = 3000)
limits = c(-0.2,0.5)

par(mfrow=c(1,2))

include = c(1,11,36,12,13)
used_color[include] = 1:8

plot(0,0, xlim=c(0,50), ylim=limits, main = "genetic gain", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
  #lines(0:50, bv_total[index,11:61] - bv_total[1,11:61], col=used_color[index], lwd=2)
  lines(0:50, bv_total[index,c(11,62:111)] - bv_total[1,c(11,62:111)], lty = 1, col=used_color[index], lwd=2)

}
legend("topleft", legend[include],
       col=used_color[include], lty = 1, lwd=2, cex = 0.7)

plot(-10000,0, xlim=c(0,50), ylim=c(0,1), main = "inbreeding", ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(0:50, inb_total[index,11:61], col=used_color[index], lwd=2)
}

dev.off()

png(filename = "EAAP_7a.png", res = 300,
    height = 1500, width = 3000)
limits = c(-0.2,0.5)

par(mfrow=c(1,2))

include = c(1,15,37,16,17)
used_color[include] = 1:8

plot(0,0, xlim=c(0,50), ylim=limits, main = "genetic gain", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
 # lines(0:50, bv_total[index,11:61] - bv_total[1,11:61], col=used_color[index], lwd=2)
  lines(0:50, bv_total[index,c(11,62:111)] - bv_total[1,c(11,62:111)], lty = 1, col=used_color[index], lwd=2)

}
legend("topleft", legend[include],
       col=used_color[include], lty = 1, lwd=2, cex = 0.7)

plot(-10000,0, xlim=c(0,50), ylim=c(0,1), main = "inbreeding", ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(0:50, inb_total[index,11:61], col=used_color[index], lwd=2)
}

dev.off()


png(filename = "EAAP_7b.png", res = 300,
    height = 1500, width = 3000)
limits = c(-0.4,1.1)

par(mfrow=c(1,2))

include = c(1,15,37,16,17)
used_color[include] = 1:8

plot(0,0, xlim=c(0,50), ylim=limits, main = "genetic gain", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
  # lines(0:50, bv_total[index,11:61] - bv_total[1,11:61], col=used_color[index], lwd=2)
  lines(0:50, bv_total[index,c(11,62:111)] - bv_total[1,c(11,62:111)], lty = 1, col=used_color[index], lwd=2)

}
legend("topleft", legend[include],
       col=used_color[include], lty = 1, lwd=2, cex = 0.7)

plot(-10000,0, xlim=c(0,50), ylim=c(0,1), main = "inbreeding", ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(0:50, inb_total[index,11:61], col=used_color[index], lwd=2)
}

dev.off()

points(((bv_total[5,c(11,62:111)] - bv_total[1,c(11,62:111)])) / (bv_total[1,c(11,62:111)] - 100)*100,
       xlab = "generation", col = 2, ylab = "in percent", main = "Relative gains to BVE selection", ylim = c(-12,10))

png(filename = "EAAP_8.png", res = 300,
    height = 1500, width = 3000)
limits = c(-0.5,2.5)
par(mfrow=c(1,2))
limits = c(-0.5,1)

include = c(1,7,15,33)
used_color[include] = 1:8

plot(0,0, xlim=c(0,50), ylim=limits, main = "genetic gain", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
#  lines(0:50, bv_total[index,11:61] - bv_total[1,11:61], col=used_color[index], lwd=2)
  lines(0:50, bv_total[index,c(11,62:111)] - bv_total[1,c(11,62:111)], lty = 1, col=used_color[index], lwd=2)

}
legend("topleft", legend[include],
       col=used_color[include], lty = 1, lwd=2, cex = 0.7)

plot(-10000,0, xlim=c(0,50), ylim=c(0,1), main = "inbreeding", ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(0:50, inb_total[index,11:61], col=used_color[index], lwd=2)
}

dev.off()


par(mfrow=c(1,2))
include = c(38,45:48)

used_color = numeric(12)
used_color[include] = 1:8


plot(0,0, xlim=c(0,50), ylim=limits, main = "genetic gain", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
  lines(0:50, bv_total[index,11:61] - bv_total[38,11:61], col=used_color[index], lwd=2)

}

legend("topleft", legend[include],
       col=used_color[include], lty = 1, lwd=2)



plot(0,0, xlim=c(10,60), ylim=c(0,1), main = "inbreeding", ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(inb_total[index,], col=used_color[index], lwd=2)
}


par(mfrow=c(1,2))
include = c(7,20)
used_color = numeric(12)
used_color[include] = 1:8


plot(0,0, xlim=c(10,60), ylim=c(-1.5,2), main = "genetic gain", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
  lines(bv_total[index,] - bv_total[7,], col=used_color[index], lwd=2)
}


plot(0,0, xlim=c(10,60), ylim=c(0,1), main = "inbreeding", ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(inb_total[index,], col=used_color[index], lwd=2)
}

## weighting bv

par(mfrow=c(1,2))
include = 7:14
used_color = numeric(12)
used_color[include] = 1:8


plot(0,0, xlim=c(10,60), ylim=c(-3,5), main = "genetic gain", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
  lines(bv_total[index,] - bv_total[7,], col=used_color[index], lwd=2)
}

legend("topleft",
       c("phaeno", "EBV", "true BV", "real good rare", "est_5", "all_5", "PC_avg_l1", "PC_top_l1",
         "1/freq_bv_wrong", "1/freq_ebv_wrong", "PC_avg_l2", "PC_top_l2", "1/freq_bv", "1/freq_ebv", "1/freq^2_bv", "1/freq^2_ebv",
         "1-Fst","1-Fst 1-maf","1 - Fst range01","1 - Fst range01 1-maf","1 - Fstunscaled","1 - Fstunsclaed 1-maf",
         "1/freq^3_bv", "1/freq^3_ebv", "1-freq^2_bv", "1-freq_ebv", "(1-freq)^2_bv", "(1-freq)^2_ebv",
         "est_10", "all_10", "est_20", "all_20", "est_30", "all_30",
         "min RelTop", "min RelTop Diag", "MinInb", "MinInbGeno", "Min RelAll")[include],
       col=used_color[include], lty = 1, lwd=2)



plot(0,0, xlim=c(10,60), ylim=c(0,1), main = "inbreeding", ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(inb_total[index,], col=used_color[index], lwd=2)
}


## Fst based things
par(mfrow=c(1,2))
include = c(1,2,17:22)
used_color = numeric(12)
used_color[include] = 1:8


plot(0,0, xlim=c(10,60), ylim=c(-1.5,2), main = "genetic gain", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
  lines(bv_total[index,] - bv_total[2,], col=used_color[index], lwd=2)
}

legend("topleft",
       c("phaeno", "EBV", "true BV", "real good rare", "est_5", "all_5", "PC_avg_l1", "PC_top_l1",
         "1/freq_bv_wrong", "1/freq_ebv_wrong", "PC_avg_l2", "PC_top_l2", "1/freq_bv", "1/freq_ebv", "1/freq^2_bv", "1/freq^2_ebv",
         "1-Fst","1-Fst 1-maf","1 - Fst range01","1 - Fst range01 1-maf","1 - Fstunscaled","1 - Fstunsclaed 1-maf",
         "1/freq^3_bv", "1/freq^3_ebv", "1-freq^2_bv", "1-freq_ebv", "(1-freq)^2_bv", "(1-freq)^2_ebv",
         "est_10", "all_10", "est_20", "all_20", "est_30", "all_30",
         "min RelTop", "min RelTop Diag", "MinInb", "MinInbGeno", "Min RelAll")[include],
       col=used_color[include], lty = 1, lwd=2)



plot(0,0, xlim=c(10,60), ylim=c(0,1), main = "inbreeding", ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(inb_total[index,], col=used_color[index], lwd=2)
}

## rare allele based things
par(mfrow=c(1,2))
include = c(2,5,6,29:34)
used_color = numeric(12)
used_color[include] = c(2:8,1)


plot(0,0, xlim=c(10,60), ylim=c(-1.5,2), main = "genetic gain", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
  lines(bv_total[index,] - bv_total[2,], col=used_color[index], lwd=2)
}

legend("topleft",
       c("phaeno", "EBV", "true BV", "real good rare", "est_5", "all_5", "PC_avg_l1", "PC_top_l1",
         "1/freq_bv_wrong", "1/freq_ebv_wrong", "PC_avg_l2", "PC_top_l2", "1/freq_bv", "1/freq_ebv", "1/freq^2_bv", "1/freq^2_ebv",
         "1-Fst","1-Fst 1-maf","1 - Fst range01","1 - Fst range01 1-maf","1 - Fstunscaled","1 - Fstunsclaed 1-maf",
         "1/freq^3_bv", "1/freq^3_ebv", "1-freq^2_bv", "1-freq_ebv", "(1-freq)^2_bv", "(1-freq)^2_ebv",
         "est_10", "all_10", "est_20", "all_20", "est_30", "all_30",
         "min RelTop", "min RelTop Diag", "MinInb", "MinInbGeno", "Min RelAll")[include],
       col=used_color[include], lty = 1, lwd=2)



plot(0,0, xlim=c(10,60), ylim=c(0,1), main = "inbreeding", ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(inb_total[index,], col=used_color[index], lwd=2)
}

par(mfrow=c(1,2))
include = c(1,2,35:39)
used_color = numeric(12)
used_color[include] = c(1:8)


plot(0,0, xlim=c(10,60), ylim=c(-1.5,2), main = "genetic gain", ylab = "genetic gain in gSD", xlab="generation")
for(index in include){
  lines(bv_total[index,] - bv_total[2,], col=used_color[index], lwd=2)
}

legend("topleft",
       c("phaeno", "EBV", "true BV", "real good rare", "est_5", "all_5", "PC_avg_l1", "PC_top_l1",
         "1/freq_bv_wrong", "1/freq_ebv_wrong", "PC_avg_l2", "PC_top_l2", "1/freq_bv", "1/freq_ebv", "1/freq^2_bv", "1/freq^2_ebv",
         "1-Fst","1-Fst 1-maf","1 - Fst range01","1 - Fst range01 1-maf","1 - Fstunscaled","1 - Fstunsclaed 1-maf",
         "1/freq^3_bv", "1/freq^3_ebv", "1-freq^2_bv", "1-freq_ebv", "(1-freq)^2_bv", "(1-freq)^2_ebv",
         "est_10", "all_10", "est_20", "all_20", "est_30", "all_30",
         "min RelTop", "min RelTop Diag", "MinInb", "MinInbGeno", "Min RelAll")[include],
       col=used_color[include], lty = 1, lwd=2)



plot(0,0, xlim=c(10,60), ylim=c(0,1), main = "inbreeding", ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(inb_total[index,], col=used_color[index], lwd=2)
}



include = c(1:6,9:10)
used_color = numeric(12)
used_color[include] = 1:8

plot(0,0, xlim=c(10,60), ylim=c(-0.2,0.2), ylab = "inbreeding level", xlab="generation")
for(index in include){
  lines(inb_total[index,] - inb_total[2,], col=used_color[index], lwd=2)
}

legend("topright",
       c("phaeno", "EBV", "true BV", "real good rare", "estimated good rare", "everything rare", "PC_avg_l1", "PC_top_l1",
         "weight_bv", "weight_ebv", "PC_avg_l2", "PC_top_l2")[include],
       col=used_color[include], lty = 1, lwd=2)

include = c(1,2,7,8,11,12)
used_color = numeric(12)
used_color[include] = 1:8

plot(0,0, xlim=c(10,60), ylim=c(-0.2,0.2))
for(index in include){
  lines(inb_total[index,] - inb_total[2,], col=used_color[index], lwd=2)
}

legend("topright",
       c("phaeno", "EBV", "true BV", "real good rare", "estimated good rare", "everything rare", "PC_avg_l1", "PC_top_l1",
         "weight_bv", "weight_ebv", "PC_avg_l2", "PC_top_l2")[include],
       col=used_color[include], lty = 1, lwd=2)




par(mfrow=c(1,3))
include = c(1,2,7,8,11,12)
include = c(1:6,9:10)
used_color = numeric(12)
used_color[include] = 1:8

plot(0,0, xlim=c(10,60), ylim=c(-0.2,0.2), main = "Share Markers fixed", ylab = "share fixed")
for(index in include){
  lines(share_fixed_qtl_total[index,] - share_fixed_qtl_total[2,], col=used_color[index], lwd=2)
}

plot(0,0, xlim=c(10,60), ylim=c(-0.2,0.2), main = "Share Markers for good variant", ylab = "share fixed")
for(index in include){
  lines(share_fixed_good_total[index,] - share_fixed_good_total[2,], col=used_color[index], lwd=2)
}

plot(0,0, xlim=c(10,60), ylim=c(-0.2,0.2), main = "Share Markers fixed for bad variant", ylab = "share fixed")
for(index in include){
  lines(share_fixed_bad_total[index,] - share_fixed_bad_total[2,], col=used_color[index], lwd=2)
}


legend("topright",
       c("phaeno", "EBV", "true BV", "real good rare", "estimated good rare", "everything rare", "PC_avg_l1", "PC_top_l1",
         "weight_bv", "weight_ebv", "PC_avg_l2", "PC_top_l2")[include],
       col=used_color[include], lty = 1, lwd=2)

