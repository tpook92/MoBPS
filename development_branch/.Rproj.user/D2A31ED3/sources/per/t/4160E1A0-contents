##############
# Additional analyses fot plotting results:  QTL
library("MoBPS")
library("jsonlite")



path <- "./Rmodules/UserScripts/"

arg <- commandArgs(TRUE)
# trait <- fromJSON("/home/nha/Simple_Cattle_mod.json", simplifyVector=FALSE)$`Trait Info`
# arg <- c("Torsten", "Simple_Cattle")
user <- arg[1]
filename <- arg[2]
trait <- fromJSON(arg[3], simplifyVector=FALSE)


load(paste(path,user,"_",filename,".RData",sep=""))
# QTL

coh <- get.cohorts(population)
ttnames <- NULL
ttrep <- NULL
for(nn in coh){
  nn_s <- strsplit(nn, "_")[[1]]
  if(!is.na(as.numeric(nn_s[length(nn_s)]))){
    ttnames <- c(ttnames, paste(nn_s[-length(nn_s)],collapse="_" ))
    ttrep <- c(ttrep, as.numeric(nn_s[length(nn_s)]))
  }else{
    ttnames <- c(ttnames, nn)
    ttrep <- c(ttrep, 0)
  }
}


# Check if directory contains multiple simulations
filesnames <- dir(path)
check_time <- file.info(paste0("Rmodules/UserScripts/", filesnames))
original <- which(filesnames==paste0(user,"_",filename,".RData"))
newer <- as.numeric(check_time[,4])>=as.numeric(check_time[original,4])
filesnames <- filesnames[newer]
filesnames <- gsub(".RData","",filesnames)
filesnames <- gsub(paste0(user,"_",filename), "", filesnames)

avail <- suppressWarnings(unique(c(NA,as.numeric(filesnames)))[-1])

result <- list()
if(length(avail)>1){
  result1 <- list()
  for(index in avail){
    result <- list()
    cat(index)
    load(paste(path,user,"_",filename, index, ".RData",sep=""))
    for(tr in 1:length(trait)){
      result[[trait[[tr]][['Trait Name']]]] <- list()
      check <- try(length(trait[[tr]][["Trait QTL Info"]] )!= 0)
      print(check)
      if(!("try-error" %in% is(check)) & check == TRUE){
        for(qtl in 1:length(trait[[tr]][["Trait QTL Info"]])){
          snp <-  population$info$real.bv.add[[tr]][qtl,]
          ttfreq <- analyze.population(population, snp[2], snp[1], cohorts=coh)
          freq <- (ttfreq[3,]+ttfreq[2,]/2)/colSums(ttfreq)
          oH <-  ttfreq[2,]/colSums(ttfreq)
          eH <- 2*freq*(1-freq)
          result[[trait[[tr]][['Trait Name']]]][[qtl]] <- by(cbind(ttnames, ttrep, freq, oH, eH), ttnames, t)
          class(result[[trait[[tr]][['Trait Name']]]][[qtl]]) <- "list"
        }
      }
    }
    if(index==avail[1]){
      result1 <- result
    } else{
      for(tr in 1:length(trait)){
        check <- try(length(trait[[tr]][["Trait QTL Info"]] )!= 0)
        print(check)
        if(!("try-error" %in% is(check)) & check == TRUE){
          for(qtl in 1:length(trait[[tr]][["Trait QTL Info"]])){
            for(coh in unique(ttnames)){
              result1[[trait[[tr]][['Trait Name']]]][[qtl]][[coh]][-1,] <- as.numeric(result1[[trait[[tr]][['Trait Name']]]][[qtl]][[coh]][-1,]) + as.numeric(result[[trait[[tr]][['Trait Name']]]][[qtl]][[coh]][-1,])
            }
          }
        }
      }
    }

  }
  for(tr in 1:length(trait)){
    check <- try(length(trait[[tr]][["Trait QTL Info"]] )!= 0)
    print(check)
    if(!("try-error" %in% is(check)) & check == TRUE){
      for(qtl in 1:length(trait[[tr]][["Trait QTL Info"]])){
        for(coh in unique(ttnames)){
          result[[trait[[tr]][['Trait Name']]]][[qtl]][[coh]][-1,] <- as.numeric(result1[[trait[[tr]][['Trait Name']]]][[qtl]][[coh]][-1,]) / length(avail)
        }
      }
    }
  }
} else{
  for(tr in 1:length(trait)){
    result[[trait[[tr]][['Trait Name']]]] <- list()
    check <- try(length(trait[[tr]][["Trait QTL Info"]] )!= 0)
    print(check)
    if(!("try-error" %in% is(check)) & check == TRUE){
      for(qtl in 1:length(trait[[tr]][["Trait QTL Info"]])){
        snp <-  population$info$real.bv.add[[tr]][qtl,]
        ttfreq <- analyze.population(population, snp[2], snp[1], cohorts=coh)
        freq <- (ttfreq[3,]+ttfreq[2,]/2)/colSums(ttfreq)
        oH <-  ttfreq[2,]/colSums(ttfreq)
        eH <- 2*freq*(1-freq)
        result[[trait[[tr]][['Trait Name']]]][[qtl]] <- by(cbind(ttnames, ttrep, freq, oH, eH), ttnames, t)
        class(result[[trait[[tr]][['Trait Name']]]][[qtl]]) <- "list"
      }
    }
  }
}




json <- as.character(toJSON(result))
write.table(json, file=paste(path,user,"_",filename,"QTL.json",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)








