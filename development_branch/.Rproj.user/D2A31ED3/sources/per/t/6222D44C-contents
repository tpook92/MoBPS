##############
# Additional analyses fot plotting results: Relationship
library("MoBPS")
library("jsonlite")



path <- "./Rmodules/UserScripts/"

arg <- commandArgs(TRUE)
user <- arg[1]
filename <- arg[2]
#cohorts <- fromJSON(arg[3])

load(paste(path,user,"_",filename,".RData",sep=""))
# Rel

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

filesnames <- dir(path)
check_time <- file.info(paste0("Rmodules/UserScripts/", filesnames))
original <- which(filesnames==paste0(user,"_",filename,".RData"))
newer <- as.numeric(check_time[,4])>=as.numeric(check_time[original,4])
filesnames <- filesnames[newer]
filesnames <- gsub(".RData","",filesnames)
filesnames <- gsub(paste0(user,"_",filename), "", filesnames)

avail <- suppressWarnings(unique(c(NA,as.numeric(filesnames)))[-1])


ttkinship <- list()
coh12 <- unique(ttnames)

if(length(avail)>1){
  for(index in avail){
    print(index)
    load(paste(path,user,"_",filename, index, ".RData",sep=""))
    for(i1 in 1:(length(coh12)-1)){
      ttkinship[[coh12[i1]]] <- list()
      for(i2 in (i1+1):length(coh12)){
        coh1 <- coh[ttnames == coh12[i1]]
        coh2 <- coh[ttnames == coh12[i1]]
        for(i in 1:min(c(length(coh1),length(coh2)))){
          ttkinship[[coh12[i1]]][[coh12[i2]]] <- c(ttkinship[[coh12[i1]]][[coh12[i2]]], 2*kinship.emp.fast(population=population, cohorts=c(coh1[i], coh2[i]),
                                                                                                           ibd.obs = 50, hbd.obs = 0)[1])
        }
      }
    }
    if(index==avail[1]){
      ttkinship1 <- ttkinship
    } else{
      for(index2 in 1:length(ttkinship)){
        for(index3 in 1:length(ttkinship[[index2]])){
          ttkinship1[[index2]][[index3]] <- ttkinship1[[index2]][[index3]] + ttkinship[[index2]][[index3]]
        }
      }
    }
  }

  for(index2 in 1:length(ttkinship)){
    for(index3 in 1:length(ttkinship[[index2]])){
      ttkinship[[index2]][[index3]] <- ttkinship1[[index2]][[index3]] / length(avail)
    }
  }

} else{
  for(i1 in 1:(length(coh12)-1)){
    ttkinship[[coh12[i1]]] <- list()
    for(i2 in (i1+1):length(coh12)){
      coh1 <- coh[ttnames == coh12[i1]]
      coh2 <- coh[ttnames == coh12[i1]]
      for(i in 1:min(c(length(coh1),length(coh2)))){
        ttkinship[[coh12[i1]]][[coh12[i2]]] <- c(ttkinship[[coh12[i1]]][[coh12[i2]]], 2*kinship.emp.fast(population=population, cohorts=c(coh1[i], coh2[i]),
                                                                                                         ibd.obs = 50, hbd.obs = 0)[1])
      }
    }
  }

}



result <- ttkinship

json <- as.character(toJSON(result))
write.table(json, file=paste(path,user,"_",filename,"RelbetweenC.json",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)




