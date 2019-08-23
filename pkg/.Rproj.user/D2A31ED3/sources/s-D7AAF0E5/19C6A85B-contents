##############
# Additional analyses fot plotting results: Relationship
library("MoBPS")
library("jsonlite")


path <- "./Rmodules/UserScripts/"

arg <- commandArgs(TRUE)
user <- arg[1]
filename <- arg[2]
sindex <- fromJSON(arg[3])

# Need way to figure out if rep was used + Efficiency stuff
if(FALSE && length(dat$'Genomic Info'$'multi-mode')>0 && dat$'Genomic Info'$'multi-mode' =="Yes"){

  for(rep in 1:as.numeric(dat$'Genomic Info'$'number-simulations')){
    cat(paste0("Analyze rep ", rep, ".\n"))
    load(paste(path,user,"_",filename,rep,".RData",sep=""))
    coh <- get.cohorts(population, extended=TRUE)
    ttnames <- NULL
    ttrep <- NULL
    for(nn in coh[,1]){
      nn_s <- strsplit(nn, "_")[[1]]
      if(suppressWarnings(!is.na(as.numeric(nn_s[length(nn_s)])))){
        ttnames <- c(ttnames, paste(nn_s[-length(nn_s)],collapse="_" ))
        ttrep <- c(ttrep, as.numeric(nn_s[length(nn_s)]))
      }else{
        ttnames <- c(ttnames, nn)
        ttrep <- c(ttrep, 0)
      }
    }

    for(i in 1:nrow(coh)){
      ani1 <- NULL
      if(coh[i,3] != 0){
        ani1 <- cbind(ani1, population$breeding[[as.numeric(coh[i,2])]][[3]][,as.numeric(coh[i,6]):(as.numeric(coh[i,6])+as.numeric(coh[i,3])-1), drop=FALSE])
      }
      if(coh[i,4] != 0){
        ani1 <- cbind(ani1, population$breeding[[as.numeric(coh[i,2])]][[4]][,as.numeric(coh[i,7]):(as.numeric(coh[i,7])+as.numeric(coh[i,4])-1), drop=FALSE])
      }

      ani2 <- NULL
      if(coh[i,3] != 0){
        ani2 <- cbind(ani2, population$breeding[[as.numeric(coh[i,2])]][[7]][,as.numeric(coh[i,6]):(as.numeric(coh[i,6])+as.numeric(coh[i,3])-1), drop=FALSE])
      }
      if(coh[i,4] != 0){
        ani2 <- cbind(ani2, population$breeding[[as.numeric(coh[i,2])]][[8]][,as.numeric(coh[i,7]):(as.numeric(coh[i,7])+as.numeric(coh[i,4])-1), drop=FALSE])
      }

      corr1 <- NULL
      for(ind in 1:nrow(sindex)){
        corr1 <- c(corr1, cor(colSums(as.numeric(unlist(sindex[ind,2:ncol(sindex)]))*ani1), colSums(as.numeric(unlist(sindex[ind,2:ncol(sindex)]))*ani2)))
      }
      if(rep==1){
        corr <- corr1 / as.numeric(dat$'Genomic Info'$'number-simulations')
      } else{
        corr <- corr + corr1 / as.numeric(dat$'Genomic Info'$'number-simulations')
      }

      Acc <- list()
      Acc[[ttnames[i]]][[as.character(ttrep[i])]] <- list(ttime=coh[i,"time point"],tval=corr)
    }
  }

  result <- Acc

  #dat <- by(result, result[,1], t)
  #class(dat) <- "list"

  json <- as.character(toJSON(result))
  write.table(json, file=paste(path,user,"_",filename,"AccBVE.json",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)

} else{

load(paste(path,user,"_",filename,".RData",sep=""))
# Rel

coh <- get.cohorts(population, extended=TRUE)
ttnames <- NULL
ttrep <- NULL
for(nn in coh[,1]){
  nn_s <- strsplit(nn, "_")[[1]]
  if(!is.na(as.numeric(nn_s[length(nn_s)]))){
    ttnames <- c(ttnames, paste(nn_s[-length(nn_s)],collapse="_" ))
    ttrep <- c(ttrep, as.numeric(nn_s[length(nn_s)]))
  }else{
    ttnames <- c(ttnames, nn)
    ttrep <- c(ttrep, 0)
  }
}

#keep <- ttnames %in% cohorts
#coh <- coh[keep,, drop=FALSE]
#ttrep <- ttrep[keep]
#ttnames <- ttnames[keep]

Acc <- list()
for(i in 1:nrow(coh)){
  ani1 <- NULL
  if(coh[i,3] != 0){
    ani1 <- cbind(ani1, population$breeding[[as.numeric(coh[i,2])]][[3]][,as.numeric(coh[i,6]):(as.numeric(coh[i,6])+as.numeric(coh[i,3])-1), drop=FALSE])
  }
  if(coh[i,4] != 0){
    ani1 <- cbind(ani1, population$breeding[[as.numeric(coh[i,2])]][[4]][,as.numeric(coh[i,7]):(as.numeric(coh[i,7])+as.numeric(coh[i,4])-1), drop=FALSE])
  }

  ani2 <- NULL
  if(coh[i,3] != 0){
    ani2 <- cbind(ani2, population$breeding[[as.numeric(coh[i,2])]][[7]][,as.numeric(coh[i,6]):(as.numeric(coh[i,6])+as.numeric(coh[i,3])-1), drop=FALSE])
  }
  if(coh[i,4] != 0){
    ani2 <- cbind(ani2, population$breeding[[as.numeric(coh[i,2])]][[8]][,as.numeric(coh[i,7]):(as.numeric(coh[i,7])+as.numeric(coh[i,4])-1), drop=FALSE])
  }

  corr <- NULL
  for(ind in 1:nrow(sindex)){
    corr <- c(corr, cor(colSums(as.numeric(unlist(sindex[ind,2:ncol(sindex)]))*ani1), colSums(as.numeric(unlist(sindex[ind,2:ncol(sindex)]))*ani2)))
  }
  Acc[[ttnames[i]]][[as.character(ttrep[i])]] <- list(ttime=coh[i,"time point"],tval=corr)
}




result <- Acc

#dat <- by(result, result[,1], t)
#class(dat) <- "list"

json <- as.character(toJSON(result))
write.table(json, file=paste(path,user,"_",filename,"AccBVE.json",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)

}


