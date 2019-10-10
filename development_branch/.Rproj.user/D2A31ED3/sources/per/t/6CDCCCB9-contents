
library("MoBPS")
library("jsonlite")



path <- "./Rmodules/UserScripts/"

arg <- commandArgs(TRUE)
user <- arg[1]
filename <- arg[2]

a<- try(load(paste(path,user,"_",filename,".RData",sep="")))
if("try-error" %in% is(a)){
  stop("<span style='color: red; font-size:40px'>No old results of this breeding program to reload! This breeding program has never been successfully processed.</span>")
}
ct <- file.info(paste(path,user,"_",filename,".RData",sep=""))$ctime

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

result <- as.list(table(ttnames))
for(rr in names(result)){
     result[[rr]] <- list(tfounder=(coh[rr,"creating.type"] =="0"),trep=result[[rr]])                            
}

json <- as.character(toJSON(result))     
write.table(json, file=paste(path,user,"_",filename,"Summary.json",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(paste("<span style='color: red; font-size:40px'> Attention: You have loaded an old result of the breeding program '",filename,"' that has been produced on ",ct,".</span>", sep=""), file=paste(path,user,"_",filename,"ReloadSim.json",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)
################################################################################
#
#                    Data has been successfull reloaded
#
################################################################################






