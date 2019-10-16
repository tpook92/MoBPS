
library("MoBPS")
library("jsonlite")



path <- "./Rmodules/UserScripts/"

arg <- commandArgs(TRUE)


dat <- try(fromJSON(arg[2], simplifyVector=FALSE))
# dat <- jsonlite::read_json(path="/home/nha/PhänoDuinIndex.json")
if("try-error" %in% is(dat)){
  print(dat)
  stop("Cannot Read JSON file!")
}

fname <- try(dat$'Genomic Info'$'Project Name')
if("try-error" %in% is(fname)){
  print(fname)
  stop("Cannot Read File Name!")
}

fname <- paste(path, arg[1],"_",fname, sep="")

fileSum  <- paste(fname,"Summary.json",sep="")
if(file.exists(fileSum)){
  file.remove(fileSum)
}
if(arg[1]==paste0("EAAP","guest")){
  stop(paste0("EAAP","guest is not allowed to run simulations."))
}

if(length(dat$'Genomic Info'$'multi-mode')>0 && dat$'Genomic Info'$'multi-mode' =="Yes"){

  ## Perform phasing
  required_phasing <- file <- NULL
  for(index in 1:length(dat$Nodes)){
    if(length(dat$Nodes[[index]]$'phasing_required')==1 && dat$Nodes[[index]]$'phasing_required' == "Yes"){
      required_phasing <- c(required_phasing, index)
      file <- c(file, dat$Nodes[[index]]$'Path')
    }
  }
  nr <- 1
  for(path in unique(file)){
    isvcf <- (substr(path, start=nchar(path)-3, stop=nchar(path))==".vcf") ||(substr(path, start=nchar(path)-3, stop=nchar(path))=="f.gz")
    pedmap.to.phasedbeaglevcf(map_path= if(!isvcf){dat$'Genomic Info'$'Own Map Path'} else{NULL},
                              ped_path = if(!isvcf){path} else{NULL},
                              vcf_path = if(isvcf){path} else{NULL},
                              db_dir=paste0("/home/nha/Plink/DB/", nr))

    for(index in required_phasing[which(file==path)]){
      dat$Nodes[[index]]$'Path' <- paste0("/home/nha/Plink/DB/", nr, "temp_vcf_phased.vcf.gz")
      dat$Nodes[[index]]$'phasing_required' <- "No"
    }
    nr <- nr + 1
  }


  if(length(dat$'Genomic Info'$'number-simulations-parallel')==1){
    ncore <- as.numeric(dat$'Genomic Info'$'number-simulations-parallel')
  } else{
    ncore <- 1
  }

  if(ncore>10 && arg[1]!="Torsten"){
    ncore <- 10
    cat("You are not allowed to run more than 10 MoBPS sessions in parallel!")
  }
  doParallel::registerDoParallel(cores=ncore)
  if(length(as.numeric(dat$'Genomic Info'$'number-simulations'))==1){
    sims <- as.numeric(dat$'Genomic Info'$'number-simulations')
  } else{
    sims <- 1
  }
  library(doParallel)
  cat("\n\n\n\n\n")
  trash <- foreach::foreach(rep=1:sims, .packages = "MoBPS") %dopar% {

    if(rep==1){
      verbose <- TRUE
    } else{
      verbose <- FALSE
    }
    t1 <- Sys.time()
    cat(paste0("Start simulation number ", rep, "\n"))
    population <- try(MoBPS::json.simulation(total=dat, verbose=verbose))
    if("try-error" %in% is(population)){
      print(population)
      stop("Cannot Simulate Project!")
    }
    t2 <- Sys.time()
    if(rep==1){
      save(population, file=paste(fname, ".RData",sep=""))
    }
    save(population, file=paste(fname, rep, ".RData",sep=""))
    cat(paste0("Finished simulation number ", rep,".\n"))
    cat(paste0("Simulation took ", round(as.numeric(t2)-as.numeric(t1), digit=3), " seconds.\n"))
  }
  load(paste(fname, ".RData",sep=""))

  doParallel::stopImplicitCluster()
} else{
  t1 <- Sys.time()
  population <- try(json.simulation(total=dat))
  if("try-error" %in% is(population)){
    print(population)
    stop("Cannot Simulate Project!")
  }

  save(population, file=paste(fname, ".RData",sep=""))
  t2 <- Sys.time()
  cat(paste0("Finished simulation.\n"))
  cat(paste0("Simulation took ", round(as.numeric(t2)-as.numeric(t1), digit=3), " seconds.\n"))
}


coh <- get.cohorts(population, extended=TRUE)
ttnames <- NULL
ttrep <- NULL
for(nn in coh[,1]){
  nn_s <- strsplit(nn, "_")[[1]]
  if(!is.na(suppressWarnings(as.numeric(nn_s[length(nn_s)])))){
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
write.table(json, file=paste(fname,"Summary.json",sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE)

################################################################################
################ Simulation success without Errors ##############################
################################################################################





