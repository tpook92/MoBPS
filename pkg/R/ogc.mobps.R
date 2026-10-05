'#
  Authors
Torsten Pook, torsten.pook@wur.nl

Copyright (C) 2017 -- 2025  Torsten Pook

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
'#


#' Breeding function
#'
#' Function to simulate a step in a breeding scheme
#' @param population Population list
#' @param animallist List of individuals to include in ocs
#' @param ocs.target Target of OSC (default: "min.sKin" - minimize inbreeding; alt: "max.BV" / "min.BV" - maximize genetic gain; both under constrains selected below)
#' @param ocs.uniform This corresponds to the uniform constrain in optiSel
#' @param ocs.lb This corresponds to the lb constrain in optiSel
#' @param ocs.ub This corresponds to the ub constrain in optiSel
#' @param ocs.ub.sKin This corresponds to the ub.sKin constrain in optiSel
#' @param ocs.lb.BV This corresponds to the lb.BV constrain in optiSel
#' @param ocs.ub.BV This corresponds to the ub.BV constrain in optiSel
#' @param ocs.eq.BV This corresponds to the eq.BV constrain in optiSel
#' @param ocs.ub.sKin.increase This corresponds to the upper bound (current sKin + ocs.ub.sKin.increase) as ub.sKin in optiSel
#' @param ocs.lb.BV.increase This corresponds to the lower bound (current BV + ocs.lb.BV.increase) as lb.BV in optiSel
#' @param relationship.matrix.ocs Method to calculate relationship matrix for OCS (Default: "pedigree", alt: "vanRaden", "CE", "non_stand", "CE2", "CM")
#' @param depth.pedigree.ocs Depth of the pedigree in generations (default: 7)
#' @param bve.pedigree.error Set to FALSE to ignore/correct for any pedigree errors
#' @param import.position.calculation XXX
#' @param decodeOriginsU XXX
#' @param store.sparse XXX
#' @param miraculix XXX
#' @param miraculix.mult XXX
#' @param bve.p_i.list XXX
#' @param verbose Set to FALSE to not display any prints
#' @examples
#' population <- creating.diploid(nsnp=1000, nindi=100, n.additive = 10)
#' population <- breeding.diploid(population, breeding.size=100, selection.size=c(25,25),
#' ocs = TRUE)
#' @return contributions of each individual in selection
#' @export

ocs.mobps = function(population, animallist,
                     relationship.matrix.ocs,
                     depth.pedigree.ocs,
                     bve.pedigree.error,
                     ocs.target = "min.sKin",
                     ocs.uniform = NULL,
                     ocs.lb = NULL,
                     ocs.ub = NULL,
                     ocs.ub.sKin = NULL,
                     ocs.lb.BV = NULL,
                     ocs.ub.BV = NULL,
                     ocs.eq.BV = NULL,
                     ocs.ub.sKin.increase = NULL,
                     ocs.lb.BV.increase = NULL,
                     bve.p_i.list = NULL,
                     miraculix = FALSE,
                     miraculix.mult = FALSE,
                     import.position.calculation = NULL,
                     decodeOriginsU = decodeOriginsR,
                     store.sparse = FALSE,
                     verbose = TRUE
                     ){


  n.animals <- nrow(animallist)

  if(sum(abs(animallist[,4]))==0){
    if(verbose) cat("No breeding values stored for OSC Use breeding value of first trait as u!\n")
    for(index in 1:nrow(animallist)){
      animallist[index,4] <- population$breeding[[animallist[index,1]]][[2+animallist[index,2]]][1,animallist[index,3]]
    }
  }

  if(sum(abs(animallist[,4]))==0){
    if(verbose) cat("No breeding value estimated available! Use genomic value of first trait as u!\n")
    for(index in 1:nrow(animallist)){
      animallist[index,4] <- population$breeding[[animallist[index,1]]][[6+animallist[index,2]]][1,animallist[index,3]]
    }
  }

  BV <- animallist[,4]
  Sex <- rep("male", length(BV))
  Sex[animallist[,6]==2] <- "female"



  {


    if(length(relationship.matrix.ocs)==1 && !relationship.matrix.ocs %in% c("kinship", "pedigree",
                                                                             "emp_kinship", "emp_fast_kinship")){
      if(miraculix){

        tmp111 = check.copy(population, animallist[,1:3])

        Z.code <- miraculix::computeSNPS(population, tmp111[,1], tmp111[,2], tmp111[,3], what="geno", output_compressed = TRUE)
      } else{
        Zt <- array(0,dim=c(sum(population$info$snp), n.animals))
        for(index in 1:n.animals){
          Zt[,index] <- colSums(computing.snps(population, animallist[index,1], animallist[index,2], animallist[index,3], import.position.calculation=import.position.calculation, decodeOriginsU=decodeOriginsU))
        }
      }


    }

    animal_ids = numeric(nrow(animallist))
    for(index in 1:nrow(animallist)){
      animal_ids[index] = population$breeding[[animallist[index,1]]][[animallist[index,2] + 14]][animallist[index,3]]
    }

    if(length(relationship.matrix.ocs) > 1){

      A <- relationship.matrix.ocs

      if(prod(animal_ids %in% colnames(A))==0){
        stop("Not all individuals provided in relationship.matrix.ocs. Indicate individuals by row & column name!")
      }
      if(prod(animal_ids %in% rownames(A))==0){
        stop("Not all individuals provided in relationship.matrix.ocs. Indicate individuals by row & column name!")
      }

    } else if(relationship.matrix.ocs=="kinship" || relationship.matrix.ocs == "pedigree"){
      A <- kinship.exp(population, database = animallist[,c(1,2,3,3)], depth.pedigree = depth.pedigree.ocs,
                       verbose=verbose, mult=2, include.error = bve.pedigree.error)

      if(store.sparse){
        A <- methods::as(A, "sparseMatrix")
      }
    } else if(relationship.matrix.ocs == "emp_kinship"){
      A <- kinship.emp(population, database = animallist[,c(1,2,3,3)]) * 2
      if(store.sparse){
        A <- methods::as(A, "sparseMatrix")
      }
    } else if(relationship.matrix.ocs == "emp_fast_kinship"){
      A <- get.kinship.snp(population, database = animallist[,c(1,2,3,3)]) * 2
      if(store.sparse){
        A <- methods::as(A, "sparseMatrix")
      }
    } else if(relationship.matrix.ocs=="vanRaden"){
      if(miraculix){

        if(length(bve.p_i.list)==0){
          p_i <- rowMeans(as.matrix(Z.code))/2 # Noch nicht implementiert?
          norm = sum(p_i != 0 & p_i != 1)>0
          A <- miraculix::relationshipMatrix(Z.code, centered=TRUE, normalized=norm)
        } else{
          p_i <- bve.p_i.list
          A <- miraculix::relationshipMatrix(Z.code, centered=FALSE, normalized=FALSE)
          A <- scaling.relationship(A, Z.code, p_i)
        }

      } else if(miraculix.mult){

        Zt_miraculix <- miraculix::genomicmatrix(Zt)

        if(length(bve.p_i.list)==0){
          p_i <- rowSums(Zt)/ncol(Zt)/2
          norm = sum(p_i != 0 & p_i != 1)>0
          A <- miraculix::relationshipMatrix(Zt_miraculix, centered=TRUE, normalized=norm)
        } else{
          p_i <- bve.p_i.list
          A <- miraculix::relationshipMatrix(Zt_miraculix, centered=FALSE, normalized=FALSE)
          A <- scaling.relationship(A, Zt_miraculix, p_i)
        }

      } else{

        if(length(bve.p_i.list)==0){
          p_i <- rowSums(Zt)/ncol(Zt)/2
        } else{
          p_i <- bve.p_i.list
        }

        Ztm <- Zt - p_i * 2
        A <- crossprod(Ztm)/ (2 * sum(p_i*(1-p_i)))
      }

    } else if(relationship.matrix.ocs=="CM"){
      #CM SCHAETZER
      Ztm <- rbind(Zt==0, Zt==1, Zt==2)
      A <- crossprod(Ztm) / ncol(Zt)
    } else if(relationship.matrix.ocs=="CE"){
      Ztm <- rbind(Zt==0, Zt==1, Zt==2)
      A <- crossprod(Ztm)
      A <- (A^2 - 0.5*A)/(nrow(Zt)^2)

    } else if(relationship.matrix.ocs=="non_stand"){
      A <- crossprod(Zt) / nrow(Zt)
    }


    if(length(relationship.matrix.ocs)==1 && !relationship.matrix.ocs %in% c("kinship", "pedigree",
                                                                             "emp_kinship", "emp_fast_kinship") ){
      colnames(A) = rownames(A) = animal_ids
    } else{

      # order is already correct when not using pedigree-matrix

      id_A = numeric(length(animal_ids))
      for(index in 1:length(animal_ids)){
        id_A[index] = get.id(population, database = animallist[index,c(1,2,3,3), drop = FALSE])
      }
      A = A[as.character(id_A), as.character(id_A)]

    }

  }



  Indiv <- paste0("Indi", 1:length(BV))

  for(index in 1:nrow(animallist)){
    Indiv[index] = population$breeding[[animallist[index,1]]][[animallist[index,2] + 14]][animallist[index,3]]
  }


  if(sum(duplicated(Indiv)) > 0){

    if(verbose){
      cat("Individuals used as both first and second parent. Ignore sexes in OCS")
    }
    BV = BV[!duplicated(Indiv)]
    Indiv = Indiv[!duplicated(Indiv)]

    Sex = rep(NA, length(Indiv))
  }

  colnames(A) <- rownames(A) <- names(BV) <- Indiv
  cont <- cbind(1,0.5,0.5)
  colnames(cont) <- c("age", "male", "female")
  Born <- rep(1, length(BV))
  Breed <- rep("Breed1", length(BV))
  herd <- rep(NA, length(BV))
  isCandidate <- rep(TRUE, length(BV))
  phen <- data.frame(Indiv, Born, Breed, BV, Sex, herd, isCandidate)

  sKin <- A/2
  colnames(sKin) <- rownames(sKin) <- Indiv
  cand <- optiSel::candes(phen = phen, sKin = sKin, cont = cont)
  con <- list()
  if(length(ocs.uniform)>0){
    con$uniform = ocs.uniform
  }
  if(length(ocs.lb)>0){
    con$lb = ocs.lb
  }
  if(length(ocs.ub)>0){


    ocs.ub.vector = rep(ocs.ub[1]/2, length(Indiv))

    if(length(ocs.ub)==2){
      ocs.ub.vector[animallist[,6]==2] = ocs.ub[2]/2
    }
    names(ocs.ub.vector) = Indiv

    if(length(ocs.uniform) ==1 && ocs.uniform == "female"){
      ocs.ub.vector = ocs.ub.vector[!(animallist[,6]==2)]
    }

    if(length(ocs.uniform) ==1 && ocs.uniform == "male"){
      ocs.ub.vector = ocs.ub.vector[!(animallist[,6]==1)]
    }

    con$ub = ocs.ub.vector
  }
  if(length(ocs.ub.sKin)>0){
    con$ub.sKin = ocs.ub.sKin
  }
  if(length(ocs.lb.BV)>0){
    con$lb.BV = ocs.lb.BV
  }
  if(length(ocs.ub.BV)>0){
    con$ub.BV = ocs.ub.BV
  }
  if(length(ocs.eq.BV)>0){
    con$eq.BV = ocs.eq.BV
  }



  if(length(ocs.ub.sKin.increase)>0){
    con$ub.sKin = ocs.ub.sKin.increase + cand$current[2,4]
  }
  if(length(ocs.lb.BV.increase)>0){
    con$lb.BV = ocs.lb.BV.increase + cand$current[1,4]
  }

  Offspring <- optiSel::opticont(ocs.target, cand, con, quiet = !verbose)

  if(sum(is.na(Sex))==length(Sex)){
    contribution <- list(Offspring$parent$oc, Offspring$parent$oc)
  } else{
    contribution <- list(Offspring$parent$oc[Sex=="male"], Offspring$parent$oc[Sex=="female"])

  }

  # optiSel seem to sometimes output negative contributions.
  contribution[[1]][contribution[[1]] < 0] = 0
  contribution[[2]][contribution[[2]] < 0] = 0

  if(verbose){
    cat(paste0(sum(contribution[[1]]>0), " male individuals with positive contribution ((", sum(contribution[[1]]>(0.001 * max(contribution[[1]]))), " with major contribution).\n"))
    cat(paste0(sum(contribution[[2]]>0), " female individuals with positive contribution ((", sum(contribution[[2]]>(0.001 * max(contribution[[2]]))), " with major contribution)\n."))
  }

  return(contribution)
}
