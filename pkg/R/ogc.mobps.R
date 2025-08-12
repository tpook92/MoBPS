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
#' @param animallist List of individuals to include in ogc
#' @param ogc.target Target of OGC (default: "min.sKin" - minimize inbreeding; alt: "max.BV" / "min.BV" - maximize genetic gain; both under constrains selected below)
#' @param ogc.uniform This corresponds to the uniform constrain in optiSel
#' @param ogc.lb This corresponds to the lb constrain in optiSel
#' @param ogc.ub This corresponds to the ub constrain in optiSel
#' @param ogc.ub.sKin This corresponds to the ub.sKin constrain in optiSel
#' @param ogc.lb.BV This corresponds to the lb.BV constrain in optiSel
#' @param ogc.ub.BV This corresponds to the ub.BV constrain in optiSel
#' @param ogc.eq.BV This corresponds to the eq.BV constrain in optiSel
#' @param ogc.ub.sKin.increase This corresponds to the upper bound (current sKin + ogc.ub.sKin.increase) as ub.sKin in optiSel
#' @param ogc.lb.BV.increase This corresponds to the lower bound (current BV + ogc.lb.BV.increase) as lb.BV in optiSel
#' @param relationship.matrix.ogc Method to calculate relationship matrix for OGC (Default: "pedigree", alt: "vanRaden", "CE", "non_stand", "CE2", "CM")
#' @param depth.pedigree.ogc Depth of the pedigree in generations (default: 7)
#' @param bve.pedigree.error Set to FALSE to ignore/correct for any pedigree errors
#' @param import.position.calculation XXX
#' @param decodeOriginsU XXX
#' @param nbits XXX
#' @param store.sparse XXX
#' @param miraculix XXX
#' @param miraculix.mult XXX
#' @param bve.p_i.list XXX
#' @param verbose Set to FALSE to not display any prints
#' @param bit.storing Set to TRUE if the MoBPS (not-miraculix! bit-storing is used)
#' @examples
#' population <- creating.diploid(nsnp=1000, nindi=100, n.additive = 10)
#' population <- breeding.diploid(population, breeding.size=100, selection.size=c(25,25),
#' ogc = TRUE)
#' @return contributions of each individual in selection
#' @export

ogc.mobps = function(population, animallist,
                     relationship.matrix.ogc,
                     depth.pedigree.ogc,
                     bve.pedigree.error,
                     ogc.target = "min.sKin",
                     ogc.uniform = NULL,
                     ogc.lb = NULL,
                     ogc.ub = NULL,
                     ogc.ub.sKin = NULL,
                     ogc.lb.BV = NULL,
                     ogc.ub.BV = NULL,
                     ogc.eq.BV = NULL,
                     ogc.ub.sKin.increase = NULL,
                     ogc.lb.BV.increase = NULL,
                     bve.p_i.list = NULL,
                     miraculix = FALSE,
                     miraculix.mult = FALSE,
                     import.position.calculation = NULL,
                     decodeOriginsU = decodeOriginsR,
                     nbits = NULL,
                     store.sparse = FALSE,
                     verbose = TRUE,
                     bit.storing = FALSE
                     ){


  n.animals <- nrow(animallist)

  if(sum(abs(animallist[,4]))==0){
    if(verbose) cat("No breeding values stored for OGC. Use breeding value of first trait as u!\n")
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

  if(relationship.matrix.ogc != "kinship" || relationship.matrix.ogc != "pedigree"){
    if(miraculix){
      Z.code <- miraculix::computeSNPS(population, animallist[,1], animallist[,2], animallist[,3], what="geno", output_compressed = TRUE)
    } else{
      Zt <- array(0,dim=c(sum(population$info$snp), n.animals))
      for(index in 1:n.animals){
        Zt[,index] <- colSums(computing.snps(population, animallist[index,1], animallist[index,2], animallist[index,3], import.position.calculation=import.position.calculation, decodeOriginsU=decodeOriginsU, bit.storing=bit.storing, nbits=nbits, output_compressed=FALSE))
      }
    }
  }

  # Verwandtschaftsmatrix:
  if(relationship.matrix.ogc=="kinship" || relationship.matrix.ogc == "pedigree"){
    A <- kinship.exp(population, database = animallist[,c(1,2,3,3)], depth.pedigree = depth.pedigree.ogc,
                     verbose=verbose, mult=2, include.error = bve.pedigree.error)

    if(store.sparse){
      A <- methods::as(A, "sparseMatrix")
    }
  } else if(relationship.matrix.ogc=="vanRaden"){
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

  } else if(relationship.matrix.ogc=="CM"){
    #CM SCHAETZER
    Ztm <- rbind(Zt==0, Zt==1, Zt==2)
    A <- crossprod(Ztm) / ncol(Zt)
  } else if(relationship.matrix.ogc=="CE"){
    Ztm <- rbind(Zt==0, Zt==1, Zt==2)
    A <- crossprod(Ztm)
    A <- (A^2 - 0.5*A)/(nrow(Zt)^2)

  } else if(relationship.matrix.ogc=="non_stand"){
    A <- crossprod(Zt) / nrow(Zt)
  }

  id_A = numeric(nrow(A))
  for(index in 1:nrow(A)){
    id_A[index] = get.id(population, database = animallist[index,c(1,2,3,3), drop = FALSE])
  }
  A = A[as.character(id_A), as.character(id_A)]


  Indiv <- paste0("Indi", 1:length(BV))

  for(index in 1:nrow(animallist)){
    Indiv[index] = population$breeding[[animallist[index,1]]][[animallist[index,2] + 14]][animallist[index,3]]
  }


  if(sum(duplicated(Indiv)) > 0){

    if(verbose){
      cat("Individuals used as both first and second parent. Ignore sexes in OGC")
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
  if(length(ogc.uniform)>0){
    con$uniform = ogc.uniform
  }
  if(length(ogc.lb)>0){
    con$lb = ogc.lb
  }
  if(length(ogc.ub)>0){


    ogc.ub.vector = rep(ogc.ub[1]/2, length(Indiv))

    if(length(ogc.ub)==2){
      ogc.ub.vector[animallist[,6]==2] = ogc.ub[2]/2
    }
    names(ogc.ub.vector) = Indiv

    if(length(ogc.uniform) ==1 && ogc.uniform == "female"){
      ogc.ub.vector = ogc.ub.vector[!(animallist[,6]==2)]
    }

    if(length(ogc.uniform) ==1 && ogc.uniform == "male"){
      ogc.ub.vector = ogc.ub.vector[!(animallist[,6]==1)]
    }

    con$ub = ogc.ub.vector
  }
  if(length(ogc.ub.sKin)>0){
    con$ub.sKin = ogc.ub.sKin
  }
  if(length(ogc.lb.BV)>0){
    con$lb.BV = ogc.lb.BV
  }
  if(length(ogc.ub.BV)>0){
    con$ub.BV = ogc.ub.BV
  }
  if(length(ogc.eq.BV)>0){
    con$eq.BV = ogc.eq.BV
  }



  if(length(ogc.ub.sKin.increase)>0){
    con$ub.sKin = ogc.ub.sKin.increase + cand$current[2,4]
  }
  if(length(ogc.lb.BV.increase)>0){
    con$lb.BV = ogc.lb.BV.increase + cand$current[1,4]
  }

  Offspring <- optiSel::opticont(ogc.target, cand, con, quiet = !verbose)

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
