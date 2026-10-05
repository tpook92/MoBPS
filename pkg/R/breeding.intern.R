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

#' Internal function to simulate one meiosis
#'
#' Internal function to simulate one meiosis
#' @param info.parent position of the parent in the dataset
#' @param parent list of information regarding the parent
#' @param population Population list
#' @param mutation.rate Mutation rate in each marker (default: 10^-5)
#' @param remutation.rate Remutation rate in each marker (default: 10^-5)
#' @param recombination.rate Average number of recombination per 1 length unit (default: 1M)
#' @param recom.f.indicator Use step function for recombination map (transform snp.positions if possible instead)
#' @param duplication.rate Share of recombination points with a duplication (default: 0 - DEACTIVATED)
#' @param duplication.length Average length of a duplication (Exponentially distributed)
#' @param duplication.recombination Average number of recombinations per 1 length uit of duplication (default: 1)
#' @param gene.editing If TRUE perform gene editing on newly generated individual
#' @param gen.architecture Used underlying genetic architecture (genome length in M)
#' @param nr.edits Number of edits to perform per individual
#' @param decodeOriginsU Used function for the decoding of genetic origins [[5]]/[[6]]
#' @param delete.same.origin If TRUE delete recombination points when genetic origin of adjacent segments is the same
#' @param recombination.function Function used to calculate position of recombination events (default: MoBPS::recombination.function.haldane())
#' @param dup_activ Internal parameter to check if duplications have to be simulated
#' @param rt_activ Internal parameter to check if RTs have to be simulated
#' @param grandsib_activ Internal parameter to check if grandsibling contributions have to be calculated
#' @param nsnp internal parameter
#' @param length_total internal parameter
#' @param length_max internal parameter
#' @param recom_prec internal parameter

#' @examples
#' data(ex_pop)
#' child_gamete <- breeding.intern(info.parent = c(1,1,1), parent = ex_pop$breeding[[1]][[1]][[1]],
#'                                 population = ex_pop)
#' @return Inherited parent gamete
#' @export
#'
breeding.intern <- function(info.parent, parent,  population , mutation.rate = 10^-5, remutation.rate = 10^-5, recombination.rate=1,
                             recom.f.indicator=NULL, duplication.rate=0, duplication.length=0.01,
                             duplication.recombination=1, delete.same.origin=FALSE,
                             gene.editing=FALSE, nr.edits= 0,
                             gen.architecture=0,
                             decodeOriginsU=MoBPS::decodeOriginsR,
                             recombination.function=MoBPS::recombination.function.haldane,
                             dup_activ = FALSE,
                             rt_activ = FALSE,
                             grandsib_activ = FALSE,
                            nsnp = 0,
                            length_total = 0,
                            length_max = 0,
                            recom_prec = 1){


  new_rec = TRUE ## setting this to FALSE will revert to use the old R based implementation instead of C++ based

  parent_rec1 = parent[[1]]
  parent_rec2 = parent[[2]]
  parent_mut1 = parent[[3]]
  parent_mut2 = parent[[4]]
  parent_ori1 = parent[[5]]
  parent_ori2 = parent[[6]]

  if(nsnp==0){
    nsnp <- sum(population$info$snp)
    length_total = population$info$length.total
    length_max = max(length_total)
    recom_prec = population$info$recombination.precision
  }
  if(gen.architecture==0){
    length.total <- length_total
  } else{
    length.total <- population$info$gen.architecture[[gen.architecture]]$length.total

    # Transform points of recombination according to position
    for(haplo in 1:2){
      for(index in 1:length(parent[[haplo]])){
        #        chromo <- find.chromo(parent[[haplo]][index], population$info$length.total)
        before <- find.snpbefore(parent[[haplo]][index], population$info$snp.position)
        if(before>0){
          p_before <- population$info$snp.position[before]
          new_p_before <- population$info$gen.architecture[[gen.architecture]]$snp.position[before]
        } else{
          p_before <- length_total[1]
          new_p_before <-length.total[1]
        }
        if(before<nsnp){
          p_after <- population$info$snp.position[before+1]
          new_p_after <- population$info$gen.architecture[[gen.architecture]]$snp.position[before+1]
        } else{
          p_after <- length_total[length(length_total)]
          new_p_after <- length.total[length(length.total)]
        }

        share <- (parent[[haplo]][index]-p_before) / (p_after-p_before)
        parent[[haplo]][index] <- new_p_before + share * (new_p_after-new_p_before)
      }
    }
  }

  n.chromosome <- length(length.total)-1

  # Ausserhalb von breeding.intern berechnen
  if(length(recom.f.indicator)!=0){
    # Fuer Polynom Numerische Bestimmung anstrengend?
    recom.f.indicator <- rbind(recom.f.indicator, c(length_max,0))
    indicator.vol <- sum((recom.f.indicator[-1,1] -recom.f.indicator[-nrow(recom.f.indicator),1])*recom.f.indicator[-nrow(recom.f.indicator),2])
    recom.vol <- indicator.vol #+ polynom.vol
  } else{
    recom.vol <- length_max*recombination.rate
  }

  valid = FALSE
  is_carrier = NULL
  rt = numeric(0)
  while(!valid){

    rt = numeric(0)
    valid = TRUE

    noc <- stats::rpois(1, recom.vol) #Anzahl Rekombinationspunkte
    if(length(recom.f.indicator)!=0){
      porc <- stats::runif(noc,0,recom.vol)
      if(length(porc)>0){
        cums <- cumsum((recom.f.indicator[-1,1] -recom.f.indicator[-nrow(recom.f.indicator),1])*recom.f.indicator[-nrow(recom.f.indicator),2])
        for(index in 1:length(porc)){
          actuel <- porc[index]
          before <- sum(actuel < cums)
          prev1 <- recom.f.indicator[nrow(recom.f.indicator)-before,]
          next1 <- recom.f.indicator[nrow(recom.f.indicator)-before+1,]

          porc[index] <- prev1[1] + (actuel-c(0,cums)[nrow(recom.f.indicator)-before]) /prev1[2]
        }
      }
    } else{
      porc <- recombination.function(noc, length_max) #Position der Rekombinationspunkte
    }


    if(recom_prec > 1){
      porc_temp = as.integer(porc * recom_prec)
    } else{
      porc_temp = porc

    }


    #start.point <- stats::rbinom(n.chromosome,1,0.5) * (1:n.chromosome) #Wechsel zu Beginn des Chromosoms
    start_tmp = length.total[which(stats::runif(n.chromosome) < 0.5)]

    if(dup_activ){
      if(duplication.rate>0){
        porc.d <- sortd(c(porc, 0, length_max))
        rpod <- c(0, (stats::rbinom(noc,1,duplication.rate) * 2:(noc+1)),0)
        pod <- porc.d[rpod]
        pod2 <- rep(0,length(pod))

        count <- 1
        add.one <- rep(0,length(pod))


        for(index in unique(c(0,rpod))[-1]){
          activ.chromosome <- sum(pod[count] > length.total)
          length.d <- stats::rexp(1, rate=(1/duplication.length)) * (-1)^(stats::rbinom(1,1,0.5))
          pod2[count] <- max(min(pod[count]+ length.d, porc.d[index+1], length.total[activ.chromosome+1]),porc.d[index-1], length.total[activ.chromosome])
          if(pod2[count]==(porc.d[index])|| sum(start_tmp==porc.d[index])){
            add.one[count] <- 1
          }
          count <- count+1
        }
        pod.start <- pod
        pod.start[pod.start>pod2] <- pod2[pod.start>pod2]
        pod.end <- pod
        pod.end[pod.end<pod2] <- pod2[pod.end<pod2]

      } else{
        pod <- pod2 <- rpod <- pod.start <- pod.end <- numeric(0)

      }
    }




    if(new_rec && !is.unsorted(porc)){

      porc <- merge_porc(
        start_tmp,
        porc,
        length_max
      )

    } else{
      porc <- sortd((c(start_tmp,porc))) # Sortieren der Rekombinationspunkte
      porc <- c(-1,porc,length_max+1)

    }

    if(recom_prec >1){
      porc = as.integer(porc * recom_prec)
    }

    if(rt_activ && (length(parent[[35]])>0 || length(parent[[36]])>0)){

      # add additional recombination events at RT positions
      activ_rt = c(parent[[35]], parent[[36]])
      activ_rt_indi = as.numeric(length(parent[[36]])>0)

      if(length(activ_rt)>6){
        stop("No matings between individuals with 2 rts defined")
      }


      if(length(is_carrier)==0){
        is_carrier = (sum(porc < activ_rt[1])+1)%%2 == (activ_rt_indi)
      } else{
        is_carrier_check = (sum(porc < activ_rt[1])+1)%%2 == (activ_rt_indi)

        if(is_carrier != is_carrier_check){
          valid = FALSE
        }
      }

      good_phase = sum(porc > activ_rt[3] & porc < activ_rt[4])%%2==1 &&       sum(porc > activ_rt[5] & porc < activ_rt[6])%%2==1


      if(!good_phase){
        valid = FALSE
      }
      # check if RTs are inherited jointly
      if(is_carrier){
        rt = activ_rt
      }

    }

  }


  if(dup_activ){
    # recombination on dupliciation sequences
    dup <- list()
    dup[[1]] <- (parent[[11]])
    dup[[2]] <- (parent[[12]])
    dup[[3]] <- "test"

    for(abc in 1:2){
      counter <- 0
      if(length(dup[[abc]])>0){
        posi <- 1:nrow(dup[[abc]])
        for(index in 1:nrow(dup[[abc]])){
          ndup <- stats::rpois(1, (dup[[abc]][index,3]- dup[[abc]][index,2])* duplication.recombination)
          if(ndup>0){
            pdup <- sort(stats::runif(ndup, min=dup[[abc]][index,2], max= dup[[abc]][index,3]))
          } else{
            pdup <- numeric(0)
          }
          if(counter==1 && dup[[abc]][index,1] == dup[[abc]][(index-1),1]){
            pdup <- c(dup[[abc]][index,2], pdup)
          }
          counter <- 0
          if(ndup%%2 && (sum(dup[[abc]][index,1]>=porc)%%2)==(abc%%2)){ # recombination on a relevant section
            porc <- sort(c(porc, dup[[abc]][index,1]))
          }
          if(ndup>0){
            if(ndup%%2==0){ #garanty odd number of elements
              pdup <- c(pdup, dup[[abc]][index,3])
              counter <- 1
            }
            dup[[abc]][index,3] <- pdup[1]
            if(ndup>1){
              for(index2 in seq(2,ndup,by=2)){
                dup[[abc]] <- rbind(dup[[abc]], dup[[abc]][index,])
                dup[[abc]][nrow(dup[[abc]]),2:3] <- pdup[index2:(index2+1)] # adding additional duplication segment
                posi <- c(posi,index)
              }
            }
          }
        }
        posi <- sort(posi,index.return=TRUE)$ix
        dup[[abc]] <- matrix(dup[[abc]][posi,], ncol=8)
        #remove duplication segements with length 0
        remove0 <- (dup[[abc]][,2] == dup[[abc]][,3])
        if(sum(remove0)>0){
          remove0 <- remove0 * 1:length(remove0)
          dup[[abc]] <- dup[[abc]][-remove0,]
        }
      }
    }
  }



  #
  new.poc <- NULL
  new.mut <- integer(0)
  new.origin <- NULL
  new.dup <- NULL
  activ <- 1

  if(length(parent_mut1)>0 || length(parent_mut2)>0 ){
    store_mut <- list(population$info$snp.position[parent_mut1], population$info$snp.position[parent_mut2])
    mut_bool <- c(length(store_mut[[1]]), length(store_mut[[2]]))==0
    if(dup_activ){
      dup_bool <- c(length(dup[[1]]), length(dup[[2]]))>0
    } else{
      dup_bool <- c(FALSE, FALSE)
    }

  } else{
    store_mut = list()
    mut_bool = c(TRUE,TRUE)
    dup_bool = c(FALSE, FALSE)
  }

  if(new_rec){

    fi <- find_intervals_joint(
      porc,
      parent_rec1,
      parent_rec2
    )

    start1 <- fi$start1
    start2 <- fi$start2
    end1   <- fi$end1
    end2   <- fi$end2

  } else{
    start1 = findInterval(porc, parent_rec1)
    start2 =  findInterval(porc, parent_rec2)
    end1 = findInterval(porc, parent_rec1, left.open = TRUE)
    end2 = findInterval(porc, parent_rec2, left.open = TRUE)
  }


  tmp1 = -2
  tmp2 = logical(0)
  if(sum(dup_bool) > 0){
    activ <- 1

    start_list = list(start1, start2)
    end_list = list(end1, end2)
    for(index in 1:(length(porc)-1)){
      if(start_list[[activ]][index] == end_list[[activ]][index+1]){
        activ.porc = NULL
      } else{
        activ.porc = (start_list[[activ]][index] + 1):end_list[[activ]][index+1]
      }
      if(dup_bool[activ]){
        activ.dup1 <- (dup[[activ]][,1] <= porc[index+1]) * (dup[[activ]][,1] >= porc[index]) * ((dup[[activ]][,8] == 1) + (dup[[activ]][,8] == 3))
        # duplication section on the start/end of a chromosome
        activ.dup2 <- (dup[[activ]][,1] < porc[index+1]) * (dup[[activ]][,1] > porc[index]) * (dup[[activ]][,8] == 2)
        # duplication section in the middle of a chromosome

        # if the adjacting porcs are because of a start.point recombination some duplication sections have to be excluded
        leftb <- sum(porc[index]==porc) - sum(porc[index] == start_tmp)
        rightb <- sum(porc[index+1]==porc) - sum(porc[index+1] == start_tmp)

        activ.dup3 <- (1-leftb) * (dup[[activ]][,1] <= porc[index+1]) * (dup[[activ]][,1] == porc[index]) * ((dup[[activ]][,8] == 1) + (dup[[activ]][,8] == 3))
        activ.dup4 <- (1-rightb) * (dup[[activ]][,1] == porc[index+1]) * (dup[[activ]][,1] >= porc[index]) * ((dup[[activ]][,8] == 1) + (dup[[activ]][,8] == 3))

        activ.dup <- activ.dup1 + activ.dup2 - ((activ.dup3+ activ.dup4)>0)
        new.dup <- rbind(new.dup, dup[[activ]][activ.dup * (1:length(activ.dup)),])
      } else{
        active.dup <- tmp2
        new.dup <- NULL
      }
      activ <- 3 - activ

    }
  }

  activ <- 1


  if(new_rec){

    res <- breeding_loop(
      porc = porc,
      start1 = start1,
      start2 = start2,
      end1 = end1,
      end2 = end2,
      parent_rec1 = parent_rec1,
      parent_rec2 = parent_rec2,
      parent_ori1 = parent_ori1,
      parent_ori2 = parent_ori2,
      mut_bool = mut_bool,
      store_mut = store_mut,
      parent_mut1 = if (is.null(parent_mut1)) integer(0) else parent_mut1,
      parent_mut2 = if (is.null(parent_mut2)) integer(0) else parent_mut2,
      recombination_precision = recom_prec,
      tmp1 = -2
    )

    new.poc    <- res$new.poc
    new.mut    <- res$new.mut
    new.origin <- res$new.origin
    tmp1       <- res$tmp1
  } else{
    for(index in 1:(length(porc)-1)){

      if(activ==1){
        tmps = start1[index]
        tmpe = end1[index+1]
        rec_activ = parent_rec1
        ori_activ = parent_ori1
      } else{
        tmps = start2[index]
        tmpe = end2[index+1]
        rec_activ = parent_rec2
        ori_activ = parent_ori2
      }

      if(tmps == tmpe){
        activ.porc = NULL
      } else{
        activ.porc = (tmps + 1):tmpe
      }

      if(mut_bool[activ]){
        activ.mut <- tmp2
      } else{
        activ.mut <- (store_mut[[activ]]<porc[index+1]) & (store_mut[[activ]]>porc[index])
      }


      tmp8 = porc[index+1]
      save1 <- tmp1!=tmp8

      if(save1){
        tmp1 = tmp8
        new.poc <- c(new.poc, rec_activ[activ.porc], tmp1)
      }

      if(length(activ.mut)>0){
        new.mut <- c(new.mut, parent[[activ+2]][activ.mut])
      }

      start <- tmps
      tmp3 = max(start, tmpe)

      if(start == 0){
        start = 1
      }

      tmp3 = min(length(ori_activ), tmp3)

      if(start <= tmp3){
        activ.origin <- start:tmp3
      } else{
        activ.origin = numeric(0)
        tmp3 = -1
      }

      if(porc[index+1]>0 && save1){
        new.origin <- c(new.origin, ori_activ[activ.origin])
      }

      activ <- 3 - activ
    }
    new.poc <- new.poc[-length(new.poc)]
  }


  # mutation process
  n_mut <- stats::rbinom(1,nsnp, mutation.rate)
  if(length(new.mut)==0){
    n_remut <- 0
  } else{
    n_remut <- stats::rbinom(1, length(new.mut), remutation.rate) ############### BOTH CHECKER
  }

  if(n_mut==0){
    mutationen <- integer(0)
  } else{
    mutationen <- sample(1:nsnp, n_mut)
  }

  if(length(new.mut)>0){
    checker <- duplicated(c(new.mut, mutationen))
    mutationen <- mutationen[(1:length(checker))[-(1:length(new.mut))]-length(new.mut)]
  }

  if(n_remut==0){
    remutationen <- integer(0)
  } else{
    remutationen <- sample(1:length(new.mut), n_remut)
  }

  if(length(remutationen)==0){
    tmp6 = c(mutationen, new.mut)
    if(length(tmp6)>0){
      new.mut <- sort(unique(tmp6))
    }

  } else{
    new.mut <- c(mutationen, new.mut[-remutationen])
    if(length(new.mut)>0){
      new.mut <- sort(new.mut)
      new.mut <- new.mut[!duplicated(new.mut)]
    }

  }



  if(dup_activ){
    new.dups <- NULL

    if(length(pod)>0){
      new.dups <- cbind(0,pod.start, pod.end, info.parent[1],info.parent[2],info.parent[3], 0, 0)
      # calculation of the position on the chromosome and the origin of the sequence
      for(index in 1:length(pod)){

        position.options <- c(length.total[sum(length.total<pod.start[index])], pod.start[index], length.total[sum(length.total<pod.start[index])+ 1])#start, mid, end
        samp<- sample(1:3,1)
        new.dups[index,1] <- position.options[samp]
        activ.chromo <- (sum(porc <= pod.start[index])+1+ add.one[index])%%2 +1
        new.dups[index,7] <- activ.chromo
        new.dups[index,8] <- samp
      }
    }

    # Duplication process

    new.dup <- rbind(new.dup, new.dups)
    if(length(new.dup)>0){
      order <- sort(new.dup[,1],index.return=TRUE)$ix
      new.dup <- new.dup[order,]
    }
  }


  new.origin_old <- new.origin
  new.poc_old <- new.poc
  if(delete.same.origin==TRUE && length(new.origin)>1){
    for(index in length(new.origin):2){
      check <- prod(new.origin[index] == new.origin[index-1])
      if(check==TRUE){
        new.origin <- new.origin[-index]
        new.poc <- new.poc[-index]
      }
    }
  }

  if(gen.architecture!=0){
    # Transform points of recombination according to position
    for(index in 1:length(new.poc)){
      #        chromo <- find.chromo(parent[[haplo]][index], population$info$length.total)
      before <- find.snpbefore(new.poc[index], population$info$gen.architecture[[gen.architecture]]$snp.position)
      if(before>0){
        new_p_before <- population$info$snp.position[before]
        p_before <- population$info$gen.architecture[[gen.architecture]]$snp.position[before]
      } else{
        new_p_before <- length_total[1]
        p_before <-length.total[1]
      }
      if(before<nsnp){
        new_p_after <- population$info$snp.position[before+1]
        p_after <- population$info$gen.architecture[[gen.architecture]]$snp.position[before+1]
      } else{
        new_p_after <- length_max
        p_after <- length.total[length(length.total)]
      }
      share <- (new.poc[index]-p_before) / (p_after-p_before)
      new.poc[index] <- new_p_before + share * (new_p_after-new_p_before)
    }
  }


  if(length(new.poc)!=(length(new.origin)+1)){

    print(porc)
    print(new.poc)
    print(new.origin)

    stop("recombination inconsistency!")

  }

  if(gene.editing){
    hap_sequence <- computing.snps_single(population, new.poc, new.mut, new.origin, decodeOriginsU=decodeOriginsU)
    ed_info <- population$info$editing_info[[length(population$info$editing_info)]]
    edits_p <- numeric(nr.edits)
    edits <- 0
    current_p <- 1
    max_p <- sum(population$info$snp)
    while(edits < nr.edits && current_p <= max_p){
      if(hap_sequence[ed_info[current_p,1]] == ed_info[current_p,2]){
        edits <- edits + 1
        edits_p[edits] <- ed_info[current_p,1]
      }
      current_p <- current_p +1
    }
    changes <- edits_p

    new.mut <- sort(c(new.mut, edits_p))
  }

  if(grandsib_activ){
    maxl <- max(length.total)
    if(length(porc)==1){
      segment_length <- diff(c(0,porc[-1], maxl))
    } else{
      segment_length <- diff(c(0,porc[-c(1,length(porc))], maxl))
    }

    share_a <- sum(segment_length[(1:length(segment_length)%%2)==1])/maxl
  } else{
    share_a = NULL
  }

  return(list(new.poc, new.mut, new.origin, info.parent, new.dup, porc, share_a, porc_temp, rt))
}





