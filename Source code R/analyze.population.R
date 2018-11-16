'#
  Authors
Torsten Pook, torsten.pook@uni-goettingen.de

Copyright (C) 2017 -- 2018  Torsten Pook

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

#' Analyze allele frequency of a single marker
#'
#' Analyze allele frequency of a single marker
#' @param population Population list
#' @param chromosome Number of the chromosome of the relevant SNP
#' @param snp Number of the relevant SNP
#' @param include include male/female (1/2) individuals in computation default( c(1,2))
#' @export

analyze.population <- function(population, chromosome, snp, include=c(1,2)){

  if (requireNamespace("miraculix", quietly = TRUE)) {
    codeOriginsU <- miraculix::codeOrigins
    decodeOriginsU <- miraculix::decodeOrigins
  } else{
    codeOriginsU <- codeOriginsR
    decodeOriginsU <- decodeOriginsR
  }
  generations <- length(population$breeding)
  state <- matrix(0, nrow=3, ncol = generations)
  p.snp <- sum(population$info$length[0:(chromosome-1)]) + population$info$position[[chromosome]][snp]
  n.snp <- sum(population$info$snp[0:(chromosome-1)]) + snp
  for(index in 1:generations){
#    active.animals <- NULL
    for(sex in include){
      n.animals <- length(population$breeding[[index]][[sex]])
      for(animal in 1:n.animals){
        position1 <- sum(population$breeding[[index]][[sex]][[animal]][[1]]< p.snp)
        position2 <- sum(population$breeding[[index]][[sex]][[animal]][[2]]< p.snp)
        ursprung1 <- decodeOriginsU(population$breeding[[index]][[sex]][[animal]][[5]],position1)
        ursprung2 <- decodeOriginsU(population$breeding[[index]][[sex]][[animal]][[6]],position2)
        mutation1 <- sum (population$breeding[[index]][[sex]][[animal]][[3]]==p.snp)
        mutation2 <- sum (population$breeding[[index]][[sex]][[animal]][[4]]==p.snp)
        gen1 <- population$breeding[[ursprung1[1]]][[ursprung1[2]]][[ursprung1[[3]]]][[ursprung1[4]+8]][n.snp] + mutation1
        gen2 <- population$breeding[[ursprung2[1]]][[ursprung2[2]]][[ursprung2[[3]]]][[ursprung2[4]+8]][n.snp] + mutation2
        gen1 <- gen1 - 2* (gen1==2)
        gen2 <- gen2 - 2* (gen2==2)
        state[(gen1+gen2+1), index] <- state[(gen1+gen2+1), index] +1
      }
    }
  }
  state.prob <- t(t(state)/colSums(state))
  maxp <- max(state.prob)
  graphics::plot(1:generations ,state.prob[1,],xlim=c(1,generations),ylim=c(0,maxp),type="l",xlab="generation", ylab="share", lwd=3,
                 main=paste("Allele frequencies in Chr", chromosome, "SNP", snp))
  graphics::lines(1:generations ,state.prob[2,],lty=2, lwd=3)
  graphics::lines(1:generations ,state.prob[3,],lty=3, lwd=3)
  graphics::legend("topleft",legend = c("hom0","hetero","hom1"),lty=c(1,2,3), lwd=c(3,3,3))
  return(state)
}
