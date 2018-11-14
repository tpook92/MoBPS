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

#' THIS NEEDs SOME REWORK BEFORE USAGE
#'
#' THIS NEEDs SOME REWORK BEFORE USAGE
#' @param A A
#' @param u u
#' @param Q Q
#' @param cAc cAc

OGC <- function(A,u,Q,cAc=NA){                # bei NA wird cAc (Increase of average relationship) minimiert

  Opt.int <- function(A, u, Q, cAc=NA){

    A1 <- solve(A)
    QAQ1 <- solve(t(Q)%*%A1%*%Q)
    minA <- round(0.25*sum(QAQ1)+0.000005,5)
    if(!is.na(cAc)){
      if(minA>cAc){stop(paste("Minimal increase possible is for cAc is ",minA,sep=""))}
    }else{
      cAc <- minA
    }
    lambda0 <- as.numeric(sqrt((t(u)%*%(A1-A1%*%Q%*%QAQ1%*%t(Q)%*%A1)%*%u)/(4*cAc-sum(QAQ1))) )
    lambda <- QAQ1%*%(t(Q)%*%A1%*%u - lambda0)

    xopt <- A1%*%(u-Q%*%lambda)/(2*lambda0)
    return(list(xopt,cAc))
  }

  rr <- Opt.int(A,u,Q,cAc)
  xopt <- rr[[1]]
  minA <- rr[[2]]
  n <- length(u)
  names(u) <- 1:n
  for(i in 1:length(xopt)){
    if((min(xopt) >= 0)){break}
    neg <- which(xopt < 0)
    A <- A[-neg,-neg]
    u <- u[-neg]
    Q <- Q[-neg,]
    if(length(u)==2){xopt=rep(0.5,2);break}
    xopt <- Opt.int(A,u,Q,cAc)[[1]]
  }
  res <- numeric(n)
  res[as.numeric(names(u))] <- xopt
  print(paste("Minimal increase is ",minA,sep=""))
  return(list("Optimal c"=res, "c'u"=t(xopt)%*%u))
}
