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

#' Calculate allele frequcenies of a miraculix Z matrix without miraculix::allele_freq()
#'
#' Calculate allele frequcenies of a miraculix Z matrix without miraculix::allele_freq()
#' @param Z.code bit-wise coded genotype matrix

allele.freq <- function(Z.code){
  if (requireNamespace("miraculix", quietly = TRUE)) {
    rowMeans(miraculix::decodeGeno(Z.code))/2
  } else{
    stop("Usage of miraculix without being installed!")
  }
}
