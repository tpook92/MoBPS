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

#' Function to analyze memory and CPU use
#'
#' Function to analyze  memory and CPU use of a generation population
#' @param population Population list
#' @examples
#' data(ex_pop)
#' get.performance.plot(ex_pop)
#' @return data frame with performance indicators of the job
#' @export
#'


get.performance.plot = function(population){

  if(population$info$tracking){

    tmp = population$info$tracking_stats

    oldpar <- graphics::par(no.readonly=TRUE)
    on.exit(graphics::par(oldpar))
    graphics::par(mfrow = c(1,2))


    minutes = 1:nrow(tmp)*population$info$tracking_interval/60
    plot(minutes, tmp[,7], main = "CPU efficiency", type = "l", ylim = c(0,max(100, unlist(tmp[,6:7])))*1.3,
         ylab = "efficency (in %)", xlab = "time (in minutes)")

    graphics::points(minutes, tmp[,6], col = 2)

    graphics::legend("topright",
           c("Average CPU",
             "Current CPU"),
           pch = 1,
           col = c(1,2),
           lwd = 3,
           lty = 0)
    graphics::plot(minutes, tmp[,2]/1024^2, main = "memory utilization",
         ylab = "memory usage (in GB)",
         xlab = "time (in minutes)",
         ylim = c(0, max(tmp[,c(2,3,4,5)])/1024^2)*1.7)

    graphics::points(minutes, tmp[,3]/1024^2, col = 2)


    graphics::points(max(minutes), utils::object.size(population)/1024^3,
           col = 3, lwd = 5)
    graphics::points(minutes, tmp[,4]/1024^2, col = 4)
    graphics::points(minutes, tmp[,5]/1024^2, col = 5)
    graphics::legend("topleft",
           c("Peak (VmHWM)",
             "Current (VmRSS)",
             "Peak (excluding MiXBLUP/child jobs)",
             "Current  (excluding MiXBLUP/child jobs)",
             "Current Population"),
           col = c(2,1,5,4,3),
           pch = 1,
           lwd =3,
           lty = 0)

    return(tmp)
  } else{
    warning("No performance tracking performed to plot")
  }

}
