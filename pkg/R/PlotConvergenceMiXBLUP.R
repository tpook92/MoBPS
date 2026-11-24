'#

This function was taken from the R-package TuttiFrutti https://github.com/bonifazi/TuttiFrutti/tree/main
that is authored and developed by Renzo Bonifazi

  Authors
  Torsten Pook, torsten.pook@wur.nl
  Renzo Bonifazi, renzo.bonifazi@wur.nl

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

#' Plot MiXBLUP model convergence
#'
#' Plot MiXBLUP model convergence
#' @param ConvFilePath (character) pathway where convergence file is located.
#' @param solver (character) Either "solver" (default) or "hpblup. "The default solver is MiXBLUP default solver,i.e. solver.exe with associated convergence file
#' named "Conlog.txt". Solver can be set to "hpblup" for hpblup.exe which has associated convergence file
#' named "convergence.dat"
#' @param title (character) title for convergence plot. Default is NA.
#' @returns ggplot2 plot object

PlotConvergeneMiXBLUP <- function(ConvFilePath, solver = "solver", title = NA) {

  if (solver == "hpblup") {
    ConvFileExt <- "/convergence.dat"

    data <- utils::read.delim(
      base::paste0(ConvFilePath, ConvFileExt),
      sep = ""
    )

    data <- data[, c(2, 5, 8:10, 12)]
    base::colnames(data)[1:4] <- c("Iteration", "tau", "CR^2", "CD^2")

    col_vec <- c(
      "log(CR^2)" = "purple",
      "log(CD^2)" = "forestgreen",
      "log(CK)" = "dodgerblue",
      "log(CM)" = "orange",
      "log(tau)" = "gold"
    )

  } else if (solver == "solver") {
    ConvFileExt <- "/Conlog.txt"

    data <- utils::read.table(
      base::paste0(ConvFilePath, ConvFileExt),
      skip = 10,
      fill = TRUE
    )

    data <- data[, 1:5]
    base::colnames(data) <- c("Iteration", "CA", "CR", "CM", "CD")

    col_vec <- c(
      "log(CR)" = "purple",
      "log(CD)" = "forestgreen",
      "log(CA)" = "dodgerblue",
      "log(CM)" = "orange"
    )

  } else {
    base::stop("Unknown solver: ", solver)
  }

  # Convert convergence criteria to log10
  criteria <- base::colnames(data)[-1]

  for (criterion in criteria) {
    data[[base::paste0("log(", criterion, ")")]] <-
      base::log10(data[[criterion]])
  }

  # Select transformed convergence criteria
  criteria_log <- base::paste0("log(", criteria, ")")
  criteria_log <- criteria_log[
    criteria_log %in% base::names(col_vec)
  ]

  # Plot convergence without interrupting execution on graphical errors
  base::tryCatch({

    # Save graphical parameters
    oldpar <- graphics::par(no.readonly = TRUE)

    base::on.exit(
      base::try(graphics::par(oldpar), silent = TRUE),
      add = TRUE
    )

    # Increase bottom margin for legend
    graphics::par(mar = c(6, 4, 4, 2) + 0.1)

    # Determine plot ranges
    xlim <- base::range(data$Iteration, finite = TRUE)

    ylim <- base::range(
      base::unlist(data[, criteria_log, drop = FALSE]),
      finite = TRUE
    )

    # Initialize plot
    graphics::plot(
      NA,
      xlim = xlim,
      ylim = ylim,
      xlab = "Iterations",
      ylab = "Convergence",
      main = if (base::is.na(title)[1]) "" else title,
      col.axis = "black",
      col.lab = "black",
      cex.axis = 0.85,
      cex.lab = 1
    )

    # Add convergence curves
    for (criterion in criteria_log) {
      graphics::lines(
        data$Iteration,
        data[[criterion]],
        col = col_vec[[criterion]],
        lwd = 1
      )
    }

    # Add legend below plot
    graphics::legend(
      "bottom",
      inset = c(0, -0.25),
      legend = criteria_log,
      col = col_vec[criteria_log],
      lty = 1,
      horiz = TRUE,
      xpd = TRUE,
      bty = "n",
      cex = 0.85
    )

  }, error = function(e) {
    base::message(
      "Convergence plot skipped: ",
      base::conditionMessage(e)
    )
  })

  # Return processed data invisibly
  base::invisible(data)
}
