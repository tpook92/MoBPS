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

#' Add a LIS-based trait
#'
#' Function to add a LIS-based trait
#' @param population population list
#' @param traits Which (existing traits to use as new traits), default: NULL --> generate new traits
#' @param meet.performer.trait Include trait for nr. meets as the performer (default: FALSE)
#' @param meet.receiver.trait Include trait for nr. meets as the receiver (default: FALSE)
#' @param interact.performer.trait Include trait for nr. interactions as the performer (default: FALSE)
#' @param interact.receiver.trait Include trait for nr. interactions as the receiver (default: FALSE)
#' @param bit.performer.trait Include trait for damage performed (default: FALSE)
#' @param bit.receiver.trait Include trait for damage received (default: FALSE)
#' @param input.traits Which (existing traits to use as auxially traits), default: NULL --> generate new traits
#' @param input.traits.qtl Corresponds to n.additive for auxillary traits
#' @param input.traits.var.target Corresponds to var.target for auxillary traits
#' @param input.traits.mean.target Corresponds to mean.target for auxiallary traits
#' @param trait.name Provide manual trait names. Defaults: c("Meet_engaged", "Meet_received","Interaction_engaged", "Interaction_received", "Damage_performed", "Damage_received")
#' @param input.traits.name Provide manual trait names. Defaults: c("performer_effect", "recipient_effect", "social_tendency", "bite_force")
#' @param lis.path Path of the LIS executable
#' @param lis.dir Directory where to write files for LIS to
#' @param step Step parameter for LIS
#' @param lis.parameter Path to LIS parameter file
#' @examples
#' data(ex_pop)
#' population <- add.array(ex_pop, marker.included = c(TRUE, FALSE), array.name="Half-density")
#' @return Population list
#' @export

add.lis = function(population, traits = NULL,
                   lis.path = "LIS",
                   lis.dir = "LIS_files/",
                   lis.parameter = "parameter.ini",
                   step = 5000,
                   meet.performer.trait = FALSE,
                   meet.receiver.trait = FALSE,
                   interact.performer.trait = FALSE,
                   interact.receiver.trait = FALSE,
                   bit.performer.trait = FALSE,
                   bit.receiver.trait = FALSE,
                   input.traits = NULL,
                   input.traits.qtl = rep(1000,4),
                   input.traits.var.target = c(0.25, 0.25, 1, 0.10),
                   input.traits.mean.target = c(0,0,0,0),
                   trait.name = NULL,
                   input.traits.name = NULL

){


  active.trait = c(meet.performer.trait, meet.receiver.trait,
                   interact.performer.trait, interact.receiver.trait,
                   bit.performer.trait, bit.receiver.trait)

  if(length(trait.name)==0){
    trait.name = c("Meet_engaged", "Meet_received",
                   "Interaction_engaged", "Interaction_received",
                   "Damage_performed", "Damage_received")[active.trait]
  }

  if(length(input.traits.name)==0){
    input.traits.name = c("performer_effect", "recipient_effect",
                   "social_tendency", "bite_force")
  }


  if(length(input.traits)==0){
    input.traits = 1:4 + population$info$bv.nr
    population <- creating.trait(population, n.traits = population$info$bv.nr + 4, verbose=FALSE,
                                 n.additive = input.traits.qtl,
                                 var.target = input.traits.var.target,
                                 mean.target = input.traits.mean.target,
                                 trait.name = input.traits.name)
  } else if(length(input.traits)!=4){
    stop("Either provide 4 input traits or none!")
  }

  if(length(traits)==0){
    traits = population$info$bv.nr + 1:sum(active.trait)
  }

  if(traits[1] > population$info$bv.nr){
    if(max(traits) == (population$info$bv.nr+sum(active.trait))){
      population <- creating.trait(population, n.traits = max(traits), verbose=FALSE,
                                   trait.name = trait.name,
                                   polygenic.variance = 0.0001)
    } else{
      stop("Illegal traits to enter")
    }
  }

  population$info$is.combi[traits] <- FALSE
  population$info$is.maternal[traits] <- FALSE
  population$info$is.paternal[traits] <- FALSE
  population$info$is.lis[traits] <- TRUE


  if(length(population$info$lis_info)==0){
    population$info$lis_info = list()
  }
  population$info$lis_info[[traits[1]]] = list(
    traits = traits,
    active.trait = active.trait,
    input.traits = input.traits,
    lis.path = lis.path,
    lis.dir = lis.dir,
    step = step,
    lis.parameter = lis.parameter
  )

  return(population)
}
