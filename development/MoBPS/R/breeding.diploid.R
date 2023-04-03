'#
  Authors
Torsten Pook, torsten.pook@uni-goettingen.de

Copyright (C) 2017 -- 2020  Torsten Pook

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
#### Selection of individuals
#' @param selection.size Number of selected individuals as parents (default: all individuals in selection.m/f.gen/database/gen - alt: positive numbers)
#' @param selection.m.gen,selection.m.cohorts,selection.m.database Generations/cohorts/groups available for selection of first/paternal parent
#' @param selection.f.gen,selection.f.cohorts,selection.f.database Generations available for selection of maternal parent
#' @param selection.criteria What to use in the selection proces (default: "bve", alt: "bv", "pheno", "random", "offpheno")
#' @param class.m,class.f For selection only individuals from this class (included in selection.m/f.gen/database/cohorts) will be considered for selection (default: 0 - which is all individuals if never used class elsewhere)
#' @param add.class.cohorts Inital classes of cohorts used in selection.m/f.cohorts are automatically added to class.m/f (default: TRUE)
#' @param multiple.bve Way to handle multiple traits in selection (default: "add" - use values directly in an index, alt: "ranking" - ignore values but only use ranking per trait)
#' @param multiple.bve.weights.m,multiple.bve.weights.f Weighting between traits (default: 1)
#' @param multiple.bve.scale.m,multiple.bve.scale.f Default: "bv_sd"; Set to "pheno_sd" when using gains per phenotypic SD, "unit" when using gains per unit, "bve" when using estimated breeding values
#' @param selection.highest If FALSE to select individuals with lowest value for the selection criterium (default c(TRUE,TRUE) - (m,w))
#' @param ignore.best Not consider the top individuals of the selected individuals (e.g. to use 2-10 best individuals)
#' @param best.selection.ratio.m,best.selection.ratio.f Ratio of the frequency of the selection of the best best individual and the worst best individual (default=1)
#' @param best.selection.criteria.m,best.selection.criteria.f Criteria to calculate this ratio (default: "bv", alt: "bve", "pheno")
#' @param best.selection.manual.ratio.m,best.selection.manual.ratio.f vector containing probability to draw from for every individual (e.g. c(0.1,0.2,0.7))
#' @param best.selection.manual.reorder Set to FALSE to not use the order from best to worst selected individual but plain order based on database-order
#' @param selection.m.random.prob,selection.f.random.prob Use this parameter to control the probablity of each individual to be selected when doing random selection
#' @param reduced.selection.panel.m,reduced.selection.panel.f Use only a subset of individuals of the potential selected ones ("Split in user-interface")
#' @param threshold.selection.index Selection index on which to access (matrix which one index per row)
#' @param threshold.selection.value Minimum value in the selection index selected individuals have to have
#' @param threshold.selection.sign Pick all individuals above (">") the threshold. Alt: ("<", "=", "<=", ">=")
#' @param threshold.selection.criteria Criterium on which to evalute the index (default: "bve", alt: "bv", "pheno")
#' @param remove.duplicates Set to FALSE to select the same individual multiple times when the gen/database/cohorts for selection contains it multiple times
#' @param selection.m.miesenberger,selection.f.miesenberger Use Weighted selection index according to Miesenberger 1997 for paternal/maternal selection
#' @param miesenberger.trafo Ignore all eigenvalues below this threshold and apply dimension reduction (default: 0 - use all)
#' @param selection.miesenberger.reliability.est If available reliability estimated are used. If not use default: "derived" (cor(BVE,BV)^2) , alt: "heritability", "estimated" (SD BVE / SD Pheno) as replacement
#' @param sort.selected.pos Set to TRUE to arrange selected individuals according to position in the database (not by breeding value)
#' @param ogc If TRUE use optimal genetic contribution theory to perform selection ( This requires the use of the R-package optiSel)
#' @param relationship.matrix.ogc Method to calculate relationship matrix for OGC (Default: "pedigree", alt: "vanRaden", "CE", "non_stand", "CE2", "CM")
#' @param depth.pedigree.ogc Depth of the pedigree in generations (default: 7)
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
#### Generation of new individuals
#' @param breeding.size Number of individuals to generate (default: 0, use vector with two entries to control offspring per sex)
#' @param breeding.size.litter Number of litters to generate (default: NULL - use breeding.size; only single positive number input allowed)
#' @param name.cohort Name of the newly added cohort
#' @param breeding.sex Share of female individuals (if single value is used for breeding size; default: 0.5)
#' @param breeding.sex.random If TRUE randomly chose sex of new individuals (default: FALSE - use expected values)
#' @param sex.s Specify which newly added individuals are male (1) or female (2)
#' @param add.gen Generation you want to add the new individuals to (default: New generation)
#' @param share.genotyped Share of individuals newly generated individuals that are genotyped (Default: 0)
#' @param phenotyping.child Starting phenotypes of newly generated individuals (default: "mean" of both parents, "obs" - regular observation, "zero" - 0)
#' @param fixed.effects.p Parametrization for the fixed effects (default: c(0,0..,0), if multiple different parametrizations are possible use a matrix with one parametrization per row)
#' @param fixed.effects.freq Frequency of each different parametrization of the fixed effects
#' @param new.class Migration level of newly generated individuals (default: 0 / use vector for different classes for different sexes)
#' @param max.offspring Maximum number of offspring per individual (default: c(Inf,Inf) - (m,w))
#' @param max.litter Maximum number of litters per individual (default: c(Inf,Inf) - (m,w))
#' @param max.mating.pair Maximum number of matings between two specific individuals (default: Inf)
#' @param avoid.mating.fullsib Set to TRUE to not generate offspring of full siblings
#' @param avoid.mating.halfsib Set to TRUE to not generate offspring from half or full siblings
#' @param fixed.breeding Set of targeted matings to perform (matrix with 7 columns: database position first parent (gen, sex, nr), database position second parent (gen,sex,nr), likelihood to be female (optional))
#' @param fixed.breeding.best Perform targeted matings in the group of selected individuals (matrix with 5 columns:  position first parent (male/female pool of selected individuals, ranking in selected animals), position second parent (male/female pool of selected individuals, ranking in selected animals), likelihood to be female (optional))
#' @param fixed.assignment Set to "bestbest" / TRUE for targeted mating of best-best individual till worst-worst (of selected). set to "bestworst" for best-worst mating
#' @param breeding.all.combination Set to TRUE to automatically perform each mating combination possible exactly ones.
#' @param repeat.mating Generate multiple mating from the same dam/sire combination (first column: number of offspring; second column: probability)
#' @param repeat.mating.copy Generate multiple copies from a copy action (combine / copy.individuals.m/f) (first column: number of offspring; second column: probability)
#' @param repeat.mating.fixed Vector containing number of times each mating is repeated. This will overwrite sampling from repeat.mating / repeat.mating.copy (default: NULL)
#' @param repeat.mating.overwrite Set to FALSE to not use the current repeat.mating / repeat.mating.copy input as the new standard values (default: TRUE)
#' @param repeat.mating.trait Trait that should be linked to the litter size
#' @param repeat.mating.max Maximum number of individuals in a litter
#' @param repeat.mating.s Use this parameter to manually provide the size of each litter generated
#' @param same.sex.activ If TRUE allow matings of individuals of same sex (Sex here is a general term with the first sex referring to the first parent, second sex second parent)
#' @param same.sex.sex Probability to use female individuals as parents (default: 0.5)
#' @param same.sex.selfing Set to TRUE to allow for selfing when using same.sex matings (default: FALSE)
#' @param selfing.mating If TRUE generate new individuals via selfing
#' @param selfing.sex Share of female individuals used for selfing (default: 0.5)
#' @param dh.mating If TRUE generate a DH-line in mating process
#' @param dh.sex Share of DH-lines generated from selected female individuals
#' @param combine Copy existing individuals (e.g. to merge individuals from different groups in a joined cohort). Individuals to use are used as the first parent
#' @param copy.individual If TRUE a copy of the first parent will be generated instead of meosis between two parents
#' @param copy.individual.m,copy.individual.f If TRUE generate exactly one copy of all selected male/female in a new cohort (or more by setting breeding.size)
#' @param copy.individual.keep.pheno Set to FALSE to not keep phenotypes in case of use of copying individuals instead of regular meoisis
#' @param copy.individual.keep.bve Set to FALSE to not keep estimated breeding value in case of use of copying individuals instead of regular meoisis
#' @param added.genotyped Share of individuals that is additionally genotyped (only for copy.individuals, default: 0)
#' @param bv.ignore.traits Vector of traits to ignore in the calculation of the genomic value (default: NULL; Only recommended for high number of traits and experienced users!)
#' @param generation.cores Number of cores used for the generation of new individuals (This will only be active when generating more than 500 individuals)
#### Genotyping (already existing individuals)
#' @param genotyped.gen,genotyped.cohorts,genotyped.database Generations/cohorts/groups to generate genotype data (that can be used in a BVE)
#' @param genotyped.share Share of individuals in genotyped.gen/database/cohort to generate genotype data from (default: 1)
#' @param genotyped.array Genotyping array used
#' @param genotyped.remove.gen,genotyped.remove.database,genotyped.remove.cohorts Generations/cohorts/groups from which to remove genotyping information (this will affect all copies of an individual unless genotyped.remove.all.copy is set to FALSE)
#' @param genotyped.remove.all.copy Set to FALSE to only change the genotyping state of this particular copy of an individual (default: TRUE)
#### Phenotyping (already existing individuals)
#' @param phenotyping Quick acces to phenotyping for (all: "all", non-phenotyped: "non_obs", non-phenotyped male: "non_obs_m", non-phenotyped female: "non_obs_f")
#' @param phenotyping.gen,phenotyping.cohorts,phenotyping.database Generations/cohorts/groups from which to generate additional phenotypes
#' @param n.observation Number of phenotypic observations generated per trait and per individuals (use repeatability to control correlation between observations)
#' @param phenotyping.class Classes of individuals for which to generate phenotypes (default: NULL --> all classes)
#' @param heritability Use sigma.e to obtain a certain heritability (default: NULL)
#' @param repeatability Set this to control the share of the residual variance (sigma.e) that is permanent (there for each observation)
#' @param multiple.observation If an already phenotyped trait is phenotyped again this will on NOT lead to an additional phenotyped observation unless this is set to TRUE
#' @param share.phenotyped Share of the individuals to phenotype (use vector for different probablities for different traits)
#' @param offpheno.parents.gen,offpheno.parents.database,offpheno.parents.cohorts Generations/groups/cohorts to consider to derive phenotype from offspring phenotypes
#' @param offpheno.offspring.gen,offpheno.offspring.cohorts,offpheno.offspring.database Active generations/cohorts/groups for import of offspring phenotypes
#' @param sigma.e Enviromental standard deviation (default: use sigma.e from last run / usually fit by use of heritability; if never provided: 10; used in BVE for variance components if manually set)
#' @param sigma.e.gen,sigma.e.cohorts,sigma.e.database Generations/cohorts/groups to consider when estimating sigma.e when using heritability
#' @param new.residual.correlation Correlation of the simulated residual variance
#' @param new.breeding.correlation Correlation of the simulated genetic variance (only impacts non-QTL based traits. Needs to be fit in creating.diploid/trait for QTL-based traits)
#### Breeding value estimation
#' @param bve If TRUE perform a breeding value estimation (default: FALSE)
#' @param bve.gen,bve.cohorts,bve.database Generations/Groups/Cohorts of individuals to consider in breeding value estimation (default: NULL)
#' @param relationship.matrix Method to calculate relationship matrix for the breeding value estimation (Default: "vanRaden", alt: "pedigree", "CE", "non_stand", "CE2", "CM")
#' @param depth.pedigree Depth of the pedigree in generations (default: 7)
#' @param singlestep.active Set FALSE remove all individuals without genomic data from the breeding value estimation
#' @param bve.ignore.traits Vector of traits to ignore in the breeding value estimation (default: NULL, use: "zero" to not consider traits with 0 index weight in multiple.bve.weights.m/.w)
#' @param bve.array Array to use in the breeding value estimation (default: NULL; chose largest possible based on used individuals in BVE)
#' @param bve.imputation Set to FALSE to not perform imputation up to the highest marker density of genotyping data that is available
#' @param bve.imputation.errorrate Share of errors in the imputation procedure (default: 0)
#' @param bve.all.genotyped Set to TRUE to act as if every individual in the breeding value estimation has been genotyped
#' @param bve.insert.gen,bve.insert.cohorts,bve.insert.database Generations/Groups/Cohorts of individuals to compute breeding values for (default: all groups in bve.database)
#' @param variance.correction Correct for "parental.mean" or "generation.mean" in the estimation of  sigma.g for BVE / sigma.e estimation (default: "none")
#' @param bve.class Consider only individuals of those class classes in breeding value estimation (default: NULL - use all)
#' @param sigma.g Genetic standard deviation (default: calculated based on individuals in BVE ; used in BVE for variance components if manually set; mostly recommended to be used for non-QTL based traits)
#' @param sigma.g.gen,sigma.g.cohorts,sigma.g.database Generations/cohorts/groups to consider when estimating sigma.g
#' @param forecast.sigma.g Set FALSE to not estimate sigma.g (Default: TRUE // in case sigma.g is set this is automatically set to FALSE)
#' @param remove.effect.position If TRUE remove real QTLs in breeding value estimation
#' @param estimate.add.gen.var If TRUE estimate additive genetic variance and heritability based on parent model
#' @param estimate.pheno.var If TRUE estimate total variance in breeding value estimation
#' @param bve.avoid.duplicates If set to FALSE multiple generatations of the same individual can be used in the bve (only possible by using copy.individual to generate individuals)
#' @param calculate.reliability Set TRUE to calculate a reliability when performing Direct-Mixed-Model BVE
#' @param estimate.reliability Set TRUE to estimate the reliablity in the BVE by calculating the correlation between estimated and real breeding values
#' @param bve.input.phenotype Select what to use in BVE (default: own phenotype ("own"), offspring phenotype ("off"), their average ("mean") or a weighted average ("weighted"))
#' @param mas.bve If TRUE use marker assisted selection in the breeding value estimation
#' @param mas.markers Vector containing markers to be used in marker assisted selection
#' @param mas.number If no markers are provided this nr of markers is selected (if single marker QTL are present highest effect markers are prioritized)
#' @param mas.effects Effects assigned to the MAS markers (Default: estimated via lm())
#' @param mas.geno Genotype dataset used in MAS (default: NULL, automatic internal calculation)
#' @param bve.parent.mean Set to TRUE to use the average parental performance as the breeding value estimate
#' @param bve.grandparent.mean Set to TRUE to use the average grandparental performance as the breeding value estimate
#' @param bve.mean.between Select if you want to use the "bve", "bv", "pheno" or "bvepheno" to form the mean (default: "bvepheno" - if available bve, else pheno)
#' @param bve.exclude.fixed.effects Vector of fixed effects to ignore in the BVE (default: NULL)
#' @param bve.beta.hat.approx Set to FALSE to use the true underlying value for beta_hat for the fixed effect in the direct BVE model. rrBLUP, BGLR, sommer will always estimate beta_hat.
#' @param bve.per.sample.sigma.e Set to FALSE to deactivate the use of a heritablity based on the number of observations generated per sample
#### Software for breeding value estimation
#' @param mobps.bve If TRUE predict BVEs in direct estimation with assumed known heritablity (default: TRUE; activating use of any other BVE method to TRUE will overwrite this)
#' @param mixblup.bve Set to TRUE to activate breeding value estimation via MiXBLUP (requires MiXBLUP license!)
#' @param emmreml.bve If TRUE use REML estimator from R-package EMMREML in breeding value estimation
#' @param rrblup.bve If TRUE use REML estimator from R-package rrBLUP in breeding value estimation
#' @param sommer.bve If TRUE use REML estimator from R-package sommer in breeding value estimation
#' @param sommer.multi.bve Set TRUE to use a mulit-trait model in the R-package sommer for BVE
#' @param BGLR.bve If TRUE use BGLR to perform breeding value estimation
#' @param pseudo.bve If set to TRUE the breeding value estimation will be simulated with resulting accuracy pseudo.bve.accuracy (default: 1)
#' @param pseudo.bve.accuracy The accuracy to be obtained in the "pseudo" - breeding value estimation
#' @param bve.solve Provide solver to be used in BVE (default: "exact" solution via inversion, alt: "pcg", function with inputs A, b and output y_hat)
#' @param mixblup.pedfile Set to FALSE to manually generate your MiXBLUP pedfile
#' @param mixblup.parfile Set to FALSE to manually generate your MiXBLUP parfile
#' @param mixblup.datafile Set to FALSE to manually write your MiXBLUP datafile
#' @param mixblup.inputfile Set to FALSE to manually write your MiXBLUP inputfile
#' @param mixblup.path Provide path to MiXBLUP.exe (default is your working directory: Windows: MixBLUP; Linux ./MixBLUP.exe)
#' @param mixblup.files Directory to generate all files generated when using MiXBLUP (default: MiXBLUP_files/ )
#' @param mixblup.genofile Set to FALSE to manually write the MiXBLUP genotypefile
#' @param mixblup.path.pedfile Path from where to import the MiXBLUP pedfile
#' @param mixblup.path.parfile Path from where to import the MiXBLUP parfile
#' @param mixblup.path.datafile Path from where to import the MiXBLUP datafile
#' @param mixblup.path.inputfile Path from where to import the MiXBLUP inputfile
#' @param mixblup.path.genofile Path from where to import the MiXBLUP genofile
#' @param mixblup.lambda Lambda parameter in MiXBLUP
#' @param mixblup.omega Omega parameter in MiXBLUP
#' @param mixblup.alpha Alpha parameter in MiXBLUP (default: 1, with alpha + beta = 1)
#' @param mixblup.beta Beta parameter in MiXBLUP (default: 0, with alpha + beta = 1)
#' @param mixblup.verbose Set to TRUE to display MiXBLUP prints
#' @param mixblup.genetic.cov Provide genetic covariance matrix to be used in MiXBLUP (lower-triangle is sufficent)
#' @param mixblup.residual.cov Provide residual covariance matrix to be used in MiXBLUP (lower-triangle is sufficent)
#' @param mixblup.apy Set to TRUE to use APY inverse in MiXBLUP (default: FALSE)
#' @param mixblup.apy.core Number of core individuals in the APY algorithm (default: 5000)
#' @param BGLR.model Select which BGLR model to use (default: "RKHS", alt: "BRR", "BL", "BayesA", "BayesB", "BayesC")
#' @param BGLR.burnin Number of burn-in steps in BGLR (default: 1000)
#' @param BGLR.iteration Number of iterations in BGLR (default: 5000)
#' @param BGLR.print If TRUE set verbose to TRUE in BGLR
#' @param BGLR.save Method to use in BGLR (default: "RKHS" - alt: NON currently)
#' @param BGLR.save.random Add random number to store location of internal BGLR computations (only needed when simulating a lot in parallel!)
#' @param miraculix If TRUE use miraculix to perform computations (ideally already generate population in creating.diploid with this; default: automatic detection from population list)
#' @param miraculix.mult If TRUE use miraculix for matrix multiplications even if miraculix is not used for storage
#' @param miraculix.chol Set to FALSE to deactive miraculix based Cholesky-decomposition (default: TRUE)
#' @param miraculix.cores Number of cores used in miraculix applications (default: 1)
#' @param miraculix.destroyA If FALSE A will not be destroyed in the process of inversion (less computing / more memory)
#### Estimation of SNP-effects, GWAS, genome editing
#' @param estimate.u If TRUE estimate u in breeding value estimation (Y = Xb + Zu + e)
#' @param fast.uhat Set to FALSE to  derive inverse of A in rrBLUP (only required when this becomes numerical unstable otherwise)
#' @param gwas.u If TRUE estimate u via GWAS (relevant for gene editing)
#' @param approx.residuals If FALSE calculate the variance for each marker separatly instead of using a set variance (doesnt change order - only p-values)
#' @param gwas.gen,gwas.cohorts,gwas.database Generations/cohorts/groups to consider in GWAS analysis
#' @param gwas.group.standard If TRUE standardize phenotypes by group mean
#' @param y.gwas.used What y value to use in GWAS study (Default: "pheno", alt: "bv", "bve")
#' @param gene.editing.offspring If TRUE perform gene editing on newly generated individuals
#' @param gene.editing.best If TRUE perform gene editing on selected individuals
#' @param gene.editing.offspring.sex Which sex to perform editing on (Default c(TRUE,TRUE), mw)
#' @param gene.editing.best.sex Which sex to perform editing on (Default c(TRUE,TRUE), mw)
#' @param nr.edits Number of edits to perform per individual
#### Culling
#' @param culling.gen,culling.cohorts,culling.database Generations/cohorst/groups to consider to culling
#' @param culling.time Age of the individuals at culling
#' @param culling.name Name of the culling action (user-interface stuff)
#' @param culling.bv1 Reference Breeding value
#' @param culling.share1 Probability of death for individuals with bv1
#' @param culling.bv2 Alternative breeding value (linear extended for other bvs)
#' @param culling.share2 Probability of death for individuals with bv2
#' @param culling.index Genomic index (default:0 - no genomic impact, use: "lastindex" to use the last selection index applied in selection)
#' @param culling.single Set to FALSE to not apply the culling module on all individuals of the cohort
#' @param culling.all.copy Set to FALSE to not kill copies of the same individual in the culling module
#### Meiosis Parameter
#' @param mutation.rate Mutation rate in each marker (default: 10^-8)
#' @param remutation.rate Remutation rate in each marker (default: 10^-8)
#' @param recombination.rate Average number of recombination per 1 length unit (default: 1M)
#' @param recombination.function Function used to calculate position of recombination events (default: MoBPS::recombination.function.haldane())
#' @param recombination.minimum.distance Minimum distance between two points of recombination (default: 0)
#' @param recombination.distance.penalty Reduced probability for recombination events closer than this value - linear penalty (default: 0)
#' @param recombination.distance.penalty.2 Reduced probability for recombination events closer than this value - quadratic penalty (default: 0)
#' @param recom.f.indicator Use step function for recombination map (transform snp.positions if possible instead)
#' @param import.position.calculation Function to calculate recombination point into adjacent/following SNP
#' @param duplication.rate Share of recombination points with a duplication (default: 0 - DEACTIVATED)
#' @param duplication.length Average length of a duplication (Exponentially distributed)
#' @param duplication.recombination Average number of recombinations per 1 length uit of duplication (default: 1)
#' @param gen.architecture.m,gen.architecture.f Genetic architecture for male/female individuals (default: 0 - no transformation)
#' @param add.architecture List with two vectors containing (A: length of chromosomes, B: position in cM of SNPs)
#' @param intern.func Chose which function will be used for simulation of meosis (default: 0, alt: 1,2) - can be faster for specific cases
#### Advanced memory savings
#' @param delete.haplotypes Generations for with haplotypes of founders can be deleted from population list for memory reduction (default: NULL)
#' @param delete.individuals Generations for with individuals are completely deleted from population list for memory reduction (default: NULL)
#' @param delete.gen Generations to entirely deleted fro population list for memory reduction (default: NULL)
#' @param delete.sex Remove all individuals from these sex from generation delete.individuals (default: 1:2 ; note:delete individuals=NULL)
#' @param delete.same.origin If TRUE delete recombination points when genetic origin of adjacent segments is the same
#' @param save.recombination.history If TRUE store the time point of each recombination event
#' @param store.sparse Set to TRUE to store the pedigree relationship matrix as a sparse matrix
#' @param storage.save Lower numbers will lead to less memory but slightly higher computing time for calculation of the pedigree relationship matrix (default: 1.5, min: 1)
#### Tracking/Reporting of breeding actions & computing time
#' @param verbose Set to FALSE to not display any prints
#' @param report.accuracy Report the accuracy of the breeding value estimation
#' @param store.breeding.totals If TRUE store information on selected individuals in $info$breeding.totals (default: FALSE)
#' @param store.bve.data If TRUE store information of bve in $info$bve.data
#' @param store.comp.times If TRUE store computation times in $info$comp.times.general (default: TRUE)
#' @param store.comp.times.bve If TRUE store computation times of breeding value estimation in $info$comp.times.bve (default: TRUE)
#' @param store.comp.times.generation If TRUE store computation times of mating simulations in $info$comp.times.generation (default: TRUE)
#' @param store.effect.freq If TRUE store the allele frequency of effect markers per generation
#' @param Rprof Store computation times of each function
#' @param randomSeed Set random seed of the process
#' @param display.progress Set FALSE to not display progress bars. Setting verbose to FALSE will automatically deactive progress bars
#' @param time.point Time point at which the new individuals are generated
#' @param creating.type Technique to generate new individuals (use mostly intended for web-based application)
#### Import / Export
#' @param import.relationship.matrix Input the wanted relationship matrix with this parameter (default: NULL - relationship matrix will be calculated from other sources)
#' @param export.selected Set to TRUE to export the list of selected individuals
#' @param export.relationship.matrix Export the relationship matrix used in the breeding value estimation
#### Still in development
#' @param pen.assignments This is a placeholder to deactivate this module for now
#' @param pen.size Pen size. When different types of pen are used: use a matrix with two columns coding Number of individuals per pen, Probability for each pen size
#' @param pen.by.sex Only individuals of the same sex are put in the same pen (default: TRUE)
#' @param pen.by.litter Only individuals of the same litter are put in the same pen (default: FALSE)
#' @param pen.size.overwrite Set to FALSE to not use the input for pen.size for down-stream use of breeding.diploid (default: TRUE)
#### Old Parameters - these are mostly keep to ensure that old scripts still work but are not recommended for use for new scripts
#' @param selection.m,selection.f (OLD! use selection criteria) Selection criteria for male/female individuals (Set to "random" to randomly select individuals - default: "function"  based on selection.criteria ((usually breeding values)))
#' @param new.bv.observation.gen,new.bv.observation.cohorts,new.bv.observation.database (OLD! use phenotyping.gen/cohorts/database) Vector of generation from which to generate additional phenotypes
#' @param best1.from.group,best1.from.cohort (OLD!- use selection.m.database/cohorts) Groups of individuals to consider as First Parent / Father (also female individuals are possible)
#' @param best2.from.group,best2.from.cohort (OLD!- use selection.f.database/cohorts) Groups of individuals to consider as Second Parent / Mother (also male individuals are possible)
#' @param new.bv.observation (OLD! - use phenotyping) Quick acces to phenotyping for (all: "all", non-phenotyped: "non_obs", non-phenotyped male: "non_obs_m", non-phenotyped female: "non_obs_f")
#' @param reduce.group (OLD! - use culling modules) Groups of individuals for reduce to a new size (by changing class to -1)
#' @param reduce.group.selection (OLD! - use culling modules) Selection criteria for reduction of groups (cf. selection.m / selection.f - default: "random")
#' @param new.bv.child (OLD! - use phenotyping.child) Starting phenotypes of newly generated individuals (default: "mean" of both parents, "obs" - regular observation, "zero" - 0)
#' @param computation.A (OLD! - use relationship.matrix) Method to calculate relationship matrix for the breeding value estimation (Default: "vanRaden", alt: "pedigree", "CE", "non_stand", "CE2", "CM")
#' @param computation.A.ogc (OLD! use relationship.matrix.ogc) Method to calculate pedigree matrix in OGC (Default: "pedigree", alt: "vanRaden", "CE", "non_stand", "CE2", "CM")
#' @param new.phenotype.correlation (OLD! - use new.residual.correlation!) Correlation of the simulated enviromental variance
#' @param offspring.bve.parents.gen,offspring.bve.parents.cohorts,offspring.bve.parents.database (OLD! use offpheno.parents.gen/database/cohorts) Generations/cohorts/groups to consider to derive phenotype from offspring phenotypes
#' @param offspring.bve.offspring.gen,offspring.bve.offspring.cohorts,offspring.bve.offspring.database (OLD! use offpheno.offspring.gen/database/cohorts) Active generations/cohorts/groups for import of offspring phenotypes
#' @param input.phenotype (OLD! use bve.input.phenotype) Select what to use in BVE (default: own phenotype ("own"), offspring phenotype ("off"), their average ("mean") or a weighted average ("weighted"))
#' @param threshold.selection Minimum value in the selection index selected individuals have to have
#' @param threshold.sign Pick all individuals above (">") the threshold. Alt: ("<", "=", "<=", ">=")
#### Other
#' @param use.recalculate.manual Set to TRUE to use recalculate.manual to calculate genomic values (all individuals and traits jointly, default: FALSE)
#' @param size.scaling Set to value to scale all input for breeding.size / selection.size (This will not work for all breeding programs / less general than json.simulation)
#' @param parallel.internal Internal parameter for the parallelization
#' @examples
#' population <- creating.diploid(nsnp=1000, nindi=100)
#' population <- breeding.diploid(population, breeding.size=100, selection.size=c(25,25))
#' @return Population-list
#### @usage breeding.diploid(population, ...)
#' @export


breeding.diploid <- function(population,
                             #### Selection of individuals
                             selection.size = 0,
                             selection.criteria = NULL,
                             selection.m.gen = NULL,
                             selection.f.gen = NULL,
                             selection.m.database = NULL,
                             selection.f.database = NULL,
                             selection.m.cohorts=NULL,
                             selection.f.cohorts=NULL,
                             class.m = 0,
                             class.f = 0,
                             add.class.cohorts = TRUE,
                             multiple.bve = "add",
                             multiple.bve.weights.m = 1,
                             multiple.bve.weights.f = NULL,
                             multiple.bve.scale.m = "bv_sd",
                             multiple.bve.scale.f = NULL,
                             selection.highest = c(TRUE,TRUE),
                             ignore.best = 0,
                             best.selection.ratio.m = 1,
                             best.selection.ratio.f = NULL,
                             best.selection.criteria.m = "bv",
                             best.selection.criteria.f = NULL,
                             best.selection.manual.ratio.m = NULL,
                             best.selection.manual.ratio.f = NULL,
                             best.selection.manual.reorder = TRUE,
                             selection.m.random.prob = NULL,
                             selection.f.random.prob = NULL,
                             reduced.selection.panel.m = NULL,
                             reduced.selection.panel.f = NULL,
                             threshold.selection.index = NULL,
                             threshold.selection.value = NULL,
                             threshold.selection.sign = ">",
                             threshold.selection.criteria = "bve",
                             threshold.selection=NULL,
                             threshold.sign=">",
                             remove.duplicates = TRUE,
                             selection.m.miesenberger=FALSE,
                             selection.f.miesenberger=NULL,
                             selection.miesenberger.reliability.est="derived",
                             miesenberger.trafo = 0,
                             sort.selected.pos = FALSE,
                             ogc = FALSE,
                             relationship.matrix.ogc = "pedigree",
                             depth.pedigree.ogc = 7,
                             ogc.target = "min.sKin",
                             ogc.uniform=NULL,
                             ogc.ub = NULL,
                             ogc.lb = NULL,
                             ogc.ub.sKin = NULL,
                             ogc.lb.BV = NULL,
                             ogc.ub.BV = NULL,
                             ogc.eq.BV = NULL,
                             ogc.ub.sKin.increase = NULL,
                             ogc.lb.BV.increase = NULL,


                             #### Generation of new individuals
                             breeding.size = 0,
                             breeding.size.litter = NULL,
                             name.cohort = NULL,
                             breeding.sex = NULL,
                             breeding.sex.random = FALSE,
                             sex.s = NULL,
                             add.gen = 0,
                             share.genotyped = 0,
                             phenotyping.child = NULL,
                             fixed.effects.p = NULL,
                             fixed.effects.freq = NULL,
                             new.class = 0L,
                             max.offspring = Inf,
                             max.litter = Inf,
                             max.mating.pair = Inf,
                             avoid.mating.fullsib=FALSE,
                             avoid.mating.halfsib=FALSE,
                             fixed.breeding = NULL,
                             fixed.breeding.best = NULL,
                             fixed.assignment = FALSE,
                             breeding.all.combination = FALSE,
                             repeat.mating = NULL,
                             repeat.mating.copy = NULL,
                             repeat.mating.fixed = NULL,
                             repeat.mating.overwrite = TRUE,
                             repeat.mating.trait = 1,
                             repeat.mating.max = NULL,
                             repeat.mating.s = NULL,
                             same.sex.activ = FALSE,
                             same.sex.sex = 0.5,
                             same.sex.selfing = FALSE,
                             selfing.mating = FALSE,
                             selfing.sex = 0.5,
                             dh.mating = FALSE,
                             dh.sex = 0.5,
                             combine = FALSE,
                             copy.individual = FALSE,
                             copy.individual.m = FALSE,
                             copy.individual.f = FALSE,
                             copy.individual.keep.bve = TRUE,
                             copy.individual.keep.pheno = TRUE,
                             added.genotyped = 0,
                             bv.ignore.traits=NULL,
                             generation.cores = NULL,

                             #### Genotying (already existing individuals)
                             genotyped.database = NULL,
                             genotyped.gen = NULL,
                             genotyped.cohorts = NULL,
                             genotyped.share = 1,
                             genotyped.array = 1,
                             genotyped.remove.gen = NULL,
                             genotyped.remove.database = NULL,
                             genotyped.remove.cohorts = NULL,
                             genotyped.remove.all.copy = TRUE,

                             #### Phenotyping (already existing individuals)
                             phenotyping = NULL,
                             phenotyping.gen = NULL,
                             phenotyping.cohorts = NULL,
                             phenotyping.database = NULL,
                             n.observation = NULL,
                             phenotyping.class = NULL,
                             heritability = NULL,
                             repeatability = NULL,
                             multiple.observation = FALSE,
                             share.phenotyped=1,
                             offpheno.parents.gen = NULL,
                             offpheno.parents.database = NULL,
                             offpheno.parents.cohorts = NULL,
                             offpheno.offspring.gen = NULL,
                             offpheno.offspring.database = NULL,
                             offpheno.offspring.cohorts = NULL,
                             sigma.e = NULL,
                             sigma.e.gen = NULL,
                             sigma.e.cohorts = NULL,
                             sigma.e.database = NULL,
                             new.residual.correlation = NULL,
                             new.breeding.correlation = NULL,

                             #### Breeding value estimation
                             bve = FALSE,
                             bve.gen = NULL,
                             bve.cohorts = NULL,
                             bve.database = NULL,
                             relationship.matrix = "vanRaden",
                             depth.pedigree = 7,
                             singlestep.active = TRUE,
                             bve.ignore.traits=NULL,
                             bve.array = NULL,
                             bve.imputation = TRUE,
                             bve.imputation.errorrate = 0,
                             bve.all.genotyped = FALSE,
                             bve.insert.gen = NULL,
                             bve.insert.cohorts = NULL,
                             bve.insert.database = NULL,
                             variance.correction = "none",
                             bve.class = NULL,
                             sigma.g = NULL,
                             sigma.g.gen = NULL,
                             sigma.g.cohorts = NULL,
                             sigma.g.database = NULL,
                             forecast.sigma.g = NULL,
                             remove.effect.position = FALSE,
                             estimate.add.gen.var = FALSE,
                             estimate.pheno.var = FALSE,
                             bve.avoid.duplicates = TRUE,
                             calculate.reliability=FALSE,
                             estimate.reliability=FALSE,
                             bve.input.phenotype="own",
                             mas.bve=FALSE,
                             mas.markers=NULL,
                             mas.number=5,
                             mas.effects=NULL,
                             mas.geno = NULL,
                             bve.parent.mean=FALSE,
                             bve.grandparent.mean=FALSE,
                             bve.mean.between="bvepheno",
                             bve.exclude.fixed.effects = NULL,
                             bve.beta.hat.approx = TRUE,
                             bve.per.sample.sigma.e=TRUE,

                             #### Software for breeding value estimation
                             mobps.bve=TRUE,
                             mixblup.bve = FALSE,
                             emmreml.bve = FALSE,
                             rrblup.bve = FALSE,
                             sommer.bve = FALSE,
                             sommer.multi.bve=FALSE,
                             BGLR.bve = FALSE,
                             pseudo.bve=FALSE,
                             pseudo.bve.accuracy=1,
                             bve.solve = "exact",
                             mixblup.pedfile=TRUE,
                             mixblup.parfile=TRUE,
                             mixblup.datafile=TRUE,
                             mixblup.inputfile=TRUE,
                             mixblup.genofile=TRUE,
                             mixblup.path=NULL,
                             mixblup.path.pedfile=NULL,
                             mixblup.path.parfile=NULL,
                             mixblup.path.datafile=NULL,
                             mixblup.path.inputfile=NULL,
                             mixblup.path.genofile=NULL,
                             mixblup.files="MiXBLUP_files",
                             mixblup.verbose = TRUE,
                             mixblup.genetic.cov = NULL,
                             mixblup.residual.cov = NULL,
                             mixblup.lambda = 1,
                             mixblup.alpha = NULL,
                             mixblup.beta = NULL,
                             mixblup.omega = NULL,
                             mixblup.apy = FALSE,
                             mixblup.apy.core = NULL,
                             BGLR.model = "RKHS",
                             BGLR.burnin = 500,
                             BGLR.iteration = 5000,
                             BGLR.print = TRUE,
                             BGLR.save = "RKHS",
                             BGLR.save.random = FALSE,
                             miraculix = NULL,
                             miraculix.cores = 1,
                             miraculix.mult = NULL,
                             miraculix.chol = TRUE,
                             miraculix.destroyA=TRUE,

                             #### Estimation of SNP-effects, GWAS, genome editing
                             estimate.u = FALSE,
                             fast.uhat = TRUE,
                             gwas.u = FALSE,
                             approx.residuals = TRUE,
                             gwas.gen = NULL,
                             gwas.cohorts = NULL,
                             gwas.database = NULL,
                             gwas.group.standard = FALSE,
                             y.gwas.used = "pheno",
                             gene.editing.offspring = FALSE,
                             gene.editing.best = FALSE,
                             gene.editing.offspring.sex = c(TRUE,TRUE),
                             gene.editing.best.sex = c(TRUE,TRUE),
                             nr.edits = 0,

                             #### Culling
                             culling.gen=NULL,
                             culling.database=NULL,
                             culling.cohorts=NULL,
                             culling.time = Inf,
                             culling.name = "Not_named",
                             culling.bv1 = 0,
                             culling.share1 = 0,
                             culling.bv2 = NULL,
                             culling.share2 = NULL,
                             culling.index = 0,
                             culling.single = TRUE,
                             culling.all.copy = TRUE,

                             #### Meiosis Parameters
                             mutation.rate = 10^-8,
                             remutation.rate = 10^-8,
                             recombination.rate = 1,
                             recombination.function = NULL,
                             recombination.minimum.distance = NULL,
                             recombination.distance.penalty = NULL,
                             recombination.distance.penalty.2 = NULL,
                             recom.f.indicator = NULL,
                             import.position.calculation = NULL,
                             duplication.rate = 0,
                             duplication.length = 0.01,
                             duplication.recombination = 1,
                             gen.architecture.m = 0,
                             gen.architecture.f = NULL,
                             add.architecture = NULL,
                             intern.func = 0,

                             #### Advanced memory savings
                             delete.haplotypes = NULL,
                             delete.individuals = NULL,
                             delete.gen = NULL,
                             delete.sex = 1:2,
                             delete.same.origin = FALSE,
                             save.recombination.history = FALSE,
                             store.sparse = FALSE,
                             storage.save=1.5,

                             #### Tracking/Reporting of breeding actions & computing time
                             verbose=TRUE,
                             report.accuracy = TRUE,
                             store.breeding.totals = FALSE,
                             store.bve.data = FALSE,
                             store.comp.times = TRUE,
                             store.comp.times.bve = TRUE,
                             store.comp.times.generation = TRUE,
                             store.effect.freq = FALSE,
                             Rprof = FALSE,
                             randomSeed = NULL,
                             display.progress = TRUE,
                             time.point = 0,
                             creating.type = 0,

                             #### Import / Export
                             import.relationship.matrix = NULL,
                             export.selected = FALSE,
                             export.relationship.matrix = FALSE,

                             #### Still development
                             pen.assignments = NULL,
                             pen.size = NULL,
                             pen.by.sex = TRUE,
                             pen.by.litter = FALSE,
                             pen.size.overwrite = TRUE,

                             #### Old Parameter
                             selection.m = NULL,
                             selection.f = NULL,
                             new.bv.observation.gen = NULL,
                             new.bv.observation.cohorts = NULL,
                             new.bv.observation.database = NULL,
                             best1.from.group = NULL,
                             best2.from.group = NULL,
                             best1.from.cohort = NULL,
                             best2.from.cohort = NULL,
                             new.bv.observation = NULL,
                             reduce.group = NULL,
                             reduce.group.selection = "random",
                             new.bv.child = NULL,
                             computation.A = NULL,
                             computation.A.ogc = NULL,
                             new.phenotype.correlation = NULL,
                             offspring.bve.parents.gen = NULL,
                             offspring.bve.parents.database = NULL,
                             offspring.bve.parents.cohorts = NULL,
                             offspring.bve.offspring.gen = NULL,
                             offspring.bve.offspring.database = NULL,
                             offspring.bve.offspring.cohorts = NULL,
                             input.phenotype=NULL,

                             #### Other
                             use.recalculate.manual = FALSE,
                             size.scaling = NULL,
                             parallel.internal = FALSE){

  # Implement ?
  # @param min.coanc.mating Set TRUE to activate the use of min coancestry mating via optiSel
  # @param min.coanc.mating.solver Function with the solver for optiSel (default: use of lpsymphony::lpsymphony_solve_LP if available, else "default" optiSel solver)
  # @param relationship.matrix.min.coanc  Method to calculate relationship matrix for min coancestry (Default: "pedigree")

  if(use.recalculate.manual){
    bv.ignore.traits = 1:population$info$bv.nr
  }

  if(length(selection.criteria)>0){
    selection.criteria[selection.criteria=="ebv"] = "bve"
    selection.criteria[selection.criteria=="gv"] = "bv"
  }
  if(length(multiple.bve.scale.m)>0){
    multiple.bve.scale.m[multiple.bve.scale.m=="ebv_sd"] = "bve_sd"
    multiple.bve.scale.m[multiple.bve.scale.m=="gv_sd"] = "bv_sd"
  }
  if(length(multiple.bve.scale.f)>0){
    multiple.bve.scale.f[multiple.bve.scale.f=="ebv_sd"] = "bve_sd"
    multiple.bve.scale.f[multiple.bve.scale.f=="gv_sd"] = "bv_sd"
  }




  if(length(best.selection.criteria.m)>0){
    best.selection.criteria.m[best.selection.criteria.m=="ebv"] = "bve"
    best.selection.criteria.m[best.selection.criteria.m=="gv"] = "bv"
  }

  if(length(best.selection.criteria.f)>0){
    best.selection.criteria.f[best.selection.criteria.f=="ebv"] = "bve"
    best.selection.criteria.f[best.selection.criteria.f=="gv"] = "bv"
  }

  if(length(threshold.selection.criteria)>0){
    threshold.selection.criteria[threshold.selection.criteria=="ebv"] = "bve"
    threshold.selection.criteria[threshold.selection.criteria=="gv"] = "bv"
  }

  if(length(bve.mean.between)>0){
    bve.mean.between[bve.mean.between=="ebv"] = "bve"
    bve.mean.between[bve.mean.between=="ebvpheno"] = "bvepheno"
    bve.mean.between[bve.mean.between=="gv"] = "bv"
  }

  if(length(y.gwas.used)>0){
    y.gwas.used[y.gwas.used=="ebv"] = "bve"
    y.gwas.used[y.gwas.used=="gv"] = "bv"
  }

  if(length(selection.criteria)>0){
    selection.criteria[selection.criteria=="ebv"] = "bve"
    selection.criteria[selection.criteria=="gv"] = "bv"
  }
  if(length(selection.criteria)>0){
    selection.criteria[selection.criteria=="ebv"] = "bve"
    selection.criteria[selection.criteria=="gv"] = "bv"
  }












  if(length(recombination.minimum.distance)>0){
    recombination.function = function(noc, length.genome){
      if(noc==0){
        return(numeric(0))
      }
      min_distance = recombination.minimum.distance
      rec = rep(-Inf, noc)
      index = 1

      rep = 0
      while(index <= noc){
        rep = rep+1
        new_rec = stats::runif(1, 0, length.genome)

        if(min(abs(new_rec - rec[1:index]))> min_distance){
          rec[index] = new_rec
          index = index + 1
        }

        if(rep > noc*10){
          cat(paste0("number of recombination points was reduced to ", index - 1))
          rec = rec[1:(index-1)]
          index = Inf
        }
      }

      return(rec)
    }
  } else if(length(recombination.distance.penalty)>0){
    recombination.function = function(noc, length.genome){

      if(noc==0){
        return(numeric(0))
      }

      min_distance = recombination.distance.penalty
      rec = rep(-Inf, noc)
      index = 1

      rep = 0
      while(index <= noc){
        rep = rep+1
        new_rec = stats::runif(1, 0, length.genome)

        dist = min(abs(new_rec - rec[1:index]))
        if(dist > min_distance || stats::rbinom(1,1,dist/min_distance) == 1){
          rec[index] = new_rec
          index = index + 1
        }

        if(rep > noc*10){
          cat(paste0("number of recombination points was reduced to ", index - 1))
          rec = rec[1:(index-1)]
          index = Inf
        }
      }


      return(rec)
    }
  } else if(length(recombination.distance.penalty.2)>0){
    recombination.function = function(noc, length.genome){

      if(noc==0){
        return(numeric(0))
      }

      min_distance = recombination.distance.penalty.2
      rec = rep(-Inf, noc)
      index = 1

      rep = 0
      while(index <= noc){
        rep = rep+1
        new_rec = stats::runif(1, 0, length.genome)

        dist = min(abs(new_rec - rec[1:index]))
        if(dist > min_distance || stats::rbinom(1,1,(dist/min_distance)^2) == 1){
          rec[index] = new_rec
          index = index + 1
        }

        if(rep > noc*10){
          cat(paste0("number of recombination points was reduced to ", index - 1))
          rec = rec[1:(index-1)]
          index = Inf
        }
      }


      return(rec)
    }
  } else if(length(recombination.function)==0){
    recombination.function = recombination.function.haldane
  }

  if(length(population$info$pool_effects) == 0){
    population$info$pool_effects = FALSE
  }
  if(mixblup.bve){
    if(length(mixblup.omega)==0){
      mixblup.omega = mixblup.lambda
    }

    if(length(mixblup.beta)==0 && length(mixblup.alpha)==1){
      mixblup.beta = 1 - mixblup.alpha
    }
    if(length(mixblup.alpha)==0 && length(mixblup.beta)==1){
      mixblup.alpha = 1 - mixblup.beta
    }
    if(length(mixblup.alpha)==0){
      mixblup.alpha = 1
    }
    if(length(mixblup.beta)==0){
      mixblup.beta = 0
    }


    if(length(mixblup.residual.cov)>0 && !is.matrix(mixblup.residual.cov)){
      if(length(mixblup.residual.cov)==population$info$bv.nr){
        mixblup.residual.cov = diag(mixblup.residual.cov)
      } else{
        mixblup.residual.cov = matrix(mixblup.residual.cov, ncol= sqrt(length(mixblup.residual.cov)))
      }

    }

    if(length(mixblup.genetic.cov)>0 && !is.matrix(mixblup.genetic.cov)){
      if(length(mixblup.genetic.cov)==population$info$bv.nr){
        mixblup.genetic.cov = diag(mixblup.genetic.cov)
      } else{
        mixblup.genetic.cov = matrix(mixblup.genetic.cov, ncol= sqrt(length(mixblup.genetic.cov)))
      }
    }

  }

  if(length(threshold.selection.index)>0){
    if(!is.matrix(threshold.selection.index)){
      threshold.selection.index = matrix(threshold.selection.index, ncol= population$info$bv.nr)
    }
  }
  if(length(threshold.selection.index)>0){
    if(length(threshold.selection.value) != nrow(threshold.selection.index)){
      stop("Please provide threshold values for all threshold indices provided!")
    }
  }

  if(length(threshold.selection.index)>0){
    if(length(threshold.selection.sign) != nrow(threshold.selection.index)){
      threshold.selection.sign = rep(threshold.selection.sign, length.out = nrow(threshold.selection.index))
    }
    if(length(threshold.selection.criteria) != nrow(threshold.selection.index)){
      threshold.selection.criteria = rep(threshold.selection.criteria, length.out = nrow(threshold.selection.index))
    }
  }


{


  selection.highest = as.logical(selection.highest)

  if(length(population$info$founder_multi)==0){
    population$info$founder_pools = 1
    population$info$founder_multi = FALSE
  }
  if(length(new.class)==1){
    new.class = rep(new.class,2)
  }
  if(length(size.scaling)>0){
    population$info$size.scaling = size.scaling
  }
  if(length(population$info$size.scaling)==0){
    population$info$size.scaling = 1
  }

  if(population$info$size.scaling != 1){
    breeding.size = ceiling(breeding.size * population$info$size.scaling)
    selection.size = ceiling(selection.size * population$info$size.scaling)
  }

  ###
  # preparations for MiXBLUP
  ###

  if(mixblup.bve){
    if(length(mixblup.path)==0){
      if(Sys.info()[['sysname']]=="Windows"){
        mixblup.path = "MiXBLUP.exe"
      } else{
        mixblup.path = "./MiXBLUP.exe"
      }
    }

    if(sum(dir() == mixblup.files)==0){
      if(verbose){ cat(paste0("Create temporary folder for MiXBLUP files: ",   mixblup.files, "/n"))}
      dir.create(mixblup.files)
    }

    if(length(mixblup.path.pedfile)>0){
      mixblup.pedfile = FALSE
    }
    if(length(mixblup.path.parfile)>0){
      mixblup.parfile = FALSE
    }
    if(length(mixblup.path.datafile)>0){
      mixblup.datafile = FALSE
    }
    if(length(mixblup.path.inputfile)>0){
      mixblup.inputfile = FALSE
    }
    if(length(mixblup.path.genofile)>0){
      mixblup.genofile = FALSE
    }

    if(mixblup.inputfile){
      mixblup.path.inputfile = paste0(mixblup.files, "/", "InpMiXBLUP.txt")
    }
    if(mixblup.datafile){
      mixblup.path.datafile = paste0(mixblup.files, "/", "data.txt")
    }
    if(mixblup.pedfile){
      mixblup.path.pedfile = paste0(mixblup.files, "/", "pedigree.txt")
    }
    if(mixblup.parfile){
      mixblup.path.parfile = paste0(mixblup.files, "/", "parfile.txt")
    }
    if(mixblup.genofile){
      mixblup.path.genofile = paste0(mixblup.files, "/", "genofile.txt")
    }

  }

}

  #######################################################################
  ## Pre-work to allow for flexiblity when inputing parameter values ####
  # Initialisize parameters that were not initialized in early versions #
  #######################################################################
  {


    if(parallel.internal){
      generation.cores <- 1
    } else if(length(generation.cores)==1 && is.numeric(generation.cores)){
      population$info$generation.cores <- generation.cores
      if(verbose) cat("New default for number of cores for generation of individuals set.\n")

    } else if(length(population$info$generation.cores)>0 ){
      generation.cores <- population$info$generation.cores
    } else{
      generation.cores <- 1
    }

    if(copy.individual || copy.individual.f || copy.individual.m){
      generation.cores <- 1
    }


    if(generation.cores > 1){
      population_parallel <- population
    }

    if(generation.cores>1){

      if (requireNamespace("miraculix", quietly = TRUE)) {

        # This is not how CRAN wants it but not sure how else to do it with RandomFieldsUtils...
        library(RandomFieldsUtils)
        if(RandomFieldsUtils::RFoptions()$basic$cores>1){
          warning("RFoptions cores has been set to 1")
          RandomFieldsUtils::RFoptions(cores=1)
        }
      }


    }


    if(intern.func==0){
      breeding.intern.activ <- breeding.intern
    } else if(intern.func==1){
      breeding.intern.activ <- breeding.intern1
    } else if(intern.func==2){
      breeding.intern.activ <- breeding.intern2
    }



    if(length(pen.assignments)==0 && length(population$info$pen.effect.covariance)>0 && sum(abs(population$info$pen.effect.covariance))>0){
      pen.assignments <- TRUE
    } else{
      pen.assignments <- FALSE
    }

    if(length(pen.size)==0){
      pen.size <- population$info$pen.size
    } else{
      if(length(pen.size)==1){
        pen.size <- cbind(pen.size, 1)
      }
      if(!is.matrix(pen.size)){
        pen.size <- cbind(pen.size, 1/length(pen.size))
      }

      if(verbose & pen.size.overwrite){
        warning("New standard for pen size set. This will be the new default for all downstream generation of offspring.")
        population$info$pen.size <- pen.size
      }

    }


    if(length(forecast.sigma.g)==0){
      if(length(sigma.g)!=0){
        forecast.sigma.g <- FALSE
      } else{
        forecast.sigma.g <- TRUE
      }
    }

    if(length(sigma.g)==0){
      sigma.g <- 10
    }


    if(length(population$info$one.sex.mode)>0 && population$info$one.sex.mode){
      breeding.sex <- 0
    } else{
      population$info$one.sex.mode <- FALSE
    }

    if(export.relationship.matrix){
      bve <- TRUE
    }

    if(length(population$info$litter.effect.active)==0){
      if(length( population$info$litter.effect.covariance)==0 || sum( abs(population$info$litter.effect.covariance))==0){
        population$info$litter.effect.active <- FALSE
      } else{
        population$info$litter.effect.active <- TRUE
      }
    }

    if(length(population$info$pen.effect.active)==0){
      if(length( population$info$pen.effect.covariance)==0 || sum( abs(population$info$pen.effect.covariance))==0){
        population$info$pen.effect.active <- FALSE
      } else{
        population$info$pen.effect.active <- TRUE
      }
    }


    if(length(bve.array)>0){
      temp1 <- which(population$info$array.name == bve.array)
      if(length(temp1)==1){
        bve.array <- temp1
      }
      if(!is.numeric(bve.array)){
        stop("Input for bve.array can not be assigned! Check your input!")
      }
    }

    if((genotyped.array)!=1){
      temp1 <- which(population$info$array.name == genotyped.array)
      if(length(temp1)==1){
        genotyped.array <- temp1
      }
      if(!is.numeric(genotyped.array)){
        stop("Input for genotyped.array can not be assigned! Check your input!")
      }
    }

    if(length(input.phenotype)==0){
      input.phenotype <- bve.input.phenotype
    }

    if(length(offpheno.parents.gen)>0){
      offspring.bve.parents.gen <- offpheno.parents.gen
    }
    if(length(offpheno.parents.database)>0){
      offspring.bve.parents.database <- offpheno.parents.database
    }
    if(length(offpheno.parents.cohorts)>0){
      offspring.bve.parents.cohorts <- offpheno.parents.cohorts
    }
    if(length(offpheno.offspring.gen)>0){
      offspring.bve.offspring.gen <- offpheno.offspring.gen
    }
    if(length(offpheno.offspring.database)>0){
      offspring.bve.offspring.database <- offpheno.offspring.database
    }
    if(length(offpheno.offspring.cohorts)>0){
      offspring.bve.offspring.cohorts <- offpheno.offspring.cohorts
    }

    if(is.function(bve.solve) || bve.solve != "exact"){

      if(!is.function(bve.solve) && bve.solve =="pcg"){
        if (requireNamespace("cPCG", quietly = TRUE)) {
          bve.solve <- cPCG::pcgsolve
        } else{
          stop("Selected solver not available!")
        }

      }

    }

    if(length(population$info$default.parameter.name)>0){
      for(index in 1:length(population$info$default.parameter.name)){
        assign(population$info$default.parameter.name[index], value = population$info$default.parameter.value[[index]])
      }
    }


    if(length(bv.ignore.traits)>0){
      temp123 <- setdiff(population$info$bv.random.activ , bv.ignore.traits)
    } else{
      temp123 <- population$info$bv.random.activ
    }

    {

      if(avoid.mating.halfsib){
        max_rel = 0
      } else if(avoid.mating.fullsib){
        max_rel = 1
      } else{
        max_rel = 2
      }


      repeat.mating.store <- population$info$repeat.mating
      repeat.mating.copy.store <- population$info$repeat.mating.copy

      if(length(repeat.mating)>0){

        temp1 <- population$info$repeat.mating

        if(length(repeat.mating)==1){

          if(repeat.mating != "genetic"){
            population$info$repeat.mating <- cbind(repeat.mating, 1)
          } else{
            population$info$repeat.mating <- repeat.mating
            population$info$repeat.mating.trait <- repeat.mating.trait
            if(length(repeat.mating.max)==0){
              repeat.mating.max = 100
              warning("Assume a maximum litter size of 100. Provide a realistic value in repeat.mating.max for efficiency.")
            }
            population$info$repeat.mating.max <- repeat.mating.max
          }

        } else{
          population$info$repeat.mating <- repeat.mating
        }
        if(verbose & repeat.mating.overwrite & ((length(temp1)==0) || (length(temp1)!=length(population$info$repeat.mating)) || prod(temp1==population$info$repeat.mating)!=1) ) warning("New standard for litter size / repeat.mating set. This will be the new default for all downstream generation of offspring via breeding")
      }

      if(length(population$info$repeat.mating)==0){
        population$info$repeat.mating <- cbind(1,1)
      }

      if(length(repeat.mating.copy)>0){
        if(length(repeat.mating.copy)==1){
          population$info$repeat.mating.copy <- cbind(repeat.mating.copy, 1)
        } else{
          population$info$repeat.mating.copy <- repeat.mating.copy
        }
        if(verbose & repeat.mating.overwrite) warning("New standard for litter size / repeat.mating.copy set. This will be the new default for all downstream generation of offspring via copy/combine.")
      }
      if(length(population$info$repeat.mating.copy)==0){
        population$info$repeat.mating.copy <- cbind(1,1)
      }
    }

    if(copy.individual){
      repeat.mating.activ <- population$info$repeat.mating.copy
    } else{
      repeat.mating.activ <- population$info$repeat.mating
    }


    if(length(breeding.size.litter)>0 && sum(breeding.size)==0){

      if(length(repeat.mating.activ)==1 && repeat.mating.activ=="genetic"){
        repeat.mating.s = rep(population$info$repeat.mating.max, breeding.size.litter)
        breeding.size = sum(repeat.mating.s)
        if(verbose){cat(paste0("Litter size has underlying genetics. Prepare space for a maximum of  ", breeding.size, " individuals from ", breeding.size.litter, " litters.\n"))}

      } else{
        repeat.mating.s = sample(repeat.mating.activ[,1], prob = repeat.mating.activ[,2], breeding.size.litter, replace = TRUE)
        breeding.size = sum(repeat.mating.s)
        if(verbose){cat(paste0("Generate ", breeding.size, " individuals from ", breeding.size.litter, " litters.\n"))}
      }

    }

    population$info$neff <- list()
    if(length(population$info$real.bv.add)>1){
      for(index in 1:(length(population$info$real.bv.add)-1)){
        if(length(population$info$real.bv.add[[index]])>0){
          population$info$neff[[index]] <- 1:nrow(population$info$real.bv.add[[index]])
        }
      }
    }



    if(Rprof){
      Rprof()
    }

    {
      # foreach variables
      indexb <- NULL
      keeps <- NULL
      pheno <- NULL
      if (requireNamespace("foreach", quietly = TRUE)) {
        `%dopar%` <- foreach::`%dopar%`
      }


    }


    if(sum(population$info$is.combi)>0){
      for(index in 1:length(population$info$combi.weights)){
        population$info$combi.weights[[index]] <- c(population$info$combi.weights[[index]], rep(0, population$info$bv.nr - length(population$info$combi.weights[[index]])))
      }
    }
    if(length(population$info$array.name)==0){
      population$info$array.name = "Full_Array"
      population$info$array.markers = list(rep(TRUE,sum(population$info$snp)))
      population$info$array.is_subset = FALSE
    }

    if(dh.mating){
      selfing.mating = TRUE
    }

    reduced.selection.panel <- list(reduced.selection.panel.m, reduced.selection.panel.f)

    if(length(selection.criteria)==0){
      selection.criteria <- c("bve", "bve")
      if(length(selection.m)==0){
        if(population$info$bv.nr>0){
          selection.m <- "function"
        } else{
          selection.m <- "random"
        }

      }
    } else if(sum(selection.criteria=="random")>0){
      if(selection.criteria[1]=="random"){
        selection.m = "random"
      } else{
        selection.m = "function"
      }
      if(length(selection.criteria)==2 && selection.criteria[2]=="random" && length(selection.m)==0){
        selection.f = "random"
      }
      if(length(selection.criteria)==1 && selection.criteria[1]=="random" && length(selection.f)==0){
        selection.f = "random"
      }
    } else {
      if(length(selection.m)==0){
        selection.m <- "function"
      }
    }
    if(culling.share1>0 || (length(culling.share2)>0 && culling.share2>0)){
      culling <- TRUE
    } else{
      culling <- FALSE
    }
    if(length(new.bv.child)>0){
      phenotyping.child <- new.bv.child
    }

    if(length(computation.A)>0){
      relationship.matrix <- computation.A
    }
    if(length(computation.A.ogc)){
      relationship.matrix.ogc <- computation.A.ogc
    }


    if(length(randomSeed)>0){
      set.seed(randomSeed)
    }

    if(length(bve.database)==0 && length(bve.gen)==0 && length(bve.cohorts)==0){
      bve.gen <- nrow(population$info$size)
    }

    if(length(population$info$phenotypic.transform)==0){
      population$info$phenotypic.transform <- rep(FALSE, population$info$bv.nr)
    }
    activ.trafo <- which(population$info$phenotypic.transform)


    # Fill databases

    genotyped.database <- get.database(population, genotyped.gen, genotyped.database, genotyped.cohorts)

    genotyped.remove.database <- get.database(population, genotyped.remove.gen, genotyped.remove.database, genotyped.remove.cohorts)

    bve.gen.input <- bve.gen
    bve.database.input <- bve.database
    bve.cohorts.input <- bve.cohorts

    bve.database <- get.database(population, bve.gen, bve.database, bve.cohorts) # NOT DONE

    if(length(bve.insert.gen)==0 && length(bve.insert.cohorts)==0 && length(bve.insert.database)==0){
      bve.insert.gen <- bve.gen.input
      bve.insert.cohorts <- bve.cohorts.input
      bve.insert.database <- bve.database.input
    }
    bve.insert.database <- get.database(population, bve.insert.gen, bve.insert.database, bve.insert.cohorts) # NOT DONE

    # Add Individuals in bve.insert that are not in bve.database
    for(index in 1:nrow(bve.insert.database)){
      included <- (bve.insert.database[index,1] == bve.database[,1]) & (bve.insert.database[index,2] == bve.database[,2]) & (bve.insert.database[index,3] >= bve.database[,3]) & (bve.insert.database[index,4] <= bve.database[,4])
      if(sum(included)==0){
        if(verbose) cat("Not all individuals with breeding values to insert included in breeding value estimation.\n")
        if(verbose) cat("Missing individuals are automatically added in BVE.\n")
        bve.database <- rbind(bve.database, bve.insert.database[index,])
        bve.database <- get.database(population,  database= bve.database)
      }
    }


    if(length(sigma.g.gen)==0 && length(sigma.g.cohorts)==0 && length(sigma.g.database)==0){
      sigma.g.gen <- bve.gen.input
      sigma.g.cohorts <- bve.cohorts.input
      sigma.g.database <- bve.database.input
    }
    sigma.g.database <- get.database(population, sigma.g.gen, sigma.g.database, sigma.g.cohorts) # NOT DONE

    if(length(gwas.gen)==0 && length(gwas.cohorts)==0 && length(gwas.database)==0){
      gwas.gen <- bve.gen.input
      gwas.cohorts <- bve.cohorts.input
      gwas.database <- bve.database.input
    }
    gwas.database <- get.database(population, gwas.gen, gwas.database, gwas.cohorts) # NOT DONE

    if(length(new.bv.observation)>0){
      phenotyping <- new.bv.observation
    }
    if(length(new.bv.observation.gen)>0){
      phenotyping.gen <- new.bv.observation.gen
    }
    if(length(new.bv.observation.cohorts)>0){
      phenotyping.cohorts <- new.bv.observation.cohorts
    }
    if(length(new.bv.observation.database)>0){
      phenotyping.database <- new.bv.observation.database
    }
    if(length(phenotyping)==1 && (phenotyping=="all" || phenotyping=="non_obs" )){
      phenotyping.gen <- 1:length(population$breeding)
      if(phenotyping=="non_obs"){
        multiple.observation <- FALSE
      }
    }
    if(length(phenotyping)==1 && phenotyping=="non_obs_m"){
      phenotyping.database <- cbind(1:length(population$breeding),1)
    }
    if(length(phenotyping)==1 && phenotyping=="non_obs_f"){
      phenotyping.database <- cbind(1:length(population$breeding),2)
    }

    phenotyping.database <- get.database(population, phenotyping.gen, phenotyping.database, phenotyping.cohorts)

    if(length(sigma.e.gen)==0 && length(sigma.e.cohorts)==0 && length(sigma.e.database)==0){

      if(length(phenotyping.database)>0){
        sigma.e.gen <- phenotyping.gen
        sigma.e.cohorts <- phenotyping.cohorts
        sigma.e.database <- phenotyping.database
      } else{
        sigma.e.gen <- bve.gen.input
        sigma.e.cohorts <- bve.cohorts.input
        sigma.e.database <- bve.database.input
      }

    } else{

      if(length(heritability)>0){

      } else if(length(population$info$last.sigma.e.heritability)>0){
        warning("No heritability given for sigma.e.gen/database/cohort. Automatically use last available")
        heritability <- population$info$last.sigma.e.heritability
      } else{
        if(population$info$bv.nr>0){
          stop("No heritability given for sigma.e.gen/database/cohort.")
        }

      }

    }

    if(length(heritability)>population$info$bv.nr){
      stop(paste0("There are only ", population$info$bv.nr), " traits! Check your heritability input!")
    }
    sigma.e.database <- get.database(population, sigma.e.gen, sigma.e.database, sigma.e.cohorts) # NOT DONE

    if((copy.individual.m + copy.individual.f + combine)>1){
      stop("Use of multiple copy parameter at the same time is forbidden!")
    }

    selection.size.calc <- FALSE
    if(length(selection.size)==1){
      selection.size <- rep(selection.size,2)
    }






    if(length(best1.from.cohort)){
      if(length(selection.m.cohorts)==0){
        selection.m.cohorts <- best1.from.cohort
      } else{
        selection.m.cohorts <- c(selection.m.cohorts,best1.from.cohort)
      }
    }
    if(length(best2.from.cohort)){
      if(length(selection.m.cohorts)==0){
        selection.f.cohorts <- best2.from.cohort
      } else{
        selection.f.cohorts <- c(selection.f.cohorts,best2.from.cohort)
      }
    }
    if(length(best1.from.group)){
      if(length(selection.m.database)==0){
        selection.m.database <- best1.from.group
      } else{
        selection.m.database <- c(selection.m.database,best1.from.group)
      }
    }
    if(length(best2.from.group)){
      if(length(selection.f.database)==0){
        selection.f.database <- best2.from.group
      } else{
        selection.f.database <- c(selection.f.database,best2.from.group)
      }
    }

    if(population$info$bve && (!population$info$bv.calculated || length(  population$info$bv.random.activ)==0 || sum(population$info$is.combi)>0)){
      population$info$bv.random.activ <- which(population$info$bv.random[1:population$info$bv.calc]==FALSE & population$info$is.combi[1:population$info$bv.calc]==FALSE)
    }

    add.selection <- sum(breeding.size)>0 & sum(selection.size)==0
    if(length(selection.m.gen)==0 && length(selection.m.database)==0 && length(selection.m.cohorts)==0 && (selection.size[1]>0 || add.selection) &&
       length(fixed.breeding)==0 ){

      if(sum(population$info$size[,1])>0){
        if(verbose) cat("No individuals for selection provided (male side). Use last available.\n")
        selection.m.database <- cbind(max(which(population$info$size[,1]>0)),1)
      } else{
        if(verbose) cat("No individuals for selection provided (male side). Non available.\n")
      }

    }
    if(length(selection.f.gen)==0 && length(selection.f.database)==0 && length(selection.f.cohorts)==0 && (selection.size[2]>0 ||add.selection)&&
       length(fixed.breeding)==0){
      if(sum(population$info$size[,2])>0){
        if(verbose) cat("No individuals for selection provided (female side). Use last available.\n")
        selection.f.database <- cbind(max(which(population$info$size[,2]>0)),2)
      } else{
        if(population$info$one.sex.mode){
          if(verbose) cat("No individuals as second parent were provided! \n")
        } else{
          if(verbose) cat("No individuals for selection provided (female side). Non available.\n")
        }

      }

    }

    selection.m.database <- get.database(population, selection.m.gen, selection.m.database, selection.m.cohorts)
    selection.f.database <- get.database(population, selection.f.gen, selection.f.database, selection.f.cohorts)



    if(copy.individual.m){
      copy.individual <- TRUE
      selfing.mating <- TRUE
      selfing.sex <- 0

      if(sum(selection.size)==0){
        selection.size.calc <- TRUE
        selection.size <- c(sum(selection.m.database[,4] -  selection.m.database[,3] +1),0)
      }
      if(sum(breeding.size)==0){
        breeding.size <- c(selection.size[1],0)
      }
      if(selection.size[1]>=breeding.size[1]){
        max.offspring <- c(1,1)
      }
    }

    if(copy.individual.f){
      copy.individual <- TRUE
      selfing.mating <- TRUE
      selfing.sex <- 1
      if(sum(selection.size)==0){
        selection.size.calc <- TRUE
        selection.size <- c(0, sum(selection.f.database[,4] -  selection.f.database[,3] +1))
      }
      if(sum(breeding.size)==0){
        breeding.size <- c(0,selection.size[2])
      }
      if(selection.size[2]>=breeding.size[2]){
        max.offspring <- c(1,1)
      }
    }


    if(combine==TRUE){

      if(sum(breeding.size)==0){
        breeding.size <- c(0,0)
        breeding.size[1] <- sum(selection.m.database[,4] -  selection.m.database[,3] +1)
        breeding.size[2] <- sum(selection.f.database[,4] -  selection.f.database[,3] +1)
      }
      if(sum(selection.size)==0){
        selection.size.calc <- TRUE
        selection.size <- c(0,0)
        selection.size[1] <- breeding.size[1]
        selection.size[2] <- breeding.size[2]
      }
      # combine is modelled via cloning with no recombination
      copy.individual <- TRUE
      selfing.mating <- TRUE
      if(selection.size[2]>0 & selection.size[1]==0){
        selfing.sex <- 1
      } else{
        selfing.sex <- 0
      }

      class.m <- unique(c(class.m, class.f))
      best1.from.cohort <- c(best1.from.cohort, best2.from.cohort)
      best2.from.cohort <- NULL
      best1.from.group <- c(best1.from.group, best2.from.group)
      best2.from.group <- NULL
      max.offspring = c(1,1)
    }


    if(length(phenotyping.child)==0){
      if(copy.individual){
        if(length(n.observation)>0 & sum(n.observation)>0){
          phenotyping.child <- "addobs"
        } else{
          phenotyping.child <- "keep"
        }

      } else{
        phenotyping.child <- "zero"
      }

    }

    if(length(n.observation)==0){
      if(copy.individual){
        n.observation <- 0L
      } else{
        n.observation <- 1L
      }

    }
    if(length(n.observation)>0){
      n.observation <- as.integer(n.observation)
    }
    if(length(n.observation)<population$info$bv.nr){
      n.observation <- rep(n.observation, length.out=population$info$bv.nr)

    }

    if(add.gen==0){
      current.gen <- length(population$breeding)
      add.gen <- current.gen + 1
    } else{
      current.gen <- add.gen - 1
    }

    if(!copy.individual && (length(selection.m.database)>0 && sum(selection.m.database[,1]>= add.gen)>0 && (sum(breeding.size)+ sum(selection.size))>0)){
      stop("Parental individuals must be in an earlier generation than offspring! Check add.gen / selection.m.gen/database/cohort!")
    }
    if(!copy.individual && (length(selection.f.database)>0 && sum(selection.f.database[,1]>= add.gen)>0 && (sum(breeding.size)+ sum(selection.size))>0)){
      stop("Parental individuals must be in an earlier generation than offspring! Check add.gen / selection.m.gen/database/cohort!")
    }
    if(length(population$info$cumsnp)==0){
      population$info$cumsnp <- cumsum(population$info$snp)
    }
    if(length(gen.architecture.f)==0){
      gen.architecture.f=gen.architecture.m
    }
    if(length(selection.criteria)==1){
      selection.criteria <- rep(selection.criteria,2)
    }
    if(selection.criteria[1]=="random"){
      selection.m = "random"
    }
    if(selection.criteria[2]=="random"){
      selection.f = "random"
    }
    if(length(max.offspring)==1){
      max.offspring <- rep(max.offspring,2)
    }
    if(length(max.litter)==1){
      max.litter <- rep(max.litter,2)
    }
    best.selection.ratio <- list()
    best.selection.ratio[[1]] <- best.selection.ratio.m
    if(length(best.selection.ratio.f)==0){
      best.selection.ratio.f <- best.selection.ratio.m
    }
    best.selection.ratio[[2]] <- best.selection.ratio.f

    best.selection.criteria <- list()
    best.selection.criteria[[1]] <- best.selection.criteria.m
    if(length(best.selection.criteria.f)==0){
      best.selection.criteria.f <- best.selection.criteria.m
    }
    best.selection.criteria[[2]] <- best.selection.criteria.f

    best.selection.manual.ratio <- list()
    best.selection.manual.ratio[[1]] <- best.selection.manual.ratio.m
    if(length(best.selection.manual.ratio.f)==0 && length(best.selection.manual.ratio.f)!=0 ){
      best.selection.manual.ratio.f <- best.selection.manual.ratio.m
    }
    best.selection.manual.ratio[[2]] <- best.selection.manual.ratio.f
    best.selection.manual.ratio[[3]] <- "placeholder"

    if(length(population$info$origin.gen)>0){
      population$info$origin.gen <- base::as.integer(population$info$origin.gen)
    } else if(population$info$miraculix){
      population$info$origin.gen <- 1:32L
    } else{
      population$info$origin.gen <- 1:64L
    }

    if(pseudo.bve && length(pseudo.bve.accuracy) < population$info$bv.nr){
      pseudo.bve.accuracy <- rep(pseudo.bve.accuracy, length.out=population$info$bv.nr)
    }



    if(length(population$info$origin)==0){
      population$info$origin <- 1:64
    }

    if(length(miraculix)==0){
      miraculix <- population$info$miraculix
      if(length(miraculix)==0){
        stop("No input on use of miraculix provided!")
      }
    }
    if(miraculix && length(miraculix.mult)==0){
      miraculix.mult <- TRUE
    } else if(length(miraculix.mult)==0){
      miraculix.mult <- FALSE
    }
    if(store.comp.times){
      comp.times <- numeric(7)
      comp.times[1] <- as.numeric(Sys.time())
    }

    if(gene.editing.offspring || gene.editing.best){
      gene.editing <- TRUE
    } else{
      gene.editing <- FALSE
    }
    if(length(population$info$bitstoring)>0){
      bit.storing <- TRUE
      nbits <- population$info$bitstoring
    } else{
      nbits <- 30
      bit.storing <- FALSE
    }


    if(length(add.architecture)>0){
      population$info$gen.architecture[[length(population$info$gen.architecture)+1]] <- list()
      population$info$gen.architecture[[length(population$info$gen.architecture)]]$length.total <- cumsum(c(0,add.architecture[[1]]))
      population$info$gen.architecture[[length(population$info$gen.architecture)]]$snp.position <- add.architecture[[2]]

    }
    if (requireNamespace("miraculix", quietly = TRUE)) {
      if(miraculix.cores>1){
        RandomFieldsUtils::RFoptions(cores = miraculix.cores)
      }
      codeOriginsU <- miraculix::codeOrigins
      decodeOriginsU <- miraculix::decodeOrigins
    } else{
      codeOriginsU <- codeOriginsR
      decodeOriginsU <- decodeOriginsR
    }

    if(length(new.phenotype.correlation)>0){
      new.residual.correlation <- new.phenotype.correlation
    }



    if(length(new.residual.correlation)>0){
      population$info$pheno.correlation <- t(chol(new.residual.correlation))
    }
    if(length(new.breeding.correlation)>0){
      population$info$bv.correlation <- new.breeding.correlation
    }


    class <- list()
    class[[1]] <- class.m
    class[[2]] <- class.f

    if(length(sigma.e)==0 && length(population$info$last.sigma.e.value)>0){
      sigma.e <- population$info$last.sigma.e.value
    } else if(length(sigma.e)>0){
      if(length(sigma.e)==1){
        sigma.e <- rep(sigma.e, length.out = population$info$bv.nr)
      }
      if(length(sigma.e)!=population$info$bv.nr){
        stop("Invalid input for sigma.e")
      }
    } else {
      sigma.e <- 10

    }

    if(length(sigma.e)==1){
      sigma.e <- rep(sigma.e, population$info$bv.nr)
    }
    if(length(sigma.g)==1){
      sigma.g <- rep(sigma.g, population$info$bv.nr)
    }

    if(store.comp.times){
      comp.times[2] <- as.numeric(Sys.time())
    }


    if(length(selection.f)==0){
      selection.f <- selection.m
    }
    if(length(multiple.bve.weights.m)< population$info$bv.nr){
      if(length(multiple.bve.weights.m)!=1){
        warning("Number of traits does not match with number of index weights. Check your input in multiple.bve.weights.m")
      }
      multiple.bve.weights.m <- rep(multiple.bve.weights.m, length.out = population$info$bv.nr)
    }
    if(length(multiple.bve.weights.f)==0){
      multiple.bve.weights.f <- multiple.bve.weights.m
    }
    if(length(bve.ignore.traits)==1 && bve.ignore.traits=="zero"){
      bve.ignore.traits <- which(multiple.bve.weights.m==0 & multiple.bve.weights.f==0)
    }
    if(length(bve.ignore.traits)==0){
      bve.keeps <- 1:population$info$bv.nr
    } else{
      bve.keeps <- (1:population$info$bv.nr)[-bve.ignore.traits]
    }
    if(length(multiple.bve.weights.f)< population$info$bv.nr){
      multiple.bve.weights.f <- rep(multiple.bve.weights.f, length.out = population$info$bv.nr)
    }
    if(length(multiple.bve.scale.f)==0){
      multiple.bve.scale.f <- multiple.bve.scale.m
    }

    if(length(selection.f.miesenberger)==0){
      selection.f.miesenberger <- selection.m.miesenberger
    }


    if(length(ignore.best)==1){
      ignore.best <- rep(ignore.best,2)
    }
    if(length(breeding.sex)==0 && length(breeding.size)==2 && sum(breeding.size)>0){
      breeding.sex <- breeding.size[2] / sum(breeding.size)
    } else if(length(breeding.sex)==0){
      breeding.sex <- 0.5
    }
    if(length(breeding.size)==1){
      breeding.size.total <- breeding.size
      breeding.size.m <- round(breeding.size*(1-breeding.sex))
      if(breeding.sex.random){
        breeding.size.m <- sum(stats::rbinom(breeding.size,1,(1-breeding.sex)))
      }
      breeding.size <- c(breeding.size.m, breeding.size - breeding.size.m)
    } else{
      breeding.size.total <- sum(breeding.size)
    }


    if(length(sex.s)==0){
      sex.animal <- rep(2, breeding.size.total)
      if(breeding.size[1] > 0){
        sex.animalr <- sample(1:breeding.size.total, breeding.size[1])
        sex.animal[sex.animalr] <- 1
      }
    } else{
      sex.animal <- rep(sex.s, length.out = breeding.size.total)
      if(sum(sex.animal==2) != breeding.size[2]){
        stop("Chosen sex ratio in sex.s does not match with breeding.size")
      }
    }


    if(breeding.all.combination){
      if(sum(selection.size)==0){
        selection.size[1] <- sum(selection.m.database[,4]-selection.m.database[,3]+1)
        selection.size[2] <- sum(selection.f.database[,4]-selection.f.database[,3]+1)
      }
      if(selection.size[1]>0 && selection.size[2]>0){
        fixed.breeding.best <- matrix(0, nrow= selection.size[1] * selection.size[2], ncol=4)
        for(index in 1:selection.size[1]){
          fixed.breeding.best[1:selection.size[2] + (index-1) * selection.size[2],] <- cbind(1,index, 2, 1:selection.size[2])
        }
      } else{
        nindi <- max(selection.size[1:2])
        sex.fixed <- 1 + as.numeric(selection.size[2]>0)
        fixed.breeding.best <- matrix(0, nrow= nindi * (nindi-1) / 2, ncol=4)
        prior <- 0
        for(index in 1:(nindi-1)){
          fixed.breeding.best[1:(nindi-index)+ prior,] <- cbind(sex.fixed,index, sex.fixed, (index+1):nindi)
          prior <- prior +nindi-index
        }
      }
    }

    if(length(fixed.breeding)>0 && ncol(fixed.breeding)==3){
      fixed.breeding <- cbind(fixed.breeding, fixed.breeding, breeding.sex)
    }
    if(length(fixed.breeding) >0 && ncol(fixed.breeding)==6){
      fixed.breeding <- cbind(fixed.breeding,breeding.sex)
    }
    if(length(fixed.breeding.best) >0 && ncol(fixed.breeding.best)==4){
      fixed.breeding.best <- cbind(fixed.breeding.best,breeding.sex)
    }




    if(nrow(repeat.mating.activ)==1 && length(fixed.breeding)>0 && repeat.mating.activ[1,1]>1 ){
      fixed.breeding <- matrix(rep(t(fixed.breeding), repeat.mating.activ[1,1]), ncol=7, byrow=TRUE)
    }
    if(nrow(repeat.mating.activ)==1 && length(fixed.breeding.best)>0 && repeat.mating.activ[1,1]>1 ){
      fixed.breeding.best <- matrix(rep(t(fixed.breeding.best), repeat.mating.activ[1,1]), ncol=5, byrow=TRUE)
    }


    if(length(fixed.effects.p)==0 && length(population$info$fixed.effects)>0){
      fixed.effects.p <- rep(0, ncol(population$info$fixed.effects))
    } else if(length(population$info$fixed.effects)==0){
      fixed.effects.p <- numeric(0)
    }

    if(length(fixed.effects.p)>0 & !is.matrix(fixed.effects.p)){
      if((length(fixed.effects.p) %% ncol(population$info$fixed.effects))!=0){
        stop("Invalid fixed.effects.p input!")
      } else{
        fixed.effects.p <- matrix(fixed.effects.p, byrow=TRUE, ncol = ncol(population$info$fixed.effects))
      }

    }

    if(length(fixed.effects.freq)==0 & length(fixed.effects.p)>0){
      fixed.effects.freq <- rep(1, nrow(fixed.effects.p))
    }
    if(length(fixed.effects.p)>0 && (length(fixed.effects.freq)!= nrow(fixed.effects.p))){
      stop("Number of different fixed fixed effect realisations does not match with probabilities of realisations! Check fixed.effects.freq!")
    }

  }

  #######################################################################
  ############## Start of the actual Simulation #########################
  #######################################################################
  {

    # Check if underlying true genomic values need to be calculated ((usually only on first use of breeding.diploid() or when traits are added))

    if(population$info$bv.calculated==FALSE || length(population$info$effect.p)==0 || (length(population$info$effect.p.add)==0 && length(population$info$effect.p.mult1)==0 && length(population$info$effect.p.dice)==0)){
      excludes <- NULL
      snp.before <- cumsum(c(0,population$info$snp))
      population$info$effect.p.add <- list()
      if(population$info$real.bv.length[1]>0){
        for(index in 1:population$info$real.bv.length[1]){
          if(length(population$info$real.bv.add[[index]])>0 && nrow(population$info$real.bv.add[[index]])>0){
            population$info$effect.p.add[[index]] <- as.integer(population$info$real.bv.add[[index]][,1]+ snp.before[population$info$real.bv.add[[index]][,2]])
            excludes <- unique(c(excludes, population$info$effect.p.add))

          }
        }
      }

      if(population$info$real.bv.length[2]>0){
        population$info$effect.p.mult1 <- list()
        population$info$effect.p.mult2 <- list()
        for(index in 1:population$info$real.bv.length[2]){
          if(length(population$info$real.bv.mult[[index]])>0 &&nrow(population$info$real.bv.mult[[index]])>0){
            population$info$effect.p.mult1[[index]] <- as.integer(population$info$real.bv.mult[[index]][,1]+ snp.before[population$info$real.bv.mult[[index]][,2]])
            population$info$effect.p.mult2[[index]] <- as.integer(population$info$real.bv.mult[[index]][,3]+ snp.before[population$info$real.bv.mult[[index]][,4]])
            excludes <- unique(c(excludes, population$info$effect.p.mult1,population$info$effect.p.mult2))
          }
        }
      }
      if(population$info$real.bv.length[3]>0){

        for(index in 1:population$info$real.bv.length[3]){
          if(length(population$info$real.bv.dice[[index]])>0){
            for(index2 in 1:length(population$info$real.bv.dice[[index]][[1]])){
              population$info$effect.p.dice <- as.integer(c(population$info$effect.p.dice, population$info$real.bv.dice[[index]][[1]][[index2]][,1]+ snp.before[population$info$real.bv.dice[[index]][[1]][[index2]][,2]]))
            }
            excludes <- unique(c(excludes, population$info$effect.p.dice))
          }
        }
      }

      population$info$effect.p <- as.integer(unlist(excludes))

      population$info$effect.p.add.same <- rep(FALSE, population$info$bv.nr)
      if(population$info$real.bv.length[1]>1){
        for(index in 2:population$info$real.bv.length[1]){
          if(length(population$info$effect.p.add)>=index && (length(population$info$effect.p.add[[index]]) == length(population$info$effect.p.add[[index-1]]) &&
                                                             length(population$info$effect.p.add[[index]]) > 0 &&
                                                             prod(population$info$effect.p.add[[index]] == population$info$effect.p.add[[index-1]]) == 1)){
            population$info$effect.p.add.same[index] <- TRUE
          }
        }
      }

    }

    temp1 <- population$info$bv.calculated
    if(population$info$bve && population$info$bv.calculated==FALSE){
      if(verbose) cat("Derive genomic values of founders. \n")


      for(index in 1:length(population$breeding)){
        for(sex in 1:2){
          nanimals <- length(population$breeding[[index]][[sex]])
          if(nanimals >0){
            for(nr.animal in 1:nanimals){
              activ_bv <- population$info$bv.random.activ
              if(length(population$info$bv.calculated.partly)>0){
                activ_bv <- setdiff(activ_bv, population$info$bv.calculated.partly)
              }
              if(length(activ_bv)>0){
                temp_out <- calculate.bv(population, index, sex, nr.animal,
                                         activ_bv, import.position.calculation=import.position.calculation,
                                         decodeOriginsU=decodeOriginsU,
                                         store.effect.freq=store.effect.freq,
                                         bit.storing=bit.storing, nbits=nbits, output_compressed=FALSE,
                                         bv.ignore.traits=bv.ignore.traits)
                population$breeding[[index]][[6+sex]][activ_bv,nr.animal] <- temp_out[[1]]
                population$breeding[[index]][[sex]][[nr.animal]][[25]] <- length(bv.ignore.traits)==0

                if(length(temp123)>0){
                  population$breeding[[index]][[sex]][[nr.animal]][[26]] <- temp123
                }

                if(store.effect.freq){
                  if(length(population$info$store.effect.freq) < index || length(population$info$store.effect.freq[[index]])==0){
                    colnames(temp_out[[2]]) <- c("Homo0", "Hetero", "Homo1")
                    rownames(temp_out[[2]]) <- population$info$snp.name[population$info$effect.p]
                    population$info$store.effect.freq[[index]] <- temp_out[[2]]
                  } else{
                    population$info$store.effect.freq[[index]] <- population$info$store.effect.freq[[index]] + temp_out[[2]]
                  }
                }
              }
            }
          }
        }
      }
      population$info$bv.calculated <- TRUE
    }

    # Calculate breeding values for non-QTL-based traits
    if(population$info$bv.calc>0 && population$info$bv.random[population$info$bv.calc]){
      mu1 <- numeric((population$info$bv.calc-1))
      # Estimate sigma.g for QTL-based traits ((correlation))
      if(population$info$bv.calc>1){

        n.animals <- 0
        for(index in 1:nrow(sigma.g.database)){
          n.animals <- n.animals + diff(sigma.g.database[index,3:4]) + 1
        }
        y_real <- array(0, dim=c(n.animals,(population$info$bv.nr))) # schaetzung sigma.g
        cindex <- 1
        temp2 <- 1:(population$info$bv.nr)
        temp3 <- which(population$info$bv.random==FALSE)
        for(index in 1:nrow(sigma.g.database)){
          k.database <- sigma.g.database[index,]
          if(diff(k.database[3:4])>=0){
            for(kindex in k.database[3]:k.database[4]){
              y_real[cindex,temp2] <- population$breeding[[k.database[[1]]]][[6+k.database[[2]]]][temp2, kindex]
              cindex <- cindex +1
            }
          }
        }
        for(bven in temp3){
          population$info$bv.random.variance[bven] <- stats::var(y_real[,bven])
          mu1[bven] <- mean(y_real[,bven])
        }

      }
      population$info$current.bv.correlation <- population$info$bv.correlation
      if(sum(is.na(population$info$bv.correlation))>0){
        emp_cor <- which(is.na(population$info$bv.correlation), arr.ind=TRUE)
        for(index in 1:nrow(emp_cor)){
          population$info$current.bv.correlation[emp_cor[index,1], emp_cor[index,2]] <-
            stats::cor(y_real[,emp_cor[index,1]], y_real[,emp_cor[index,2]])
        }
      }
      if(temp1==FALSE){
        population$info$current.bv.random.variance  <- population$info$bv.random.variance
        if(population$info$bv.calc==1){
          bv.var <- diag.mobps(sqrt(population$info$current.bv.random.variance)) %*%population$info$current.bv.correlation %*% diag.mobps(sqrt(population$info$current.bv.random.variance))
        } else{
          AA <- diag.mobps(sqrt(population$info$current.bv.random.variance)[1:(population$info$bv.calc-1)]) %*% population$info$current.bv.correlation[1:(population$info$bv.calc-1), 1:(population$info$bv.calc-1)]%*% diag.mobps(sqrt(population$info$current.bv.random.variance)[(1:(population$info$bv.calc-1))])
          BB <- diag.mobps(sqrt(population$info$current.bv.random.variance)[1:(population$info$bv.calc-1)]) %*%population$info$current.bv.correlation[1:(population$info$bv.calc-1), -(1:(population$info$bv.calc-1))]%*% diag.mobps(sqrt(population$info$current.bv.random.variance)[-(1:(population$info$bv.calc-1))])
          CC <- diag.mobps(sqrt(population$info$current.bv.random.variance)[-(1:(population$info$bv.calc-1))]) %*%population$info$current.bv.correlation[-(1:(population$info$bv.calc-1)), -(1:(population$info$bv.calc-1))] %*% diag.mobps(sqrt(population$info$current.bv.random.variance)[-(1:(population$info$bv.calc-1))])
          if (requireNamespace("MASS", quietly = TRUE)) {
            bv.var <- CC - t(BB) %*% MASS::ginv(AA) %*% BB
            mu_mult <- t(BB) %*% MASS::ginv(AA)
          } else{
            bv.var <- CC - t(BB) %*% solve(AA) %*% BB
            mu_mult <- t(BB) %*% solve(AA)
          }
        }
        bv.var_chol <- t(chol(bv.var))
        population$info$bv.correlation_col <- bv.var_chol
        for(index in 1:length(population$breeding)){
          for(sex in 1:2){
            nanimals <- length(population$breeding[[index]][[sex]])
            if(nanimals >0){
              for(nr.animal in 1:nanimals){
                if(population$info$bv.calc==1){
                  bv.mean <- population$info$base.bv
                } else{
                  bv.mean <- population$info$base.bv[-(1:(population$info$bv.calc-1))] +  mu_mult %*% (population$breeding[[index]][[sex+6]][1:(population$info$bv.calc-1), nr.animal] - mu1)
                }
                population$breeding[[index]][[6+sex]][population$info$bv.calc:population$info$bv.nr, nr.animal] <- bv.mean + bv.var_chol %*% stats::rnorm(population$info$bv.nr-population$info$bv.calc+1,0,1)
              }
            }
          }
        }
      }

    }

    if(temp1==FALSE && sum(population$info$is.combi)>0){
      combis <- which(population$info$is.combi)
      for(index in 1:length(population$breeding)){
        for(sex in 1:2){
          nanimals <- length(population$breeding[[index]][[sex]])
          if(nanimals > 0){
            for(combi in combis){
              population$breeding[[index]][[6+sex]][combi,] <- colSums(population$info$combi.weights[[combi]] * population$breeding[[index]][[6+sex]])
            }
          }
        }
      }

    }
    if(store.comp.times){
      comp.times[3] <- as.numeric(Sys.time())
    }


    # Compute sigma.e to fullfil a target heritablity in the reference population
    if(length(heritability)>0){

      if(length(population$info$last.sigma.e.database)==0 || !(nrow(population$info$last.sigma.e.database)==nrow(sigma.e.database) && prod(population$info$last.sigma.e.database==sigma.e.database)==1) || length(heritability)>0 && (length(population$info$last.sigma.e.heritability)==0 || prod(population$info$last.sigma.e.heritability==heritability)==0)){
        if(verbose) cat("Start deriving enviromental variance (according to given heritability).\n")
        if(length(heritability)!=population$info$bv.nr){
          heritability <- rep(heritability, length.out = population$info$bv.nr)
        }
        n.animals <- 0
        for(index in 1:nrow(sigma.e.database)){
          n.animals <- n.animals + diff(sigma.e.database[index,3:4]) + 1
        }
        y_real <- array(0, dim=c(n.animals,(population$info$bv.nr))) # schaetzung sigma.g

        y_gen <- numeric(n.animals)
        if(variance.correction=="parental.mean" || variance.correction=="generation.mean"){
          generation_mean <- matrix(0, ncol=max(sigma.e.database[,1]), nrow=population$info$bv.nr)
          for(index in unique(sigma.e.database[,1])){
            generation_mean[,index] <- rowMeans(get.bv(population, gen=index), na.rm = TRUE)
          }
        }

        if(variance.correction=="parental.mean"){
          y_p1 <- y_p2 <- array(0,dim=c(n.animals,population$info$bv.nr))
        }


        cindex <- 1
        temp1 <- 1:(population$info$bv.nr)
        for(index in 1:nrow(sigma.e.database)){
          k.database <- sigma.e.database[index,]
          if(diff(k.database[3:4])>=0){

            y_gen[cindex:(cindex+diff(k.database[3:4]))] <- k.database[1]
            for(kindex in k.database[3]:k.database[4]){
              y_real[cindex,temp1] <- population$breeding[[k.database[1]]][[6+k.database[2]]][temp1,kindex]

              if(variance.correction=="parental.mean"){

                p1 <- population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[7]]
                y_p1[cindex,] <- population$breeding[[p1[1]]][[6+p1[2]]][,p1[3]]
                p2 <- population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[8]]
                y_p2[cindex,] <- population$breeding[[p2[1]]][[6+p2[2]]][,p2[3]]

                if(p1[1]==k.database[1]){
                  y_p1[cindex,] <- generation_mean[,p1[1]]
                }

                if(p2[1]==k.database[1]){
                  y_p2[cindex,] <- generation_mean[,p2[1]]
                }


              }
              cindex <- cindex +1
            }
          }
        }



        for(bven in 1:population$info$bv.nr){
          if(forecast.sigma.g){

            if(variance.correction=="parental.mean"){
              sigma.g2.temp <- stats::var(y_real[,bven]-y_p1[,bven]/2 - y_p2[,bven]/2, na.rm = TRUE)

            } else if(variance.correction =="generation.mean"){
              sigma.g2.temp <- stats::var(y_real[,bven]-generation_mean[bven,y_gen], na.rm = TRUE)
            } else {
              sigma.g2.temp <- stats::var(y_real[,bven], na.rm = TRUE)

              if(nrow(sigma.e.database)>=2){

                n1 <- sum(sigma.e.database[1:(nrow(sigma.e.database)/2),4] - sigma.e.database[1:(nrow(sigma.e.database)/2),3 ] +1)
                n2 <- nrow(y_real) - n1
                if(n1>1 & n2>1){
                  test <- stats::t.test(y_real[1:n1,bven], y_real[-(1:n1),bven])
                  if(test$p.value<1e-10){
                    warning("Fitting of sigma.e does not account for population structure when estimating sigma.g. Consider using variance.correction for this fitting estimation")
                  }
                }


              }
            }

          }
          sigma.e[bven] <- sqrt(((1- heritability[bven]) * sigma.g2.temp)/ heritability[bven])
        }

        if(verbose){
          cat("Estimated residual variances:", round(sigma.e^2, digits=4), "\n")
        }

      } else{
        sigma.e <- population$info$last.sigma.e.value
      }

    }

    if(length(heritability)==0){
      heritability <- population$info$last.sigma.e.heritability
    }

    if(length(repeatability)>0){
      if(length(repeatability)!= population$info$bv.nr){
        repeatability <- rep(repeatability, length.out = population$info$bv.nr)
      }

      repeatability[repeatability<heritability] <-  heritability[repeatability<heritability]
      population$info$repeatability <- repeatability
    } else{

      if(length(population$info$repeatability)==0){
        population$info$repeatability <- heritability
      }

      repeatability <- population$info$repeatability

    }

    sigma.e2 <- sigma.e ^ 2

    if(length(heritability)>0 & length(repeatability)>0){

      sigma.g2.temp1 <- (sigma.e2 * heritability) /  (1 - heritability)
      sigma.total.temp1 <-  sigma.g2.temp1 + sigma.e2

      sigma.e2.perm <- repeatability * sigma.total.temp1 - sigma.g2.temp1
      sigma.e2.rest <- sigma.e2 - sigma.e2.perm

      sigma.e2.perm[sigma.e2.perm<0] <- 0
      sigma.e2.rest[sigma.e2.rest<0] <- 0

      sigma.e2.perm[is.na(sigma.e2.perm)] <- 0
      sigma.e2.rest[is.na(sigma.e2.rest)] <- 0



    } else{
      sigma.e2.perm <- rep(0, population$info$bv.nr)
      sigma.e2.rest <- sigma.e2
    }

    # Genotyping

    if(length(genotyped.database)>0){
      for(index in 1:nrow(genotyped.database)){
        if((genotyped.database[index,4]-genotyped.database[index,3])>=0){
          for(index2 in genotyped.database[index,3]:genotyped.database[index,4]){
            temp1 <- stats::rbinom(1, 1, genotyped.share)
            population$breeding[[genotyped.database[index,1]]][[genotyped.database[index,2]]][[index2]][[16]] <- max(population$breeding[[genotyped.database[index,1]]][[genotyped.database[index,2]]][[index2]][[16]], stats::rbinom(1, 1, genotyped.share))
            if(temp1==1){
              population$breeding[[genotyped.database[index,1]]][[genotyped.database[index,2]]][[index2]][[22]] <-
                c(population$breeding[[genotyped.database[index,1]]][[genotyped.database[index,2]]][[index2]][[22]], genotyped.array)
            }

          }
        }
      }
    }

    if(length(genotyped.remove.database)>0){
      for(index in 1:nrow(genotyped.remove.database)){
        if((genotyped.remove.database[index,4]-genotyped.remove.database[index,3])>=0){
          for(index2 in genotyped.remove.database[index,3]:genotyped.remove.database[index,4]){

            activ <- c(genotyped.remove.database[index,1:2], index2)
            if(genotyped.remove.all.copy){
              activ <- population$breeding[[activ[1]]][[activ[2]]][[activ[3]]][[21]]
            } else{
              activ <- t(activ)
            }

            for(index3 in 1:nrow(activ)){

              population$breeding[[activ[index3,1]]][[activ[index3,2]]][[activ[index3,3]]][[16]] <- 0
              population$breeding[[activ[index3,1]]][[activ[index3,2]]][[activ[index3,3]]][[22]] <- numeric(0)
            }

          }
        }
      }
    }

    # Phenotyping

    if(length(phenotyping.database)>0 && population$info$bve && sum(n.observation)>0){

      if(verbose) cat("Start simulating phenotypes.\n")

      for(index in 1:nrow(phenotyping.database)){
        gen <- phenotyping.database[index,1]
        sex <- phenotyping.database[index,2]
        if(diff(phenotyping.database[index,3:4])>=0){

          if(share.phenotyped<1){
            to_phenotype <- sample(phenotyping.database[index,3]:phenotyping.database[index,4], share.phenotyped *
                                     length(phenotyping.database[index,3]:phenotyping.database[index,4]))
          } else{
            to_phenotype <- phenotyping.database[index,3]:phenotyping.database[index,4]
          }

          if(length(phenotyping.class)>0){
            class_pheno <- get.class(population, database = phenotyping.database[index,])
            to_phenotype <- intersect( to_phenotype, which(class_pheno %in% phenotyping.class) + phenotyping.database[index,3] -1)
          }
          for(nr.animal in to_phenotype){

            prior_pheno <- population$breeding[[gen]][[sex]][[nr.animal]][[27]]
            if(length(fixed.effects.p)>0){
              population$breeding[[gen]][[sex]][[nr.animal]][[28]] <- fixed.effects.p[sample(1:nrow(fixed.effects.p), 1, prob= fixed.effects.freq),]
            } else{
              population$breeding[[gen]][[sex]][[nr.animal]][[28]] <- numeric(0)
            }


            if( length(population$breeding[[gen]][[sex]][[nr.animal]][[23]]) == 0){
              population$breeding[[gen]][[sex]][[nr.animal]][[23]] <- stats::rnorm(population$info$bv.nr,0,1)
            }
            multi_check <- (population$breeding[[gen]][[sex]][[nr.animal]][[15]]==0)

            n.observation_temp <- n.observation * (multi_check | multiple.observation)
            population$breeding[[gen]][[sex]][[nr.animal]][[15]] <- population$breeding[[gen]][[sex]][[nr.animal]][[15]] + n.observation_temp

            obsmax <- max(population$breeding[[gen]][[sex]][[nr.animal]][[15]])
            if(length(population$breeding[[gen]][[sex]][[nr.animal]][[24]])==0 ||
               obsmax > ncol(population$breeding[[gen]][[sex]][[nr.animal]][[24]]) &&
               obsmax > 0){

              population$breeding[[gen]][[sex]][[nr.animal]][[24]] <- cbind(population$breeding[[gen]][[sex]][[nr.animal]][[24]],
                                                                            matrix(stats::rnorm(obsmax * population$info$bv.nr - length(population$breeding[[gen]][[sex]][[nr.animal]][[24]]),0,1),
                                                                                   nrow = population$info$bv.nr))
            }

            if(!population$breeding[[gen]][[sex]][[nr.animal]][[25]]){

              if(length(population$breeding[[gen]][[sex]][[nr.animal]][[26]])==0 || sum(n.observation_temp[-population$breeding[[gen]][[sex]][[nr.animal]][[26]]]>0)>0){

                to_ignore_temp <- which(n.observation_temp==0)

                activ_bv <- population$info$bv.random.activ


                temp1234 <- setdiff(population$info$bv.random.activ , to_ignore_temp)

                temp_out <- calculate.bv(population, gen, sex, nr.animal,
                                         activ_bv, import.position.calculation=import.position.calculation,
                                         decodeOriginsU=decodeOriginsU,
                                         store.effect.freq=store.effect.freq,
                                         bit.storing=bit.storing, nbits=nbits, output_compressed=FALSE,
                                         bv.ignore.traits=to_ignore_temp)

                if(length(population$breeding[[gen]][[sex]][[nr.animal]][[26]])>0){
                  activ_replace <- !duplicated(c(population$breeding[[gen]][[sex]][[nr.animal]][[26]], activ_bv))[-(1:length(population$breeding[[gen]][[sex]][[nr.animal]][[26]]))]
                } else{
                  activ_replace <- rep(TRUE, length(activ_bv))
                }

                population$breeding[[gen]][[6+sex]][activ_bv[activ_replace],nr.animal] <- temp_out[[1]][activ_replace]

                if(length(c(population$breeding[[gen]][[sex]][[nr.animal]][[26]], temp1234))>0){
                  population$breeding[[gen]][[sex]][[nr.animal]][[26]] <- c(population$breeding[[gen]][[sex]][[nr.animal]][[26]], temp1234)
                }

                if(length(population$breeding[[gen]][[sex]][[nr.animal]][[26]])==population$info$bv.nr){
                  population$breeding[[gen]][[sex]][[nr.animal]][[25]] <- TRUE
                }


              }



            }
            for(bven in intersect(which(n.observation_temp>0), setdiff(1:population$info$bv.nr, activ.trafo))){
              if(population$breeding[[gen]][[sex]][[nr.animal]][[15]][bven]>=1){

                temp1 <- (sqrt(sigma.e2.rest) * (population$info$pheno.correlation %*% population$breeding[[gen]][[sex]][[nr.animal]][[24]][,1:population$breeding[[gen]][[sex]][[nr.animal]][[15]][bven]]) +
                            as.numeric(sqrt(sigma.e2.perm) * population$info$pheno.correlation %*% population$breeding[[gen]][[sex]][[nr.animal]][[23]]) +
                            population$breeding[[gen]][[6+sex]][, nr.animal])[bven,]  + sum(population$info$fixed.effects[bven,] * population$breeding[[gen]][[sex]][[nr.animal]][[28]])

                if(length(temp1)>0){
                  prior_pheno[[bven]] <- temp1
                }


                population$breeding[[gen]][[8+sex]][bven, nr.animal] <- mean(temp1)
              }

            }
            for(bven in intersect( which(n.observation_temp>0), intersect(1:population$info$bv.nr, activ.trafo))){
              if(population$breeding[[gen]][[sex]][[nr.animal]][[15]][bven]>=1){
                new_pheno <- (sqrt(sigma.e2.rest) * population$info$pheno.correlation %*% population$breeding[[gen]][[sex]][[nr.animal]][[24]][,1:population$breeding[[gen]][[sex]][[nr.animal]][[15]][bven]])[bven,] +
                  (sqrt(sigma.e2.perm) * population$info$pheno.correlation %*% population$breeding[[gen]][[sex]][[nr.animal]][[23]])[bven] +
                  population$breeding[[gen]][[6+sex]][bven, nr.animal] + sum(population$info$fixed.effects[bven,] * population$breeding[[gen]][[sex]][[nr.animal]][[28]])

                temp1 <- sapply(new_pheno, FUN = population$info$phenotypic.transform.function[[bven]])

                if(length(temp1)>0){
                  prior_pheno[[bven]] <- temp1
                }

                population$breeding[[gen]][[8+sex]][bven, nr.animal] <- mean(temp1)
              }
            }

            population$breeding[[gen]][[sex]][[nr.animal]][[27]] <- prior_pheno

          }
        }
      }
    }


    # Import offspring phenotypes for parents
    offspring.bve.parents.database <- get.database(population, offspring.bve.parents.gen, offspring.bve.parents.database, offspring.bve.parents.cohorts)
    if(length(offspring.bve.parents.database)>0){
      offspring.bve <- TRUE
    } else{
      offspring.bve <- FALSE
    }

    if(length(offspring.bve.offspring.gen)>0 || length(offspring.bve.offspring.database)>0 || length(offspring.bve.offspring.cohorts)>0){
      offspring.bve.offspring.database <- get.database(population, offspring.bve.offspring.gen, offspring.bve.offspring.database, offspring.bve.offspring.cohorts)

      first_gen <- min(offspring.bve.parents.database[,1])
      list_of_copy <- list()
      for(gen_check in 1:nrow(offspring.bve.parents.database)){
        for(gen_check2 in offspring.bve.parents.database[gen_check,3]:offspring.bve.parents.database[gen_check,4]){
          first_gen <- min(first_gen, population$breeding[[offspring.bve.parents.database[gen_check,1]]][[offspring.bve.parents.database[gen_check,2]]][[gen_check2]][[21]][,1])
          list_of_copy[[length(list_of_copy)+1]] <- rbind(offspring.bve.parents.database[gen_check,1], offspring.bve.parents.database[gen_check,2], gen_check2 ,
                                                          t(population$breeding[[offspring.bve.parents.database[gen_check,1]]][[offspring.bve.parents.database[gen_check,2]]][[gen_check2]][[21]]), deparse.level = 0)
        }
      }
      list_of_copy <- matrix(unlist(list_of_copy), ncol=6, byrow=TRUE)
      list_of_copy <- list_of_copy[!(list_of_copy[,1]==list_of_copy[,4] & list_of_copy[,2]==list_of_copy[,5] & list_of_copy[,3]==list_of_copy[,6]),]

    } else if(offspring.bve){
      if(verbose) cat("No potential offspring for phenotype import given. Consider all potential individuals.\n")
      first_gen <- min(offspring.bve.parents.database[,1])
      list_of_copy <- list()
      for(gen_check in 1:nrow(offspring.bve.parents.database)){
        for(gen_check2 in offspring.bve.parents.database[gen_check,3]:offspring.bve.parents.database[gen_check,4]){
          first_gen <- min(first_gen, population$breeding[[offspring.bve.parents.database[gen_check,1]]][[offspring.bve.parents.database[gen_check,2]]][[gen_check2]][[21]][,1])
          list_of_copy[[length(list_of_copy)+1]] <- rbind(offspring.bve.parents.database[gen_check,1], offspring.bve.parents.database[gen_check,2], gen_check2 ,
                                                          t(population$breeding[[offspring.bve.parents.database[gen_check,1]]][[offspring.bve.parents.database[gen_check,2]]][[gen_check2]][[21]]), deparse.level = 0)
        }
      }
      list_of_copy <- matrix(unlist(list_of_copy), ncol=6, byrow=TRUE)
      list_of_copy <- list_of_copy[!(list_of_copy[,1]==list_of_copy[,4] & list_of_copy[,2]==list_of_copy[,5] & list_of_copy[,3]==list_of_copy[,6]),]
      ped_off <- get.pedigree(population, gen = first_gen:nrow(population$info$size), raw=TRUE)
      candidates <- rep(TRUE,nrow(ped_off))

      # Remove Individuals themselves
      for(index in 1:nrow(offspring.bve.parents.database)){
        act <- offspring.bve.parents.database[index,]
        candidates[ped_off[,1]==act[1] & ped_off[,2] == act[2] & ped_off[,3]>= act[3] & ped_off[,6]<=act[4]] <- FALSE
      }

      # Remove founder
      candidates[ped_off[,1]==ped_off[,4] & ped_off[,1]==ped_off[,7] &
                   ped_off[,2]==ped_off[,5] & ped_off[,2]==ped_off[,8] &
                   ped_off[,3]==ped_off[,6] & ped_off[,3]==ped_off[,9]] <- FALSE

      ped_off <- ped_off[candidates,]

      #candidates <- rep(FALSE,nrow(ped_off))
      #for(index in 1:nrow(offspring.bve.parents.database)){
      #  act <- offspring.bve.parents.database[index,]
      #  candidates[ped_off[,4]==act[1] & ped_off[,5] == act[2] & ped_off[,6]>= act[3] & ped_off[,6]<=act[4]] <- TRUE
      #  candidates[ped_off[,7]==act[1] & ped_off[,8] == act[2] & ped_off[,9]>= act[3] & ped_off[,9]<=act[4]] <- TRUE
      #}
      candidates <- TRUE
      male_candidates <- (ped_off[candidates & ped_off[,2]==1 ,])
      female_candidates <- (ped_off[candidates & ped_off[,2]==2 ,])

      male_gen <- unique(male_candidates[,1])
      female_gen  <- unique(female_candidates[,1])

      offspring.bve.offspring.database <- NULL
      if(length(male_gen)>0){
        for(index in male_gen){
          offspring.bve.offspring.database <- rbind(offspring.bve.offspring.database,
                                                    c(index,1,min(male_candidates[male_candidates[,1]==index,3]),
                                                      max(male_candidates[male_candidates[,1]==index,3])))
        }
      }
      if(length(female_gen)>0){
        for(index in female_gen){
          offspring.bve.offspring.database <- rbind(offspring.bve.offspring.database,
                                                    c(index,2,min(female_candidates[female_candidates[,1]==index,3]),
                                                      max(female_candidates[female_candidates[,1]==index,3])))
        }
      }

    }

    if(offspring.bve){
      offspring_id <- get.id(population, offspring.bve.offspring.database)
      is_double <- duplicated(offspring_id[length(offspring_id):1])[length(offspring_id):1]
      if(verbose) cat("Import phenotypes of offspring.\n")
      for(index in 1:nrow(offspring.bve.parents.database)){
        activ.parents <- offspring.bve.parents.database[index,]
        new.bv <- counter <- matrix(0, nrow=population$info$bv.nr, ncol=activ.parents[4]-activ.parents[3]+1)

        next_indi <- 1
        for(index2 in 1:nrow(offspring.bve.offspring.database)){
          activ.offspring <- offspring.bve.offspring.database[index2,]
          n.animals <- activ.offspring[4]-activ.offspring[3] +1
          if(n.animals>0){
            for(index3 in (activ.offspring[3]:activ.offspring[4])[!is_double[next_indi:(next_indi+n.animals-1)]]){
              parent1 <- population$breeding[[activ.offspring[1]]][[activ.offspring[2]]][[index3]][[7]]
              parent2 <- population$breeding[[activ.offspring[1]]][[activ.offspring[2]]][[index3]][[8]]
              if(length(list_of_copy)>0){
                replace1 <- which(parent1[1]==list_of_copy[,4] & parent1[2]==list_of_copy[,5] & parent1[3]==list_of_copy[,6])
                replace2 <- which(parent2[1]==list_of_copy[,4] & parent2[2]==list_of_copy[,5] & parent2[3]==list_of_copy[,6])
                if(length(replace1)>=1){
                  parent1 <- list_of_copy[replace1[1],1:3]
                }
                if(length(replace2)>=1){
                  parent2 <- list_of_copy[replace2[1],1:3]
                }
              }
              if((nrow(population$breeding[[activ.offspring[1]]][[activ.offspring[2]]][[index3]][[21]]) ==
                  nrow(population$breeding[[parent1[1]]][[parent1[2]]][[parent1[3]]][[21]])) &&
                 (prod(population$breeding[[activ.offspring[1]]][[activ.offspring[2]]][[index3]][[21]][1,] ==
                       population$breeding[[parent1[1]]][[parent1[2]]][[parent1[3]]][[21]][1,])==1)){
                parent1 <- c(-1,-1,-1)
              }
              if((nrow(population$breeding[[activ.offspring[1]]][[activ.offspring[2]]][[index3]][[21]]) ==
                  nrow(population$breeding[[parent2[1]]][[parent2[2]]][[parent2[3]]][[21]])) &&
                 (prod(population$breeding[[activ.offspring[1]]][[activ.offspring[2]]][[index3]][[21]][1,] ==
                       population$breeding[[parent2[1]]][[parent2[2]]][[parent2[3]]][[21]][1,])==1)){
                parent2 <- c(-1,-1,-1)
              }
              if(parent1[1]==activ.parents[1] && parent1[2]==activ.parents[2] && parent1[3]>= activ.parents[3] && parent1[3]<= activ.parents[4]){
                activ.take <- population$breeding[[activ.offspring[1]]][[activ.offspring[2]]][[index3]][[15]]>0 & !is.na(population$breeding[[activ.offspring[1]]][[activ.offspring[2]+8]][,index3])
                new.bv[activ.take,parent1[3] - activ.parents[3]+1] <- (new.bv[,parent1[3]- activ.parents[3]+1] + population$breeding[[activ.offspring[1]]][[activ.offspring[2]+8]][,index3] * population$breeding[[activ.offspring[1]]][[activ.offspring[2]]][[index3]][[15]])[activ.take]
                counter[activ.take,parent1[3]- activ.parents[3]+1] <- (counter[,parent1[3]- activ.parents[3]+1] + population$breeding[[activ.offspring[1]]][[activ.offspring[2]]][[index3]][[15]])[activ.take]
              }
              if(parent2[1]==activ.parents[1] && parent2[2]==activ.parents[2] && parent2[3]>= activ.parents[3] && parent2[3]<= activ.parents[4]){
                activ.take <-  population$breeding[[activ.offspring[1]]][[activ.offspring[2]]][[index3]][[15]]>0 & !is.na(population$breeding[[activ.offspring[1]]][[activ.offspring[2]+8]][,index3])
                new.bv[activ.take,parent2[3]- activ.parents[3]+1] <- (new.bv[,parent2[3]- activ.parents[3]+1] + population$breeding[[activ.offspring[1]]][[activ.offspring[2]+8]][,index3] * population$breeding[[activ.offspring[1]]][[activ.offspring[2]]][[index3]][[15]])[activ.take]
                counter[activ.take,parent2[3]- activ.parents[3]+1] <- (counter[,parent2[3]- activ.parents[3]+1] + population$breeding[[activ.offspring[1]]][[activ.offspring[2]]][[index3]][[15]])[activ.take]
              }
              if(sum(is.na(new.bv))>0){
                stop()
              }
            }
          }
          next_indi <- next_indi + n.animals


        }
        population$breeding[[activ.parents[1]]][[activ.parents[2]+26]][,activ.parents[3]:activ.parents[4]] <- new.bv / counter
        population$breeding[[activ.parents[1]]][[activ.parents[2]+28]][,activ.parents[3]:activ.parents[4]] <- counter
        if(sum(counter==0)>0){
          if(verbose) cat(paste0(sum(counter==0), " phenotype entries without valid offspring for phenotype import from offspring! Set offspring phenotype to NA.\n"))
          population$breeding[[activ.parents[1]]][[activ.parents[2]+26]][,activ.parents[3]:activ.parents[4]][counter==0] <- NA
        }

      }

      if(sum(counter)==0){
        if(input.phenotype!="own"){
          input.phenotype <- "own"
          if(verbose) cat("No phenotypes to import. Automatically set input.phenotype to own")
        }
      }
    }


    ## Culling Module
    if(culling){
      culling.database <- get.database(population, gen=culling.gen, database=culling.database, cohorts=culling.cohorts)

      if(length(culling.index)==1 && culling.index=="lastindex"){
        culling.index <- get.selectionindex(population, cohorts=culling.database)
      }
      if(population$info$bv.nr>0){
        culling.bv <- colSums(get.bv(population, database = culling.database) * culling.index)
      } else{
        culling.bv <- numeric(length(get.id(population, database= culling.database)))
      }
      n_animals <- length(culling.bv)

      if(sum(abs(culling.index))==0){
        culling.prob <- rep(culling.share1, n_animals)
      } else{
        culling.prob <- (culling.bv - culling.bv1) * (culling.share1 - culling.share2)/ (culling.bv1 - culling.bv2) + culling.share1
        culling.prob[culling.prob>1] <- 1
        culling.prob[culling.prob<0] <- 0
      }
      if(length(culling.prob)>length(culling.single)){
        culling.single <- rep(culling.single, length.out=length(culling.prob))
      }
      culling.prob <- culling.prob * culling.single
      culling.action <- stats::rbinom(n_animals,1,culling.prob)==1


      k <- 1
      for(index in 1:nrow(culling.database)){
        active.database <- culling.database[index,]
        n_animals_group <- active.database[4] - active.database[3] +1
        culling.action.group <- culling.action[k:(k + n_animals_group -1)]
        k <- k + n_animals_group

        store <- population$breeding[[active.database[1]]][[active.database[2]+4]][active.database[3]:active.database[4]]
        population$breeding[[active.database[1]]][[active.database[2]+4]][active.database[3]:active.database[4]][culling.action.group] <- (-1)
        population$breeding[[active.database[1]]][[active.database[2]+24]][active.database[3]:active.database[4]][culling.action.group] <-
          population$breeding[[active.database[1]]][[active.database[2]+22]][active.database[3]:active.database[4]][culling.action.group] + culling.time
        new_death <- population$breeding[[active.database[1]]][[active.database[2]+4]][active.database[3]:active.database[4]] != store

        if(culling.all.copy){
          for(kindex in (active.database[3]:active.database[4])[new_death]){
            active_indi <- c(active.database[1:2], kindex)
            if(nrow(population$breeding[[active_indi[1]]][[active_indi[2]]][[kindex]][[21]])){
              for(switch in 1:nrow(population$breeding[[active_indi[1]]][[active_indi[2]]][[kindex]][[21]])){
                active_copy <- population$breeding[[active_indi[1]]][[active_indi[2]]][[kindex]][[21]][switch,]

                population$breeding[[active_copy[1]]][[active_copy[2]+4]][active_copy[3]] <-
                  population$breeding[[active_indi[1]]][[active_indi[2]+4]][active_indi[3]]
                population$breeding[[active_copy[1]]][[active_copy[2]+24]][active_copy[3]] <-
                  population$breeding[[active_indi[1]]][[active_indi[2]+24]][active_indi[3]]

              }
            }
          }
        }
        n_death <- sum(new_death)

        if(n_death>0 && length(culling.cohorts)==1 && prod(active.database == get.database(population, cohorts = culling.cohorts))==1){
          population$breeding[[active.database[1]]][[active.database[2]+16]][active.database[3]:active.database[4]][which(new_death)] <- population$breeding[[active.database[1]]][[active.database[2]+10]][active.database[3]:active.database[4]][which(new_death)]
          active_cohort <- which(population$info$cohorts[,1]==culling.cohorts)
          if(length(active_cohort)>0){
            if(length(population$info$culling.stats)<active_cohort){
              population$info$culling.stats[[active_cohort]] <- cbind(culling.name, n_death, deparse.level = 0)
            } else{
              population$info$culling.stats[[active_cohort]] <- rbind(population$info$culling.stats[[active_cohort]],
                                                                      c(culling.name, n_death), deparse.level = 0)

            }
          }


        }

      }



    }
  }

  #######################################################################
  ################ Breeding value estimation ############################
  #######################################################################
  {
    ## Track computing times
    if(store.comp.times){
      comp.times[4] <- as.numeric(Sys.time())
    }

    if(store.comp.times.bve){
      comp.times.bve <- numeric(5)
      z_chol <- 0 # Cholesky decomposition
      z_uhat <- 0 # rrblup
      zcalc <- 0 # calculation of Z
      z_ped <- 0 # calculation of pedigree relationship matrix
      z_h <- 0 # calculation single-step H matrix
      comp.times.bve[1] <- as.numeric(Sys.time())
    }


    u_hat <- NULL
    gwas_hat <- NULL
    u_hat_possible <- FALSE

    ## BVE according to average performance of parents / grandparents
    if(bve.parent.mean || bve.grandparent.mean){
      addpheno <- FALSE
      if(bve.mean.between=="pheno"){
        addmean <- 8
      } else if(bve.mean.between=="bv"){
        addmean <- 6
      } else if(bve.mean.between=="bve"){
        addmean <- 2
      } else{
        addmean <- 2
        addpheno <- TRUE
      }
      for(index in 1:nrow(bve.insert.database)){
        activ.base <- bve.insert.database[index,]
        if(diff(activ.base[3:4])>=0){
          for(index2 in activ.base[3]:activ.base[4]){
            reference <- rbind(population$breeding[[activ.base[1]]][[activ.base[2]]][[index2]][[7]],population$breeding[[activ.base[1]]][[activ.base[2]]][[index2]][[8]])
            if(bve.grandparent.mean){
              reference <- rbind(population$breeding[[reference[1,1]]][[reference[1,2]]][[reference[1,3]]][[7]],
                                 population$breeding[[reference[1,1]]][[reference[1,2]]][[reference[1,3]]][[8]],
                                 population$breeding[[reference[2,1]]][[reference[2,2]]][[reference[2,3]]][[7]],
                                 population$breeding[[reference[2,1]]][[reference[2,2]]][[reference[2,3]]][[8]])
            }
            bve_import <- matrix(0, nrow=nrow(reference), ncol=population$info$bv.nr)
            for(index3 in 1:nrow(reference)){
              bve_import[index3,] <- population$breeding[[reference[index3,1]]][[reference[index3,2]+addmean]][,reference[index3,3]]
            }
            if(addpheno){
              for(index3 in 1:nrow(reference)){
                bve_import[index3,bve_import[index3,]==0] <- population$breeding[[reference[index3,1]]][[reference[index3,2]+8]][bve_import[index3,]==0,reference[index3,3]]
              }
            }
            bve_import[bve_import==0] <- NA
            population$breeding[[activ.base[1]]][[activ.base[2]+2]][,index2] <- colMeans(bve_import, na.rm=TRUE)
            if(sum(is.na(population$breeding[[activ.base[1]]][[activ.base[2]+2]][,index2]))>0){
              population$breeding[[activ.base[1]]][[activ.base[2]+2]][is.na(population$breeding[[activ.base[1]]][[activ.base[2]+2]][,index2]),index2] <- 0
            }

          }
        }
        if(report.accuracy){
          y_real_report <- NULL
          y_hat_report <- NULL
          for(index in 1:nrow(bve.insert.database)){
            activ.base <- bve.insert.database[index,]
            if(diff(activ.base[3:4])>=0){
              y_real_report <- cbind(y_real_report, population$breeding[[activ.base[1]]][[activ.base[2]+6]][,activ.base[3]:activ.base[4], drop=FALSE], deparse.level = 0)
              y_hat_report <- cbind(y_hat_report, population$breeding[[activ.base[1]]][[activ.base[2]+2]][,activ.base[3]:activ.base[4], drop=FALSE], deparse.level = 0)
            }
          }

          y_hat_report[y_hat_report==0] <- NA
          if(verbose) cat("Correlation between genetic values and BVE (parent-mean BVE):\n")
          acc <- suppressWarnings(stats::cor(t(y_real_report), t(y_hat_report), use="pairwise.complete.obs"))


          if(sum(is.na(acc))>0){
            acc[is.na(acc)] <- 0
          }
          if(verbose) cat(diag(acc))
          if(verbose) cat("\n")
        }

      }
    } else if(pseudo.bve){

      ## Simulate a breeding value estimation
      ## This lacks structure (e.g. related individuals are usally all overestimated or all underestimated)
      ## but its super fast!

      n.animals <- sum(bve.insert.database[,4]-bve.insert.database[,3]+1)
      y_real <- y_hat <- array(0, dim=c(n.animals,population$info$bv.nr)) # schaetzung sigma.g
      cindex <- 1
      for(index in 1:nrow(bve.insert.database)){
        k.database <- bve.insert.database[index,]
        if(diff(k.database[3:4])>=0){
          # t() not needed just to be safe when using multiple individuals at once later
          y_real[cindex:(cindex+k.database[4]-k.database[3]),] <- t(population$breeding[[k.database[[1]]]][[6+k.database[[2]]]][,k.database[3]:k.database[4]])
          cindex <- cindex + k.database[4] - k.database[3] +1
        }
      }
      for(index in 1:population$info$bv.nr){
        var_trait <- stats::var(y_real[,index])
        if(pseudo.bve.accuracy[index]==0){
          var_residual <- (var_trait - (pseudo.bve.accuracy[index])^2 * var_trait) / (pseudo.bve.accuracy[index])^2
          y_real[,index] <- y_hat[,index] <- 0
        } else{
          var_residual <- (var_trait - (pseudo.bve.accuracy[index])^2 * var_trait) / (pseudo.bve.accuracy[index])^2
          y_hat[,index] <- (y_real[,index] + stats::rnorm(n.animals, sd= sqrt(var_residual)))

          if(pseudo.bve.accuracy[index]<0){
            temp1 <-  -y_hat[,index]
            temp1 <- temp1 + max(y_hat[,index]) -max(temp1)
            y_hat[,index] <- temp1
          }
        }
        if(!is.matrix(y_hat)){
          y_hat <- matrix(y_hat, ncol=1)
          y_real <- matrix(y_real, ncol=1)
        }

      }

      # Enter BVE
      cindex <- 1
      for(index in 1:nrow(bve.insert.database)){
        k.database <- bve.insert.database[index,]
        if(diff(k.database[3:4])>=0){
          # t() not needed just to be safe when using multiple individuals at once later
          population$breeding[[k.database[[1]]]][[2+k.database[[2]]]][,k.database[3]:k.database[4]] <- t(y_hat[cindex:(cindex+k.database[4]-k.database[3]),])
          cindex <- cindex + k.database[4] - k.database[3] +1
        }
      }

      if(report.accuracy){
        y_hat_report <- y_hat
        y_real_report <- y_real
        y_hat_report[y_hat_report==0] <- NA
        if(verbose) cat("Correlation between genetic values and simulated BVE:\n")
        acc <- suppressWarnings(stats::cor(y_real_report, y_hat_report, use="pairwise.complete.obs"))


        if(sum(is.na(acc))>0){
          acc[is.na(acc)] <- 0
        }
        if(verbose) cat(diag(acc))
        if(verbose) cat("\n")
      }

    } else if(bve){
      ## "proper" breeding value estimation


      dense = NULL # File used for MixBLUP

      if(verbose && (relationship.matrix!="kinship" && relationship.matrix !="pedigree")) cat("Start genomic BVE.\n")
      if(verbose && (relationship.matrix=="kinship" || relationship.matrix =="pedigree")) cat("Start pedigree BVE.\n")

      ## Derive if individuals are using multiple times in the BVE
      ## Import Phenotypes from all copies of an individual (potentially use individual entry [[21]] for efficiency?)
      loop_elements_list <- derive.loop.elements(population=population, bve.database=bve.database,
                                                 bve.class=bve.class, bve.avoid.duplicates=bve.avoid.duplicates,
                                                 store.which.adding = TRUE, list.of.copys = TRUE)
      loop_elements <- loop_elements_list[[1]]
      loop_elements_copy <- loop_elements_list[[3]]

      n.animals <- nrow(loop_elements)
      genotyped <- numeric(n.animals)

      stay.loop.elements <- NULL
      n.rep <- 0
      if(length(loop_elements_copy)>0){
        n.rep <- nrow(loop_elements_copy)
      }
      # remove non-genotyped samples in case no pedigree-based estimation // single-step
      if(singlestep.active==FALSE && (relationship.matrix!="kinship" && relationship.matrix !="pedigree")){
        for(index in 1:n.animals){

          k.database <- bve.database[loop_elements[index,3],]
          kindex <- loop_elements[index,2]
          genotyped[index] <- population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[16]]
        }

        if(n.rep>0){
          genotyped.copy <- numeric(n.rep)
          for(index in 1:n.rep){
            kindex <- loop_elements_copy[index,2]
            k.database <- bve.database[loop_elements_copy[index,3],]
            if(population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[16]]==1){
              genotyped[loop_elements_copy[index,6]] <-  genotyped.copy[index] <- population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[16]]

            }
          }
          remove.loop.elements <- which(genotyped.copy==0)

          if(length(remove.loop.elements)>0 && FALSE){
            loop_elements_list[[3]] <- loop_elements_list[[3]][-remove.loop.elements,]
          }
        }
        stay.loop.elements <- which(genotyped==1)
        remove.loop.elements <- which(genotyped==0)
        if(length(remove.loop.elements)>0){
          loop_elements_list[[1]] <- loop_elements_list[[1]][-remove.loop.elements,]
          loop_elements_list[[2]] <- loop_elements_list[[2]][-remove.loop.elements]
        }

        if(length(remove.loop.elements)>0){
          warning("Use of GBLUP although some individuals are not genotyped. Non-genotyped individuals have been automatically removed from the BVE.\nConsider use ssGBLUP (singlestep.active = TRUE) oder pBLUP (relationship.matrix ='pedigree') ")
        }


        loop_elements <- loop_elements_list[[1]]
        n.animals <- nrow(loop_elements)
        loop_elements_list[[1]][,1] <- 1:n.animals
        genotyped <- numeric(n.animals)

      }

      if((relationship.matrix!="kinship" && relationship.matrix !="pedigree")){
        Zt <- array(0L,dim=c(sum(population$info$snp), n.animals))
      }
      y <- y_real <- y_real2 <- y_hat <- y_reli <- y_parent <- array(0,dim=c(n.animals,population$info$bv.nr))


      if(variance.correction=="parental.mean" || variance.correction=="generation.mean"){
        generation_mean <- matrix(0, ncol=max(bve.database[,1]), nrow=population$info$bv.nr)
        for(index in unique(bve.database[,1])){
          generation_mean[,index] <- rowMeans(get.bv(population, gen=index), na.rm = TRUE)
        }
      }

      if(variance.correction=="parental.mean"){
        y_p1 <- y_p2 <- array(0,dim=c(n.animals,population$info$bv.nr))
      }

      X_fixed <- matrix(0, nrow=n.animals,ncol=ncol(population$info$fixed.effects))


      grid.position <- numeric(n.animals)
      cindex <- 1
      size <- cumsum(c(0,as.vector(t(population$info$size))))
      y_obs <- matrix(0, nrow=n.animals, ncol=population$info$bv.nr)
      genotyped_array <- matrix(FALSE, ncol=length(population$info$array.name), nrow=n.animals)


      # Import phenotypes / genomic values
      for(index in 1:n.animals){

        k.database <- bve.database[loop_elements[index,3],]
        kindex <- loop_elements[index,2]
        y[index,] <- population$breeding[[k.database[[1]]]][[8+k.database[[2]]]][,kindex]
        y_real[index,] <- y_real2[index,] <-population$breeding[[k.database[[1]]]][[6+k.database[[2]]]][,kindex]

        if(variance.correction=="parental.mean"){

          p1 <- population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[7]]
          y_p1[index,] <- population$breeding[[p1[1]]][[6+p1[2]]][,p1[3]]
          p2 <- population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[8]]
          y_p2[index,] <- population$breeding[[p2[1]]][[6+p2[2]]][,p2[3]]

          if(p1[1]==k.database[1]){
            y_p1[index,] <- generation_mean[,p1[1]]
          }

          if(p2[1]==k.database[1]){
            y_p2[index,] <- generation_mean[,p2[1]]
          }

        }


        if(population$breeding[[k.database[1]]][[k.database[[2]]]][[kindex]][[25]]==FALSE){
          if(length(population$breeding[[k.database[1]]][[k.database[[2]]]][[kindex]][[26]])==0){
            y_real2[index,] <- NA
          } else{
            y_real2[index,-population$breeding[[k.database[1]]][[k.database[[2]]]][[kindex]][[26]]] <- NA
          }

        }
        X_fixed[index,] <- population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[28]]
        y_obs[index,] <- population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[15]]
        grid.position[index] <- kindex + size[sum(k.database[1:2]*c(2,1))-2] # how many individuals are in earlier generations
        genotyped[index] <- population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[16]]
        genotyped_array[index, population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[22]]] <- TRUE
        if(bve.all.genotyped){
          genotyped[index] <- 1
          genotyped_array[index, ] <- TRUE
        }

        if(estimate.add.gen.var){
          father <- population$breeding[[k.database[1]]][[k.database[2]]][[kindex]][[7]]
          mother <- population$breeding[[k.database[1]]][[k.database[2]]][[kindex]][[8]]
          if(sum(father[1:3]==c(k.database[1:2], kindex))==3 ||sum(mother[1:3]==c(k.database[1:2], kindex))==3){
            if(verbose) cat("Schaetzung der additiv genetischen Varianz extrem problematisch. Kein Elterntier fuer jedes Tier vorhanden!\n")
          }
          y_parent[index,] <- mean(population$breeding[[father[1]]][[8+father[2]]][,father[3]],population$breeding[[mother[1]]][[8+mother[2]]][,mother[3]])
        }
      }

      if(n.rep>0){
        for(index in 1:n.rep){
          kindex <- loop_elements_copy[index,2]
          k.database <- bve.database[loop_elements_copy[index,3],]
          if(length(stay.loop.elements)>0){
            non_copy <- which(stay.loop.elements==loop_elements_copy[index,6])
          } else{
            non_copy <- loop_elements_copy[index,6]
          }

          if(length(non_copy)==1 && sum(population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[15]]> y_obs[non_copy,])>=1){
            switches <- population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[15]]> y_obs[non_copy,]
            y[non_copy,switches] <- population$breeding[[k.database[[1]]]][[8+k.database[[2]]]][switches,kindex]
            y_obs[non_copy,switches] <- population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[15]][switches]

            if(length(switches)>0){
              X_fixed[non_copy,] <- population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[28]]
            }
          }
        }
      }

      # Check if some copy is genotyped
      if(nrow(loop_elements_list[[3]])>0){
        for(index in 1:nrow(loop_elements_list[[3]])){
          kindex <- loop_elements_list[[3]][index,2]
          k.database <- bve.database[loop_elements_list[[3]][index,3],]
          if(length(stay.loop.elements)>0){
            non_copy <- which(stay.loop.elements==loop_elements_copy[index,6])
          } else{
            non_copy <- loop_elements_copy[index,6]
          }
          genotyped[non_copy] <- max(genotyped[non_copy],population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[16]])
          genotyped_array[non_copy, population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[22]]] <- TRUE
        }
      }
      genotype.included <- which(genotyped==1)

      # in case of copied individuals use the copy that has been phenotyped the most


      skip_check <- FALSE

      if(input.phenotype!="own"){

        y_keep <- matrix(FALSE, nrow=nrow(y), ncol=ncol(y))

        if(length(offspring.bve.parents.database)==0){
          offspring.bve.parents.database <- loop_elements[,c(4,5,2,2,2), drop=FALSE]
          skip_check <- TRUE
          offspring.bve.parents.database2 <- loop_elements_copy[,c(4,5,2,2,6), drop=FALSE]

          is_copy <- c(rep(FALSE, nrow(offspring.bve.parents.database)), rep(TRUE, nrow(offspring.bve.parents.database2)))

          offspring.bve.parents.database <- rbind(offspring.bve.parents.database, offspring.bve.parents.database2)
        } else{
          is_copy <- rep(FALSE, sum(offspring.bve.offspring.database[,4]- offspring.bve.offspring.database[,3]+1))
        }



        for(index2 in 1:nrow(offspring.bve.parents.database)){
          activ.database <- offspring.bve.parents.database[index2,]

          if(skip_check){
            loop <- index2
            if(is_copy[index2]){
              loop <- offspring.bve.parents.database[index2,5]
            }
          } else{
            loop <- 1:nrow(loop_elements)
          }

          for(index in loop){
            if(is_copy[index2]){
              kindex <- loop_elements[offspring.bve.parents.database[index2,5],2]
              kindex2 <- offspring.bve.parents.database[index2,4]
            } else{
              kindex <- kindex2 <-  loop_elements[index,2]
            }

            k.database <- bve.database[loop_elements[index,3],]
            activ.indi <- population$breeding[[k.database[[1]]]][[k.database[[2]]]][[kindex]][[21]]

            if(skip_check){
              import <- TRUE
            } else{
              import <- activ.indi[,1]==activ.database[1] & activ.indi[,2]==activ.database[2] & activ.indi[,3]>=activ.database[3] & activ.indi[,3] <= activ.database[4]
            }

            if(sum(import)>0){
              own_pheno <- y[index,]
              n_obs <- y_obs[index,]
              off_pheno <- population$breeding[[activ.database[1]]][[activ.database[2]+26]][,kindex2]
              n_off <- population$breeding[[activ.database[1]]][[activ.database[2]+28]][,kindex2]
              off_pheno[n_off==0] <- NA ##
              if(input.phenotype=="off"){

                replace <- ceiling(n_off / 2) >=  y_obs[index,]
                y[index,replace] <- off_pheno[replace]

                y_obs[index,replace] <- ceiling(n_off / 2)[replace]
                y_keep[index,replace] <- TRUE

              } else if(input.phenotype=="mean"){

                replace <- ceiling(n_obs + n_off/2) >= y_obs[index,]

                take_both <- !is.na(own_pheno) & !is.na(off_pheno) & replace
                take_own <- !is.na(own_pheno) & is.na(off_pheno) & replace
                take_off <- is.na(own_pheno) & !is.na(off_pheno) & replace


                y[index,take_both] <- ((own_pheno + off_pheno)/((own_pheno!=0) + (off_pheno!=0)))[take_both]
                y[index,take_own] <- (own_pheno /own_pheno!=0)[take_own]
                y[index,take_off] <- (off_pheno/off_pheno!=0)[take_off]
                y_obs[index,] <- ceiling(n_obs + n_off/2)

              } else if(input.phenotype=="weighted"){

                replace <- ceiling(n_obs + n_off/2) >= y_obs[index,]

                take_both <- !is.na(own_pheno) & !is.na(off_pheno) & replace
                take_own <- !is.na(own_pheno) & is.na(off_pheno) & replace
                take_off <- is.na(own_pheno) & !is.na(off_pheno) & replace
                y[index,take_both] <- ((own_pheno*n_obs*2 + off_pheno*n_off)/(n_obs*2 + n_off))[take_both]
                y[index,take_own] <- (own_pheno /own_pheno!=0)[take_own]
                y[index,take_off] <- (off_pheno/off_pheno!=0)[take_off]
                y_obs[index,] <- ceiling(n_obs + n_off / 2)
              }
            }
          }
        }


        if(input.phenotype=="off"){
          y[!y_keep] <- NA
          y_obs[!y_keep] <- 0
        }


      }


      # Import Z
      if((relationship.matrix != "pedigree" && relationship.matrix != "kinship" && length(mas.geno) == 0)){
        if(miraculix){

          if(store.comp.times.bve){
            before <- as.numeric(Sys.time())
          }

          if (requireNamespace("miraculix", quietly = TRUE)) {
            Z.code <- miraculix::computeSNPS(population, loop_elements[,4], loop_elements[,5], loop_elements[,2], what="geno", output_compressed=TRUE)
            if(relationship.matrix!="vanRaden"){
              Zt <- miraculix::computeSNPS(population, loop_elements[,4], loop_elements[,5], loop_elements[,2], what="geno", output_compressed=FALSE)
            }
          }
          if(store.comp.times.bve){
            after <- as.numeric(Sys.time())
            zcalc <- zcalc + after - before
          }

        } else{
          if(store.comp.times.bve){
            before <- as.numeric(Sys.time())
          }

          for(index in 1:n.animals){
            k.database <- bve.database[loop_elements[index,3],]
            cindex <- index
            kindex <- loop_elements[index,2]
            Zt[,cindex] <- base::as.integer(colSums(compute.snps(population, k.database[1],k.database[2],kindex, import.position.calculation=import.position.calculation,decodeOriginsU=decodeOriginsU, bit.storing=bit.storing, nbits=nbits, output_compressed=FALSE)))
          }
          if(store.comp.times.bve){
            after <- as.numeric(Sys.time())
            zcalc <- zcalc + after - before
          }
        }

      }

      to_remove <- NULL
      genotyped_array_variants <- unique(genotyped_array)
      used_arrays <- colSums(genotyped_array)
      if(sum(used_arrays == nrow(genotyped_array)  & !population$info$array.is_subset)>0){
        to_remove <- NULL
      } else{

        if(bve.imputation){
          keeps <- rep(FALSE, sum(population$info$snp))
        } else{
          keeps <- rep(TRUE, sum(population$info$snp))
        }

        for(index in 1:nrow(genotyped_array_variants)){
          activ_keep <- rep(FALSE, sum(population$info$snp))
          for(index2 in which(genotyped_array_variants[index,])){
            activ_keep <- activ_keep | population$info$array.markers[[index2]]
          }

          if(bve.imputation){
            keeps <- keeps | activ_keep
          } else{
            keeps <- keeps & activ_keep
          }

        }
        to_remove <- c(to_remove, which(!keeps))

      }

      if(remove.effect.position){
        to_remove <- c(to_remove, population$info$effect.p)
      }
      if(length(bve.array)){
        to_remove <- c(to_remove, which(!population$info$array.markers[[bve.array]]))
      }

      if(length(to_remove)>0){
        to_remove <- unique(to_remove)
      }

      if(length(to_remove)==sum(population$info$snp)){
        to_remove <- NULL
      }


      if(length(to_remove)>0 && (relationship.matrix != "pedigree" && relationship.matrix != "kinship")){

        if(miraculix && exists("Z.code")){
          if (requireNamespace("miraculix", quietly = TRUE)) {
            Z.code <- as.matrix(Z.code)
            Z.code <- miraculix::genomicmatrix(Z.code[-to_remove,])
            #Z.code <- miraculix::zeroNthGeno(Z.code, to_remove) ### currently not implemented in miraculix
          }
        } else{
          Zt <- Zt[-to_remove,]
        }


        if(verbose) cat(paste0(sum(population$info$snp) - length(to_remove), " markers survived filtering for BVE.\n"))


      }

      if(bve.imputation.errorrate>0 && (relationship.matrix != "pedigree" && relationship.matrix != "kinship")){

        if(length(to_remove)>0){
          is_imputed <- get.genotyped.snp(population, database = cbind(loop_elements[,4], loop_elements[,5], loop_elements[,2]))[-to_remove,]==0
        } else{
          is_imputed <- get.genotyped.snp(population, database = cbind(loop_elements[,4], loop_elements[,5], loop_elements[,2]))==0
        }


        if(miraculix && exists("Z.code")){
          if (requireNamespace("miraculix", quietly = TRUE)) {
            temp1 <- as.matrix(Z.code)[is_imputed]
          }
        } else{
          temp1 <- Zt[is_imputed]
        }


        temp2 <- temp1
        zeros <- which(temp1==0)
        ones <- which(temp1==1)
        twos <- which(temp1==2)

        if(length(zeros)>0){
          temp1[zeros][stats::rbinom(length(zeros),1,bve.imputation.errorrate^2)==1] <- 2
          temp1[zeros][stats::rbinom(length(zeros),1,2*bve.imputation.errorrate * (1-bve.imputation.errorrate))==1 & temp1[zeros]==0] <- 1
        }


        if(length(ones)>0){
          temp1[ones][stats::rbinom(length(ones),1,bve.imputation.errorrate*(1-bve.imputation.errorrate))==1] <- 0
          temp1[ones][stats::rbinom(length(ones),1,bve.imputation.errorrate*(1-bve.imputation.errorrate))==1 & temp1[ones]==1] <- 2
        }


        if(length(twos)>0){
          temp1[twos][stats::rbinom(length(twos),1,bve.imputation.errorrate^2)==1] <- 0
          temp1[twos][stats::rbinom(length(twos),1,2*bve.imputation.errorrate * (1-bve.imputation.errorrate))==1 & temp1[twos]==2] <- 1
        }


        if(miraculix && exists("Z.code")){
          if (requireNamespace("miraculix", quietly = TRUE)) {

            Z.code <- as.matrix(Z.code)
            Z.code[is_imputed] <- temp1
            Z.code <- miraculix::genomicmatrix(Z.code)
          }
        } else{
          Zt[is_imputed] <- temp1
        }

        if(verbose) cat(paste0("A total of ", length(temp1), " entries in the genotype matrix were imputed.\n"))
        if(verbose) cat(paste0("Accuracy of the imputation: ", mean(temp1==temp2), ".\n"))

      }

      if(store.comp.times.bve){
        comp.times.bve[2] <- as.numeric(Sys.time())
      }

      if(mixblup.bve){
        options("scipen"=999)
      }
      if(!mas.bve){


        # Derive relationship matrix and single-step only for vanRaden based genomic relationship
        if(relationship.matrix=="kinship" || relationship.matrix == "pedigree"){
          z_ped <- z_ped - as.numeric(Sys.time())

          if(!mixblup.bve){

            A <- kinship.exp(population, database=bve.database, depth.pedigree=depth.pedigree, elements = loop_elements_list[[2]], mult=2, verbose=verbose,
                             storage.save = storage.save)

            if(store.sparse){
              A <- methods::as(A, "sparseMatrix")
            }

          } else if(mixblup.pedfile){
            write.pedigree.mixblup(population, path = mixblup.path.pedfile, database = bve.database,
                                   depth.pedigree=depth.pedigree, storage.save = storage.save, verbose=verbose)

          }


          z_ped <- z_ped + as.numeric(Sys.time())
        } else if(relationship.matrix=="vanRaden"){


          if(singlestep.active && (BGLR.model == "RKHS" || !BGLR.bve)){
            n.geno <- length(genotype.included)
            n.ped <- n.animals - n.geno
            if(n.ped==0){
              if(verbose) cat("Use Genomic BLUP (All individuals are genotyped)!\n")
              singlestep.active <- FALSE
            } else if(n.geno==0){
              z_ped <- z_ped - as.numeric(Sys.time())
              if(verbose) cat("Use Pedigree BLUP (No individuals genotyped!)\n")
              singlestep.active <- FALSE
              relationship.matrix <- "kinship"
              if(!mixblup.bve){
                A <- kinship.exp(population, database=bve.database, depth.pedigree=depth.pedigree, elements = loop_elements_list[[2]], mult = 2, verbose=verbose)
                if(store.sparse){
                  A <- methods::as(A, "sparseMatrix")
                }
              }
              z_ped <- z_ped + as.numeric(Sys.time())
            } else{
              if(verbose) cat("Use Single Step GBLUP\n")
            }
          }

          if(mixblup.bve && (TRUE || singlestep.active) && (BGLR.model == "RKHS" || !BGLR.bve)){

            z_ped <- z_ped - as.numeric(Sys.time())
            write.pedigree.mixblup(population, path = mixblup.path.pedfile, database = bve.database, depth.pedigree=depth.pedigree, storage.save = storage.save)
            z_ped <- z_ped + as.numeric(Sys.time())
          }


            if(BGLR.model == "RKHS" || !BGLR.bve){
              if(singlestep.active){

                if(!mixblup.bve){

                  if(verbose) cat("Start derive Single-step relationship matrix \n")
                  if(verbose) cat(paste0("Construct pedigree matrix for ", length(loop_elements_list[[2]]), " individuals.\n"))
                  z_ped <- z_ped - as.numeric(Sys.time())


                  A_pedigree <-  kinship.exp(population, database=bve.database, depth.pedigree=depth.pedigree, elements = loop_elements_list[[2]], mult = 2, verbose=verbose)

                  if(store.sparse){
                    A_pedigree <- methods::as(A_pedigree, "sparseMatrix")
                  }

                  z_ped <- z_ped + as.numeric(Sys.time())
                  if(verbose) cat(paste0("Derived pedigree matrix in  ", round(z_ped, digits=2), " seconds.\n"))
                  if(verbose) cat("Start deriving of H matrix for", length(genotype.included), "genotyped and", nrow(A_pedigree)-length(genotype.included), "non-genotyped individuals.\n")


                }




                if(miraculix){
                  if (requireNamespace("miraculix", quietly = TRUE)) {
                    # Avoid having the full genotype matrix in memory without bitwise storage!
                    if(bve.imputation.errorrate>0){
                      Z.code.small <- miraculix::genomicmatrix(as.matrix(Z.code)[,genotype.included])
                      if(length(to_remove)>0){
                        Z.code.small <- miraculix::genomicmatrix(as.matrix(Z.code.small)[-to_remove,])
                      }
                    } else{
                      Z.code.small <- miraculix::computeSNPS(population, loop_elements[genotype.included,4], loop_elements[genotype.included,5], loop_elements[genotype.included,2], what="geno", output_compressed=TRUE)
                      if(length(to_remove)>0){
                        Z.code.small <- miraculix::genomicmatrix(as.matrix(Z.code.small)[-to_remove,])
                      }
                    }

                    if(mixblup.bve){


                      geno = as.matrix(Z.code.small)
                      dense <- cbind(get.id(population, database = bve.database)[loop_elements_list[[2]]][genotype.included], 0)
                      for(index in 1:nrow(dense)){
                        dense[index,2] <- paste0(geno[,index], collapse = "")
                      }


                    } else{
                      #p_i <- miraculix::allele_freq(Z.code.small)
                      p_i <- rowMeans(as.matrix(Z.code.small))/2
                      A_geno <- miraculix::relationshipMatrix(Z.code.small, centered=TRUE, normalized=TRUE)
                    }

                    rm(Z.code.small)
                  }
                } else if(miraculix.mult){
                  if (requireNamespace("miraculix", quietly = TRUE)) {
                    p_i <- rowSums(Zt[,genotype.included])/ncol(Zt[,genotype.included])/2
                    Zt_miraculix <- miraculix::genomicmatrix(Zt[,genotype.included])
                    if(mixblup.bve){
                      dense <- cbind(get.id(population, database = bve.database)[loop_elements_list[[2]]][genotype.included], 0)
                      for(index in 1:nrow(dense)){
                        dense[index,2] <- paste0(Zt_miraculix[,index], collapse = "")
                      }
                    } else{
                      A_geno <- miraculix::relationshipMatrix(Zt_miraculix, centered=TRUE, normalized=TRUE)
                    }

                    rm(Zt_miraculix)
                  }
                } else{
                  if(mixblup.bve){
                    geno = Zt[,genotype.included]
                    dense <- cbind(get.id(population, database = bve.database)[loop_elements_list[[2]]][genotype.included], 0)
                    for(index in 1:nrow(dense)){
                      dense[index,2] <- paste0(geno[,index], collapse = "")
                    }
                    rm(geno)
                  } else{
                    p_i <- rowSums(Zt[,genotype.included])/ncol(Zt[,genotype.included])/2
                    Ztm <- Zt[,genotype.included] - p_i * 2
                    A_geno <- crossprod(Ztm)/ (2 * sum(p_i*(1-p_i)))
                  }

                }

                if(!mixblup.bve){
                  z_h <- z_h - as.numeric(Sys.time())

                  reorder <- c((1:n.animals)[-genotype.included], genotype.included)

                  test1 <- as.numeric(A_geno)
                  test2 <- as.numeric(A_pedigree[genotype.included, genotype.included])
                  a_step <- mean(test2) - mean(test1)
                  A_geno <- A_geno * (1-a_step/2) + a_step # Modification according to Vitezica 2011

                  A <- ssGBLUP(A11 = A_pedigree[-genotype.included, -genotype.included],
                               A12 = A_pedigree[-genotype.included, genotype.included],
                               A22 = A_pedigree[genotype.included, genotype.included], G = A_geno)

                  #rm(A_geno)
                  #rm(A_pedigree)

                  rest <- (1:n.animals)[-genotype.included]
                  A[c(genotype.included, rest), c(genotype.included, rest)] <- A[c((ncol(A)-length(genotype.included)+1):ncol(A),1:(ncol(A)-length(genotype.included)) ), c((ncol(A)-length(genotype.included)+1):ncol(A),1:(ncol(A)-length(genotype.included)))]
                  z_h <- z_h + as.numeric(Sys.time())
                  if(verbose) cat(paste0("Derived H matrix in  ", round(z_h, digits=2), " seconds.\n"))

                }

              } else if(relationship.matrix=="vanRaden"){
                if(miraculix){
                  if (requireNamespace("miraculix", quietly = TRUE)) {

                    if(mixblup.bve){

                      geno = as.matrix(Z.code)
                      dense <- cbind(get.id(population, database = bve.database), 0)
                      for(index in 1:nrow(dense)){
                        dense[index,2] <- paste0(geno[,index], collapse = "")
                      }
                      rm(geno)

                    } else{
                      #p_i <- miraculix::allele_freq(Z.code) # Noch nicht implementiert?
                      p_i <- rowMeans(as.matrix(Z.code))/2 # Noch nicht implementiert?
                      A <- miraculix::relationshipMatrix(Z.code, centered=TRUE, normalized=TRUE)
                    }

                  }
                } else if(miraculix.mult){
                  if (requireNamespace("miraculix", quietly = TRUE)) {
                    if(mixblup.bve){
                      dense <- cbind(get.id(population, database = bve.database), 0)
                      for(index in 1:nrow(dense)){
                        dense[index,2] <- paste0(Zt[,index], collapse = "")
                      }
                    } else{
                      p_i <- rowSums(Zt)/ncol(Zt)/2
                      Zt_miraculix <- miraculix::genomicmatrix(Zt)
                      A <- miraculix::relationshipMatrix(Zt_miraculix, centered=TRUE, normalized=TRUE)
                    }

                  }
                } else{

                  if(mixblup.bve){
                    dense <- cbind(get.id(population, database = bve.database), 0)
                    for(index in 1:nrow(dense)){
                      dense[index,2] <- paste0(Zt[,index], collapse = "")
                    }
                  } else{
                    p_i <- rowSums(Zt)/ncol(Zt)/2
                    Ztm <- Zt - p_i * 2
                    A <- crossprod(Ztm)/ (2 * sum(p_i*(1-p_i)))
                  }

                }
              }
            }



        } else if(relationship.matrix=="CM"){
          #CM SCHAETZER
          Ztm <- rbind(Zt==0, Zt==1, Zt==2)
          A <- crossprod(Ztm) / nrow(Zt)
        } else if(relationship.matrix=="CE"){

          Ztm <- rbind(Zt==0, Zt==1, Zt==2)
          A <- crossprod(Ztm)
          A <- (A^2 - 0.5*A)/(nrow(Zt)^2)
        } else if(relationship.matrix=="CE2"){
          A2 <- crossprod(Zt)
          A <- 0.5 * A2 * A2 + 0.5 * crossprod(Zt*Zt)
          A <- A/mean(diag(A))
        } else if(relationship.matrix=="CE3"){
          A2 <- crossprod(Zt)
          A <- 0.5 * A2 * A2 - 0.5 * crossprod(Zt*Zt)
          A <- A/mean(diag(A))
        } else if(relationship.matrix=="non_stand"){
          A <- crossprod(Zt) / nrow(Zt)
        } else if(relationship.matrix=="vanRaden2"){
          p_i <- rowSums(Zt)/ncol(Zt)/2
          Ztm <- Zt - p_i * 2
          A <- crossprod(Ztm)/ (2 * sum(p_i*(1-p_i)))
        }
      }




      sigma.a2.hat <- numeric(length(sigma.g))
      sigma.e2.hat <- sigma.e2

      prev_rest_take <- NULL

      n.rep <- 0
      if(length(bve.database)==length(bve.insert.database) && prod(bve.database==bve.insert.database)==1 && nrow(loop_elements_list[[3]])==0){
        bve.insert <- rep(TRUE, n.animals)
        bve.insert.copy <- NULL
      } else{
        bve.insert <- rep(FALSE, n.animals)
        before <- 0
        for(index in 1:nrow(bve.insert.database)){
          add_insert <- loop_elements[,2] <= bve.insert.database[index,4] & loop_elements[,2] >= bve.insert.database[index,3] & loop_elements[,4] == bve.insert.database[index,1] & loop_elements[,5] == bve.insert.database[index,2]
          bve.insert[add_insert] <- TRUE
        }

        loop_elements_copy <- loop_elements_list[[3]]
        n.rep <- nrow(loop_elements_copy)
        bve.insert.copy <- rep(FALSE, n.rep)
        if(n.rep>0){
          for(index in 1:nrow(bve.insert.database)){
            add_insert <- loop_elements_copy[,2] <= bve.insert.database[index,4] & loop_elements_copy[,2] >= bve.insert.database[index,3] & loop_elements_copy[,4] == bve.insert.database[index,1] & loop_elements_copy[,5] == bve.insert.database[index,2]
            bve.insert.copy[add_insert] <- TRUE
          }
        }

      }
      # Breeding value estimation (Solving of the mixed-model / bayesian generalized linear regression) - all single trait models except multi-sommer

      if(length(bve.ignore.traits)>0){
        text <- population$info$trait.name[bve.ignore.traits][1]
        for(write in population$info$trait.name[bve.ignore.traits][-1]){
          text <- paste0(text, ", ", write)
        }
        if(verbose) cat(paste0("BVE for: ", text, " has been skipped.\n"))
      }


      combi_trait <- NULL
      skip_trait <- NULL
      replace_trait <- NULL

      if(length(population$info$trait.combine.name)>0){

        activ_bves <- (1:population$info$bv.nr)[bve.keeps]

        for(check_combi in 1:length(population$info$trait.combine.included)){

          to_combine <- intersect(population$info$trait.combine.included[[check_combi]], activ_bves)

          if(length(to_combine)>1){
            if(verbose){
              cat(paste0("Perform a joint BVE for trait ", population$info$trait.combine.name[[check_combi]], ".\n"))

              text <- "This includes traits: "
              for(combi2 in to_combine){
                text <- paste0(text, population$info$trait.name[combi2], " ")
              }
            }

            y_temp <- y[,to_combine]
            ny_temp <- y_obs[,to_combine]
            y_temp[is.na(y_temp)] <- 0

            new_y <- rowSums(y_temp *ny_temp)
            new_ny <- rowSums(ny_temp)
            y[,to_combine] <- new_y/new_ny
            y_obs[,to_combine] <- new_ny

            combi_trait <- c(combi_trait, to_combine[1])
            skip_trait <- c(skip_trait, to_combine[-1])
            replace_trait <- c(replace_trait, rep(to_combine[1], length(to_combine)-1))
          }
        }

      }

      beta_hat <-  colMeans(y, na.rm=TRUE) # Rest faellt weg! (X' R^-1 X)^-1 X' R^-1 y
      beta_hat_check <- colMeans(y_real, na.rm=TRUE)
      beta_hat_var <- diag(stats::var(y, use = "pairwise.complete.obs"))

      if(sum(is.na(beta_hat))>0){
        skip_check <- rep(FALSE, population$info$bv.nr)
        skip_check[skip_trait] <- TRUE
        left <- population$info$trait.name[which(is.na(beta_hat) & skip_check)]
        if(length(left)>0){
          if(verbose) cat("No phenotypes available for traits:", left,"\n")
          if(verbose) cat("Set all breeding value estimates for these trait(s) to 0. \n")
          beta_hat[is.na(beta_hat)] <- 0
        }

      }

      if((rrblup.bve | sommer.bve | BGLR.bve)){

        if(length(bve.exclude.fixed.effects)>0){
          X_fixed_temp <- X_fixed[,-bve.exclude.fixed.effects, drop=TRUE]
        } else{
          X_fixed_temp <- X_fixed
        }
        X_fixed1 <- cbind(1, X_fixed_temp)

        q <- qr(X_fixed1)
        X_fixed1 <- X_fixed1[,q$pivot[seq(q$rank)], drop=FALSE]

      }


      if(export.relationship.matrix){

        sexa <- loop_elements[,5]
        sexa[sexa==1] <- "M"
        sexa[sexa==2] <- "F"
        names_temp <- paste0(sexa,loop_elements[,2], "_", loop_elements[,4])

        colnames(A) <- rownames(A) <- names_temp
        return(A)
      }

      if(length(import.relationship.matrix) > 0){
        A = import.relationship.matrix
      }

      ## MiXBLUP setup
      {
        if(mixblup.bve && mixblup.parfile){

          if(verbose) cat(paste0("Start writting parfile at ", mixblup.path.parfile,"\n"))

          if(length(mixblup.genetic.cov)>0){
            gen_cor = mixblup.genetic.cov
          } else{

            if(variance.correction=="parental.mean"){
              gen_cor = stats::cov(y_real2-y_p1/2 - y_p2/2, use = "pairwise.complete.obs")[bve.keeps,bve.keeps]
            } else if(variance.correction =="generation.mean"){

              gen_cor = stats::cov(y_real2-t(generation_mean[,loop_elements[,4], drop=FALSE]), use = "pairwise.complete.obs")[bve.keeps,bve.keeps]

            } else {
              gen_cor = stats::cov(y_real2, use = "pairwise.complete.obs")[bve.keeps,bve.keeps]
            }



            if(length(gen_cor)==1){
              gen_cor = as.matrix(gen_cor)
            }

            gen_cor = round(gen_cor, digits = 6)
            # PLACE HOLDER TO MAKE SURE MIXBLUP GETS A POSITIVE DEFINIT MATRIX
            diag(gen_cor) = diag(gen_cor) + 10^(-5)
          }

          scal1 = diag(gen_cor)

          gen_cor = diag(sqrt(1/scal1), nrow=length(scal1)) %*% gen_cor %*% diag(sqrt(1/scal1), nrow=length(scal1))
          gen_cor = round(gen_cor, digits = 6)


          k = 6
          while(eigen(gen_cor)$values[nrow(gen_cor)]<= (10^(-5))){
            diag(gen_cor) = diag(gen_cor) + 10^(-k)
            gen_cor = gen_cor / gen_cor[1,1]
            gen_cor = round(gen_cor, digits = 6)
            if(k<=3){if(verbose) cat("Diagonal of Genetic covariance matrix was stretched to ensure definitness (1)!")}
            k = k-1
          }

          gen_cor = gen_cor * diag(sqrt(scal1), nrow = length(scal1)) %*% gen_cor %*% diag(sqrt(scal1), nrow = length(scal1))

          k = 6
          while(eigen(gen_cor)$values[nrow(gen_cor)]<= (10^(-5))){
            diag(gen_cor) = diag(gen_cor) + 10^(-k)
            if(k<=3){if(verbose) cat("Diagonal of Genetic covariance matrix was stretched to ensure definitness (2)!")}
            k = k-1
          }

          rownames(gen_cor) <- paste0("phen", 1:length(bve.keeps),"(animal)")
          utils::write.table(file=mixblup.path.parfile, "G", row.names = FALSE, col.names = FALSE , quote=FALSE)
          for(index in 1:nrow(gen_cor)){
            utils::write.table(file=mixblup.path.parfile, gen_cor[index,1:index, drop=FALSE], row.names = TRUE, col.names = FALSE , quote=FALSE, append=TRUE)
          }
          utils::write.table(file=mixblup.path.parfile, "", row.names = FALSE, col.names = FALSE , quote=FALSE, append=TRUE)

          if(length(mixblup.residual.cov)>0){
            res_cor = mixblup.residual.cov
          } else{
            res_cor <- stats::cov(y-y_real2, use = "pairwise.complete.obs")[bve.keeps,bve.keeps]
            if(length(res_cor)==1){
              res_cor = as.matrix(res_cor)
            }

            if(sum(is.na(res_cor))>0){
              if(verbose){
                cat("Residual correlation for MiXBLUP highly questionable!")
                print(res_cor)
                res_cor[is.na(res_cor)] = 0
              }
            }

            # PLACE HOLDER TO MAKE SURE MIXBLUP GETS A POSITIVE DEFINIT MATRIX
            diag(res_cor) = diag(res_cor) + 10^(-5)
          }



          scal2 = diag(res_cor)

          res_cor = diag(sqrt(1/scal2), nrow = length(scal2)) %*% res_cor %*% diag(sqrt(1/scal2), nrow = length(scal2))

          res_cor = round(res_cor, digits = 6)

          k = 6
          while(eigen(res_cor)$values[nrow(res_cor)]<=(10^(-5))){
            diag(res_cor) = diag(res_cor) + 10^(-k)
            res_cor = res_cor / res_cor[1,1]
            res_cor = round(res_cor, digits = 6)
            if(k<=3){if(verbose) cat("Diagonal of Residual covariance matrix was stretched to ensure definitness (1)!")}
            k = k-1
          }

          res_cor = diag(sqrt(scal2), nrow = length(scal2)) %*% res_cor %*% diag(sqrt(scal2), nrow = length(scal2))

          k = 6
          while(eigen(res_cor)$values[nrow(res_cor)]<=(10^(-5))){
            diag(res_cor) = diag(res_cor) + 10^(-k)
            if(k<=3){if(verbose) cat("Diagonal of Residual covariance matrix was stretched to ensure definitness (2)!")}
            k = k-1
          }
          rownames(res_cor) <- paste0("phen", 1:length(bve.keeps),"(animal)")
          utils::write.table(file=mixblup.path.parfile, "Res", row.names = FALSE, col.names = FALSE , quote=FALSE, append = TRUE)
          for(index in 1:nrow(res_cor)){
            utils::write.table(file=mixblup.path.parfile, res_cor[index,1:index, drop=FALSE], row.names = TRUE, col.names = FALSE , quote=FALSE, append=TRUE)
          }

          for(index in 1:nrow(res_cor)){
            if(verbose){cat(paste0("Variance components in BVE: sigma_g^2 = ", round(gen_cor[index,index], digits=4), "; sigma_e^2 = ", round(res_cor[index,index], digits=4), " ; h^2 = ", round(gen_cor[index,index]/ (gen_cor[index,index] + res_cor[index,index]), digits=3), "(Trait: ", population$info$trait.name[bve.keeps][index] ,")\n"))}
          }



        }

        if(mixblup.bve && mixblup.datafile){
          if(verbose) cat(paste0("Start writting datafile at ", mixblup.path.datafile,"\n"))
          y_temp = y
          y_temp[is.na(y_temp)] = -99

          pheno_table = cbind(get.id(population, database = bve.database)[loop_elements_list[[2]]], 1, y_temp)
          utils::write.table(file=mixblup.path.datafile, pheno_table, row.names = FALSE, col.names = FALSE , quote=FALSE)
        }

        if(mixblup.bve && mixblup.genofile && length(dense)>0){



          if (requireNamespace("data.table", quietly = TRUE)) {
            data.table::fwrite(file=mixblup.path.genofile, dense, sep = " ", col.names = FALSE)
          } else{
            utils::write.table(file=mixblup.path.genofile, dense, col.names = FALSE, row.names = FALSE, quote = FALSE)
          }

        }

        if(mixblup.bve && mixblup.inputfile){
          if(verbose) cat(paste0("Start writting input file at ", mixblup.path.inputfile,"\n"))

          utils::write.table(file=mixblup.path.inputfile, paste0("TITLE MoBPS v", utils::sessionInfo()$otherPkgs$MoBPS$Version, " at " ,Sys.time()),  row.names = FALSE, col.names = FALSE, quote=FALSE)
          utils::write.table(file=mixblup.path.inputfile, "", append = TRUE, row.names = FALSE, col.names = FALSE , quote=FALSE)
# , c("!MISSING -99", rep("",2+length(bve.keeps)))
          datafile_temp <- cbind(c("DATAFILE", "animal", "mean", paste0("phen", 1:length(bve.keeps))), c(mixblup.path.datafile, "I" ,"I", rep("T", length(bve.keeps))))
          datafile_temp <- cbind(datafile_temp, "")
          datafile_temp[1,3] = "!MISSING -99"
          utils::write.table(file=mixblup.path.inputfile, datafile_temp, append = TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)
          utils::write.table(file=mixblup.path.inputfile, "", append = TRUE, row.names = FALSE, col.names = FALSE , quote=FALSE)

          if(length(dense)>0){

            if(mixblup.apy){

              if(length(mixblup.apy.core)==0  || nrow(dense) < mixblup.apy.core){
                mixblup.apy.core = nrow(dense)
              }

              if(length(mixblup.apy.core)==1 && nrow(dense) <= mixblup.apy.core){
                mixblup.apy.core = nrow(dense)
                if(verbose){ cat("Number of genotyped individuals is smaller than the core!\nNo use of APY")}

                utils::write.table(file = mixblup.path.inputfile, cbind("ERMFILE", mixblup.path.genofile, "!CONSTRUCT", "ssmat", "!SINGLESTEP", "!DENSE", 2, "!MAF", 0.01),
                                   append = TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)
              } else{
                utils::write.table(file = mixblup.path.inputfile, cbind("ERMFILE", mixblup.path.genofile, "!CONSTRUCT", "ssmat", "!SINGLESTEP", "!DENSE", 2, "!MAF", 0.01, "!APY", "!APYCoreRan", mixblup.apy.core),
                                   append = TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)

              }





            } else{
              utils::write.table(file = mixblup.path.inputfile, cbind("ERMFILE", mixblup.path.genofile, "!CONSTRUCT", "ssmat", "!SINGLESTEP", "!DENSE", 2, "!MAF", 0.01),
                                 append = TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)
            }

            utils::write.table(file = mixblup.path.inputfile, cbind("animal", "I"),
                        append = TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)
            utils::write.table(file = mixblup.path.inputfile, cbind("!METHOD", "VanRaden"),
                        append = TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)
            utils::write.table(file = mixblup.path.inputfile, cbind("!Lambda", mixblup.lambda),
                        append = TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)
            utils::write.table(file = mixblup.path.inputfile, cbind("!ALPHA", mixblup.alpha),
                        append = TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)
            utils::write.table(file = mixblup.path.inputfile, cbind("!Beta", mixblup.beta),
                               append = TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)
            utils::write.table(file = mixblup.path.inputfile, cbind("!OMEGA", mixblup.omega),
                               append = TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)
            utils::write.table(file=mixblup.path.inputfile, "", append = TRUE, row.names = FALSE, col.names = FALSE , quote=FALSE)
          }

          if(TRUE || length(dense)==0 || singlestep.active==TRUE){
            utils::write.table(file=mixblup.path.inputfile, cbind("PEDFILE", mixblup.path.pedfile,  "!CalcInbr",  "!Groups", "1.0"), append = TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)
            utils::write.table(file=mixblup.path.inputfile, cbind(c("animal", "sire", "dam"), c("I", "I", "I")), append = TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)
            utils::write.table(file=mixblup.path.inputfile, "", append = TRUE, row.names = FALSE, col.names = FALSE , quote=FALSE)
          }

          utils::write.table(file=mixblup.path.inputfile, cbind("PARFILE", mixblup.path.parfile), append = TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE)
          utils::write.table(file=mixblup.path.inputfile, "", append = TRUE, row.names = FALSE, col.names = FALSE , quote=FALSE)

          utils::write.table(file=mixblup.path.inputfile, "MODEL", append = TRUE, row.names = FALSE, col.names = FALSE, quote=FALSE )
          for(index in 1:length(bve.keeps)){
            utils::write.table(file=mixblup.path.inputfile, paste0("phen", index, " ~ mean !RANDOM G(animal)"), append = TRUE, row.names = FALSE, col.names = FALSE , quote=FALSE)
          }

          utils::write.table(file=mixblup.path.inputfile, "", append = TRUE, row.names = FALSE, col.names = FALSE , quote=FALSE)

          utils::write.table(file=mixblup.path.inputfile, "SOLVING", append = TRUE, row.names = FALSE, col.names = FALSE , quote=FALSE)
          utils::write.table(file=mixblup.path.inputfile, "!MAXIT 50", append = TRUE, row.names = FALSE, col.names = FALSE , quote=FALSE)

        }
      }

      if(store.comp.times.bve){
        comp.times.bve[3] <- as.numeric(Sys.time())
      }

      for(bven in (1:population$info$bv.nr)[bve.keeps]){

        if(sum(bven==skip_trait)==0){

          if(forecast.sigma.g){

            if(variance.correction=="parental.mean"){
              sigma.g[bven] <- sqrt(stats::var(y_real2[,bven]-y_p1[,bven]/2 - y_p2[,bven]/2, na.rm = TRUE))

            } else if(variance.correction =="generation.mean"){
              sigma.g[bven] <- sqrt(stats::var(y_real2[,bven]-generation_mean[loop_elements[,4]], na.rm = TRUE))
            } else {
              sigma.g[bven] <- sqrt(stats::var(y_real2[,bven], na.rm = TRUE))
            }



          }
          if(estimate.pheno.var){
            sigma.e2.hat[bven] <- max(0, stats::var(y[,bven], na.rm=TRUE) - sigma.g[bven]^2)
          }
          sigma.a2.hat[bven] <- sigma.g[bven]^2
          if(estimate.add.gen.var){
            sigma.a2.hat[bven] <- max(min(stats::lm(y[,bven]~y_parent[,bven])$coefficients[2],1),0.001) * sigma.e2.hat[bven]
          }
          if(mas.bve){
            if(length(mas.markers)==0){
              if(length(population$info$real.bv.add[[bven]])>24){
                if(verbose) cat("No markers for MAS provided. Use 5 true effect markers with highest effect.\n")
                temp1 <- sort(rowSums(abs(population$info$real.bv.add[[bven]][,3:5])), index.return=TRUE, decreasing=TRUE)$ix[1:mas.number]
                active.snp <- population$info$real.bv.add[[bven]][temp1,]
                mas.markers.temp <- c(0,population$info$cumsnp)[active.snp[,2]]+ active.snp[,1]
              } else{
                mas.markers.temp <- sample(population$info$effect.p, min(mas.number, length(population$info$effect.p)))
              }

            } else{
              mas.markers.temp <- mas.markers
            }

            if(length(mas.geno)>0){
              mas_geno <- mas.geno[mas.markers.temp,]
            } else if(miraculix){
              mas_geno <- as.matrix(Z.code)[mas.markers.temp,]
            } else{
              mas_geno <- Zt[mas.markers.temp,]
            }

            if(length(mas.effects)==0){
              model <- stats::lm(y~t(mas_geno))
              y_hat[,bven] <- model$fitted.values
            } else{
              mas.effects <- rep(mas.effects, length.out=length(mas.markers.temp))
              y_hat[,bven] <- mas.effects %*% mas_geno
            }

          } else if(BGLR.bve){

            if(miraculix){
              Zt <- as.matrix(Z.code)
            }
            Zt <- t(scale(t(Zt), center=TRUE, scale=FALSE))
            fixed <- which(is.na(Zt[,1]))
            if(length(fixed)>0){
              Zt <- Zt[-fixed,]
            }
            if(BGLR.model=="RKHS"){
              ETA <- list(list(X = X_fixed1, model="FIXED"), list(K=A, model='RKHS'))
            } else if(BGLR.model=="BayesA"){
              ETA <- list(list(X = X_fixed1, model="FIXED"), list(X=t(Zt), model='BayesA'))
            } else if(BGLR.model=="BayesB"){
              ETA <- list(list(X = X_fixed1, model="FIXED"), list(X=t(Zt), model='BayesB'))
            } else if(BGLR.model=="BayesC"){
              ETA <- list(list(X = X_fixed1, model="FIXED"), list(X=t(Zt), model='BayesC'))
            } else if(BGLR.model=="BL"){
              ETA <- list(list(X = X_fixed1, model="FIXED"), list(X=t(Zt), model='BL'))
            } else if(BGLR.model=="BRR"){
              ETA <- list(list(X = X_fixed1, model="FIXED"), list(X=t(Zt), model='BRR'))
            }

            if(min(y[!is.na(y[,bven]),bven])==max(y[!is.na(y[,bven]),bven])){
              y_hat[,bven] <- rep(max(y[!is.na(y[,bven]),bven]), length(y[,bven]))
            } else{
              if(BGLR.save.random){
                BGLR.save <- paste0(BGLR.save, sample(1:10000,1))
              }
              if(requireNamespace("BGLR", quietly = TRUE)){
                fm <- BGLR::BGLR(y=y[,bven], ETA=ETA, nIter=BGLR.iteration, burnIn=BGLR.burnin, saveAt=BGLR.save, verbose=BGLR.print)
              } else{
                stop("Use of BGLR without being installed!")
              }

              if(ncol(X_fixed1)==1){
                temp1 <- 0
              } else{
                temp1 <- X_fixed1[,-1, drop=FALSE] %*% fm$ETA[[1]]$b[-1]
              }

              if(estimate.u){
                u_hat_possible <- TRUE
                if(BGLR.model=="RKHS") stop("RKHS does not provide marker estimates in BGLR. Please select BGLR.model = BayesA or similar")
                u_hat <- cbind(u_hat, fm$ETA[[2]]$b)
              }

              y_hat[,bven] <- fm$yHat - temp1
            }

          } else if(rrblup.bve){

            check <- sum(is.na(y[,bven]))
            if(verbose) cat(paste0(length(y[,bven]) - check, " phenotyped individuals in BVE (Trait: ", population$info$trait.name[bven],").\n"))
            if(verbose) cat(paste0(length(y[,bven]), " individuals considered in BVE.\n"))
            if(check >= (length(y[,bven])-1)){
              if(verbose) cat(paste0("No phenotyped individuals for trait ", population$info$trait.name[bven], "\n"))

              if(verbose) cat(paste0("Skip this BVE.\n"))
              next
            }

            t1 <- as.numeric(Sys.time())

            if(estimate.u){
              u_hat_possible <- TRUE
              if(miraculix){
                Zt <- as.matrix(Z.code)
              }
              test <- rrBLUP::mixed.solve(y[,bven], Z = t(Zt), X = X_fixed1, method="REML", bounds = c(1e-9,1e9))
              y_hat[,bven] <- as.numeric(test$beta[1]) + as.numeric(test$u %*% Zt)
              t2 <- as.numeric(Sys.time())

              u_hat <- cbind(u_hat, test$u)
            } else{
              if (requireNamespace("rrBLUP", quietly = TRUE)) {
                test <- rrBLUP::mixed.solve(y[,bven], K = A, X = X_fixed1, method="REML", bounds = c(1e-9,1e9))
              } else{
                stop("Use of rrBLUP without being installed!")
              }

              t2 <- as.numeric(Sys.time())


              y_hat[,bven] <- as.numeric(test$beta[1]) + test$u

              if(verbose){
                cat(paste0("Variance components in BVE: sigma_g^2 = ", round(test$Vu, digits=4), "; sigma_e^2 = ", round(test$Ve, digits=4), " ; h^2 = ", round(test$Vu / (test$Vu + test$Ve), digits=3), "\n"))
              }

            }
            if(verbose) cat(paste0(round(t2-t1, digits=2), " seconds for BVE.\n"))


          } else if(emmreml.bve){
            check <- sum(is.na(y[,bven]))
            if(verbose) cat(paste0(length(y[,bven]) - check, " phenotyped individuals in BVE (Trait: ", population$info$trait.name[bven],").\n"))
            if(verbose) cat(paste0(length(y[,bven]), " individuals considered in BVE.\n"))
            if(check >= (length(y[,bven])-1)){
              if(verbose) cat(paste0("No phenotyped individuals for trait ", population$info$trait.name[bven], "\n"))
              if(verbose) cat(paste0("Skip this BVE."))
              next
            }
            if(check>0){
              if(verbose) cat(paste0("Breeding value estimation with ", check, " NA phenotypes! EMMREML does not support this!\n"))
              if(verbose) cat(paste0("No estimation is performed to NA individuals. \n"))
              take <- which(!is.na(y[,bven]))
            } else{
              take <- 1:length(y[,bven])
            }
            n <- length(y[take,bven])
            p <- nrow(Zt)
            stopifnot(n == ncol(Zt[,take]))

            if(requireNamespace("EMMREML", quietly = TRUE)){
              fm <- EMMREML::emmreml(
                y[take,bven],
                matrix(1,nrow=n),
                diag(n),
                A[take,take])
            } else{
              stop("Usage of EMMREML without being installed!")
            }



            y_hat[take,bven] <- as.numeric(fm$uhat) + as.numeric(fm$betahat)
            if(estimate.u){
              u_hat_possible <- TRUE
              if(miraculix){
                Zt <- as.matrix(Z.code)
              }
              u_hat <- cbind(u_hat, alpha_to_beta(drop(fm$uhat),A[take,take],t(Zt[,take])), deparse.level = 0)
            }

          } else if(sommer.bve){

            check <- sum(is.na(y[,bven]))
            if(verbose) cat(paste0(length(y[,bven]) - check, " phenotyped individuals in BVE (Trait: ", population$info$trait.name[bven],").\n"))
            if(verbose) cat(paste0(length(y[,bven]), " individuals considered in BVE.\n"))
            if(check >= (length(y[,bven])-1)){
              if(verbose) cat(paste0("No phenotyped individuals for trait ", population$info$trait.name[bven], "\n"))
              if(verbose) cat(paste0("Skip this BVE.\n"))
              next
            }

            traitnames <- (paste0("name", 1:ncol(y)))
            traitnames[bven] <- "name"
            traitnames <- as.factor(traitnames)
            colnames(y) <- traitnames
            id <- as.factor(paste0("P", 1:nrow(y)))
            y_som <- data.frame(y, id)
            rownames(y_som) <- id
            colnames(A) <- rownames(A) <- id

            if (requireNamespace("sommer", quietly = TRUE)) {
              test <- sommer::mmer(name ~ X_fixed1, random=~sommer::vs(id, Gu=A), rcov = ~units, data=y_som)
            } else{
              stop("Use of sommer without being installed!")
            }

            cat(paste0("Variance components in BVE: sigma_g^2 = ", round(test$sigma$`u:id`, digits=4), "; sigma_e^2 = ", round(test$sigma$units, digits=4), " ; h^2 = ", round(test$sigma$`u:id` / (test$sigma$`u:id` + test$sigma$units), digits=3), "\n"))


            y_hat[sort(as.character(id), index.return=TRUE)$ix,bven] <- test$U[[1]][[1]] + as.numeric((test$Beta[[3]][1]))

          } else if(sommer.multi.bve){

            check <- sum(is.na(y))

            if(verbose) cat(paste0(length(y[,bven]) - check, " phenotyped individuals in BVE (Trait: ", population$info$trait.name[bven],").\n"))
            if(verbose) cat(paste0(length(y[,bven]), " individuals considered in BVE.\n"))
            if(check >= (length(y[,bven])-1)){
              if(verbose) cat(paste0("No phenotyped individuals for multi-trait mixed model\n"))
              if(verbose) cat(paste0("Skip this BVE.\n"))
              next
            }


            if(bven==max((1:population$info$bv.nr)[bve.keeps])){

              traitnames <- paste0("name", 1:ncol(y))
              colnames(y) <- as.factor(traitnames)
              id <- as.factor(paste0("P", 1:nrow(y)))
              y_som <- data.frame(y, id)
              rownames(y_som) <- id
              colnames(A) <- rownames(A) <- id
              text <- "cbind("
              for(index in 1:length(traitnames)){
                if(index==length(traitnames)){
                  text <- paste0(text, traitnames[index], ")")
                } else{
                  text <- paste0(text, traitnames[index], ",")
                }

              }
              text <- paste0("sommer::mmer(",text,"~1, random=~sommer::vs(id, Gu=A, Gtc=sommer::unsm(bven)), rcov = ~sommer::vs(units, Gtc=diag(bven)), data=y_som)")
              if (requireNamespace("sommer", quietly = TRUE)) {
                test <- eval(parse(text=text))
              } else{
                stop("Use of sommer without being installed!")
              }

              for(bven1 in 1:population$info$bv.nr){
                y_hat[sort(as.character(id), index.return=TRUE)$ix,bven1] <- test$U[[1]][[bven1]]
              }


              if(estimate.u){
                if(verbose) cat("U estimation not available in sommer")
              }
            }
          } else if(mixblup.bve){


            if(store.comp.times.bve){
              before <- as.numeric(Sys.time())
            }

            check <- sum(is.na(y[,bven]))

            if(verbose) cat(paste0(length(y[,bven]) - check, " phenotyped individuals in BVE (Trait: ", population$info$trait.name[bven],").\n"))
            if(verbose) cat(paste0(length(y[,bven]), " individuals considered in BVE.\n"))
            if(FALSE && check >= (length(y[,bven])-1)){
              if(verbose) cat(paste0("No phenotyped individuals for multi-trait mixed model\n"))
              if(verbose) cat(paste0("Skip this BVE.\n"))
              next
            }


            if(bven==max((1:population$info$bv.nr)[bve.keeps])){

              if(verbose) cat(paste0("Start MiXBLUP BVE\n"))

              if(sum(dir()=="Solani.txt")>0){
                cat("Remove previous Solani.txt file!\n")
                file.remove("Solani.txt")
              }
              tick_mix = Sys.time()
              system(paste0(mixblup.path, " ", mixblup.path.inputfile), intern=verbose&&mixblup.verbose,
                     ignore.stdout=!(verbose&&mixblup.verbose), ignore.stderr=!(verbose&&mixblup.verbose))

              tack_mix = Sys.time()
              if(verbose) cat(paste0("Finished MiXBLUP BVE\nEstimate took ", round(tack_mix-tick_mix, digits = 2), " seconds.\n"))

              estimated_bve <- as.matrix(utils::read.table("Solani.txt"))
              fixed_eff <- as.matrix(utils::read.table("Solfix.txt", header=TRUE)[,c(2,5)])

              if(nrow(fixed_eff) < length(bve.keeps)){
                fixed_eff = rbind(fixed_eff, cbind((1:length(bve.keeps))[-fixed_eff[,1]], 0))
                fixed_eff = fixed_eff[sort(fixed_eff[,1], index.return=TRUE)$ix,]
              }

              fixed_eff = fixed_eff[,2]

              for(index in 1:length(bve.keeps)){
                estimated_bve[,index+3] <- estimated_bve[,index+3] + fixed_eff[index]
              }

              sorting <- sort(estimated_bve[,1], index.return=TRUE)
              estimated_bve <- estimated_bve[sorting$ix,,drop=FALSE]

              #### THIS IS EXTREMLY INEFFICIENT AND NEEDS REWORK!
              matching = get.id(population, database = bve.database)[loop_elements_list[[2]]]

              if(min(matching)==matching[1] && prod(matching==estimated_bve[1:length(matching) + which(matching[1]==estimated_bve[,1])-1,1])){
                takes = 1:length(matching) + which(matching[1]==estimated_bve[,1])-1
              } else{
                takes = numeric(length(matching))
                for(index5 in 1:length(matching)){
                  if(index5 > 1 && matching[index5]==estimated_bve[takes[index5-1]+1,1]){
                    takes[index5] = takes[index5-1]+1
                  } else{
                    takes[index5] = which(estimated_bve[,1]==matching[index5])

                  }
                }
              }

              y_hat = estimated_bve[takes,4:ncol(estimated_bve),drop=FALSE]

              if(store.comp.times.bve){
                after <- as.numeric(Sys.time())
                z_chol <- after - before
              }

            }
          } else if(sigma.e2[bven]>0 || sum(is.na(y[,bven]>0))){
            # sigma.a2.hat / sigma.e2.hat are modified to avoid h=0, h=1 cases (numeric instability)


            if(input.phenotype!="own"){
              warning("Direct BVE method is not really designed to be used when using offspring phentypes!\n Check variance components / accuracies! \n Consider using rrblup.bve=TRUE")
            }

            if(population$info$phenotypic.transform[bven]){
              warning("Direct BVE method is not really designed when using a phenotypic transformation!\n Check variance components / accuracies! \n Consider using rrblup.bve=TRUE")
            }
            if(ncol(X_fixed)>=1 && prod(X_fixed==X_fixed[1,1])!=1){
              warning("Direct BVE method does not support the use of fixed effects!\nAssume all fixed effects to be known!")

              beta_fixed_store <- as.numeric(X_fixed %*% t(population$info$fixed.effects[bven,,drop=FALSE]) )

              X_fixed1 <- matrix(1, nrow=nrow(X_fixed), ncol=1)
            } else{
              beta_fixed_store <- rep(0, nrow(X_fixed))
              X_fixed1 <- X_fixed
            }

            if(population$info$phenotypic.transform[bven] && sum(abs(beta_fixed_store))>0){
              stop("Direct BVE method with fixed effects does not allow for a phenotypic transformation!")
              beta_fixed_store <- rep(0, nrow(X_fixed))
            }

            if(!bve.beta.hat.approx){
              if(population$info$phenotypic.transform[bven]==FALSE || min(y_obs[,bven])>0){
                beta_hat[bven] <- beta_hat_check[bven]
              } else{
                pop1 <- breeding.diploid(population, phenotyping.database = loop_elements[,c(4,5,2)], verbose=TRUE)
                y_temp <- get.pheno(pop1, database = loop_elements[,c(4,5,2)])[bven,]
                beta_hat_check <- mean(y_temp)
                beta_hat[bven] <- beta_hat_check[bven]
                rm(pop1)
              }

            }

            if(population$info$phenotypic.transform[bven]==FALSE & !is.na(beta_hat_var)[bven] & (abs(beta_hat - beta_hat_check) / sqrt(beta_hat_var))[bven]>3){
              warning("Direct BVE method does not estimate intercept in the fixed effect, but just uses the mean phenotype!\n This estimate here seems to be not good. Consider using a different BVE method or setting bve.beta.hat.approx = FALSE to use the true underlying beta value")
            }



            bve.direct.est.now <- mobps.bve
            check <- sum(is.na(y[,bven]))

            u_hat_possible <- TRUE
            # This will be more complicated in case other fixed-effects are added to the model !
            if(ncol(X_fixed1)>0){
              multi <- y[,bven] - X_fixed1 %*% population$info$fixed.effect[bven,] - beta_hat[bven]
            } else{
              multi <- y[,bven] - beta_hat[bven]
            }

            if(!population$info$phenotypic.transform[bven]){
              multi <- multi - beta_fixed_store
            }


            rrblup.required <- FALSE

            if(verbose) cat(paste0(length(y[,bven]) - check, " phenotyped individuals in BVE (Trait: ", population$info$trait.name[bven],").\n"))
            if(verbose) cat(paste0(length(y[,bven]), " individuals considered in BVE.\n"))
            if(check >= (length(y[,bven])-1)){
              if(verbose) cat(paste0("No phenotyped individuals for trait ", population$info$trait.name[bven], "\n"))
              if(verbose) cat(paste0("Skip this BVE.\n"))
              next
            }

            # take Individuals used to trait the mixed model on
            # take2 individuals for which to enter a estimated breeding value
            # take3 individuals to estimate via rrblup
            if(check>0){
              if(sum(genotyped==1 & is.na(y[,bven]))==sum(genotyped==1)){
                if(bve.direct.est.now==FALSE){
                  if(verbose) cat("No genotyped and phenotyped individuals. Application of rrBLUP not possible!\n")
                  if(verbose) cat("Assume non phenotyped individuals to have average phenotype.\n")
                  take <- take2 <- 1:length(y[,bven])
                } else{
                  take <- which(!is.na(y[,bven]))
                  take2 <- 1:length(y[,bven])
                }


              } else if(sum(genotyped==1 & is.na(y[,bven]))>0){
                if(verbose) cat("Some genotyped individuals without phenotype.\n")
                take <- which(!is.na(y[,bven])) # individuals for GBLUP
                if(!bve.direct.est.now){
                  rrblup.required <- TRUE
                  take2 <- which(!is.na(y[,bven]) & genotyped==1) # individuals for rrBLUP
                } else if(FALSE && estimate.u){
                  take2 <- which(!is.na(y[,bven]) & genotyped==1) # individuals for rrBLUP
                } else{
                  take2 <- 1:length(y[,bven])
                }


                bve.insert.full <- bve.insert
                if(length(bve.insert.copy)>0){

                  bve.insert.full[loop_elements_copy[bve.insert.copy,6]] <- TRUE
                }
                take3 <- which(is.na(y[,bven]) & genotyped==1 & bve.insert.full) # individuals to estimate via rrBLUP
                if(length(take3)>0 && !bve.direct.est.now){
                  if(verbose) cat("Use rrBLUP to estimate breeding value for those individuals!\n")
                }
              } else{
                take <- take2 <- which(!is.na(y[,bven]))
                if(bve.direct.est.now){
                  take2 <- 1:length(y[,bven])
                }
              }
            } else{
              take <- take2 <- 1:length(y[,bven])
            }

            if(length(take)==length(take2)){
              bve.direct.est.now <- FALSE
            }


            if(length(take2) > length(take)){
              if(forecast.sigma.g){
                sigma.g[bven] <- sqrt(stats::var(y_real2[take,bven], na.rm = TRUE))
              }
              if(estimate.pheno.var){
                sigma.e2.hat[bven] <- max(0, stats::var(y[take,bven], na.rm=TRUE) - sigma.g[bven]^2)
              }
              sigma.a2.hat[bven] <- sigma.g[bven]^2
              if(estimate.add.gen.var){
                sigma.a2.hat[bven] <- max(min(stats::lm(y[take,bven]~y_parent[take,bven])$coefficients[2],1),0.001) * sigma.e2.hat[bven]
              }
            }

            if(sigma.a2.hat[bven]==0){
              if(sigma.e2.hat[bven]>0){
                sigma.a2.hat[bven] <- sigma.e2.hat[bven] * 0.001
              } else{
                sigma.e2.hat[bven] <- 1
                sigma.a2.hat[bven] <- 1000
              }
            }
            if(sigma.e2.hat[bven]==0){
              sigma.e2.hat[bven] <- sigma.a2.hat[bven] * 0.001
            }


            if(verbose){
              cat(paste0("Variance components in BVE: sigma_g^2 = ", round(sigma.a2.hat[bven], digits=4), "; sigma_e^2 = ", round(sigma.e2.hat[bven], digits=4), " ; h^2 = ", round(sigma.a2.hat[bven] / (sigma.e2.hat[bven] + sigma.a2.hat[bven]), digits=3), "\n"))
            }

            if(store.comp.times.bve){
              before <- as.numeric(Sys.time())
            }
            # Solve mixed model which assumed known heritability
            # This is not 100% accurate but massively reduces computing time and should be ok for large scale breeding programs
            t1 <- as.numeric(Sys.time())

            if(length(take)==nrow(A)){
              skip.copy = TRUE
            } else{
              skip.copy = FALSE
            }

            if(bve.per.sample.sigma.e & max(y_obs, na.rm = TRUE)>1){
              if(length(repeatability) < bven || length(heritability) <bven || repeatability[bven]==heritability[bven]){
                heri_factor <- 1 / y_obs[take,bven]
              } else{
                heri_factor <- (sigma.e2.perm[bven] + sigma.e2.rest[bven] / y_obs[take,bven]) / (sigma.e2.perm[bven] + sigma.e2.rest[bven])
                if(sum(heri_factor< 0.01)>0 || sum(heri_factor>1)){
                  warning("Breeding value estimation with multiple observations with extrem residual variance scaling!")
                }
              }


            } else{
              heri_factor <- 1
            }

            if(prod(heri_factor==1)!=1 || is.function(bve.solve)){
              miraculix_possible <- FALSE
            } else{
              miraculix_possible <- TRUE
            }



            if(calculate.reliability){
              if(skip.copy){
                GR1 <- chol2inv(chol(add.diag(A, sigma.e2.hat[bven] / sigma.a2.hat[bven] * heri_factor)))
              } else{
                GR1 <- chol2inv(chol(add.diag(A[take,take], sigma.e2.hat[bven] / sigma.a2.hat[bven] * heri_factor)))
              }


              Rest_term <- GR1 %*% multi[take]

              if(bve.direct.est.now){
                y_hat[take2,bven] <- A[take2,take] %*% (Rest_term)  + beta_hat[bven]
                y_reli[take2,bven] <- diag( A[take2,take] %*% GR1 %*% A[take,take2])
              } else{
                y_hat[take,bven] <- A[take,take] %*% (Rest_term)  + beta_hat[bven]
                y_reli[take,bven] <- diag( A[take,take] %*% GR1 %*% A[take,take])
              }
            } else if(miraculix && miraculix.chol && miraculix_possible){
              if (requireNamespace("miraculix", quietly = TRUE)) {

                 if(bve.direct.est.now){
                  temp1 = miraculix::solveRelMat(A[take,take], sigma.e2.hat[bven] / sigma.a2.hat[bven], multi[take],betahat = beta_hat[bven], destroy_A = miraculix.destroyA)
                  y_hat[take2,bven] <-  A[take2,take] %*% temp1[[1]] + beta_hat[bven]
                  Rest_term = temp1[[1]]
                } else{
                  if(skip.copy){
                    temp1 = miraculix::solveRelMat(A, sigma.e2.hat[bven] / sigma.a2.hat[bven], multi[take],beta_hat[bven], destroy_A = miraculix.destroyA)
                    y_hat[,bven] <- temp1[[2]]
                    Rest_term = temp1[[1]]
                  } else{

                    temp1 = miraculix::solveRelMat(A[take,take], sigma.e2.hat[bven] / sigma.a2.hat[bven], multi[take],beta_hat[bven], destroy_A = miraculix.destroyA)
                    y_hat[take,bven] <- temp1[[2]]
                    Rest_term = temp1[[1]]
                  }
                }

              }
            } else{
              if(bve.direct.est.now){
                if(is.function(bve.solve)){
                  Rest_term = bve.solve(add.diag(A[take,take], sigma.e2.hat[bven] / sigma.a2.hat[bven] * heri_factor), b= multi[take])
                } else{
                  Rest_term = (chol2inv(chol(add.diag(A[take,take], sigma.e2.hat[bven] / sigma.a2.hat[bven] * heri_factor))) %*% multi[take])
                }
                y_hat[take2,bven] <- A[take2,take] %*% Rest_term + beta_hat[bven]
              } else{
                if(skip.copy){
                  if(is.function(bve.solve)){
                    Rest_term = bve.solve(A=add.diag(A,sigma.e2.hat[bven] / sigma.a2.hat[bven] * heri_factor), b=multi[take])
                  } else{
                    Rest_term = (chol2inv(chol(add.diag(A,sigma.e2.hat[bven] / sigma.a2.hat[bven] * heri_factor))) %*% multi[take])
                  }
                  y_hat[,bven] <- A %*%  Rest_term + beta_hat[bven]
                } else{
                  if(is.function(bve.solve)){
                    Rest_term = bve.solve(A=add.diag(A[take,take],sigma.e2.hat[bven] / sigma.a2.hat[bven] * heri_factor), b=multi[take])
                  } else{
                    Rest_term = (chol2inv(chol(add.diag(A[take,take],sigma.e2.hat[bven] / sigma.a2.hat[bven] * heri_factor))) %*% multi[take])
                  }
                  y_hat[take,bven] <- A[take,take] %*%  Rest_term + beta_hat[bven]
                }
              }
            }


            y_hat[take2, bven] <- y_hat[take2,bven] # + beta_fixed_store

            t2 <- as.numeric(Sys.time())
            if(verbose) cat(paste0(round(t2-t1, digits=2), " seconds for BVE.\n"))

            if(store.comp.times.bve){
              after <- as.numeric(Sys.time())
              z_chol <- after - before
            }

            # Estimation of marker effects via rrBLUP

            if(estimate.u || rrblup.required){

              nsnp_bve <- sum(population$info$snp) - length(to_remove)


              while(bven>1 && (length(u_hat)==0 || ncol(u_hat)<(bven-1))){
                u_hat <- cbind(u_hat,rep(0,nsnp_bve), deparse.level = 0)
              }

              rest_take <- which(duplicated(c(take,take2))[-(1:length(take))])

              if(store.comp.times.bve){
                before <- as.numeric(Sys.time())
              }
              rest_take <- which(duplicated(c(take,take2))[-(1:length(take))])
              if(!(length(rest_take)==length(prev_rest_take) && prod(rest_take==prev_rest_take)==1)){
                if(miraculix && length(take2[rest_take])!=nrow(loop_elements)){
                  if (requireNamespace("miraculix", quietly = TRUE)) {
                    Z.code2 <- miraculix::computeSNPS(population, loop_elements[take2[rest_take],4], loop_elements[take2[rest_take],5], loop_elements[take2[rest_take],2], what="geno",
                                                      output_compressed = TRUE)

                    if(length(to_remove)>0){
                      Z.code2 <- as.matrix(Z.code2)
                      Z.code2 <- miraculix::genomicmatrix(Z.code2[-to_remove,])
                    }

                  }
                } else if(length(take2)!=nrow(loop_elements)){
                  Zt2 <- Zt[,take2[rest_take]]
                } else if(miraculix){
                  Z.code2 <- Z.code
                } else {
                  Zt2 <- Zt
                }
                if(!fast.uhat){
                  if (requireNamespace("MASS", quietly = TRUE)) {
                    A1 <- MASS::ginv(A[take2[rest_take], take2[rest_take]])
                  } else{
                    stop("Use of MASS without being installed!")
                  }
                }

                prev_rest_take <- rest_take
              } else if(length(prev_rest_take)==0){
                if(miraculix){
                  Z.code2 <- Z.code
                } else{
                  Zt2 <- Zt
                }

              }

              if(relationship.matrix!="vanRaden"){
                if(miraculix){
                  if (requireNamespace("miraculix", quietly = TRUE)) {
                    #p_i <- miraculix::allele_freq(Z.code2)
                    p_i <- rowMeans(as.matrix(Z.code2))/2
                  }
                } else{
                  p_i <- rowSums(Zt[,take2[rest_take]])/2
                }
              }

              if(miraculix){
                if (requireNamespace("miraculix", quietly = TRUE)) {
                  if(fast.uhat){
                    u_hat <- cbind(u_hat, 1/ 2 / sum(p_i*(1-p_i))* ((as.matrix(Z.code2) -  2* p_i) %*%  Rest_term), deparse.level = 0)
                  } else{
                    u_hat <- cbind(u_hat, 1/ 2 / sum(p_i*(1-p_i)) * ((as.matrix(Z.code2) -  2* p_i) %*%  A1 %*% (y_hat[take2[rest_take],bven] - beta_hat[bven])), deparse.level = 0)
                  }
                }

              } else {
                if(fast.uhat){
                  u_hat <- cbind(u_hat, 1/ 2 / sum(p_i*(1-p_i))*(Ztm[,take2[rest_take]] %*% Rest_term), deparse.level = 0)

                } else{
                  u_hat <- cbind(u_hat, 1/ 2 / sum(p_i*(1-p_i))*(Ztm[,take2[rest_take]] %*% (A1 %*% (y_hat[take2[rest_take],bven] - beta_hat[bven]))), deparse.level = 0)
                }
              }


              if(store.comp.times.bve){
                after <- as.numeric(Sys.time())
                z_uhat <- z_uhat + after - before
              }

              if(rrblup.required){
                if(miraculix){
                  if (requireNamespace("miraculix", quietly = TRUE)) {
                    y_hat[take3,bven] <- u_hat[,bven] %*% (as.matrix(Z.code)[,take3]-2*p_i) + beta_hat[bven]
                  }
                } else {
                  y_hat[take3,bven] <- u_hat[,bven] %*% Ztm[,take3] + beta_hat[bven]
                }
              }
            }
          } else{
            y_hat[,bven] <- y[,bven]
          }

        } else{

          which_skip <- which(skip_trait==bven)

          y_hat[,bven] <- y_hat[, replace_trait[which_skip]]
        }


        if(store.comp.times.bve){
          comp.times.bve[4] <- as.numeric(Sys.time())
        }




        if(estimate.reliability){
          suppressWarnings(temp1 <- stats::cor(y_hat[bve.insert,bven], y_real[bve.insert, bven])^2)
          if(is.na(temp1)){
            if(stats::var(y_hat[bve.insert,bven])==0){
              warning("Reliability estimation = 0 because all individuals have the same breeding value")
            }
            if(stats::var(y_real[bve.insert,bven])==0){
              warning("Reliability estimation = 0 because all individuals have the same genomic value")
            }
            temp1 <- 0
          }
          y_reli[,bven] <- temp1

        }

        if(n.rep>0){
          for(index in (1:nrow(loop_elements_copy))[bve.insert.copy]){
            if(length(stay.loop.elements)>0){
              non_copy <- which(stay.loop.elements==loop_elements_copy[index,6])
            } else{
              non_copy <- loop_elements_copy[index,6]
            }
            if(length(non_copy)==1){
              population$breeding[[loop_elements_copy[index,4]]][[loop_elements_copy[index,5]+2]][bve.keeps, loop_elements_copy[index,2]] <- y_hat[non_copy,bve.keeps]

            }
          }
        }



        #  GWAS CODE NOT WRITEN FOR PARALLEL COMPUTING
        if(gwas.u){

          if(y.gwas.used=="pheno"){
            y_gwas <- y
          } else if(y.gwas.used=="bv"){
            y_gwas <- y_real
          } else if(y.gwas.used=="bve"){
            y_gwas <- y_hat
          }

          if(nrow(gwas.database)!=nrow(bve.database) || prod(gwas.database==bve.database)==0){

            loop_elements_gwas_list <- derive.loop.elements(population=population, bve.database=bve.database,
                                                            bve.class=bve.class, bve.avoid.duplicates=bve.avoid.duplicates,
                                                            store.adding=TRUE)
            loop_elements_gwas <- loop_elements_gwas_list[[1]]
            n.animals.gwas <- nrow(loop_elements_gwas)
          }
          if(gwas.group.standard){
            gwas_start <- loop_elements_gwas_list[[2]]
          }

          if(nrow(gwas.database)!=nrow(bve.database) || prod(gwas.database==bve.database)==0){
            Zt <- array(0L,dim=c(sum(population$info$snp), n.animals.gwas))
            y_gwas <- array(0, dim=c(n.animals.gwas, population$info$bv.nr))
            cindex <- 1
            for(index in 1:n.animals.gwas){
              k.database <- gwas.database[loop_elements_gwas[index,3],]
              if(miraculix){
                if (requireNamespace("miraculix", quietly = TRUE)) {
                  Zt[,cindex] <- miraculix::computeSNPS(population, k.database[1],k.database[2],kindex, what="geno", output_compressed=FALSE)
                }
              } else{
                Zt[,cindex] <- base::as.integer(colSums(compute.snps(population, k.database[1],k.database[2],kindex, import.position.calculation=import.position.calculation, decodeOriginsU=decodeOriginsU, bit.storing=bit.storing, nbits=nbits, output_compressed=FALSE)))
              }
              for(bven in 1:population$info$bv.nr){
                # HIER EVENTUELL SNPs auslesen!
                if(y.gwas.used=="pheno"){
                  y_gwas[cindex,bven] <- population$breeding[[k.database[1]]][[8+k.database[2]]][bven,kindex]
                } else if(y.gwas.used=="bv"){
                  y_gwas[cindex,bven] <- population$breeding[[k.database[1]]][[6+k.database[2]]][bven,kindex]
                } else if(y.gwas.used=="bve"){
                  y_gwas[cindex,bven] <- population$breeding[[k.database[1]]][[2+k.database[2]]][bven,kindex]
                }
              }
              cindex <- cindex +1
            }
          }

          if(gwas.group.standard){

            for(indexg in 1:(length(gwas_start)-1)){
              if(gwas_start[indexg]>=(gwas_start[indexg+1]-1)){
                y_gwas[gwas_start[indexg]:(gwas_start[indexg+1]-1), bven] <- y_gwas[gwas_start[indexg]:(gwas_start[indexg+1]-1), bven] - mean(y_gwas[gwas_start[indexg]:(gwas_start[indexg+1]-1), bven])
              }
            }
          }
          x_mean <- rowMeans(Zt)
          x2_mean <- rowMeans(Zt^2)
          #          xy_mean <- colMeans(Z*y_gwas[,bven])
          xy_mean <- colMeans(t(Zt)*y_gwas[,bven])


          n <- length(y_gwas[,bven])
          y_mean <- mean(y_gwas[,bven])

          b1 <- (n *xy_mean - n *x_mean * y_mean) / (x2_mean*n - n *x_mean^2)


          if(approx.residuals==FALSE){
            sigma1 <- 1/(n * (x2_mean-(x_mean)^2))
            b0 <- y_mean - b1 * x_mean
            var1 <- numeric(length(sigma1))
            for(index in 1:length(sigma1)){
              var1[index] <- sigma1[index] * stats::var(y_gwas[,bven] - b1[index] * Zt[index,] - b0[index]) * (n-1)/(n-2)
            }

          } else{
            var1 <- 1/(n * (x2_mean-(x_mean)^2)) * stats::var(y_gwas[,bven])
          }
          test <- b1/sqrt(var1)
          gwas_hat <- cbind(gwas_hat, test, deparse.level = 0)
          #sorted <- sort(abs(test), index.return=TRUE)
        }

      }

      for(index in (1:nrow(loop_elements))[bve.insert]){
        population$breeding[[loop_elements[index,4]]][[loop_elements[index,5]+2]][bve.keeps, loop_elements[index,2]] <- y_hat[index,bve.keeps]
      }
      if(calculate.reliability || estimate.reliability){
        for(index in (1:nrow(loop_elements))[bve.insert]){
          population$breeding[[loop_elements[index,4]]][[loop_elements[index,5]+18]][bve.keeps, loop_elements[index,2]] <- y_reli[index,bve.keeps]
        }
      }

      if(report.accuracy){
        if(verbose) cat("Correlation between genetic values and BVE:\n")
        if(n.rep==0){
          y_hat_temp <- y_hat
          y_hat_temp[y_hat_temp==0] <- NA
          if(length(bve.ignore.traits)>0){
            acc <- suppressWarnings(stats::cor(y_real2[bve.insert,-bve.ignore.traits], y_hat_temp[bve.insert,-bve.ignore.traits], use="pairwise.complete.obs"))
          } else{
            acc <- suppressWarnings(stats::cor(y_real2[bve.insert,], y_hat_temp[bve.insert,], use="pairwise.complete.obs"))
          }

        } else{
          insert.temp <- numeric(length(bve.insert.copy))

          if(length(stay.loop.elements)>0){
            for(index in (1:nrow(loop_elements_copy))[bve.insert.copy]){
              inserter <- which(stay.loop.elements==loop_elements_copy[index,6])
              insert.temp[index] <- if(length(inserter)==1){ inserter} else{NA}
            }
          } else{
            for(index in (1:nrow(loop_elements_copy))[bve.insert.copy]){
              insert.temp[index] <- loop_elements_copy[index,6]
            }

          }
          y_hat_temp <- rbind(y_hat[bve.insert,,drop=FALSE], y_hat[insert.temp,,drop=FALSE])
          y_hat_temp[y_hat_temp==0] <- NA

          if(length(bve.ignore.traits)>0){
            acc <- suppressWarnings(stats::cor(rbind(y_real2[bve.insert,-bve.ignore.traits,drop=FALSE], y_real2[insert.temp,-bve.ignore.traits, drop=FALSE]),
                                               y_hat_temp[,-bve.ignore.traits,drop=FALSE], use="pairwise.complete.obs"))
          } else{
            acc <- suppressWarnings(stats::cor(rbind(y_real2[bve.insert,,drop=FALSE], y_real2[insert.temp,, drop=FALSE]),
                                               y_hat_temp[,,drop=FALSE], use="pairwise.complete.obs"))
          }

        }
        if(length(acc)==1){
          acc <- matrix(acc,nrow=1)
        }

        if(sum(is.na(acc))>0){
          acc[is.na(acc)] <- 0
        }
        if(verbose) cat(diag(acc))
        if(verbose) cat("\n")
      }

    }
    if(u_hat_possible && bve && estimate.u && relationship.matrix=="vanRaden"){

      if(length(to_remove)==0){
        rownames(u_hat) = population$info$snp.name
      } else{
        rownames(u_hat) = population$info$snp.name[-to_remove]
      }


      population$info$u_hat[[length(population$info$u_hat)+1]] <- u_hat
      population$info$u_hat_single[[length(population$info$u_hat)]] <- list()
      for(bven in 1:ncol(u_hat)){
        population$info$u_hat_single[[length(population$info$u_hat)]][[bven]] <- cbind((-2*p_i) *u_hat[,bven],(-2*p_i+1) *u_hat[,bven],(-2*p_i+2) *u_hat[,bven], deparse.level = 0)
      }
    } else if(u_hat_possible && bve && estimate.u && relationship.matrix=="CM"){
      population$info$u_hat[[length(population$info$u_hat)+1]] <- u_hat
      population$info$u_hat_single[[length(population$info$u_hat)]] <- list()
      for(bven in 1:ncol(u_hat)){
        population$info$u_hat_single[[length(population$info$u_hat)]][[bven]] <- cbind(u_hat[1:nrow(Zt),bven],u_hat[1:nrow(Zt)+ nrow(Zt),bven],u_hat[1:nrow(Zt)+2*nrow(Zt),bven], deparse.level = 0)
      }
    }
    if(gwas.u){
      if(sum(is.na(gwas_hat)>0)){
        gwas_hat[is.na(gwas_hat)] <- 0
      }
      population$info$gwas_hat[[length(population$info$gwas_hat)+1]] <- gwas_hat
    }

    if(store.comp.times.bve){
      comp.times.bve[5] <- as.numeric(Sys.time())
    }

    if(store.comp.times){
      comp.times[5] <- as.numeric(Sys.time())
    }
  }

  #######################################################################
  ############## Selection of top individuals (sire/dams) ###############
  #######################################################################
  {

    if(length(selection.m.database)>0 & selection.size[1]==0){
      selection.size[1] <- sum(selection.m.database[,4] - selection.m.database[,3] + 1)
      if(verbose) cat("No selection.size provided. Use all available selected individuals.\n")
      if(verbose) cat(paste0(selection.size[1], " male individuals selected.\n"))
      selection.size.calc <- TRUE
    }
    if(length(selection.f.database)>0 & selection.size[2]==0){
      selection.size[2] <- sum(selection.f.database[,4] - selection.f.database[,3] + 1)
      if(verbose) cat("No selection.size provided. Use all available selected individuals.\n")
      if(verbose) cat(paste0(selection.size[2], " female individuals selected.\n"))
      selection.size.calc <- TRUE
    }



    if(length(threshold.selection)==1 && is.na(threshold.selection)){
      threshold.selection <- NULL
    }

    if(length(threshold.selection)>0){
      if(threshold.sign=="<"){
        selection.highest <- c(FALSE, FALSE)
      }
    }

    {
      best <- list(matrix(0, nrow=selection.size[1],ncol=5),
                   matrix(0, nrow=selection.size[2],ncol=5))


      addsel <- c(2,2)
      if(selection.criteria[1]=="bv"){
        addsel[1] = 6
      }
      if(selection.criteria[2]=="bv"){
        addsel[2] = 6
      }
      if(selection.criteria[1]=="pheno"){
        addsel[1] = 8
      }
      if(selection.criteria[2]=="pheno"){
        addsel[2] = 8
      }

      if(selection.criteria[1]=="offpheno"){
        addsel[1] = 26
      }
      if(selection.criteria[2]=="offpheno"){
        addsel[2] = 26
      }

      if(add.class.cohorts && length(selection.m.cohorts)>0){
        for(index in 1:length(selection.m.cohorts)){
          take <- which(population$info$cohorts[,1]==selection.m.cohorts[index])
          class.m <- c(class.m, as.numeric(population$info$cohorts[take,5]))
        }
        class[[1]] <- class.m
      }
      if(add.class.cohorts && length(selection.f.cohorts)>0){
        for(index in 1:length(selection.f.cohorts)){
          take <- which(population$info$cohorts[,1]==selection.f.cohorts[index])
          class.f <- c(class.f, as.numeric(population$info$cohorts[take,5]))
        }
        class[[2]] <- class.f
      }




    }



    activ_groups_list <- list(selection.m.database, selection.f.database)
    chosen.animals.list <- list()
    selection.sex <- c(selection.m, selection.f)
    selection.miesenberger <- c(selection.m.miesenberger, selection.f.miesenberger)

    selection.random.prob = list(selection.m.random.prob, selection.f.random.prob)
    multiple.bve.weights <- list(multiple.bve.weights.m, multiple.bve.weights.f)
    multiple.bve.scale <- c(multiple.bve.scale.m, multiple.bve.scale.f)

    multiple.bve.scale[multiple.bve.scale=="bve"] <- "bve_sd"
    multiple.bve.scale[multiple.bve.scale=="pheno"] <- "pheno_sd"
    multiple.bve.scale[multiple.bve.scale=="bv"] <- "bv_sd"


    sd_scaling <- rep(1, population$info$bv.nr)
    if(length(fixed.breeding)==0 || length(fixed.breeding.best)>0){
      if(sum(selection.size)>0){
        if(verbose) cat("Start selection procedure.\n")
        for(sex in (1:2)[selection.size>0]){

          if(sex==1){
            if(verbose){cat("Selection male size:\n")}
          } else{
            if(verbose){cat("Selection female size:\n")}
          }
          possible_animals <- NULL
          activ_groups <- activ_groups_list[[sex]]
          for(index5 in 1:nrow(activ_groups)){
            possible_animals <- rbind(possible_animals, cbind(activ_groups[index5,1], activ_groups[index5,2], activ_groups[index5,3]:activ_groups[index5,4], population$breeding[[activ_groups[index5,1]]][[4+activ_groups[index5,2]]][activ_groups[index5,3]:activ_groups[index5,4]]))
          }

          if(length(reduced.selection.panel[[sex]])>0 && length(possible_animals)>0){
            possible_animals <- possible_animals[reduced.selection.panel[[sex]],,drop=FALSE]
          }
          relevant.animals <- rep(possible_animals[,4], length(class[[sex]])) == rep(class[[sex]],each=nrow(possible_animals))
          relevant.animals <- unique(c(0,relevant.animals * 1:nrow(possible_animals)))[-1]
          possible_animals <- rbind(NULL, possible_animals[unique(relevant.animals),])

          ids = get.id(population, database = possible_animals[,1:3])

          duplicates = (duplicated(ids[length(ids):1]))[length(ids):1]

          if(remove.duplicates && sum(duplicates)>0){
            possible_animals = possible_animals[!duplicates,]
          }




          if(length(threshold.selection.index)>0){

            keeps = rep(TRUE, nrow(possible_animals))
            for(index in 1:nrow(threshold.selection.index)){
              if(threshold.selection.criteria[index]=="bve"){
                check_value =get.bve(population, database = possible_animals[,1:3])
              }

              index_value  = colSums(check_value * threshold.selection.index[index,])

              keeps_partial = NULL

              eval(parse(text=paste0("keeps_partial <- index_value",threshold.selection.sign[index],"threshold.selection.value[index]")))

              keeps = keeps & keeps_partial
              if(verbose){
                cat(paste0("Threshold selection ", index, " is fullfilled by ", sum(keeps_partial), " out of ", nrow(possible_animals), ".\n"))
              }
            }

            if(verbose){
              cat(paste0("After threshold selection ", sum(keeps), " individuals remain to select ", selection.size[sex], ".\n"))
            }

            possible_animals = possible_animals[keeps,]
          } else{
            if(verbose){
              cat(paste0("Select ", selection.size[sex] , " individuals out of  ", nrow(possible_animals), ".\n"))
            }
          }

          n.animals <- nrow(possible_animals)

          is_selection <- selection.size[sex] > n.animals

          if(n.animals==0){
            stop("No available individuals for selection provided - check gen/database/cohorts and classes of your input!")
          } else if(n.animals<selection.size[sex]){

            if(!selection.size.calc){
              warning(paste0("Less individuals available for selection than given in selection.size.\n Automatically reduce the number of selected individuals to " ,n.animals))

            }
            selection.size[sex] <- n.animals
            best[[sex]] <- matrix(0, nrow=selection.size[sex],ncol=5)

            if(copy.individual==TRUE && sum(max.offspring * selection.size < breeding.size)>0){
              breeding.size[breeding.size > (selection.size * max.offspring)] <- (selection.size * max.offspring)[breeding.size > (selection.size * max.offspring)]
              breeding.size.total <- sum(breeding.size)
            }
          }

          if(selection.sex[sex]=="random"){
            if(population$info$bv.nr==0){
              import.bv <- rep(0, n.animals)
            } else if(population$info$bv.nr==1){
              import.bv <- rep(0, n.animals)
              for(index5 in 1:nrow(possible_animals)){
                import.bv[index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+addsel[possible_animals[index5,2]]]][,possible_animals[index5,3]]
              }
            } else{
              if(multiple.bve=="add"){
                breeding.values <- matrix(0, ncol=nrow(possible_animals), nrow= population$info$bv.nr)
                for(index5 in 1:nrow(possible_animals)){
                  breeding.values[,index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+addsel[sex]]][,possible_animals[index5,3]]
                }


                if(multiple.bve.scale[sex]=="bve_sd" || multiple.bve.scale[sex]=="pheno_sd" || multiple.bve.scale[sex]=="bv_sd"){

                  sd_store <- numeric(population$info$bv.nr)

                  if(multiple.bve.scale[sex]=="pheno_sd" || multiple.bve.scale[sex] =="bv_sd"){
                    pheno.values <- genomic.values <-  matrix(0, ncol=nrow(possible_animals), nrow= population$info$bv.nr)
                    for(index5 in 1:nrow(possible_animals)){
                      pheno.values[,index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+8]][,possible_animals[index5,3]]
                      genomic.values[,index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+6]][,possible_animals[index5,3]]
                    }
                  }

                  for(bven in 1:population$info$bv.nr){
                    breeding.values[bven,] <- breeding.values[bven,] - mean(breeding.values[bven,])
                    if(ncol(breeding.values)!=1){
                      if(multiple.bve.scale[sex]=="bve_sd"){
                        sd_store[bven] <- stats::sd(breeding.values[bven,])
                      } else if(multiple.bve.scale[sex]=="pheno_sd"){
                        sd_store[bven] <- stats::sd(pheno.values[bven,])
                        if(is.na(sd_store[bven])){
                          sd_store[bven] <- 0
                        }
                        if(sd_store[bven]==0){
                          if(verbose) cat("No observed phenotypes in the group of selected individuals. Sure you want to scale according to phenotypes?\n")
                          if(verbose) cat("Use residual variance for scaling!\n")
                          sd_store[bven] <- sqrt(sigma.e2[bven])
                          if(is.na(sd_store[bven])){
                            sd_store[bven] <- 0
                          }
                        }
                      } else {
                        sd_store[bven] <- stats::sd(genomic.values[bven,])
                        if(sd_store[bven]==0){
                          if(verbose & multiple.bve.weights[[sex]][bven]!=0 & is_selection) cat("No variation in true genomic values. Sure you want to scale according to true genomic values?\n")
                          if(verbose & multiple.bve.weights[[sex]][bven]!=0 & is_selection) cat("Use residual variance for scaling!\n")
                          sd_store[bven] <- sqrt(sigma.e2[bven])
                          if(is.na(sd_store[bven])){
                            sd_store[bven] <- 0
                          }
                        }
                      }

                    } else{
                      sd_store[bven] <- 1
                    }

                    if(sd_store[bven]>0){
                      breeding.values[bven,] <- breeding.values[bven,] / sd_store[bven]
                    }
                  }
                }

                bve.sum <- colSums(breeding.values * multiple.bve.weights[[sex]])

              } else if(multiple.bve=="ranking"){
                breeding.values <- matrix(0, ncol=nrow(possible_animals), nrow= population$info$bv.nr)
                for(index5 in 1:nrow(possible_animals)){
                  breeding.values[,index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+addsel[sex]]][,possible_animals[index5,3]]
                }
                ranking <- matrix(0, ncol=nrow(possible_animals), nrow= population$info$bv.nr)
                for(bven in 1:population$info$bv.nr){
                  order <- sort(breeding.values[bven,], index.return=TRUE, decreasing=selection.highest[sex])$ix
                  ranking[bven,order] <- length(order):1
                }
                bve.sum <- colSums(ranking*multiple.bve.weights[[sex]])
              }
              import.bv <- bve.sum
            }

            if(nrow(possible_animals)>0 & sum(abs(import.bv))>0){
              for(entry in 1:nrow(possible_animals)){
                population$breeding[[possible_animals[entry,1]]][[possible_animals[entry,2]+30]][possible_animals[entry,3]] <- import.bv[entry]
              }
            }


            chosen.animals <- sample(1:n.animals, selection.size[sex], prob = selection.random.prob[[sex]])
            best[[sex]][,1:3] <- possible_animals[chosen.animals,1:3]
            best[[sex]][,4] <- import.bv[chosen.animals]
          } else if(selection.sex[sex]=="function"){
            if(population$info$bv.nr==1){
              import.bv <- numeric(nrow(possible_animals))
              for(index5 in 1:nrow(possible_animals)){
                import.bv[index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+addsel[sex]]][,possible_animals[index5,3]]
              }
              import.bv[is.na(import.bv)] <- -Inf
              chosen.animals <- sort(import.bv, index.return=TRUE, decreasing=selection.highest[sex])$ix[1:sum(selection.size[[sex]])] # Diese 3er werden zu 4 in Weiblich
              if(sum(import.bv[chosen.animals[length(chosen.animals)]]==import.bv)>1){
                cutoff1 <- import.bv[chosen.animals[length(chosen.animals)]]
                chosen.animals <- NULL
                if(selection.highest[sex]){
                  chosen.animals <- which(import.bv>cutoff1)
                } else{
                  chosen.animals <- which(import.bv<cutoff1)
                }
                chosen.animals <- c(chosen.animals, sample(which(import.bv==cutoff1),  sum(selection.size[[sex]]) - length(chosen.animals)))
              }
              best[[sex]][,1:4] <- cbind(possible_animals[chosen.animals, 1], possible_animals[chosen.animals, 2], possible_animals[chosen.animals, 3], import.bv[chosen.animals])
            } else{
              if(multiple.bve=="add"){
                breeding.values <- matrix(0, ncol=nrow(possible_animals), nrow= population$info$bv.nr)
                for(index5 in 1:nrow(possible_animals)){
                  breeding.values[,index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+addsel[sex]]][,possible_animals[index5,3]]
                }
                breeding.values[is.na(breeding.values)] <- -Inf
                if((multiple.bve.scale[sex]=="bve_sd" || multiple.bve.scale[sex]=="pheno_sd" || multiple.bve.scale[sex] =="bv_sd") && selection.miesenberger[sex]==FALSE){

                  sd_scaling <- numeric(population$info$bv.nr)

                  if(multiple.bve.scale[sex]=="pheno_sd" || multiple.bve.scale[sex] =="bv_sd"){
                    pheno.values <- genomic.values <-  matrix(0, ncol=nrow(possible_animals), nrow= population$info$bv.nr)
                    for(index5 in 1:nrow(possible_animals)){
                      pheno.values[,index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+8]][,possible_animals[index5,3]]
                      genomic.values[,index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+6]][,possible_animals[index5,3]]
                    }
                  }

                  if(sum(abs(breeding.values))==0){
                    if(verbose) cat(paste0("No selection criteria or breeding values for selection provided. Continue with random selection (", if(sex==1){"male side"}else{"female side"} ,").\n"))
                    sd_scaling <- rep(1, population$info$bv.nr)
                  } else{
                    for(bven in 1:population$info$bv.nr){
                      breeding.values[bven,] <- breeding.values[bven,] - mean(breeding.values[bven,])
                      if(ncol(breeding.values)!=1){
                        if(multiple.bve.scale[sex]=="bve_sd"){
                          sd_scaling[bven] <- stats::sd(breeding.values[bven,], na.rm=TRUE)
                          if(sd_scaling[bven]==0 || is.na(sd_scaling[bven])){
                            if(verbose & multiple.bve.weights[[sex]][bven]!=0) cat(paste0("No estimated breeding values entered in the group of selected individuals (Trait: ", population$info$trait.name[bven], " - ", if(sex==1){"male side"}else{"female side"}, "). Please check your input variables!?\n"))
                            if(verbose & multiple.bve.weights[[sex]][bven]!=0) cat("Use residual variance!\n")
                            sd_scaling[bven] <- sqrt(sigma.e2[bven])
                          }
                        } else if(multiple.bve.scale[sex]=="pheno_sd") {
                          sd_scaling[bven] <- stats::sd(pheno.values[bven,], na.rm=TRUE)

                          if(is.na(sd_scaling[bven])){
                            sd_scaling[bven] <- 0
                          }
                          if(sd_scaling[bven]==0){
                            if(verbose & multiple.bve.weights[[sex]][bven]!=0) cat(paste0("No observed phenotypes in the group of selected individuals (Trait: ", population$info$trait.name[bven], " - ", if(sex==1){"male side"}else{"female side"},"). Sure you want to scale according to phenotypes?\n"))
                            if(verbose & multiple.bve.weights[[sex]][bven]!=0) cat("Expected phenotypic sd based on one observation was used!\n")
                            sd_scaling[bven] <- sqrt(stats::var(genomic.values[bven,]) + sigma.e2[bven])

                            if(is.na(sd_scaling[bven])){
                              sd_scaling[bven] <- 0
                            }

                          }
                        } else{
                          sd_scaling[bven] <- stats::sd(genomic.values[bven,], na.rm=TRUE)
                          if(sd_scaling[bven]==0){
                            if(verbose & multiple.bve.weights[[sex]][bven]!=0 & is_selection) cat(paste0("No variation in true genomic values in the group of selected individuals (Trait: ", population$info$trait.name[bven], " - ", if(sex==1){"male side"}else{"female side"},"). Sure you want to scale according to true genomic values?\n"))
                            if(verbose & multiple.bve.weights[[sex]][bven]!=0 & is_selection) cat("Expected phenotypic sd based on one observation was used!\n")
                            sd_scaling[bven] <- sqrt(stats::var(genomic.values[bven,]) + sigma.e2[bven])
                            if(is.na(sd_scaling[bven])){
                              sd_scaling[bven] <- 0
                            }
                          }
                        }

                      } else{
                        sd_scaling[bven] <- 1
                      }

                      if(sd_scaling[bven]>0){
                        breeding.values[bven,] <- breeding.values[bven,] / sd_scaling[bven]
                      }
                    }
                  }

                }


                if(selection.miesenberger[sex]){
                  genomic.values <- matrix(0, ncol=nrow(possible_animals), nrow= population$info$bv.nr)
                  bve.values <- matrix(0, ncol=nrow(possible_animals), nrow= population$info$bv.nr)
                  pheno.values <- matrix(0, ncol=nrow(possible_animals), nrow= population$info$bv.nr)
                  reliability.values <- matrix(0, ncol=nrow(possible_animals), nrow= population$info$bv.nr)
                  bve.index <- matrix(0, ncol=nrow(possible_animals), nrow= population$info$bv.nr)

                  for(index5 in 1:nrow(possible_animals)){
                    genomic.values[,index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+6]][,possible_animals[index5,3]]
                    bve.values[,index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+addsel[sex]]][,possible_animals[index5,3]]
                    pheno.values[,index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+8]][,possible_animals[index5,3]]
                    reliability.values[,index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+18]][,possible_animals[index5,3]]
                  }
                  if(sum(bve.values==0)==length(bve.values)){
                    if(verbose) cat("No breeding values for selected individuals! Miesenberger selection was skip.\n")

                    index.weights <- multiple.bve.weights[[sex]]
                    for(index5 in 1:nrow(possible_animals)){
                      bve.index[,index5] <- index.weights
                      population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+20]][,possible_animals[index5,3]] <- bve.index[,index5]
                    }
                  } else{
                    active_traits <- which(rowSums(bve.values==0)<nrow(bve.values))
                    if(length(active_traits)<population$info$bv.nr){
                      if(verbose) cat("Only use traits", population$info$trait.name[active_traits], "for Miesenberger.\n")
                      if(verbose) cat("Traits", population$info$trait.name[-active_traits], "with no BVE.\n")
                    }
                    genomic.cov <- stats::var(t(genomic.values[active_traits,,drop=FALSE]))
                    bve.cov <- stats::var(t(bve.values[active_traits,,drop=FALSE]))
                    pheno.cov <- stats::var(t(pheno.values[active_traits,colSums(is.na(pheno.values))==0,drop=FALSE]))

                    #                  diag(genomic.cov) <- diag(genomic.cov) + 10^-8
                    #                  diag(bve.cov) <- diag(bve.cov) + 10^-8
                    #                  diag(pheno.cov) <- diag(pheno.cov) + 10^-8

                    heri.emp <- diag(genomic.cov) / diag(pheno.cov)
                    derived.reliability <- diag(stats::cor(t(genomic.values[active_traits,,drop=FALSE]), t(bve.values[active_traits,,drop=FALSE])))^2
                    if(addsel[sex]==2){
                      estimated.reliability <- sqrt(diag(bve.cov) / diag(pheno.cov))
                    } else{
                      estimated.reliability <- numeric(nrow(pheno.cov))
                    }

                    estimated.reliability[diag(pheno.cov)==0 | is.na(diag(pheno.cov))] <- derived.reliability[diag(pheno.cov)==0 | is.na(diag(pheno.cov))]
                    estimated.reliability[estimated.reliability>1] <- 1


                    for(bven in 1:length(active_traits)){
                      trait_index <- active_traits[bven]
                      if(selection.miesenberger.reliability.est=="estimated"){
                        reliability.values[trait_index, reliability.values[trait_index,]==0] <- estimated.reliability[bven]
                      } else if(selection.miesenberger.reliability.est=="heritability"){
                        reliability.values[trait_index, reliability.values[trait_index,]==0] <- heri.emp[bven]
                      } else if(selection.miesenberger.reliability.est=="derived"){
                        reliability.values[trait_index, reliability.values[trait_index,]==0] <- derived.reliability[bven]
                      }

                    }
                    reliability.values[is.na(reliability.values)] <- 0

                    V <- bve.cov
                    V1 <- MASS::ginv(V)

                    # V1 <- chol2inv(chol(V))

                    G_cov <- genomic.cov
                    RG <- sqrt(diag(1/diag(G_cov))) %*% G_cov %*% sqrt(diag(1/diag(G_cov)))

                    if(multiple.bve.scale[sex]=="pheno_sd"){
                      index.weights <- multiple.bve.weights[[sex]][active_traits] / sqrt(diag(pheno.cov))
                      if(sum(is.na(diag(pheno.cov))) >0 || sum(diag(pheno.cov)==0)>0){
                        if(verbose) cat("Are you sure you want to scale Miesenberger w by phenotypic sd? No phenotypes for some traits!\n")
                        if(verbose) cat("Expected phenotypic sd was used!\n")
                        index.weights <- multiple.bve.weights[[sex]][active_traits] / sqrt(diag(genomic.cov) + sigma.e2[active_traits])
                      }
                    } else if(multiple.bve.scale[sex]=="unit"){
                      index.weights <- multiple.bve.weights[[sex]][active_traits]
                    } else if(multiple.bve.scale[sex]=="bve_sd"){
                      index.weights <- multiple.bve.weights[[sex]][active_traits] / sqrt(diag(bve.cov))
                      index.weights[diag(bve.cov)==0] <- 0
                    } else if(multiple.bve.scale[sex]=="bv_sd"){
                      index.weights <- multiple.bve.weights[[sex]][active_traits] / sqrt(diag(genomic.cov))
                      index.weights[diag(genomic.cov)==0] <- 0
                    }
                    if(miesenberger.trafo>0){

                      eigenV <- eigen(V)
                      keeps <- which(eigenV$values > (miesenberger.trafo*sum(eigenV$values)))

                      P <- eigenV$vectors[,keeps]

                      V_trafo <- t(P) %*% V %*% P
                      G_trafo <- t(P) %*% G_cov %*% P

                      V1_trafo <- MASS::ginv(V_trafo)
                      RG_trafo <- sqrt(diag(1/diag(G_trafo))) %*% G_trafo %*% sqrt(diag(1/diag(G_trafo)))

                      w_trafo <- t(P) %*% index.weights

                      breeding.values1 <- t(P) %*% breeding.values[active_traits, ]
                      bve.index1 <- matrix(0, nrow=nrow(breeding.values1), ncol=ncol(breeding.values1))
                      for(index5 in 1:nrow(possible_animals)){
                        r_trafo <- sqrt(colSums(reliability.values[active_traits,index5] * (eigenV$vectors[,keeps]^2)))
                        bve.index1[,index5] <- miesenberger.index(V1 = V1_trafo, V = V_trafo, G_trafo, RG = RG_trafo, r = r_trafo, w = w_trafo)
                      }


                    } else{
                      for(index5 in 1:nrow(possible_animals)){
                        bve.index[active_traits,index5] <- miesenberger.index(V1=V1, V= V, G = G_cov, RG = RG, r = sqrt(reliability.values[active_traits,index5]), w = index.weights)
                        population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+20]][,possible_animals[index5,3]] <- bve.index[,index5]
                      }

                    }

                    if(verbose) cat("Average reliablity estimated per trait:\n")
                    if(verbose) cat(rowMeans(reliability.values))

                    if(verbose) cat("\nAverage index used for in selection according to Miesenberger:\n")

                    if(miesenberger.trafo>0){
                      temp1 <- P %*% bve.index1
                      bve.index2 <- matrix(0, nrow=population$info$bv.nr, ncol=ncol(bve.index1))
                      bve.index2[active_traits,] <-  P %*% bve.index1
                    } else{
                      bve.index2 <- bve.index
                    }
                    if(multiple.bve.scale[sex]=="pheno_sd"){
                      if(sum(is.na(diag(pheno.cov)))>0 || sum(diag(pheno.cov)==0)>0){
                        if(verbose) cat(paste0(population$info$trait.name[active_traits], ": ", round(rowMeans(bve.index2)[active_traits] * sqrt(diag(genomic.cov) + sigma.e2[active_traits]), digits=2), " avg. index weightings per pSD.\n"))
                      } else{
                        if(verbose) cat(paste0(population$info$trait.name[active_traits], ": ",round(rowMeans(bve.index2)[active_traits] * sqrt(diag(pheno.cov)), digits=2), " avg. index weightings per pSD.\n"))
                      }
                    } else if(multiple.bve.scale[sex]=="unit"){
                      if(verbose) cat(paste0(population$info$trait.name[active_traits], ": ",round(rowMeans(bve.index2)[active_traits], digits=2), " avg. index weightings per Unit.\n"))
                    } else if(multiple.bve.scale[sex]=="bve_sd"){
                      if(verbose) cat(paste0(population$info$trait.name[active_traits], ": ",round(rowMeans(bve.index2)[active_traits] * sqrt(diag(bve.cov)), digits=2), " avg. index weightings per bveSD.\n"))
                    } else if(multiple.bve.scale[sex]=="bv_sd"){
                      if(verbose) cat(paste0(population$info$trait.name[active_traits], ": ",round(rowMeans(bve.index2)[active_traits] * sqrt(diag(genomic.cov)), digits=2), " avg. index weightings per gSD.\n"))
                    }


                  }

                  if(miesenberger.trafo>0){
                    bve.sum <- colSums(bve.index1 * (breeding.values1-rowMeans(breeding.values1)), na.rm=TRUE)
                  } else{
                    bve.sum <- colSums(bve.index * (breeding.values-rowMeans(breeding.values)), na.rm=TRUE)
                  }


                } else{
                  bve.sum <- colSums(breeding.values * multiple.bve.weights[[sex]], na.rm=TRUE)
                }
              } else if(multiple.bve=="ranking"){
                breeding.values <- matrix(0, ncol=nrow(possible_animals), nrow= population$info$bv.nr)
                for(index5 in 1:nrow(possible_animals)){
                  breeding.values[,index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+addsel[sex]]][,possible_animals[index5,3]]
                }
                breeding.values[is.na(breeding.values)] <- -Inf
                ranking <- matrix(0, ncol=nrow(possible_animals), nrow= population$info$bv.nr)
                for(bven in 1:population$info$bv.nr){
                  order <- sort(breeding.values[bven,], index.return=TRUE, decreasing=selection.highest[sex])$ix
                  ranking[bven,order] <- length(order):1
                }
                bve.sum <- colSums(ranking*multiple.bve.weights[[sex]], na.rm=TRUE)
              }

              if(nrow(possible_animals)>0 & sum(abs(bve.sum))>0){
                for(entry in 1:nrow(possible_animals)){
                  population$breeding[[possible_animals[entry,1]]][[possible_animals[entry,2]+30]][possible_animals[entry,3]] <- bve.sum[entry]
                }
              }
              chosen.animals <- sort(bve.sum, index.return=TRUE, decreasing=selection.highest[sex])$ix[1:sum(selection.size[sex])] # Diese 3er werden zu 4 in Weiblich

              if(sum(bve.sum[chosen.animals[length(chosen.animals)]]==bve.sum)>1){
                cutoff1 <- bve.sum[chosen.animals[length(chosen.animals)]]
                chosen.animals <- NULL
                if(selection.highest[sex]){
                  chosen.animals <- which(bve.sum>cutoff1)
                } else{
                  chosen.animals <- which(bve.sum<cutoff1)
                }
                chosen.animals <- c(chosen.animals, sample(which(bve.sum==cutoff1),  sum(selection.size[[sex]]) - length(chosen.animals)))
              }
              best[[sex]][,1:4] <- cbind(possible_animals[chosen.animals,1:3, drop=FALSE], bve.sum[chosen.animals])
            }

          }

          if(selection.miesenberger[sex]==FALSE){
            for(index5 in 1:nrow(possible_animals)){
              population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+20]][,possible_animals[index5,3]] <- multiple.bve.weights[[sex]]
            }
          }

          ## Just some follow up computation in case individuasl are not picked with equal share
          if(population$info$bv.nr==1){
            breeding.values <- numeric(length(chosen.animals))
            for(running in 1:length(chosen.animals)){
              if(best.selection.criteria[[sex]]=="bv"){
                breeding.values[running] <- population$breeding[[possible_animals[chosen.animals[running],1]]][[possible_animals[chosen.animals[running],2]+6]][1,possible_animals[chosen.animals[running],3]]
              } else if(best.selection.criteria[[sex]]=="pheno"){
                breeding.values[running] <- population$breeding[[possible_animals[chosen.animals[running],1]]][[possible_animals[chosen.animals[running],2]+8]][1,possible_animals[chosen.animals[running],3]]
              } else{
                breeding.values[running] <- population$breeding[[possible_animals[chosen.animals[running],1]]][[possible_animals[chosen.animals[running],2]+2]][1,possible_animals[chosen.animals[running],3]]
              }
            }
            breeding.values[is.na(breeding.values)] <- -Inf
            best[[sex]][,5] <- breeding.values
          } else if(population$info$bv.nr>1){
            if(multiple.bve=="add"){
              breeding.values <- matrix(0, nrow=population$info$bv.nr, ncol=length(chosen.animals))
              breeding.values.scaling <- matrix(0, nrow=population$info$bv.nr, ncol=length(chosen.animals))
              for(running in 1:length(chosen.animals)){
                if(best.selection.criteria[[sex]]=="bv"){
                  breeding.values[,running] <- population$breeding[[possible_animals[chosen.animals[running],1]]][[possible_animals[chosen.animals[running],2]+6]][,possible_animals[chosen.animals[running],3]]
                } else if(best.selection.criteria[[sex]]=="pheno"){
                  breeding.values[,running] <- population$breeding[[possible_animals[chosen.animals[running],1]]][[possible_animals[chosen.animals[running],2]+8]][,possible_animals[chosen.animals[running],3]]
                } else{
                  breeding.values[,running] <- population$breeding[[possible_animals[chosen.animals[running],1]]][[possible_animals[chosen.animals[running],2]+2]][,possible_animals[chosen.animals[running],3]]
                }
                breeding.values.scaling[,running] <- population$breeding[[possible_animals[chosen.animals[running],1]]][[possible_animals[chosen.animals[running],2]+20]][,possible_animals[chosen.animals[running],3]]

              }
              breeding.values[is.na(breeding.values)] <- -Inf

              if((multiple.bve.scale[sex]=="bve_sd" || multiple.bve.scale[sex]=="pheno_sd" || multiple.bve.scale[sex]=="bv_sd") && selection.miesenberger[sex]==FALSE){

                sd_scaling <- numeric(population$info$bv.nr)
                if(ncol(breeding.values)>1){
                  P <- stats::var(t(breeding.values))
                } else{
                  P <- diag(nrow(breeding.values))
                }

                if(multiple.bve.scale[sex]=="pheno_sd" || multiple.bve.scale[sex]=="bv_sd"){
                  pheno.values <- genomic.values <-  matrix(0, ncol=nrow(possible_animals), nrow= population$info$bv.nr)
                  for(index5 in 1:nrow(possible_animals)){
                    pheno.values[,index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+8]][,possible_animals[index5,3]]
                    genomic.values[,index5] <- population$breeding[[possible_animals[index5,1]]][[possible_animals[index5,2]+6]][,possible_animals[index5,3]]
                  }
                }
                if(sum(abs(breeding.values))==0){
                  sd_scaling <- rep(1, population$info$bv.nr)
                } else{
                  for(bven in 1:population$info$bv.nr){
                    breeding.values[bven,] <- breeding.values[bven,] - mean(breeding.values[bven,])
                    if(ncol(breeding.values)!=1){
                      if(multiple.bve.scale[sex]=="bve_sd"){
                        sd_scaling[bven] <- stats::sd(breeding.values[bven,])
                        if(sd_scaling[bven]==0 || is.na(sd_scaling[bven])){
                          if(verbose & multiple.bve.weights[[sex]][bven]!=0 & selection.sex[sex]!= "random") cat("No estimated breeding values entered in the group of selected individuals. Please check your input variables!?\n")
                          if(verbose & multiple.bve.weights[[sex]][bven]!=0 & selection.sex[sex]!= "random") cat("Use residual variance!\n")
                          sd_scaling[bven] <- sqrt(sigma.e2[bven])
                        }
                      } else if(multiple.bve.scale[sex] =="pheno_sd"){
                        sd_scaling[bven] <- stats::sd(pheno.values[bven,])
                        if(is.na(sd_scaling[bven]) || sd_scaling[bven]==0){
                          if(verbose & multiple.bve.weights[[sex]][bven]!=0 & selection.sex[sex]!= "random") cat("No observed phenotypes in the group of selected individuals. Sure you want to scale according to phenotypes?\n")
                          if(verbose & multiple.bve.weights[[sex]][bven]!=0 & selection.sex[sex]!= "random") cat("Expected phenotypic sd based on one observation was used!\n")
                          sd_scaling[bven] <- sqrt(stats::var(genomic.values[bven,]) + sigma.e2[bven])
                        }
                      } else{
                        sd_scaling[bven] <- stats::sd(genomic.values[bven,])
                        if(sd_scaling[bven]==0){
                          if(verbose & multiple.bve.weights[[sex]][bven]!=0 & selection.sex[sex]!= "random") cat("No variation in true genomic values in the group of selected individuals. Sure you want to scale according to true genomic values?\n")
                          if(verbose & multiple.bve.weights[[sex]][bven]!=0 & selection.sex[sex]!= "random") cat("Expected phenotypic sd based on one observation was used!\n")
                          sd_scaling[bven] <- sqrt(stats::var(genomic.values[bven,]) + sigma.e2[bven])
                          if(is.na(sd_scaling[bven])){
                            sd_scaling[bven] <- 0
                          }
                        }
                      }

                    } else{
                      sd_scaling[bven] <- 1
                    }

                    if(sd_scaling[bven]>0){
                      breeding.values[bven,] <- breeding.values[bven,] / sd_scaling[bven]
                    }
                  }
                }

              }
              bve.sum <- colSums(breeding.values * breeding.values.scaling, na.rm=TRUE)
            } else if(multiple.bve=="ranking"){
              ranking <- matrix(0, nrow=population$info$bv.nr, ncol=length(chosen.animals))
              for(running in 1:length(chosen.animals)){
                if(best.selection.criteria[[sex]]=="bv"){
                  ranking[,running] <- population$breeding[[possible_animals[chosen.animals[running],1]]][[possible_animals[chosen.animals[running],2]+6]][,possible_animals[chosen.animals[running],3]]
                } else if(best.selection.criteria[[sex]]=="pheno"){
                  ranking[,running] <- population$breeding[[possible_animals[chosen.animals[running],1]]][[possible_animals[chosen.animals[running],2]+8]][,possible_animals[chosen.animals[running],3]]
                } else{
                  ranking[,running] <- population$breeding[[possible_animals[chosen.animals[running],1]]][[possible_animals[chosen.animals[running],2]+2]][,possible_animals[chosen.animals[running],3]]
                }
              }
              ranking[is.na(ranking)] <- -Inf
              for(bven in 1:population$info$bv.nr){
                order <- sort(ranking[bven,], index.return=TRUE, decreasing=selection.highest[sex])$ix
                ranking[bven,order] <- length(order):1
              }
              bve.sum <- colSums(ranking*breeding.values.scaling, na.rm=TRUE)
            }
            best[[sex]][,5] <- bve.sum
          }

          reorder <- sort(best[[sex]][,4], index.return=TRUE, decreasing = selection.highest[sex])$ix
          best[[sex]] <- best[[sex]][reorder,,drop=FALSE]
          chosen.animals <- chosen.animals[reorder]
          if(ignore.best[sex]>0){
            best[[sex]] <- best[[sex]][-(1:ignore.best[sex]),,drop=FALSE]
            chosen.animals <- chosen.animals[-(1:ignore.best[sex])]
            selection.size[sex] <- selection.size[sex] - ignore.best[sex]

          }

          chosen.animals.list[[sex]] <- chosen.animals

          if(!best.selection.manual.reorder){
            new_order <- numeric(nrow(best[[sex]]))
            temp1 <- t(best[[sex]][,1:3])
            for(index in 1:nrow(best[[sex]])){
              new_order[index] <- which(colSums(temp1==possible_animals[index,1:3])==3)
            }
            best.selection.manual.ratio[[sex]][new_order] <- best.selection.manual.ratio[[sex]]
          }
        }

        if(fixed.assignment!=FALSE){
          nm <- nrow(best[[1]])
          nw <- nrow(best[[2]])
          male.rounds <- breeding.size.total / nm
          female.rounds <- breeding.size.total / nw
          full.m <- nm*(male.rounds - round(male.rounds - 0.5))
          full.f <- nw*(female.rounds - round(female.rounds -0.5))
          fixed.assignment.m <- c(rep(1:full.m,each=male.rounds-full.m/nm+1), rep((full.m+1):nm, each =male.rounds-full.m/nm))
          fixed.assignment.f <- c(rep(1:full.f,each=female.rounds - full.f/nw +1), rep((full.f+1):nw, each =female.rounds-full.f/nw))
          fixed.assignment.m <- fixed.assignment.m[1:breeding.size.total]
          fixed.assignment.f <- fixed.assignment.f[1:breeding.size.total]
          if(fixed.assignment=="bestworst" || fixed.assignment=="worstbest"){
            fixed.assignment.f <- fixed.assignment.f[breeding.size.total:1]
          }
          fixed.breeding <- cbind(best[[1]][fixed.assignment.m,1:3], best[[2]][fixed.assignment.f,1:3], sex.animal)
        }
        if(ogc){
          animallist <- rbind(cbind(best[[1]],1), cbind(best[[2]],2))
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
                Zt[,index] <- colSums(compute.snps(population, animallist[index,1], animallist[index,2], animallist[index,3], import.position.calculation=import.position.calculation, decodeOriginsU=decodeOriginsU, bit.storing=bit.storing, nbits=nbits, output_compressed=FALSE))
              }
            }
          }


          # Verwandtschaftsmatrix:
          if(relationship.matrix.ogc=="kinship" || relationship.matrix.ogc == "pedigree"){
            A <- kinship.exp(population, database = animallist[,c(1,2,3,3)], depth.pedigree = depth.pedigree.ogc, verbose=verbose, mult=2)

            if(store.sparse){
              A <- methods::as(A, "sparseMatrix")
            }
          } else if(relationship.matrix.ogc=="vanRaden"){
            if(miraculix){
              #p_i <- miraculix::allele_freq(Z.code) # Noch nicht implementiert?
              p_i <- rowMeans(as.matrix(Z.code))/2 # Noch nicht implementiert?
              A <- miraculix::relationshipMatrix(Z.code, centered=TRUE, normalized=TRUE)

            } else if(miraculix.mult){
              p_i <- rowSums(Zt)/ncol(Zt)/2
              Zt_miraculix <- miraculix::genomicmatrix(Zt)
              A <- miraculix::relationshipMatrix(Zt_miraculix, centered=TRUE, normalized=TRUE)
            } else{
              p_i <- rowSums(Zt)/ncol(Zt)/2
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

          Indiv <- paste0("Indi", 1:length(BV))
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
            con$ub = ogc.ub
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

          contribution <- list(Offspring$parent$oc[Sex=="male"], Offspring$parent$oc[Sex=="female"])

          if(verbose){
            cat(paste0(sum(contribution[[1]]>0), " male individuals with positive contribution ((", sum(contribution[[1]]>(0.001 * max(contribution[[1]]))), " with major contribution).\n"))
            cat(paste0(sum(contribution[[2]]>0), " female individuals with positive contribution ((", sum(contribution[[2]]>(0.001 * max(contribution[[2]]))), " with major contribution)\n."))
          }
        }
      }

    } else{
      breeding.size.total <- nrow(fixed.breeding)
      sex.animal <- fixed.breeding[,7] <- stats::rbinom(breeding.size.total, 1, fixed.breeding[,7]) +1
      breeding.size <- c(sum(fixed.breeding[,7]==1), sum(fixed.breeding[,7]==2))

      if(sum(is.na(breeding.size))>0){
        stop("Sex assignment for fixed breeding went wrong!\n 7th column in fixed breeding should contain probability to be female (value between 0 and 1).")
      }

    }

    if(length(threshold.selection)>0){
      if(length(best[[1]])>0){
        eval(parse(text=paste0("keeps <- which(best[[1]][,4]",threshold.sign,"threshold.selection)")))

        best[[1]] <- best[[1]][keeps,]

        if(verbose) cat(paste0(length(keeps), " sires fullfil threshold selection!\n"))
        selection.size[1] <- nrow(best[[1]])
        if(max.offspring[1]*selection.size[1] < breeding.size[1]){
          breeding.size[1] <- max.offspring[1]*selection.size[1]
          if(verbose) cat("Breeding size reduced based on threshold selection.\n")
        }
      }
      if(length(best[[2]])>0){
        eval(parse(text=paste0("keeps <- which(best[[2]][,4]",threshold.sign,"threshold.selection)")))
        best[[2]] <- best[[2]][keeps,]
        if(verbose) cat(paste0(length(keeps), " dams fullfil threshold selection!\n"))
        selection.size[2] <- nrow(best[[2]])
        if(max.offspring[2]*selection.size[2] < breeding.size[2]){
          breeding.size[2] <- max.offspring[2]*selection.size[2]
          if(verbose) cat("Breeding size reduced based on threshold selection")
        }
      }

      breeding.size.total <- sum(breeding.size)
    }

    sample_prob <- list()
    for(sex in 1:2){
      if(length(best[[sex]])>0){
        if(ogc){
          sample_prob[[sex]] <- contribution[[sex]]
        } else if(length(best.selection.manual.ratio[[sex]])==0){
          individual_bv <- best[[sex]][,5]
          min_bv <- min(individual_bv)
          max_bv <- max(individual_bv)
          if(min_bv==max_bv){
            sample_prob[[sex]] <- rep(1, length(individual_bv))
          } else{
            sample_prob[[sex]] <- (individual_bv- min_bv) / (max_bv-min_bv) * (best.selection.ratio[[sex]]-1) + 1
          }
        } else{
          sample_prob[[sex]] <- best.selection.manual.ratio[[sex]]
        }

        if(length(unique(sample_prob[[sex]]))==1){
          sample_prob[[sex]] <- NULL
        }
      }
    }
    sample_prob[[3]] <- "placeholder"



    if(length(fixed.breeding.best)>0){
      fixed.breeding <- matrix(0, nrow=nrow(fixed.breeding.best), ncol=7)

      for(index in 1:nrow(fixed.breeding.best)){
        fixed.breeding[index,1:3] <- best[[fixed.breeding.best[index,1]]][fixed.breeding.best[index,2], 1:3]
        fixed.breeding[index,4:6] <- best[[fixed.breeding.best[index,3]]][fixed.breeding.best[index,4], 1:3]
      }

      fixed.breeding[,7] <- fixed.breeding.best[,5]
      if(breeding.size.total==0){
        breeding.size.total <- nrow(fixed.breeding)
      }
      sex.animal <- fixed.breeding[,7] <- stats::rbinom(breeding.size.total, 1, fixed.breeding[,7]) +1
      breeding.size <- c(sum(fixed.breeding[,7]==1), sum(fixed.breeding[,7]==2))
    }



    if(sort.selected.pos){
      for(sex in 1:2){
        if(length(best[[sex]])>0){
          pos <- best[[sex]][,1] * (max(population$info$size)+5) * (length(population$breeding)+5) + best[[sex]][,2] * (max(population$info$size)+5) + best[[sex]][,3]

          best[[sex]] <- best[[sex]][sort(pos, index.return=TRUE)$ix,]
        }

      }
    }
    if(export.selected){
      return(best)
    }

    if(store.comp.times){
      comp.times[6] <- as.numeric(Sys.time())
    }
  }

  #######################################################################
  ############# Gene editing based on Simianer et al. 2018 ##############
  #######################################################################
  {
    if(gene.editing){
      if(estimate.u){
        effect_order<- sort(abs(population$info$u_hat[[length(population$info$u_hat)]]), index.return=TRUE, decreasing=TRUE)$ix
        direction <- (population$info$u_hat[[length(population$info$u_hat)]] < 0)
        distance_actual <- integer(length(effect_order))
        for(index in 1:length(effect_order)){
          distance_actual[index] <- min(abs(effect_order[index]-population$info$effect.p))
        }
        if(length(population$info$editing_info)==0){
          population$info$editing_info <- list()
        }
        population$info$editing_info[[length(population$info$editing_info)+1]] <- cbind(effect_order, direction[effect_order], distance_actual, deparse.level = 0)
      }
      if(gwas.u){
        effect_order <- sort(abs(population$info$gwas_hat[[length(population$info$gwas_hat)]]), index.return=TRUE, decreasing=TRUE)$ix
        direction <- (population$info$gwas_hat[[length(population$info$gwas_hat)]] <0 )
        distance_actual <- integer(length(effect_order))
        for(index in 1:length(effect_order)){
          distance_actual[index] <- min(abs(effect_order[index]-population$info$effect.p))
        }
        if(length(population$info$editing_info)==0){
          population$info$editing_info <- list()
        }
        population$info$editing_info[[length(population$info$editing_info)+1]] <- cbind(effect_order, direction[effect_order], distance_actual, deparse.level = 0)
      }
    }



    if(gene.editing.best){
      for(sex in (1:2)[gene.editing.best.sex]){
        if(length(best[[sex]])>0){
          for(index in 1:nrow(best[[sex]])){
            activ <- best[[sex]][index,]
            population$breeding[[activ[1]]][[activ[2]]][[activ[3]]] <- edit_animal(population, activ[1], activ[2], activ[3], nr.edits, decodeOriginsU=decodeOriginsU, bit.storing=bit.storing, nbits=nbits)
          }
          for(index in 1:nrow(best[[sex]])){
            activ <- best[[sex]][index,]
            if(length(population$breeding[[activ[1]]][[activ[2]]][[activ[3]]])>=17){
              population$breeding[[activ[1]]][[activ[2]]][[activ[3]]][[17]] <- c(population$breeding[[activ[1]]][[activ[2]]][[activ[3]]][[17]],population$breeding[[activ[1]]][[6+activ[2]]][,activ[3]])
            } else{
              population$breeding[[activ[1]]][[activ[2]]][[activ[3]]][[17]] <- population$breeding[[activ[1]]][[6+activ[2]]][,activ[3]]
            }
            population$breeding[[activ[1]]][[activ[2]]][[activ[3]]][[15]] <- rep(0L, population$info$bv.nr)
            activ_bv <- population$info$bv.random.activ
            if(length(activ_bv)>0){
              temp_out <- calculate.bv(population, activ[1], activ[2], activ[3], activ_bv, import.position.calculation=import.position.calculation, decodeOriginsU=decodeOriginsU, bit.storing=bit.storing, nbits=nbits, output_compressed=FALSE, bv.ignore.traits=bv.ignore.traits)
              population$breeding[[activ[1]]][[6+activ[2]]][activ_bv,activ[3]] <- temp_out[[1]]
            }
            population$breeding[[activ[1]]][[activ[2]]][[activ[3]]][[17]] <- c(population$breeding[[activ[1]]][[activ[2]]][[activ[3]]][[17]] ,population$breeding[[activ[1]]][[6+activ[2]]][,activ[3]])
          }
        }
      }
    }
  }

  #######################################################################
  ############## Old Culling model - make ready for delete ##############
  #######################################################################
  {
    if(is.matrix(reduce.group)==0 && length(reduce.group)>0){
      reduce.group <- t(reduce.group)
    }
    if(length(reduce.group)>0){
      reduce.group <- reduce.group[reduce.group[,1]!=0,]
    }
    if(length(reduce.group)>0){
      for(index in 1:nrow(reduce.group)){
        activ.reduce <- reduce.group[index,]
        group.animals <- (population$breeding[[activ.reduce[1]]][[activ.reduce[2]+4]] == activ.reduce[4])
        n.animals <- sum(group.animals)
        to.kill <- n.animals - activ.reduce[3]
        animal.position <- unique(c(0,group.animals * 1:length(group.animals)))[-1]
        if(to.kill>0){
          if(reduce.group.selection=="random"){
            death <- sample(animal.position, to.kill)

          } else if(reduce.group.selection=="selection"){
            if(multiple.bve=="add"){
              if(verbose) cat(population$breeding[[activ.reduce[1]]][[activ.reduce[2]+2]])
              bve.sum <- colSums(rbind(population$breeding[[activ.reduce[1]]][[activ.reduce[2]+2]][,animal.position]*multiple.bve.weights[[sex]],0))
            } else if(multiple.bve=="ranking"){
              ranking <- population$breeding[[index]][[sex+2]][,relevant.animals]
              for(bven in 1:population$info$bv.nr){
                order <- sort(ranking[bven,], index.return=TRUE, decreasing=selection.highest[activ.reduce[2]])$ix
                ranking[bven,order] <- length(order):1
              }
              bve.sum <- colSums(ranking*multiple.bve.weights[[sex]])
            }
            chosen.animals <- sort(bve.sum, index.return=TRUE, decreasing=selection.highest[activ.reduce[2]])$ix
            death <- chosen.animals[(length(chosen.animals)-to.kill+1):length(chosen.animals)]
          }
          for(modanimal in death){
            population$breeding[[activ.reduce[1]]][[activ.reduce[2]]][[modanimal]][[18]] <- c(population$breeding[[activ.reduce[1]]][[activ.reduce[2]+4]][modanimal], current.gen)
            population$breeding[[activ.reduce[1]]][[activ.reduce[2]+4]][modanimal] <- -1
          }
        }


      }
    }

  }

  #######################################################################
  #################### Generation of new individuals ####################
  #######################################################################
  {

    # preparation of the population - list // counting for max.offspring etc.
    {

      if(length(population$breeding)<=current.gen && sum(breeding.size)>0 ){
        population$breeding[[current.gen+1]] <- list()
        population$info$size <- rbind(population$info$size, 0)
      }

      selection.rate <- list(numeric(selection.size[1]), numeric(selection.size[2]))
      selection.rate.litter <- list(numeric(selection.size[1]), numeric(selection.size[2]))
      activ.selection.size <- selection.size
      availables.m <- 1:selection.size[1]
      availables.f <- 1:selection.size[2]
      availables <- list(availables.m, availables.f)

      current.size <- start.size <- c(1,1)
      for(sex in 1:2){
        if(sum(breeding.size)>0){
          if(length(population$breeding[[current.gen+1]])<=2 || length(population$breeding[[current.gen+1]][[sex]])==0){

            temp0 = matrix(0L, nrow=population$info$bv.nr, ncol=breeding.size[sex]) # this matrix will be enter in multiple positions
            population$breeding[[current.gen+1]][[2+sex]] <- temp0
            population$breeding[[current.gen+1]][[4+sex]] <- rep(new.class[sex], breeding.size[sex])
            population$breeding[[current.gen+1]][[6+sex]] <- temp0
            population$breeding[[current.gen+1]][[8+sex]] <- matrix(NA, nrow=population$info$bv.nr, ncol=breeding.size[sex])
            population$breeding[[current.gen+1]][[10+sex]] <- rep(time.point, breeding.size[sex])
            population$breeding[[current.gen+1]][[12+sex]] <- rep(creating.type, breeding.size[sex])
            if(copy.individual){
              population$breeding[[current.gen+1]][[14+sex]] <- rep(0L, breeding.size[sex])
            } else{
              population$breeding[[current.gen+1]][[14+sex]] <- seq(population$info$next.animal, population$info$next.animal + breeding.size[sex] -1, length.out= breeding.size[sex])
              population$info$next.animal <- population$info$next.animal + breeding.size[sex]
            }
            population$breeding[[current.gen+1]][[16+sex]] <- rep(NA, breeding.size[sex])
            population$breeding[[current.gen+1]][[18+sex]] <- temp0
            population$breeding[[current.gen+1]][[20+sex]] <- temp0
            if(copy.individual){
              #placeholder
              population$breeding[[current.gen+1]][[22+sex]] <- rep(-1, breeding.size[sex])
            } else{
              population$breeding[[current.gen+1]][[22+sex]] <- rep(time.point, breeding.size[sex])
            }
            population$breeding[[current.gen+1]][[24+sex]] <- rep(NA, breeding.size[sex])
            population$breeding[[current.gen+1]][[26+sex]] <- temp0
            population$breeding[[current.gen+1]][[28+sex]] <- temp0

            population$breeding[[current.gen+1]][[30+sex]] <- rep(0L, ncol=breeding.size[sex])

            population$breeding[[current.gen+1]][[32+sex]] <- rep(0L, breeding.size[sex])
            population$breeding[[current.gen+1]][[34+sex]] <- rep(0L, breeding.size[sex])
            population$breeding[[current.gen+1]][[36+sex]] <- rep(NA, breeding.size[sex])

          } else{

            temp0 = matrix(0L, nrow=population$info$bv.nr, ncol=breeding.size[sex]) # this matrix will be enter in multiple positions

            current.size[sex] <- start.size[sex] <- length(population$breeding[[current.gen+1]][[4+sex]]) + 1
            population$breeding[[current.gen+1]][[2+sex]] <- cbind(population$breeding[[current.gen+1]][[sex+2]], temp0)
            population$breeding[[current.gen+1]][[4+sex]] <- c(population$breeding[[current.gen+1]][[sex+4]], rep(new.class[sex], breeding.size[sex]))
            population$breeding[[current.gen+1]][[6+sex]] <- cbind(population$breeding[[current.gen+1]][[6+sex]], temp0)
            population$breeding[[current.gen+1]][[8+sex]] <- cbind(population$breeding[[current.gen+1]][[8+sex]], matrix(NA, nrow= population$info$bv.nr, ncol=breeding.size[sex]))
            population$breeding[[current.gen+1]][[10+sex]] <- c(population$breeding[[current.gen+1]][[sex+10]], rep(time.point, breeding.size[sex]))
            population$breeding[[current.gen+1]][[12+sex]] <- c(population$breeding[[current.gen+1]][[sex+12]], rep(creating.type, breeding.size[sex]))

            if(copy.individual){
              population$breeding[[current.gen+1]][[14+sex]] <- c(population$breeding[[current.gen+1]][[14+sex]], rep(0L,breeding.size[sex]))
            } else{
              population$breeding[[current.gen+1]][[14+sex]] <- c(population$breeding[[current.gen+1]][[14+sex]], seq(population$info$next.animal, population$info$next.animal + breeding.size[sex] -1, length.out= breeding.size[sex]))
              population$info$next.animal <- population$info$next.animal + breeding.size[sex]
            }
            population$breeding[[current.gen+1]][[16+sex]] <- c(population$breeding[[current.gen+1]][[sex+16]], rep(NA, breeding.size[sex]))
            population$breeding[[current.gen+1]][[18+sex]] <- cbind(population$breeding[[current.gen+1]][[18+sex]], temp0)
            population$breeding[[current.gen+1]][[20+sex]] <- cbind(population$breeding[[current.gen+1]][[20+sex]], temp0)
            if(copy.individual){
              #placeholder
              population$breeding[[current.gen+1]][[22+sex]] <- c(population$breeding[[current.gen+1]][[sex+22]], rep(-1, breeding.size[sex]))

            } else{
              population$breeding[[current.gen+1]][[22+sex]] <- c(population$breeding[[current.gen+1]][[sex+22]], rep(time.point, breeding.size[sex]))

            }#
            population$breeding[[current.gen+1]][[24+sex]] <- c(population$breeding[[current.gen+1]][[sex+24]], rep(NA, breeding.size[sex]))
            population$breeding[[current.gen+1]][[26+sex]] <- cbind(population$breeding[[current.gen+1]][[26+sex]], temp0)
            population$breeding[[current.gen+1]][[28+sex]] <- cbind(population$breeding[[current.gen+1]][[28+sex]], temp0)
            population$breeding[[current.gen+1]][[30+sex]] <- c(population$breeding[[current.gen+1]][[30+sex]], rep(0L, breeding.size[sex]))

            population$breeding[[current.gen+1]][[32+sex]] <- c(population$breeding[[current.gen+1]][[sex+32]], rep(0L, breeding.size[sex]))
            population$breeding[[current.gen+1]][[34+sex]] <- c(population$breeding[[current.gen+1]][[sex+34]], rep(0L, breeding.size[sex]))
            population$breeding[[current.gen+1]][[36+sex]] <- c(population$breeding[[current.gen+1]][[sex+36]], rep(NA, breeding.size[sex]))

          }
          if(length(population$breeding[[current.gen+1]][[sex]])==0){
            population$breeding[[current.gen+1]][[sex]] <- list()
          }

        }
      }

      ##  store.effect.freq and multiple correlated bvs deactivated
      if(store.comp.times.generation){
        pre_stuff <- 0
        generation_stuff <- 0
        bv_stuff <- 0
        parallel_gen <- 0
        parallel_joint <- 0

      }

      store.comp.times.generation.temp <- store.comp.times.generation
      if(generation.cores>1){
        store.comp.times.generation.temp2 <- store.comp.times.generation.temp
        store.comp.times.generation.temp <- FALSE
      } else{
        store.comp.times.generation.temp2 <- FALSE
      }
      if(length(name.cohort)==0 && breeding.size.total>0){
        name.cohort <- paste0("Cohort_", population$info$cohort.index)
        population$info$cohort.index <- population$info$cohort.index + 1
      }

      if(copy.individual==FALSE && breeding.size.total>0){
        if(nrow(best[[1]])==0){
          if(same.sex.activ==FALSE || same.sex.sex > 0){
            if(verbose){
              if(!dh.mating && !selfing.mating) {cat("No male / first parents (selection.m.gen/database/cohorts) provided for reproduction. Automatically allow female X female (second parent x second parent) matings.\n")}
            }
          }
          same.sex.activ <- TRUE
          same.sex.sex <- selfing.sex <- 1
        } else if(nrow(best[[2]])==0){
          if(same.sex.activ==FALSE || same.sex.sex < 1){
            if(verbose){
              if(!dh.mating && !selfing.mating) cat("No females / second parents (selection.f.gen/database/cohorts) provided for reproduction. Automatically allow male X male (first parent x first parent) matings.\n")
            }
          }
          same.sex.activ <- TRUE
          same.sex.sex <- selfing.sex <-  0
        }
      }

    }

    if(breeding.size.total>0){
      if(length(name.cohort)>0){
        if(verbose) cat(paste0("Start generation of new individuals (cohort: ", name.cohort,").\n"))
      } else{
        if(verbose) cat("Start generation of new individuals.\n")
      }

      pb_temp <- round(breeding.size.total/100)

      if(generation.cores>1){
        generation.cores <- min(ceiling(breeding.size.total/500), generation.cores)
        # Not worth to run something in parallel with low number of individuals individuals
      }

      if(generation.cores > 1){
        fixed.breeding_parallel <- matrix(0, nrow=breeding.size.total, ncol=6)
      }

      if(display.progress & verbose & generation.cores==1){
        pb <- utils::txtProgressBar(min = 0, max = breeding.size.total, style = 3)
      }

      runs <- 0

      if(max.mating.pair < Inf){
        info_father_list <- info_mother_list <- matrix(0, nrow=breeding.size.total, ncol=5)
      }

      if(length(repeat.mating.fixed)>0){
        repeat.mating.fixed <- rep(repeat.mating.fixed, length.out = breeding.size.total)
      }
      rmf <- 1 # repeat.mating temp-variable


      if(pen.assignments){


        assignment <- numeric(breeding.size.total)


        if(pen.by.sex == FALSE){

          for(sex in 1:2){
            assigned <- 0
            breeding.size.total_temp <- sum(sex.animal==sex)
            while(assigned < breeding.size.total_temp){
              next_pen <- min(sample(pen.size[,1],1, prob = pen.size[,2]), breeding.size.total_temp-assigned)

              assignment[sample((1:breeding.size.total)[assignment==0 & sex.animal==sex], next_pen)] <- population$info$next.pen
              population$info$next.pen <- population$info$next.pen + 1
              assigned <- assigned + next_pen
              print(c(assigned, sum(assignment[sex.animal==sex]>0)))
            }
          }

        } else{
          assigned <- 0
          while(assigned < breeding.size.total){
            next_pen <- min(sample(pen.size[,1],1, prob = pen.size[,2]), breeding.size.total-assigned)
            if(pen.by.litter){
              assignment[1:next_pen + assigned] <- population$info$next.pen
            } else{
              assignment[sample((1:breeding.size.total)[assignment==0], next_pen)] <- population$info$next.pen
            }

            population$info$next.pen <- population$info$next.pen + 1
            assigned <- assigned + next_pen
          }
        }




      }



      max_check <- sum(c(max.litter<Inf , max.offspring<Inf, selfing.mating==TRUE, same.sex.activ==TRUE))>0

      if(!max_check && length(sample_prob[[1]])==0 && length(sample_prob[[2]])==0) {
        fast_mode <- TRUE
        sex1 <- 1
        sex2 <- 2
        number1_vector <- availables[[sex1]][sample(activ.selection.size[sex1],breeding.size.total, replace=TRUE, prob=sample_prob[[sex1]][availables[[sex1]]])]
        number2_vector <- availables[[sex2]][sample(activ.selection.size[sex2],breeding.size.total, replace=TRUE, prob=sample_prob[[sex2]][availables[[sex2]]])]

      } else{
        fast_mode <- FALSE
      }

      if(store.comp.times.generation.temp2){
        tick <- as.numeric(Sys.time())
        tack <- as.numeric(Sys.time())
      }

      if(length(repeat.mating.activ)==1 && repeat.mating.activ=="genetic"){
        repeat.mating.genetic = TRUE
      } else{
        repeat.mating.genetic = FALSE
      }

      litter_counter = 1
      for(animal.nr in 1:breeding.size.total){
        if(store.comp.times.generation.temp){
          tick <- as.numeric(Sys.time())
        }

        if(repeat.mating.genetic && litter_counter > breeding.size.litter){
          break()
        }
        sex <- sex.animal[animal.nr]
        new.bv <-  new.bve <- new.reli <- individual.id <- rep(0L, population$info$bv.nr)
        new.bv_approx <- rep(NA, population$info$bv.nr)
        calc_litter = FALSE
        if(length(fixed.breeding)>0){
          info.father <- fixed.breeding[animal.nr,1:3]
          info.mother <- fixed.breeding[animal.nr,4:6]
        } else{
          if(runs>0){
            runs <- runs - 1

          } else{
            if(length(repeat.mating.s)>0){

              if(repeat.mating.genetic){
                calc_litter = TRUE
              }

              repeat.mating.temp = repeat.mating.s[litter_counter]
              litter_counter = litter_counter + 1
            } else{

              if(copy.individual){

                tt <- nrow(population$info$repeat.mating.copy)
                if(tt>1){
                  mating.temp <- sample(1:tt, prob = population$info$repeat.mating.copy[,2], 1)
                } else{
                  mating.temp <- 1
                }
                repeat.mating.temp <- population$info$repeat.mating.copy[mating.temp,1]

              } else{
                tt <- nrow(population$info$repeat.mating)
                if(tt>1){
                  mating.temp <- sample(1:tt, prob = population$info$repeat.mating[,2], 1)
                } else{
                  mating.temp <- 1
                }
                repeat.mating.temp <- population$info$repeat.mating[mating.temp,1]
                if(population$info$litter.effect.active){
                  litter.effect <-  as.numeric(population$info$litter.effect.covariance %*% stats::rnorm(population$info$bv.nr, 0, 1))
                  activ_litter <-  population$info$next.litter
                  population$info$next.litter <- activ_litter + 1
                }
              }
            }

            if(length(repeat.mating.fixed)>0){
              repeat.mating.temp <- repeat.mating.fixed[rmf]
              rmf <- rmf +1
            }
            runs <- repeat.mating.temp - 1




            # Determine Sire/dam combination
            check_rel <- FALSE
            accepted <- FALSE
            max_counter <- 0
            max_rel_temp <- max_rel
            while(!check_rel || accepted==FALSE){
              if(selfing.mating==FALSE){
                sex1 <- 1
                sex2 <- 2
                if(fast_mode){
                  number1 <- number1_vector[animal.nr]
                  number2 <- number2_vector[animal.nr]

                  info.father <- best[[sex1]][number1,]
                  info.mother <- best[[sex2]][number2,]
                } else if(same.sex.activ==FALSE){
                  number1 <- availables[[sex1]][sample(activ.selection.size[sex1],1, prob=sample_prob[[sex1]][availables[[sex1]]])]
                  number2 <- availables[[sex2]][sample(activ.selection.size[sex2],1, prob=sample_prob[[sex2]][availables[[sex2]]])]
                  info.father <- best[[sex1]][number1,]
                  info.mother <- best[[sex2]][number2,]
                } else{
                  sex1 <- stats::rbinom(1,1,same.sex.sex) + 1 # ungleichviele tiere erhoeht
                  sex2 <- stats::rbinom(1,1,same.sex.sex) + 1
                  number1 <- availables[[sex1]][sample(activ.selection.size[sex1],1, prob=sample_prob[[sex1]][availables[[sex1]]])]
                  number2 <- availables[[sex2]][sample(activ.selection.size[sex2],1, prob=sample_prob[[sex2]][availables[[sex2]]])]
                  test <- 1
                  while(same.sex.selfing==FALSE && number1==number2 && test < 100){
                    number2 <- availables[[sex2]][sample(1:activ.selection.size[sex2],1, prob=sample_prob[[sex2]][availables[[sex2]]])]
                    test <- test+1
                    if(test==100 && number1==number2){
                      warning("Only one remaining individual in the selected cohorts.")
                    }
                  }


                  info.father <- best[[sex1]][number1,]
                  info.mother <- best[[sex2]][number2,]



                }
              } else{
                sex1 <- sex2 <- stats::rbinom(1,1,selfing.sex)+1
                if(length(availables[[sex1]])==0){
                  if(verbose) cat(paste0("No ", if(sex1==1){"male"} else{"female"}, " individuals left. Use other sex individuals\n"))
                  sex1 <- sex2 <- 3 - sex1
                }
                number1 <- number2 <- availables[[sex1]][sample(activ.selection.size[sex1],1, prob=sample_prob[[sex1]][availables[[sex1]]])]
                info.father <- best[[sex1]][number1,]
                info.mother <- best[[sex2]][number2,]

              }


              accepted <- TRUE
              if(max.mating.pair < Inf){

                if(info.father[2]==info.mother[2]){
                  if(info.father[1]>info.mother[1] || (info.father[1]==info.mother[1] & info.father[3]> info.mother[3])){
                    info_mother_temp <- info.mother
                    info.mother <- info.father
                    info.father <- info_mother_temp

                  }
                }

                if( sum((colSums(t(info_father_list) == info.father) + colSums(t(info_mother_list) == info.mother))==10) >= max.mating.pair){
                  accepted <- FALSE
                }
              }

              if(max_rel_temp!=2){
                check_rel = check.parents(population, info.father, info.mother, max.rel=max_rel_temp)
              } else{
                check_rel <- TRUE
              }

              max_counter <- max_counter + 1
              if(max_counter>=25){
                if(max_rel<2){
                  warning("No remaining possible mating via avoid.mating. Proceed with silbing-mating.\n")
                  max_rel_temp <- max_rel + 1
                  if(max_counter >= 50){
                    check_rel <- TRUE
                    max_rel_temp <- max_rel + 2
                  }
                  if(max_counter >= 250){
                    accepted <- TRUE
                    warning("No remaining possible mating via max.mating.pair. Proceed without limitation\n")
                  }
                } else{
                  if(max_counter >= 1000){
                    accepted <- TRUE
                    warning("A mating was executed more times than allowed in max.mating.pair (check if the total number of allowed mating is sufficient!)\n")
                  }
                }

              }


              if(accepted && max.mating.pair< Inf){
                info_father_list[animal.nr,] <- info.father
                info_mother_list[animal.nr,] <- info.mother
              }



            }
            if(max_check){

              selection.rate[[sex1]][number1] <- selection.rate[[sex1]][number1] + repeat.mating.temp
              selection.rate[[sex2]][number2] <- selection.rate[[sex2]][number2] + repeat.mating.temp

              selection.rate.litter[[sex1]][number1] <- selection.rate.litter[[sex1]][number1] + 1
              selection.rate.litter[[sex2]][number2] <- selection.rate.litter[[sex2]][number2] + 1

              if(selection.rate[[sex1]][number1] >= max.offspring[sex1]){
                activ.selection.size[sex1] <-  activ.selection.size[sex1] -1
                availables[[sex1]] <- availables[[sex1]][availables[[sex1]]!=number1]
              } else if(selection.rate.litter[[sex1]][number1] >= max.litter[sex1]){
                activ.selection.size[sex1] <-  activ.selection.size[sex1] -1
                availables[[sex1]] <- availables[[sex1]][availables[[sex1]]!=number1]
              }
              if(selection.rate[[sex2]][number2] >= max.offspring[sex2]){
                if(sex1!= sex2 || number1 != number2){
                  activ.selection.size[sex2] <-  activ.selection.size[sex2] -1
                }
                availables[[sex2]] <- availables[[sex2]][availables[[sex2]]!=number2]
              } else if(selection.rate.litter[[sex2]][number2] >= max.litter[sex2]){
                if(sex1!= sex2 || number1 != number2){
                  activ.selection.size[sex2] <-  activ.selection.size[sex2] -1
                }
                availables[[sex2]] <- availables[[sex2]][availables[[sex2]]!=number2]
              }
            }


            if(calc_litter){
              mother = population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]]

              prior_pheno <- mother[[27]]



              if(length(fixed.effects.p)>0){
                mother[[28]] <- fixed.effects.p[sample(1:nrow(fixed.effects.p), 1, prob= fixed.effects.freq),]
              } else{
                mother[[28]] <- numeric(0)
              }

              if( length(mother[[23]]) == 0){
                mother[[23]] <- stats::rnorm(population$info$bv.nr,0,1)
              }
              multi_check <- (mother[[15]]==0)

              n.observation_temp <- numeric(population$info$bv.nr)
              n.observation_temp[population$info$repeat.mating.trait] = 1
              mother[[15]] <- mother[[15]] + n.observation_temp

              obsmax <- max(mother[[15]])
              if(length(mother[[24]])==0 ||
                 obsmax > ncol(mother[[24]]) &&
                 obsmax > 0){
                mother[[24]] <- cbind(mother[[24]], matrix(stats::rnorm(obsmax * population$info$bv.nr - length(mother[[24]]),0,1),                                                                               nrow = population$info$bv.nr))
              }

              if(!mother[[25]]){

                if(length(mother[[26]])==0 || sum(n.observation_temp[-mother[[26]]]>0)>0){

                  to_ignore_temp <- which(n.observation_temp==0)
                  activ_bv <- population$info$bv.random.activ

                  temp1234 <- setdiff(population$info$bv.random.activ , to_ignore_temp)

                  temp_out <- calculate.bv(population, info.mother[1], info.mother[2], info.mother[3],
                                           activ_bv, import.position.calculation=import.position.calculation,
                                           decodeOriginsU=decodeOriginsU,
                                           store.effect.freq=store.effect.freq,
                                           bit.storing=bit.storing, nbits=nbits, output_compressed=FALSE,
                                           bv.ignore.traits=to_ignore_temp)

                  if(length(mother[[26]])>0){
                    activ_replace <- !duplicated(c(mother[[26]], activ_bv))[-(1:length(mother[[26]]))]
                  } else{
                    activ_replace <- rep(TRUE, length(activ_bv))
                  }

                  population$breeding[[info.mother[1]]][[6+info.mother[2]]][activ_bv[activ_replace],info.mother[3]] <- temp_out[[1]][activ_replace]

                  if(length(c(mother[[26]], temp1234))>0){
                    mother[[26]] <- c(mother, temp1234)
                  }

                  if(length(mother[[26]])==population$info$bv.nr){
                    mother[[25]] <- TRUE
                  }
                }
              }

              bven = population$info$repeat.mating.trait

              new_pheno <- (sqrt(sigma.e2.rest) * population$info$pheno.correlation %*% mother[[24]][,1:mother[[15]][bven]])[bven,] +
                (sqrt(sigma.e2.perm) * population$info$pheno.correlation %*% mother[[23]])[bven] +
                population$breeding[[info.mother[1]]][[6+info.mother[2]]][bven, info.mother[3]] + sum(population$info$fixed.effects[bven,] * mother[[28]])

              temp1 <- sapply(new_pheno, FUN = population$info$phenotypic.transform.function[[bven]])

              if(length(temp1)>0){
                prior_pheno[[bven]] <- temp1
              }

              population$breeding[[info.mother[1]]][[8+info.mother[2]]][bven, info.mother[3]] <- mean(temp1)


              mother[[27]] <- prior_pheno

              population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]] = mother

              runs <- temp1[length(temp1)] - 1




            }

          }


        }
        if(store.comp.times.generation.temp){
          tack <- as.numeric(Sys.time())
          pre_stuff <- pre_stuff + tack -tick
        }

        if(generation.cores>1){
          fixed.breeding_parallel[animal.nr,] <- c(info.father[1:3], info.mother[1:3])
        } else{


        father <- population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]]
        mother <- population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]]

        if(dh.mating){
          if(stats::rbinom(1,1,dh.sex)==0){
            child1 <- breeding.intern.activ(info.father, father, population,
                                            mutation.rate, remutation.rate, recombination.rate,
                                            recom.f.indicator, duplication.rate, duplication.length,
                                            duplication.recombination, delete.same.origin=delete.same.origin,
                                            gene.editing=(gene.editing.offspring*gene.editing.offspring.sex[1]), nr.edits= nr.edits,
                                            gen.architecture=gen.architecture.m,
                                            decodeOriginsU=decodeOriginsU,
                                            recombination.function = recombination.function)
            child2 <- child1
          } else{
            child2 <- breeding.intern.activ(info.mother, mother, population,
                                            mutation.rate, remutation.rate, recombination.rate,
                                            recom.f.indicator, duplication.rate, duplication.length,
                                            duplication.recombination, delete.same.origin=delete.same.origin,
                                            gene.editing=(gene.editing.offspring * gene.editing.offspring.sex[1]) , nr.edits= nr.edits,
                                            gen.architecture=gen.architecture.f,
                                            decodeOriginsU=decodeOriginsU,
                                            recombination.function = recombination.function)
            child1 <- child2
          }
        } else if(copy.individual){
          info.mother <- info.father
          child1 <- list(father[[1]], father[[3]], father[[5]], father[[7]], father[[11]], 0, if(length(father)>19){father[[19]]} else{0})
          child2 <- list(father[[2]], father[[4]], father[[6]], father[[8]], father[[12]], 0, if(length(father)>19){father[[20]]} else{0})
        } else{
          child1 <- breeding.intern.activ(info.father, father, population,
                                          mutation.rate, remutation.rate, recombination.rate,
                                          recom.f.indicator, duplication.rate, duplication.length,
                                          duplication.recombination, delete.same.origin=delete.same.origin,
                                          gene.editing=(gene.editing.offspring*gene.editing.offspring.sex[1]), nr.edits= nr.edits,
                                          gen.architecture=gen.architecture.m,
                                          decodeOriginsU=decodeOriginsU,
                                          recombination.function = recombination.function)

          child2 <- breeding.intern.activ(info.mother, mother, population,
                                          mutation.rate, remutation.rate, recombination.rate,
                                          recom.f.indicator, duplication.rate, duplication.length,
                                          duplication.recombination, delete.same.origin=delete.same.origin,
                                          gene.editing=(gene.editing.offspring * gene.editing.offspring.sex[1]) , nr.edits= nr.edits,
                                          gen.architecture=gen.architecture.f,
                                          decodeOriginsU=decodeOriginsU,
                                          recombination.function = recombination.function)
        }


        # Put together offspring

        child_temp <- list(child1[[1]], child2[[1]], child1[[2]], child2[[2]], child1[[3]], child2[[3]], child1[[4]], child2[[4]])


        if(copy.individual){
          population$breeding[[current.gen+1]][[sex+22]][current.size[sex]] <- population$breeding[[info.father[1]]][[info.father[2]+22]][info.father[3]]
        }
        population$info$size[current.gen+1 ,sex] <- population$info$size[current.gen+1,sex] + 1

        if(is.vector(child1[[5]])){
          child_temp[[11]] <- t(as.matrix(child1[[5]]))
        } else{
          child_temp[[11]] <- child1[[5]]
        }
        if(is.vector(child2[[5]])){
          child_temp[[12]] <- t(as.matrix(child2[[5]]))
        } else{
          child_temp[[12]] <- child2[[5]]
        }
        if(save.recombination.history && current.gen==1){
          if(length(child1[[6]][-c(1,length(child1[[6]]))])>0){
            child_temp[[13]] <- cbind(current.gen, child1[[6]][-c(1,length(child1[[6]]))], deparse.level = 0)
          } else{
            child_temp[[13]] <- cbind(0,0, deparse.level = 0)
          }
          if(length( child2[[6]][-c(1,length(child2[[6]]))])>0){
            child_temp[[14]] <- cbind(current.gen, child2[[6]][-c(1,length(child2[[6]]))], deparse.level = 0)
          } else{
            child_temp[[14]] <- cbind(0,0, deparse.level = 0)
          }

        } else if(save.recombination.history && current.gen>1){
          if(length(child1[[6]][-c(1,length(child1[[6]]))])>0){
            child_temp[[13]] <- rbind(population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[13]], population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[14]], cbind(current.gen, child1[[6]][-c(1,length(child1[[6]]))], deparse.level = 0), deparse.level = 0)
          } else{
            child_temp[[13]] <- rbind(population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[13]], population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[14]], deparse.level = 0)

          }
          if(length( child2[[6]][-c(1,length(child2[[6]]))])>0){
            child_temp[[14]] <- rbind(population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]][[13]], population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]][[14]], cbind(current.gen, child2[[6]][-c(1,length(child2[[6]]))]))
          } else{
            child_temp[[14]] <- rbind(population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]][[13]], population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]][[14]], deparse.level = 0)

          }

        } else{
          #population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[13]] <- "test"
        }

        if(share.phenotyped==1){
          is.obs <- TRUE
        } else if(share.phenotyped==0){
          is.obs <- FALSE
        } else{
          is.obs <- stats::rbinom(1,1, share.phenotyped)==1
        }

        if(length(phenotyping.class)>0 & !(new.class[sex] %in% phenotyping.class)){
          is.obs <- FALSE
        }

        if(phenotyping.child=="obs"){
          child_temp[[15]] <- n.observation * is.obs
        } else if(phenotyping.child=="addobs"){
          child_temp[[15]] <-  colMeans(rbind(population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]][[15]],
                                              population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[15]]))
          switch <- (child_temp[[15]]< (n.observation*is.obs))
          child_temp[[15]][switch] <- (n.observation*is.obs)[switch]

        } else{
          child_temp[[15]] <- rep(0L, population$info$bv.nr)
        }

        if(copy.individual && !copy.individual.keep.pheno){
          child_temp[[15]] <- rep(0L, population$info$bv.nr)
        }

        child_temp[[27]] <- list()
        if(length(fixed.effects.p)==0){
          child_temp[[28]] <- numeric(0)
        } else{
          child_temp[[28]] <- fixed.effects.p[sample(1:nrow(fixed.effects.p), 1, prob= fixed.effects.freq),]
        }

        if(population$info$litter.effect.active){
          child_temp[[29]] <- litter.effect
          child_temp[[31]] <- activ_litter
        } else{
          #child_temp[[29]] <- rep(0L, population$info$bv.nr)
          #child_temp[[31]] <- 0L
        }

        if(population$info$pen.effect.active){
          child_temp[[30]] <- rep(0L, population$info$bv.nr)
          child_temp[[32]] <- 0L
        } else{
          #child_temp[[30]] <- rep(0L, population$info$bv.nr)
          #child_temp[[32]] <- 0L
        }



        child_temp[[33]] <- "placeholder"

        if(copy.individual){
          child_temp[[28]] <- population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[28]]
          child_temp[[16]] <- population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[16]]
          if(length(population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[22]])>0){
            child_temp[[22]] <- population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[22]]
          }

          if(added.genotyped>0 && child_temp[[16]]==0){
            if(stats::rbinom(1,1,added.genotyped)==1){
              child_temp[[16]] <- 1
              child_temp[[22]] <- c(child_temp[[22]], genotyped.array)
            }
          }

          first_copy <- population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[21]]
          new_copy <- rbind(first_copy,
                            c(current.gen+1L, sex, current.size[sex]), deparse.level = 0)
          storage.mode(new_copy) <- "integer"
          if(nrow(new_copy)>1){
            for(index7 in 1:(nrow(new_copy)-1)){
              if(length(population$breeding[[new_copy[index7,1]]])>1 && ((length(population$breeding[[new_copy[index7,1]]][[new_copy[index7,2]]]) > 1) || (population$breeding[[new_copy[index7,1]]][[new_copy[index7,2]]]!="removed"))){
                population$breeding[[new_copy[index7,1]]][[new_copy[index7,2]]][[new_copy[index7,3]]][[21]] <- new_copy
              }

            }
          }
          child_temp[[21]] <- new_copy

        } else{
          if(share.genotyped==1 || (share.genotyped!=0 && stats::rbinom(1,1,share.genotyped)==1)){
            child_temp[[16]] <- 1
            child_temp[[22]] <- genotyped.array
          } else{
            child_temp[[16]] <- 0
          }

          child_temp[[21]] <- cbind(current.gen+1, sex, current.size[sex], deparse.level = 0)
          storage.mode(child_temp[[21]]) <- "integer"
        }
        if(length(child1[[7]])>0){
          child_temp[[19]] <- child1[[7]]
        }
        if(length(child2[[7]])>0){
          child_temp[[20]] <- child2[[7]]
        }

        if(copy.individual){
          child_temp[[25]] <- population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[25]]

          if(length(population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[26]])>0){
            child_temp[[26]] <- population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[26]]
          }

        } else{
          child_temp[[25]] <- length(bv.ignore.traits)==0
          if(length(temp123)>0){
            child_temp[[26]] <- temp123
          }

        }



        population$breeding[[current.gen+1]][[sex]][[current.size[sex]]] <- child_temp
        if(store.comp.times.generation){
          tock <- as.numeric(Sys.time())
          generation_stuff <- generation_stuff + tock -tack
        }


        # calculate underlying true genomic value
        ## revisit computations for non-QTL based effects for efficiently

        enter <- FALSE
        if(population$info$bve){

          activ_bv <- population$info$bv.random.activ


          if(length(bv.ignore.traits)!=length(activ_bv)){
            enter <- TRUE
            if(length(activ_bv)>0){
              if(!copy.individual || store.effect.freq){
                temp_out <- calculate.bv(population, current.gen+1, sex, current.size[sex], activ_bv, import.position.calculation=import.position.calculation, decodeOriginsU=decodeOriginsU, store.effect.freq=store.effect.freq, bit.storing=bit.storing, nbits=nbits, output_compressed=FALSE, bv.ignore.traits=bv.ignore.traits)
                new.bv[activ_bv] <- temp_out[[1]]

                if(store.effect.freq){
                  if(length(population$info$store.effect.freq) < (current.gen+1) || length(population$info$store.effect.freq[[current.gen+1]])==0){
                    colnames(temp_out[[2]]) <- c("Homo0", "Hetero", "Homo1")
                    rownames(temp_out[[2]]) <- population$info$snp.name[population$info$effect.p]
                    population$info$store.effect.freq[[current.gen+1]] <- temp_out[[2]]
                  } else{
                    population$info$store.effect.freq[[current.gen+1]] <- population$info$store.effect.freq[[current.gen+1]] + temp_out[[2]]
                  }
                }

              } else{
                activ_indi <- info.father
                new.bv[activ_bv] <- population$breeding[[activ_indi[1]]][[activ_indi[2]+6]][activ_bv, activ_indi[3]]
              }


            }



            if(population$info$bv.calc > 0  && population$info$bv.random[population$info$bv.calc] && prod(population$info$is.combi[population$info$bv.calc:population$info$bv.nr])==0){

              #Means passt (Korrelation exakt wie gewuenscht)
              means <- 0.5*(population$breeding[[info.father[1]]][[6+info.father[2]]][population$info$bv.calc:population$info$bv.nr,info.father[3]] + population$breeding[[info.mother[1]]][[6+info.mother[2]]][population$info$bv.calc:population$info$bv.nr,info.mother[3]])

              # Berechnung i (0-0.5)
              varp <- kinship.emp(list(population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]], population$breeding[[info.mother[1]]][[info.mother[2]]][[info.mother[3]]]))
              varp <- (2 - (2* (varp[1,1]-0.5) + 2 * (varp[2,2]- 0.5)))/4

              if(FALSE){
                new.bv[population$info$bv.calc:population$info$bv.nr] <- stats::rnorm(1, mean=means, sd= sqrt(var*population$info$bv.random.variance[bven]))
              } else{
                if(population$info$bv.calc==1){
                  population$info$current.bv.random.variance <- varp * population$info$bv.random.variance
                  bv.var <- diag.mobps(sqrt(population$info$current.bv.random.variance)) %*%population$info$current.bv.correlation %*% diag.mobps(sqrt(population$info$current.bv.random.variance))
                  single.mean <- means
                } else{
                  #population$info$current.bv.random.variance <- c(population$info$bv.random.variance[1:(population$info$bv.calc-1)], population$info$bv.random.variance[population$info$bv.calc:population$info$bv.nr])
                  population$info$current.bv.random.variance <- c(population$info$bv.random.variance[1:(population$info$bv.calc-1)],varp * population$info$bv.random.variance[population$info$bv.calc:population$info$bv.nr])

                  AA <- diag.mobps(sqrt(population$info$current.bv.random.variance)[1:(population$info$bv.calc-1)]) %*% population$info$current.bv.correlation[1:(population$info$bv.calc-1), 1:(population$info$bv.calc-1)]%*% diag.mobps(sqrt(population$info$current.bv.random.variance)[(1:(population$info$bv.calc-1))])
                  BB <- diag.mobps(sqrt(population$info$current.bv.random.variance)[1:(population$info$bv.calc-1)]) %*%population$info$current.bv.correlation[1:(population$info$bv.calc-1), -(1:(population$info$bv.calc-1))]%*% diag.mobps(sqrt(population$info$current.bv.random.variance)[-(1:(population$info$bv.calc-1))])
                  CC <- diag.mobps(sqrt(population$info$current.bv.random.variance)[-(1:(population$info$bv.calc-1))]) %*%population$info$current.bv.correlation[-(1:(population$info$bv.calc-1)), -(1:(population$info$bv.calc-1))] %*% diag.mobps(sqrt(population$info$current.bv.random.variance)[-(1:(population$info$bv.calc-1))])
                  if (requireNamespace("MASS", quietly = TRUE)) {
                    bv.var <- CC - t(BB) %*% MASS::ginv(AA) %*% BB
                    single.mean <- means + t(BB) %*% MASS::ginv(AA) %*% ( new.bv[1:(population$info$bv.calc-1)]-mu1[1:(population$info$bv.calc-1)])
                  } else{
                    bv.var <- CC - t(BB) %*% solve(AA) %*% BB
                    single.mean <- means + t(BB) %*% solve(AA) %*% ( new.bv[1:(population$info$bv.calc-1)]-mu1[1:(population$info$bv.calc-1)])
                  }


                }

                bv.var_chol <- t(chol(bv.var))
                population$info$bv.correlation_col <- bv.var_chol

                #            random_part <- varp * bv.var_chol %*% stats:rnorm(population$info$bv.nr - population$info$bv.calc +1, 0,1 )
                random_part <- bv.var_chol %*% stats::rnorm(population$info$bv.nr - population$info$bv.calc +1, 0,1 )
                means_part <- single.mean


                new.bv[population$info$bv.calc:population$info$bv.nr] <-   single.mean  + random_part
                #new.bv[population$info$bv.calc:population$info$bv.nr] <- bv.var_chol %*% stats::rnorm(population$info$bv.nr - population$info$bv.calc +1, 0,1 ) + single.mean
              }


            }
            if(sum(population$info$is.combi)>0){
              combis <- which(population$info$is.combi)
              for(combi in combis){
                new.bv[combi] <- sum(new.bv * population$info$combi.weights[[combi]])
              }
            }

            if(copy.individual && copy.individual.keep.bve){
              new.bve <- population$breeding[[info.father[1]]][[2+info.father[2]]][,info.father[3]]
              new.reli <- population$breeding[[info.father[1]]][[18+info.father[2]]][,info.father[3]]
            }

            if(copy.individual && length(population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[23]])>0){
              population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[23]] <-
                population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[23]]
            }



            if(phenotyping.child=="mean" || phenotyping.child=="addobs"){
              if(copy.individual){
                if(length(population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[24]])>0){
                  population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[24]] <-
                    population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[24]]
                }
              } else{
                new.bv_approx <- 0.5 * (population$breeding[[info.father[1]]][[8+info.father[2]]][,info.father[3]] + population$breeding[[info.mother[1]]][[8+info.mother[2]]][,info.mother[3]])
              }

            }




            # phenotyping



            if(phenotyping.child=="obs" || phenotyping.child=="addobs" || copy.individual){
              if(sum(n.observation)>0 || (copy.individual && copy.individual.keep.pheno)){

                if( length(population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[23]]) == 0){
                  population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[23]] <- stats::rnorm(population$info$bv.nr,0,1)
                }

                obsmax <- max(population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[15]])

                if(length(population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[24]])==0 ||
                   obsmax > ncol(population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[24]]) &&
                   obsmax > 0){

                  population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[24]] <- cbind(population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[24]],
                                                                                                  matrix(stats::rnorm(obsmax * population$info$bv.nr - length(population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[24]]),0,1),
                                                                                                         nrow = population$info$bv.nr))
                }



                prior_pheno <- population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[27]]
                if(!copy.individual){

                  if(length(fixed.effects.p)==0){
                    population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[28]] <- numeric(0)
                  } else{
                    population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[28]] <- fixed.effects.p[sample(1:nrow(fixed.effects.p), 1, prob= fixed.effects.freq),]
                  }

                }


                if(sum(population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[15]])>0){
                  for(bven in setdiff((1:population$info$bv.nr)[population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[15]]>0], activ.trafo)){
                    if(population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[15]][bven]>=1){

                      replace <- TRUE
                      if(copy.individual && n.observation[bven]==0){

                        if(length(population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[27]])>=bven){
                          temp1 <- population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[27]][[bven]]
                          new.bv_approx[bven] <- mean(temp1)
                        } else if(population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[15]][bven]>0){

                          temp1 <- (sqrt(sigma.e2.rest) * (population$info$pheno.correlation %*% population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[24]][,1:population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[15]][bven]]) +
                                      as.numeric(sqrt(sigma.e2.perm) * population$info$pheno.correlation %*% population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[23]]) +
                                      new.bv)[bven,] + sum(population$info$fixed.effects[bven,] * population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[28]])


                          new.bv_approx[bven] <- mean(temp1)

                        } else {
                          replace <- FALSE
                        }

                      } else{

                        temp1 <- (sqrt(sigma.e2.rest) * (population$info$pheno.correlation %*% population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[24]][,1:population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[15]][bven]]) +
                                    as.numeric(sqrt(sigma.e2.perm) * population$info$pheno.correlation %*% population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[23]]) +
                                    new.bv)[bven,] + sum(population$info$fixed.effects[bven,] * population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[28]])


                        new.bv_approx[bven] <- mean(temp1)
                      }

                      if(length(temp1)>0 & replace){
                        prior_pheno[[bven]] <- temp1
                      }


                    }

                  }
                }


                if(sum(population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[15]])>0){
                  for(bven in intersect((1:population$info$bv.nr)[population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[15]]>0], activ.trafo)){
                    if(population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[15]][bven]>=1){
                      new_pheno <- (sqrt(sigma.e2.rest) * population$info$pheno.correlation %*% population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[24]][,1:population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[15]][bven]])[bven,] +
                        (sqrt(sigma.e2.perm) * population$info$pheno.correlation %*% population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[23]])[bven] +
                        new.bv[bven] + sum(population$info$fixed.effects[bven,] * population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[28]])

                      if(copy.individual && n.observation[bven]==0 && length(population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[27]])>bven){

                        temp1 <- population$breeding[[info.father[1]]][[info.father[2]]][[info.father[3]]][[27]][[bven]]
                        new.bv_approx[bven] <- mean(temp1)

                      } else{
                        temp1 <- sapply(new_pheno, FUN = population$info$phenotypic.transform.function[[bven]])
                        new.bv_approx[bven] <- mean(temp1)
                      }

                      if(length(temp1)>0){
                        prior_pheno[[bven]] <- temp1
                      }

                    }
                  }
                }
                population$breeding[[current.gen+1]][[sex]][[current.size[sex]]][[27]] <- prior_pheno

              }


            }



            if(phenotyping.child=="keep" && copy.individual){
              new.bv_approx <- population$breeding[[info.father[1]]][[info.father[2]+8]][,info.father[3]]
            }
          }



        }



        if(enter){
          population$breeding[[current.gen+1]][[2+sex]][,current.size[sex]] <- new.bve
          population$breeding[[current.gen+1]][[6+sex]][,current.size[sex]] <- new.bv
          population$breeding[[current.gen+1]][[8+sex]][,current.size[sex]] <- new.bv_approx
          population$breeding[[current.gen+1]][[18+sex]][,current.size[sex]] <- new.reli
        }


        if(copy.individual){
          population$breeding[[current.gen+1]][[14+sex]][current.size[sex]] <- population$breeding[[info.father[1]]][[14+info.father[2]]][info.father[3]]
          population$breeding[[current.gen+1]][[20+sex]][,current.size[sex]] <- population$breeding[[info.father[1]]][[20+info.father[2]]][,info.father[3]]
          population$breeding[[current.gen+1]][[26+sex]][,current.size[sex]] <- population$breeding[[info.father[1]]][[26+info.father[2]]][,info.father[3]]
          population$breeding[[current.gen+1]][[28+sex]][,current.size[sex]] <- population$breeding[[info.father[1]]][[28+info.father[2]]][,info.father[3]]
        }


        if(store.comp.times.generation){
          tock2 <- as.numeric(Sys.time())
          bv_stuff <- bv_stuff+tock2-tock
        }
        if(display.progress & verbose & generation.cores==1 & (breeding.size.total < 100 || animal.nr%%pb_temp==0)){
          utils::setTxtProgressBar(pb, animal.nr)
        }

        current.size[sex] <- current.size[sex] +1
        }
      }

      if(sum(breeding.size.litter>0)>0){

        breeding.size = current.size - start.size
        breeding.size.total = sum(breeding.size)
        if(verbose && generation.cores==1){
          cat(paste0("Successfully generated ", breeding.size.total, " individuals from ", litter_counter-1, " litters.\n"))
        }

        for(sex in 1:2){

          n_gen = length(population$breeding[[current.gen+1]][[4+sex]])
          if( n_gen >= (current.size[sex])){
            remove1 = current.size[sex]:n_gen
            population$breeding[[current.gen+1]][[2 + sex]] =  population$breeding[[current.gen+1]][[2 + sex]][,-remove1,drop=FALSE]
            population$breeding[[current.gen+1]][[4 + sex]] =  population$breeding[[current.gen+1]][[4 + sex]][-remove1]
            population$breeding[[current.gen+1]][[6 + sex]] =  population$breeding[[current.gen+1]][[6 + sex]][,-remove1,drop=FALSE]
            population$breeding[[current.gen+1]][[8 + sex]] =  population$breeding[[current.gen+1]][[8 + sex]][,-remove1,drop=FALSE]
            population$breeding[[current.gen+1]][[10 + sex]] =  population$breeding[[current.gen+1]][[10 + sex]][-remove1]
            population$breeding[[current.gen+1]][[12 + sex]] =  population$breeding[[current.gen+1]][[12 + sex]][-remove1]
            population$breeding[[current.gen+1]][[14 + sex]] =  population$breeding[[current.gen+1]][[14 + sex]][-remove1]
            population$breeding[[current.gen+1]][[16 + sex]] =  population$breeding[[current.gen+1]][[16 + sex]][-remove1]
            population$breeding[[current.gen+1]][[18 + sex]] =  population$breeding[[current.gen+1]][[18 + sex]][,-remove1,drop=FALSE]
            population$breeding[[current.gen+1]][[20 + sex]] =  population$breeding[[current.gen+1]][[20 + sex]][,-remove1,drop=FALSE]
            population$breeding[[current.gen+1]][[22 + sex]] =  population$breeding[[current.gen+1]][[22 + sex]][-remove1]
            population$breeding[[current.gen+1]][[24 + sex]] =  population$breeding[[current.gen+1]][[24 + sex]][-remove1]
            population$breeding[[current.gen+1]][[26 + sex]] =  population$breeding[[current.gen+1]][[26 + sex]][,-remove1,drop=FALSE]
            population$breeding[[current.gen+1]][[28 + sex]] =  population$breeding[[current.gen+1]][[28 + sex]][,-remove1,drop=FALSE]
            population$breeding[[current.gen+1]][[30 + sex]] =  population$breeding[[current.gen+1]][[30 + sex]][-remove1]
            population$breeding[[current.gen+1]][[32 + sex]] =  population$breeding[[current.gen+1]][[32 + sex]][-remove1]
            population$breeding[[current.gen+1]][[34 + sex]] =  population$breeding[[current.gen+1]][[34 + sex]][-remove1]
            population$breeding[[current.gen+1]][[36 + sex]] =  population$breeding[[current.gen+1]][[36 + sex]][-remove1]
          }

        }
      } else{
        if(verbose && generation.cores==1){
          cat(paste0("Successfully generated ", breeding.size.total, " individuals.\n"))
        }
      }



      if(store.comp.times.generation.temp2){
        tack <- as.numeric(Sys.time())
        pre_stuff <- pre_stuff + tack -tick
      }
      if(generation.cores > 1){



        if(ncol(fixed.breeding_parallel)<7){
          fixed.breeding_parallel <- cbind(fixed.breeding_parallel, sex.animal-1)
        }


        parallel_batch <- list()

        size_batch <- ceiling(breeding.size.total / generation.cores)
        for(batch in 1:generation.cores){
          parallel_batch[[batch]] <- fixed.breeding_parallel[(1+(batch-1)*size_batch):min(breeding.size.total,(batch*size_batch)) ,]
        }



        if(Sys.info()[['sysname']]=="Windows"){
          # Windows Parallel
          doParallel::registerDoParallel(cores=generation.cores)
          p1 <- foreach::foreach(indexb=1:generation.cores,
                                 .packages="MoBPS") %dopar% {

                                   pop1 <- breeding.diploid(population_parallel, fixed.breeding = parallel_batch[[indexb]],

                                                            mutation.rate = mutation.rate,
                                                            remutation.rate = remutation.rate,
                                                            recombination.rate = recombination.rate,
                                                            recom.f.indicator = recom.f.indicator,
                                                            duplication.rate = duplication.rate,
                                                            duplication.length = duplication.length,
                                                            duplication.recombination = duplication.recombination, delete.same.origin=delete.same.origin,
                                                            gene.editing.offspring = gene.editing.offspring,
                                                            gene.editing.best = gene.editing.best,
                                                            gene.editing.offspring.sex = gene.editing.offspring.sex,
                                                            gene.editing.best.sex = gene.editing.best.sex,
                                                            nr.edits= nr.edits,
                                                            gen.architecture.m=gen.architecture.m,
                                                            gen.architecture.f=gen.architecture.f,
                                                            phenotyping.child = phenotyping.child,
                                                            n.observation = n.observation,
                                                            copy.individual = copy.individual,
                                                            copy.individual.keep.pheno = copy.individual.keep.pheno,
                                                            copy.individual.keep.bve = copy.individual.keep.bve,
                                                            added.genotyped = added.genotyped,
                                                            share.genotyped = share.genotyped,
                                                            bv.ignore.traits = bv.ignore.traits,
                                                            heritability = heritability,
                                                            sigma.e =sigma.e,
                                                            dh.mating = dh.mating,
                                                            dh.sex = dh.sex,
                                                            same.sex.activ = same.sex.activ,
                                                            same.sex.sex = same.sex.sex,
                                                            same.sex.selfing = same.sex.selfing,
                                                            selfing.mating = selfing.mating,
                                                            selfing.sex = selfing.sex,
                                                            verbose=FALSE,
                                                            store.comp.times.generation = FALSE,

                                                             parallel.internal=TRUE)
                                 }

          doParallel::stopImplicitCluster()
        } else{

          p1 <- parallel::mclapply(1:generation.cores, function(x) breeding.diploid(population_parallel, fixed.breeding = parallel_batch[[x]],
                                                                                    mutation.rate = mutation.rate,
                                                                                    remutation.rate = remutation.rate,
                                                                                    recombination.rate = recombination.rate,
                                                                                    recom.f.indicator = recom.f.indicator,
                                                                                    duplication.rate = duplication.rate,
                                                                                    duplication.length = duplication.length,
                                                                                    duplication.recombination = duplication.recombination, delete.same.origin=delete.same.origin,
                                                                                    gene.editing.offspring = gene.editing.offspring,
                                                                                    gene.editing.best = gene.editing.best,
                                                                                    gene.editing.offspring.sex = gene.editing.offspring.sex,
                                                                                    gene.editing.best.sex = gene.editing.best.sex,
                                                                                    nr.edits= nr.edits,
                                                                                    gen.architecture.m=gen.architecture.m,
                                                                                    gen.architecture.f=gen.architecture.f,
                                                                                    phenotyping.child = phenotyping.child,
                                                                                    n.observation = n.observation,
                                                                                    copy.individual = copy.individual,
                                                                                    copy.individual.keep.pheno = copy.individual.keep.pheno,
                                                                                    copy.individual.keep.bve = copy.individual.keep.bve,
                                                                                    added.genotyped = added.genotyped,
                                                                                    share.genotyped = share.genotyped,
                                                                                    bv.ignore.traits = bv.ignore.traits,
                                                                                    heritability = heritability,
                                                                                    sigma.e =sigma.e,
                                                                                    dh.mating = dh.mating,
                                                                                    dh.sex = dh.sex,
                                                                                    same.sex.activ = same.sex.activ,
                                                                                    same.sex.sex = same.sex.sex,
                                                                                    same.sex.selfing = same.sex.selfing,
                                                                                    selfing.mating = selfing.mating,
                                                                                    selfing.sex = selfing.sex,
                                                                                    parallel.internal=TRUE,
                                                                                    store.comp.times.generation = FALSE,
                                                                                    verbose=FALSE),
                                   mc.cores=generation.cores)

        }

        if(store.comp.times.generation){
          tock <- as.numeric(Sys.time())
          parallel_gen  <- parallel_gen  + tock -tack
        }
        activ_male <- length(population$breeding[[current.gen+1]][[1]])
        activ_female <- length(population$breeding[[current.gen+1]][[2]])
        size_sex <- c(activ_male, activ_female)

        for(index in 1:generation.cores){
          for(sex in 1:2){

            adds <- c(length(p1[[index]][[sex]]))

            if(length(p1[[index]][[sex]])>0){

              for(tab in c(2,6,8,18,20,26,28)){
                population$breeding[[current.gen+1]][[tab+sex]][,1:adds+size_sex[sex]] <- p1[[index]][[tab+sex]]
              }
              for(vec in c(4,10,12,14,16,22,24,30,32,34,36)){
                population$breeding[[current.gen+1]][[tab+sex]][,1:adds+size_sex[sex]] <- p1[[index]][[tab+sex]]
              }
              ## 32, 34 are placeholder for pen/litter number! Need to be implemented!!

              size_sex[sex] <- size_sex[sex] + adds

              population$breeding[[current.gen+1]][[sex]] <- c(population$breeding[[current.gen+1]][[sex]], p1[[index]][[sex]])

            }


          }


        }
        population$info$size[current.gen+1,] <- size_sex
        current.size <- size_sex + 1

        if(store.comp.times.generation){
          tock2 <- as.numeric(Sys.time())
          parallel_joint  <- parallel_joint +tock2-tock
        }

        if(verbose && sum(breeding.size.litter>0)>0){
          cat(paste0("Successfully generated ", breeding.size.total, " individuals from ", litter_counter-1, " litters.\n"))
        } else if(verbose){
          cat(paste0("Successfully generated ", breeding.size.total, " individuals.\n"))
        }

      }


      if(display.progress & verbose & generation.cores==1){
        close(pb)
      }

    }

  }


  #######################################################################
  ########################## Final housekeeping #########################
  #######################################################################
  {
    delete.haplotypes <- delete.haplotypes
    if(length(delete.haplotypes)>0){
      for(gen_r in delete.haplotypes){
        for(sex in 1:2){
          n.animals <- length(population$breeding[[gen_r]][[sex]])
          if(n.animals>0){
            for(index2 in 1:n.animals){
              population$breeding[[gen_r]][[sex]][[index2]][[9]] <- "removed"
              population$breeding[[gen_r]][[sex]][[index2]][[10]] <- "removed"
            }
          }
        }
      }
    }
    delete.individuals <- delete.individuals
    if(length(delete.individuals)>0){
      for(gen_r in delete.individuals){
        for(sex in delete.sex){
          n.animals <- length(population$breeding[[gen_r]][[sex]])
          if(n.animals>0){
            population$breeding[[gen_r]][[sex]]<- "removed"

          }
        }
      }
    }

    delete.gen <- delete.gen
    if(length(delete.gen)>0){
      for(gen_r in delete.gen){
        population$breeding[[gen_r]]<- "removed"
      }
    }



    # Add delete.gen for complete removal of a generation including [[3]] - [[X]] inputs

    if(store.breeding.totals){
      cur <- length(population$info$breeding.totals) + 1L
      population$info$breeding.totals[[cur]] <- list()
      population$info$breeding.totals[[cur]][[1]] <- current.gen +1L
      population$info$breeding.totals[[cur]][[2]] <- breeding.size
      population$info$breeding.totals[[cur]][[3]] <- best[[1]]
      population$info$breeding.totals[[cur]][[4]] <- best[[2]]
      population$info$breeding.totals[[cur]][[5]] <- selection.rate[[1]]
      population$info$breeding.totals[[cur]][[6]] <- selection.rate[[2]]
      population$info$breeding.totals[[cur]][[7]] <- chosen.animals.list
    }

    if(store.bve.data && bve){
      cur <- length(population$info$bve.data) + 1
      population$info$bve.data[[cur]] <- list()
      population$info$bve.data[[cur]][[1]] <- current.gen
      population$info$bve.data[[cur]][[2]] <- sigma.e2
      population$info$bve.data[[cur]][[3]] <- sigma.e2.hat
      population$info$bve.data[[cur]][[4]] <- sigma.g^2
      population$info$bve.data[[cur]][[5]] <- sigma.a2.hat
      population$info$bve.data[[cur]][[6]] <- bve.database
      population$info$bve.data[[cur]][[7]] <- phenotyping
      population$info$bve.data[[cur]][[8]] <- y_real
      population$info$bve.data[[cur]][[9]] <- y_hat
      population$info$bve.data[[cur]][[10]] <- y
    }



    if(length(population$info$last.sigma.e.value)< population$info$bv.nr || prod(population$info$last.sigma.e.value==sigma.e)!=1){
      population$info$last.sigma.e.heritability <- heritability
      population$info$last.sigma.e.database <- sigma.e.database
      population$info$last.sigma.e.value <- sigma.e
    }


    if(!repeat.mating.overwrite){
      population$info$repeat.mating <- repeat.mating.store
      population$info$repeat.mating.copy <- repeat.mating.copy.store
    }

    if(store.comp.times){
      comp.times[7] <- as.numeric(Sys.time())
      comp.times <- c(comp.times[-1] - comp.times[-length(comp.times)], comp.times[length(comp.times)]-comp.times[1])

      comp.times[comp.times<0] <- 0
      comp.times[comp.times>10e6] <- 0

      population$info$comp.times.general<- round(rbind(population$info$comp.times.general, comp.times, deparse.level = 0), digits=4)
      if(nrow(population$info$comp.times.general)==1){
        colnames(population$info$comp.times.general) <- c("preparation", "new real BV", "phenotypes", "BVE","selection","generate new individuals","total")
      }
    }
    if(store.comp.times.bve){
      comp.times.bve <- c(comp.times.bve[-1] - comp.times.bve[-length(comp.times.bve)], zcalc, z_chol, z_uhat, z_ped, z_h, comp.times.bve[length(comp.times.bve)]-comp.times.bve[1])
      comp.times.bve[comp.times.bve<0] <- 0
      comp.times.bve[comp.times.bve>10e6] <- 0

      population$info$comp.times.bve <- round(rbind(population$info$comp.times.bve, comp.times.bve, deparse.level = 0), digits=4)
      if(nrow(population$info$comp.times.bve)==1){
        colnames(population$info$comp.times.bve) <- c("y_z_import", "A genomic", "solveMixed","Gwas_stuff", "Derive Z", "A inversion", "rrBlup", "A-Pedigree","SingleStep H", "Total")
      }
    }
    if(store.comp.times.generation){
      comp.times.generation <- c(pre_stuff, generation_stuff, bv_stuff, parallel_gen, parallel_joint, sum(pre_stuff, generation_stuff, bv_stuff, parallel_gen, parallel_joint))
      population$info$comp.times.generation <- round(rbind(population$info$comp.times.generation, comp.times.generation, deparse.level = 0), digits=4)
      if(nrow(population$info$comp.times.generation)==1){
        colnames(population$info$comp.times.generation) <- c("Preparation", "Generation", "BV-Calculation", "Parallel Generation", "Parallel Merging", "Total")
      }
    }

    if(length(name.cohort)==1 && breeding.size[1] > 0 && breeding.size[2]>0){
      name.cohort = paste0(name.cohort, c("_M", "_F"))
      if(verbose) cat("Added _M, _F to cohort names!\n")
    }

    if(length(name.cohort)>0){
      if(breeding.size[1] > 0 && breeding.size[2]>0){
        population$info$cohorts <- rbind(population$info$cohorts, c(name.cohort[1], current.gen+1, breeding.size[1],0, new.class[1], (current.size-breeding.size)[1], 0,
                                                                    time.point, creating.type),
                                         c(name.cohort[2], current.gen+1, 0, breeding.size[2], new.class[2], 0, (current.size-breeding.size)[2],time.point, creating.type))

        rownames(population$info$cohorts)[(nrow(population$info$cohorts)-1):nrow(population$info$cohorts)] <- name.cohort

        if(verbose){
          posi <- get.database(population, cohorts = name.cohort[1])
          cat(paste0("Successfully generated cohort: ",  name.cohort[1], "\n",
                     "Database position: ", posi[1], " (gen), ", posi[2], " (sex), ", posi[3], " (first), ", posi[4], " (last).\n" ))
          posi <- get.database(population, cohorts = name.cohort[2])
          cat(paste0("Successfully generated cohort: ",  name.cohort[2], "\n",
                     "Database position: ", posi[1], " (gen), ", posi[2], " (sex), ", posi[3], " (first), ", posi[4], " (last).\n" ))
        }

      } else if(sum(breeding.size)>0){
        population$info$cohorts <- rbind(population$info$cohorts, c(name.cohort, current.gen+1, breeding.size[1:2], new.class[1+ as.numeric(breeding.size[2]>0)], current.size-breeding.size,
                                                                    time.point, creating.type))
        rownames(population$info$cohorts)[nrow(population$info$cohorts)] <- paste0(name.cohort)


        if(verbose){
          posi <- get.database(population, cohorts = name.cohort)
          cat(paste0("Successfully generated cohort: ", name.cohort, "\n",
                     "Database position: ", posi[1], " (gen), ", posi[2], " (sex), ", posi[3], " (first), ", posi[4], " (last).\n" ))
        }
      } else{
        warning(paste0("No individuals generated for cohort ", name.cohort))
      }

      if(sum(breeding.size)>0){
        if(verbose & store.comp.times){cat(paste0("Generation of new individuals took ", population$info$comp.times.general[nrow(population$info$comp.times.general),6]," seconds\n"))}
      }

      if(nrow(population$info$cohorts)<=2){
        colnames(population$info$cohorts) <- c("name","generation", "male individuals", "female individuals", "class", "position first male", "position first female",
                                               "time point", "creating.type")    }
    }
    if(Rprof){
      Rprof(NULL)
      population$info$Rprof[[length(population$info$Rprof)+1]] <- utils::summaryRprof()
    }
  }

  if(use.recalculate.manual){
    population = recalculate.manual(population, cohorts = name.cohort, store.comp.times= store.comp.times)

    if(store.comp.times){
      population$info$comp.times.general[nrow(population$info$comp.times.general)-1, ] = population$info$comp.times.general[nrow(population$info$comp.times.general)-1, ] +
        population$info$comp.times.general[nrow(population$info$comp.times.general), ]
      population$info$comp.times.general = population$info$comp.times.general[-nrow(population$info$comp.times.general), ,drop = FALSE]
    }

  }

  if(parallel.internal){
    return(population$breeding[[length(population$breeding)]])
  } else{
    return(population)
  }

}

