# MoBPS
This repository contains our R-package MoBPS and the associated packages (miraculix/RandomFieldsUtils/MoBPSmaps).
A manuscript for the package is published in G3 Genes Genomes Genetics (https://www.g3journal.org/content/early/2020/03/30/g3.120.401193)
 
The package is designed in a way to allow for a maximum of flexibility and possible extensions in basically any step of the simulation. In case you feel like a specific functionally or option is missing in the program just contact me (torsten.pook@uni-goettingen.de). 
I would highly appreciate detailled explanation of the genetics/breeding you are trying to model to make it easier for me to add in the options needed in an efficient manner. 

Same goes for questions regarding the tool or how to set up your simulation. We are always happy for questions as it really helps improve the tool. For quick reply it would help to provide a small example of our problem (ideally the population-list as a .RData object)

Hopefully our extensive User Manual containing some exemplary simulation should answer most questions but it can definitly be still improved. Slide from presentations / working paper are available on request.

We are always thankful for advice, for additional things to implement, feedback and/or reports of errors.

# Web-based interface
The web-based interface is finally online at www.mobps.de.

For testing we are providing a Guest User. This user is not allowed to use our backend server simulations, but can at least try to set up his simulations and download resulting json-files that can than be simulated via json.simulation() in the MoBPS R-package.

A separte publication of the web-based framework is uploaded on bioRvix (https://www.biorxiv.org/content/10.1101/2020.07.08.193227v1) and extended documenation for the interface can be found at www.mobps.de directly.


# Update

Updates since initial release:

### Version 1.5.28 (22.07.20)

Added maternal, paternal traits (is.maternal / is.paternal in creating.diploid() / creating.trait() ) 

Added feature that traits can be combination of other traits (e.g. partly maternal traits) (add.combi() ) 

Minor hotfix for subpopulations in json.simulation()

Updated Guidelines to MoBPS (including hyperlinks)

### Version 1.5.25 (08.07.20)

Added feature to model related founders (add.founder.kinship() ) 

"Fixed" when deriving pedigree-matrix when offspring where in earlier generations than parents (Now Forbidden!)

Minor bug-fix with repeatability / offspring.phenotypes / correlated traits with epistatic QTL effects

insert.bv can now also insert NA values for phenotypes

Minor updates for json.simulation in accordance to  www.mobps.de


### Version 1.5.16 (12.06.20)

Added repeatability concept (multiple phenotyping observation are not necessaryly independent anymore)

Fixed a bug that the residual variation of traits were not correlated when phenotyping for different traits was not done at the same time

Automated Log-file generation in json.simulation()

Updated Guidelines

### Version 1.5.11 (27.05.20)

Added avoid.mating.fullsib / avoid.mating.halfsib in breeding.diploid() to not generate offspring from mating of fullsib/halfsibs

Added option to scale according to variance of true genomic values in selection (multiple.bve.scale = "bv")

Fixed bug in kinship.exp.store that lead to crashes when parents / offspring were stored in the same generation

MoBPSmaps 0.1.9 - new genetic map for salmon (Tsai et al 2016)

### Version 1.5.04 (17.05.20) - offical 1.5 release

Minor bugfixes

Improved documentation

Updated Guidelines

### Version 1.5.0 (12.05.20) - only in development branch

Non-phenotyped individuals now have the phenotype NA instead of 0 

Default of bve.0isNA is now FALSE (as phenotypes of 0 usually do not code NA anymore)

Share phenotyped of the selected cohorts can now be controlled via share.phenotyped

Introduction of multiple genotyping arrays (select array used via genotyped.array and share genotyped via genotyped.share)

Add different genotyping arrays via add.array()-function, Default array has all markers

Renaming MoBPS parameters to more intuitive names (old ones are still usable):
new.bv.observation -> phenotyping (.gen/database/cohort)
new.bv.child -> phenotyping.child
computation.A -> relationship.matrix
computation.A.ogc -> relationship.matrix.ogc
new.phenotype.correlation -> new.residual.correlation

get.pheno now can extract phenotypes from all copies of an individual (set use.all.copy = TRUE)

Added get.selectionbve() to export the estimated breeding value from the last applied selection procedure 

Added sex.s as a parameter in breeding.diploid to controll offspring sex.

Fixed a bug when partically phenotyping individuals

Fixed a bug when traits were deleted from the population-list in creating.trait()

#### Still in work: 

simulation of imputation errors for BVE

Use of genotyping arrays for GWAS


### Version 1.4.92 (18.04.09)

Added copy.individual.m / copy.individual.f in breeding.diploid() for a more convenient way to copy selected individuals

Added genotyped.share in breeding.diploid() to add genotype data for genotyped.gen/database/cohorts after initial generation

Added bve.ignore.traits in breeding.diploid() to skip breeding value estimation for selected traits

Fixed issue that bpcm.conversion in breeding.diploid() was actually expected BP/M conversion and not BP/CM

Added set.class() function to manually change the class of selected gen/database/cohorts

Added plot() for population-lists to apply common analysis function bv.development() / kinship.development(), get.pca()

Updated Guidelines to MoBPS (e.g. examples on the use of how to generate genotyping/phenotyping data for subcohorts, generate traits based on real-world genotype+phenotype data)

Fixes Bug in miraculix when more than 2 million markers were used. Miraculix 0.9.19 now automatically detects compiler settings and installs the most efficient algorith depending on the system (only for Linux)

### Version 1.4.87 (27.03.09)

CRAN RELEASE VERSION

Only minor updates to documentation to pass remaining CRAN checks

### Version 1.4.85 (24.03.09)

Adding verbose statements for all functions

Fixed issue when manually setting per QTL variants

Added export function to extract the map of a population list (get.map)

Added coloring in get.pca according to class

### Version 1.4.82 (28.02.09)

Fixed bug in fixed.breeding when no individuals were selected (selection.size=0)

Updates to Guidelines_to_MoBPS

### Version 1.4.81 (24.02.09)

Renamed selection.criteria to selection.highest

Renamed selection.criteria.type to selection.criteria

Default of selection.m / selection.f is now automatically derived based on input in selection.criteria. "function" in case selection.criteria is used and "random" if not - unless of course manually set

selection.size is automatically calculated to be all available individuals if not provided

Remove pedmap.to.phasedbeaglevcf from the list of exported function (use of BEAGLE-jar not in line with CRAN policies)

par -setting in bv.development ect. are not overwritten automatically anymore

All print can now be deactivated via use of verbose in the respective function

Added missing uses of requireNamespace for suggested packages

Added warnings() / stop() for critical issues

Minor improvements to documentation

### Version 1.4.78 (13.02.20)

BVE can be skipped for selected traits via bve.ignore.traits in breeding.diploid()

Position (gen/database/cohort) is now automatically displayed at point of generation 

Added Value/Example to all exported functions

Removed last CRAN check notes

Miraculix/RandomFieldUtils updates (0.9.19 / 0.5.17). Mostly for automatic detection of AVX2 on Linux.
On windows, MoBPS will run without AVX2 unless specified while installing miraculix.

### Version 1.4.67 (05.02.20)

Default for miraculix.chol has been set to TRUE (If miraculix is available this leads to fast version of chol2inv(chol()) without downside)

### Version 1.4.62 (22.01.20) - Only Development branch

Improvements to storage / memory usage (use integers / remove attributes)

Implemented stability of A/G in Single-Step via Vitezica 2011


### Version 1.4.56 (10.01.20) - Only Development branch

Minor updates in json.simulation()

Removed bug when generating traits in n.additive etc.

Computational improvements when working with high number of QTL 

### Version 1.4.51 (08.01.20) - Only Development branch

Added separate storage for ownspring phenotypes and own phenotypes

### Version 1.4.49 (06.01.20) - Only Development branch

Exclude non-required indivividual in pedigree calculation

Fixed bug in "Use-offspring" when parent had no phenotyped offspring

Minor updates in json.simulation

### Version 1.4.43 (03.01.20) - Only Development branch

Added selection according to a threshold (instead of fixed number of individuals) threshold.selection in breeding diploid()

Generate PCA for selected cohorts (get.pca)

Memory efficient implementation of pedigree-matrix

Commutational more efficent bv.standardization()

Removed typo in get.vcf()

Transformation function for phenotypes

### Version 1.4.28 (04.12.19) - Only Development branch

Fixed bug when deriving pedigree matrix in a population with individuals that are generated via copy.individual

New implementation of get.vcf() that is not requiring the R-package synbreed anymore

get.pedmap() now as reasonal input of family (based on cohort if available), sex, paternal IDs

Better tracking of individual ages (e.g. added time of death for each individual)

Automatically reduce cohort sizes in json.simulation() when individual exit system via culling


### Version 1.4.22 (20.11.19) - Only Development branch

Updates to json.simulation

Added Folder for IMAGE workshop

now taking of age + time the cohort is there (only different when copy.individual is used)

Culling now uses age - not time the cohort is there

### Version 1.4.19 (13.11.19) - Only Development branch

Function bv.standardization for standardization of trait mean and variance

Marker Assisted Selection via mas.bve

Updates to json.simulation for web-based application 


### Version 1.4.15 (29.10.19)

Solve issues with running MoBPS without miraculix

Added BayesA,B,C, BL, BRR to BGLR options for BVE

New function: get.genotyped() to export which individuals are genotyped

Renamed ogc_cAc to ogc.cAc for general uniformity

Updated documentation (especially for web-based application)

Minor update for miraculix (v0.9.7) - not MoBPS related

Added frozen version of the current that (Submission-version)

### Version 1.4.10 (22.10.19)

Mostly improvments to reduce memory requirement

New function: add.diag(). R-matrix is not required anymore

Clean-up of memory in ssGBLUP

Better handling of duplicated individuals (generated via copy.indiduals)

### Version 1.4.3 (16.10.19)
MoBPSmaps 0.1.6 (Including maps for Wheat and Sorghum)

Removed typos in summary.population()

Updates to creating.diploid for trait generation via n.additive/dominant etc.

### Version 1.4.2 (14.10.19)
Variety of additions to json.simulation and user-interface

Added direct-mixed-model BVE for individuals without phenotype (vanRaden 2008)

Added use of Parent/Grandparent mean as breeding values

Further updates for miraculix/RandomFieldsUtils for compiler independent computing

### Version 1.3.1 (23.08.19)
New MoBPS Web-based application www.mobps.de

Hotfixes in Single Step

Generation/Database/Cohort based selection module

Selection Index according to Hazel and Lush + Miesenberger

Additions for json.simulation and Web-based application

Culling module

Improved import for progeny phenotypes

Detection of compiler settings for miraculix/RandomFieldsUtils


### Version 1.2.9 (25.06.19)
Improved import of offspring phenotypes (offspring.bve.XXX)

Additional functionality for handling ungenotyped individuals in BVE

Massiv! speed-up of kinship.exp()

population-list is now of class population with generic function for summary()

get.pedigree now also provides raw-output (before "M"/"F" coding)

Ensembl-maps now included in R-package MoBPSmaps

Major update to miraculix (0.7.8) including:

Prepration for CRAN submission

Direct generation of genotype data via miraculix (reduced memory needs in creating.diploid() )




### Version 1.1.39 (06.06.19)
Added Implemenation of Single Step GBLUP with H according to Legarra et al 2014

Added share.genotyped / genotyped.s to control which individuals are genotyped

Hot-fixes of typos in 1.1.35.

### Version 1.1.35 (05.06.19)

New server at http://134.76.18.242/ (Now providing 20 cores, 64GB RAM) 

Faster and improved pedigree-based breeding value estimation. New parameter depth.pedigree to control the depth of the pedigree

kinship.exp now supports gen/database/cohort structure. Old version of the function still available at kinship.exp.old

Generalized get.database()

Bug-fix in analyze.bv()

Automatically report accuracies of breeding value estimation

Individuals in bve.insert are now automatically added to pool of individuals to consider in bve

Added [[15/16]] to population$breeding[[generation]] to assign each individual an ID and avoid use of duplicated in breeding value estimation (generated via copy.individual).

Improved computation speed in json.simulation() and updates allow with the web-based application.

Updated documentation

### Version 1.1.24 (23.05.19)

Minor updates to json.simulation()

### Version 1.1.23 (20.05.19)

Implemented preliminary versions of OGC & ssBLUP (not efficient and no general used recommended)

Improved functions for derivation of empirical kinship (kinship.emp)

Exemplary maps added to the package (as size is high there will be a separate package (MoBPS_maps) coming soon!

Web-interface now hosted under http://134.76.137.69/ (For access to closed-beta contact me (torsten.pook@uni-goettingen.de)

Function to analyze correlation between bv/bve/pheno added analyze.bv()

Cohort based economic evaluation added compute.costs.cohorts()

Function to generate pedmap (PLINK-format) files added (get.pedmap())

Minor changes for nicer looking plots - mostly same functionality & input parameters

### Version 1.1.12 (18.04.19)

Updated versions of miraculix and RandomFieldsUtils

Minor modification to MoBPS itself - mostly related web-interface (still close-beta)

### Version 1.1.7 (15.04.19)

Minor updates to kinship.development, improved documentation.

Package can now be installed via install_github("tpook92/MoBPS", subdir="pkg")

### Version 1.1.6 (11.04.19)
Updates to documentation and guidelines. Especially example on how to generate a base population with LD and a hard-sweep.

Added function kinship.development

### Version 1.1.4 (09.04.19):
Lot of minor hotfixes for cohort-based implementations

Possible different number of phenotypes generated for each trait

Import of map from Ensembl (ensembl.map() )

Maps as potential input in creating.diploid

Remove single individuals from the group of individuals to select from (mostly relevant for Web-based application)

Boxplot alternative to bv.development ((bv.development.box))

### Version 1.1.1 (22.03.19):

Fixed some TYPOS in 1.1.0.

Cohort-based implementation of breeding value estimation, sigma.e / sigma.g computation, gwas, bve insertion.

Renamed sigma.s to sigma.g

Removed new.bv.observation.sex

Default for randomly generated SNP-datasets is now uniformly distributed minor allele frequency instead of all 0 - use dataset <- "all0" for the old version.

### Version 1.1.0 (20.03.19):

- As there are minor changes to the the data storage we highly discourage the use of population list generated in MoBPS 1.0.X or earlier in MoBPS 1.1+

Cohort/group/generation simulation of phenotypes

Multiple matings from the same dam/sire combination

Tracking of time point / generation type (mostly relevant for web-based application)

Improved version of bv.development

Empirical (and fast) version of kinship.emp

Combining of cohorts

### Version 1.0.1 (21.02.19):

Cohort-based versions of compute.costs and analyze.population added

### Version 1.0.0 (14.02.19):

Mostly fixes of the documentation.

### Version 0.14.22 (05.02.19):

Lots of minor updates & fixes resulting from practial use, paralellization now for both windows & linux for generation of new individuals, Conversion from json to R-script for the user-interface.

### Version 0.14.7 (08.01.19):

Hot-fix for chicken template in case small chromosomes contain less than 1 marker, add ignore.best as a selection technique

### Version 0.14.6 (03.01.19):

Minor fixes in creating.diploid (bp to cM conversion) + progress bars

Version 0.14.5 (20.12.18) including miraculix 0.3.2 & RandomFieldsUtils 0.4.0 (17.12.18):

C-implementation for remove.effect.position=TRUE// unification of RandomFieldsUtils-code to HaploBlocker. Add bp to cM conversion in creating.diploid

### Version 0.14.4 (13.12.18):

Minor Updates to get.pedigree, get.pedigree2, get.pedigree3 (display "0" for founder individuals)
Uniform display of individual names in get.recombi

### Version 0.14.3 (16.11.18):

Integer-storing in get.geno & get.haplo, fixed bug in creating.diploid when generating multiple traits and chromosomes jointly

### Version 0.14.2 (14.11.18):

Improved documentation with traditional help function in R
Corrected some non-uniform parameter namings
Fixed CRAN-check Warnings

### Version 0.14.1 (13.11.18):

Renaming parameters in breeding.diploid()
Adding new version of bv.development()
