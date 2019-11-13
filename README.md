# MoBPS
This repository contains our R-package MoBPS and the associated packages (miraculix/RandomFieldsUtils/MoBPSmaps).
A manuscript for the package is finally available in preprint (https://www.biorxiv.org/content/10.1101/829333v1)
 
The package is designed in a way to allow for a maximum of flexibility and possible extensions in basically any step of the simulation. In case you feel like a specific functionally or option is missing in the program just contact me (torsten.pook@uni-goettingen.de). 
I would highly appreciate detailled explanation of the genetics/breeding you are trying to model to make it easier for me to add in the options needed in an efficient manner. 

Same goes for questions regarding the tool or how to set up your simulation. We are always happy for questions as it really helps improve the tool. For quick reply it would help to provide a small example of our problem (ideally the population-list as a .RData object)

Hopefully our extensive User Manual containing some exemplary simulation should answer most questions but it can definitly be still improved. Slide from presentations / working paper are available on request.

We are always thankful for advice, for additional things to implement, feedback and/or reports of errors.

# Web-based interface
The web-based interface is finally online at www.mobps.de (http://134.76.18.242/).

For testing we are providing a test user (EAAPguest, password: eaap2019 ). This user is not allowed to use our backend server simulations, but can at least try to set up his simulations and download resulting json-files that can than be simulated via json.simulation() in the MoBPS R-package.

Note that this web-interface is still under active development with the final structure of data storage not 100% set. Additional providing feedback on "False" input for the tool is still a major concern to guide the user in inputing resonable / expected things. This web-based interface is expilicitly NOT part of our recently submitted paper // biorvix manuscript. A seperate publication on this will follow, but do not expect this to be ready soon. Until offer close collaboration and assistants in setting up your  simulation study (contact: torsten.pook@uni-goettingen.de for that). We are interested in collaborations with partners from both industry and academia.

# Update

Updates since release:

### Version 1.4.19 (13.11.19) - Only Development branch

Function bv.standardization for standardization of trait mean and variance

Function get.genotyped to extract which individuals were genotyped

Marker Assisted Selection via mas.bve

Updates to json.simulation for web-based application 

Added BGLR estimators (BL, BRR, BayesA, BayesB, BayesC) 

Renamed ogc_cAc to ogc.cAc in breeding.diploid()


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
