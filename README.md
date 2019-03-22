# MoBPS
This repository contains our R-package MoBPS and its dependencies (miraculix/RandomFieldsUtils).
 
The package is designed in a way to allow for a maximum of flexibility and possible extensions in basically any step of the simulation. In case you feel like a specific functionally or option is missing in the program just contact me (torsten.pook@uni-goettingen.de). 
I would highly appreciate detailled explanation of the genetics/breeding you are trying to model to make it easier for me to add in the options needed in an efficient manner. 

Same goes for questions regarding the tool or how to set up your simulation. We are always happy for questions as it really helps improve the tool. For quick reply it would help to provide a small example of our problem (ideally the population-list as a .RData object)

Hopefully our extensive User Manual containing some exemplary simulation should answer most questions but it can definitly be still improved. Slide from presentations / working paper are available on request.


Note that the currently available version is still an OPEN-BETA version and all results obtained should be viewed with caution and we do not take any liability for errors nor guaranty warranty. We are always thankful for advice, for additional things to implement, feedback and/or reports of errors.

# Web-based interface
We are now able to perform simulations based on a json-file generated in our web-based application. You can find some of the code to transfer the json-file into usage R-code for the tool. Development issues currently are mostly setting up a user/login structure and server to perform smaller simulations internally.
If interested in testing the web-based application just contact me (torsten.pook@uni-goettingen.de)

# Update

Updates since release:

User-Manual according to version 1.1.0!

Version 1.1.1 (22.03.19):

Fixed a bug in the generation of phenotypes + typos in 1.1.0.

Cohort-based implementation of breeding value estimation, sigma.e / sigma.g computation, gwas, bve insertion.

Renamed sigma.s to sigma.g

Removed new.bv.observation.sex

Default for randomly generated SNP-datasets is now uniformly distributed minor allele frequency instead of all 0 - use dataset <- "all0" for the old version.

Version 1.1.0 (20.03.19):

As there are minor changes to the the data storage we highly discourage the use of population list generated in MoBPS 1.0.X or earlier in MoBPS 1.1+

Cohort/group/generation simulation of phenotypes

Multiple matings from the same dam/sire combination

Tracking of time point / generation type (mostly relevant for web-based application)

Improved version of bv.development

Empirical (and fast) version of kinship.emp

Combining of cohorts

Version 1.0.1 (21.02.19):

Cohort-based versions of compute.costs and analyze.population added

Version 1.0.0 (14.02.19):

Mostly fixes of the documentation.

Version 0.14.22 (05.02.19):

Lots of minor updates & fixes resulting from practial use, paralellization now for both windows & linux for generation of new individuals, Conversion from json to R-script for the user-interface.

Version 0.14.7 (08.01.19):

Hot-fix for chicken template in case small chromosomes contain less than 1 marker, add ignore.best as a selection technique

Version 0.14.6 (03.01.19):

Minor fixes in creating.diploid (bp to cM conversion) + progress bars

Version 0.14.5 (20.12.18) including miraculix 0.3.2 & RandomFieldsUtils 0.4.0 (17.12.18):

C-implementation for remove.effect.position=TRUE// unification of RandomFieldsUtils-code to HaploBlocker. Add bp to cM conversion in creating.diploid

Version 0.14.4 (13.12.18):

Minor Updates to get.pedigree, get.pedigree2, get.pedigree3 (display "0" for founder individuals)
Uniform display of individual names in get.recombi

Version 0.14.3 (16.11.18):

Integer-storing in get.geno & get.haplo, fixed bug in creating.diploid when generating multiple traits and chromosomes jointly

Version 0.14.2 (14.11.18):

Improved documentation with traditional help function in R
Corrected some non-uniform parameter namings
Fixed CRAN-check Warnings

Version 0.14.1 (13.11.18):

Renaming parameters in breeding.diploid()
Adding new version of bv.development()
