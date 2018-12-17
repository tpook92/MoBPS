# MoBPS
This repository contains our R-package MoBPS and its dependencies (miraculix/RandomFieldsUtils).
 
The package is designed in a way to allow for a maximum of flexibility and possible extensions in basically any step of the simulation. In case you feel like a specific functionally or option is missing in the program just contact me (torsten.pook@uni-goettingen.de). 
I would highly appreciate detailled explanation of the genetics/breeding you are trying to model to make it easier for me to add in the options needed in an efficient manner.

The current version of the tool does not contain a web-based interface but this should be openly available soon (hopefully this year but you never know!). Instead we will provide an extensive User Manual and some examplary simulations that should help you perform your simulation. Slide from presentations / working paper are available on request.

Note that the currently available version is still an OPEN-BETA version and all results obtained should be viewed with caution and we do not take any liability for errors nor guaranty warranty.

We are always thankful for advice, for additional things to implement, feedback and/or reports of errors.

# Update
Start of Open-Beta: 12nd November 2018

Current programming projects:
1. Extension of the economic cost tracking
2. Web-based userinterface + direct link to R-package
3. Cleaning up the source code and removing old functions

Updates since release:

Version 0.14.1 (13.11.18):

Renaming parameters in breeding.diploid()
Adding new version of bv.development()

Version 0.14.2 (14.11.18):

Improved documentation with traditional help function in R
Corrected some non-uniform parameter namings
Fixed CRAN-check Warnings

Version 0.14.3 (16.11.18):
Integer-storing in get.geno & get.haplo, fixed bug in creating.diploid when generating multiple traits and chromosomes jointly

Version 0.14.4 (13.12.18):
Minor Updates to get.pedigree, get.pedigree2, get.pedigree3 (display "0" for founder individuals)
Uniform display of individual names in get.recombi

miraculix 0.3.2 & RandomFieldsUtils 0.4.0 (17.12.18):
C-implementation for remove.effect.position=TRUE// unification of RandomFieldsUtils-code to HaploBlocker 

