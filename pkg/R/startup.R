.onAttach <- function(libname, pkgname) {
  mobps_info <- utils::sessionInfo()$otherPkgs$MoBPS
  miraculix_loaded <- requireNamespace("miraculix", quietly = TRUE)

  msg <- paste0(
    "#############################################################\n",
    "############ Modular Breeding Program Simulator #############\n",
    "#############################################################\n",
    "############## Version: ", mobps_info$Version, " (", mobps_info$Date, ") ###############\n",
    "#############################################################\n"
  )

  if (miraculix_loaded) {
    miraculix_info <- utils::sessionInfo()$otherPkgs$miraculix
    msg <- paste0(
      msg,
      "######## Miraculix detected and successfully loaded #########\n",
      "##################### Version: ", miraculix_info$Version, " ######################\n"
    )

    if (length(miraculix_info$Version) != 1 || miraculix_info$Version != "1.5.1.1") {
      msg <- paste0(msg, "###### Consider upgrading miraculix to version 1.5.1.1 ######\n")
    }
  } else {
    msg <- paste0(
      msg,
      "#################### Miraculix not found. ###################\n",
      "######## Consider installing to speed-up computations #######\n"
    )
  }

  msg <- paste0(
    msg,
    "#############################################################\n",
    "######## To update to the most recent stable version: #######\n",
    "## devtools::install_github('tpook92/MoBPS', subdir='pkg') ##\n",
    "#############################################################\n",
    "################ Web-interface: www.mobps.de ################\n",
    "### Extended documentation: www.github.com/tpook92/MoBPS ####\n",
    "#############################################################"
  )

  packageStartupMessage(msg)
}
