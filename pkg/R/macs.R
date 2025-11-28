
## TO BE IMPLEMENTED IN THE FUTURE!
function(x){
  if(FALSE){
    library(AlphaSimR)
    output = runMacs2(nInd = 100, nChr = 5, segSites = c(1000,500,500,200,100), Ne = 25)
    
    dataset = t(pullSegSiteHaplo(output))
    dim(dataset)
    population = creating.diploid(dataset = dataset)
    
    ld.decay(population, gen = 1, plot = TRUE)
  }

}
