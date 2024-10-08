% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/founder.simulation.R
\name{founder.simulation}
\alias{founder.simulation}
\title{Founder simulation}
\usage{
founder.simulation(
  nindi = 100,
  sex.quota = 0.5,
  nsnp = 0,
  n.gen = 100,
  nfinal = NULL,
  sex.quota.final = NULL,
  big.output = FALSE,
  plot = TRUE,
  display.progress = TRUE,
  depth.pedigree = 7,
  dataset = NULL,
  vcf = NULL,
  chr.nr = NULL,
  bp = NULL,
  snp.name = NULL,
  hom0 = NULL,
  hom1 = NULL,
  bpcm.conversion = 0,
  freq = "beta",
  sex.s = "fixed",
  chromosome.length = NULL,
  length.before = 5,
  length.behind = 5,
  snps.equidistant = NULL,
  snp.position = NULL,
  position.scaling = FALSE,
  bit.storing = FALSE,
  nbits = 30,
  randomSeed = NULL,
  miraculix = TRUE,
  miraculix.dataset = TRUE,
  template.chip = NULL,
  beta.shape1 = 1,
  beta.shape2 = 1,
  map = NULL,
  verbose = TRUE,
  vcf.maxsnp = Inf,
  plot.ld = TRUE,
  plot.allele.freq = TRUE,
  xlim = NULL,
  ylim = NULL,
  nclass = 20,
  ylim.af = NULL,
  mutation.rate = 10^-8,
  remutation.rate = 10^-8,
  estimate.ne = TRUE,
  estimate.ld = TRUE
)
}
\arguments{
\item{nindi}{number of inidividuals to generate in a random dataset}

\item{sex.quota}{Share of newly added female individuals (deterministic if sex.s="fixed", alt: sex.s="random")}

\item{nsnp}{number of markers to generate in a random dataset}

\item{n.gen}{Number of generations to simulate (default: 100)}

\item{nfinal}{Number of final individuals to include (default: nindi)}

\item{sex.quota.final}{Share of female individuals in the final generation}

\item{big.output}{Set to TRUE to export map, population list and pedigree relationship}

\item{plot}{Set to FALSE to not generate LD-decay plot and allele frequency spectrum}

\item{display.progress}{Set FALSE to not display progress bars. Setting verbose to FALSE will automatically deactive progress bars}

\item{depth.pedigree}{Depth of the pedigree in generations (default: 7)}

\item{dataset}{SNP dataset, use "random", "allhetero" "all0" when generating a dataset via nsnp,nindi}

\item{vcf}{Path to a vcf-file used as input genotypes (correct haplotype phase is assumed!)}

\item{chr.nr}{Vector containing the assosiated chromosome for each marker (default: all on the same)}

\item{bp}{Vector containing the physical position (bp) for each marker (default: 1,2,3...)}

\item{snp.name}{Vector containing the name of each marker (default ChrXSNPY - XY chosen accordingly)}

\item{hom0}{Vector containing the first allelic variant in each marker (default: 0)}

\item{hom1}{Vector containing the second allelic variant in each marker (default: 1)}

\item{bpcm.conversion}{Convert physical position (bp) into a cM position (default: 0 - not done)}

\item{freq}{frequency of allele 1 when randomly generating a dataset}

\item{sex.s}{Specify which newly added individuals are male (1) or female (2)}

\item{chromosome.length}{Length of the newly added chromosome (default: 5)}

\item{length.before}{Length before the first SNP of the dataset (default: 5)}

\item{length.behind}{Length after the last SNP of the dataset (default: 5)}

\item{snps.equidistant}{Use equidistant markers (computationally faster! ; default: TRUE)}

\item{snp.position}{Location of each marker on the genetic map}

\item{position.scaling}{Manual scaling of snp.position}

\item{bit.storing}{Set to TRUE if the MoBPS (not-miraculix! bit-storing is used)}

\item{nbits}{Bits available in MoBPS-bit-storing}

\item{randomSeed}{Set random seed of the process}

\item{miraculix}{If TRUE use miraculix package for data storage, computations and dataset generation}

\item{miraculix.dataset}{Set FALSE to deactive miraculix package for dataset generation}

\item{template.chip}{Import genetic map and chip from a species ("cattle", "chicken", "pig")}

\item{beta.shape1}{First parameter of the beta distribution for simulating allele frequencies}

\item{beta.shape2}{Second parameter of the beta distribution for simulating allele frequencies}

\item{map}{map-file that contains up to 5 colums (Chromsome, SNP-id, M-position, Bp-position, allele freq - Everything not provides it set to NA). A map can be imported via MoBPSmaps::ensembl.map()}

\item{verbose}{Set to FALSE to not display any prints}

\item{vcf.maxsnp}{Maximum number of SNPs to include in the genotype file (default: Inf)}

\item{plot.ld}{Set FALSE to not generate the LD plot (default; TRUE)}

\item{plot.allele.freq}{Set FALSE to not generate the allele frequency spectrum plot (default: TRUE)}

\item{xlim}{Axis limits for the x-axis in the LD plot (default: NULL)}

\item{ylim}{Axis limits for the y-axis in the LD plot (default: NULL)}

\item{nclass}{Number of classes to consider in the allele frequency spectrum plot (default: 20)}

\item{ylim.af}{Axis limits for the allele frequency spectrum plot (default: NULL)}

\item{mutation.rate}{Mutation rate in each marker (default: 10^-8)}

\item{remutation.rate}{Remutation rate in each marker (default: 10^-8)}

\item{estimate.ne}{Set to FALSE to not estimate the effective population size (default: TRUE)}

\item{estimate.ld}{Set to FALSE to not estimate the ld decay (default: TRUE)}
}
\description{
Function to generate founder genotypes
}
\examples{
population <- founder.simulation(nindi=100, nsnp=1000, n.gen=5)
}
