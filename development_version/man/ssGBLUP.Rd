% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssGBLUP.R
\name{ssGBLUP}
\alias{ssGBLUP}
\title{Single Step GBLUP}
\usage{
ssGBLUP(A11, A12, A22, G)
}
\arguments{
\item{A11}{pedigree relationship matrix of non-genotyped individuals}

\item{A12}{pedigree relationship matrix between non-genotyped and genotyped individuals}

\item{A22}{pedigree relationship matrix of genotyped individuals}

\item{G}{genomic relationship matrix of genotyped individuals}
}
\value{
Single step relationship matrix
}
\description{
Function to perform single step GBLUP according to Legarra 2014
}
