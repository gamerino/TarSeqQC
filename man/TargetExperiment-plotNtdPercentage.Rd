% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TargetExperiment-plotNtdPercentage.R
\docType{methods}
\name{plotNtdPercentage}
\alias{plotNtdPercentage}
\alias{plotNtdPercentage,TargetExperiment-method}
\alias{plotNtdPercentage-methods}
\title{Plot nucleotide read percentages for a particular feature.}
\usage{
plotNtdPercentage(object, featureID, BPPARAM = bpparam())

\S4method{plotNtdPercentage}{TargetExperiment}(object, featureID,
    BPPARAM = bpparam())
}
\arguments{
\item{object}{a TargetExperiment object.}

\item{featureID}{a character indicating the feature ID.}

\item{BPPARAM}{An optional BiocParallelParam instance defining the parallel
back-end to be used during evaluation.

returned by the function.}
}
\value{
ggplot2 graphics
}
\description{
\code{plotNtdPercentage} plots the percentages of the ocurrence of each
nucleotide in each position for a selected feature.
}
\note{
see full example in \code{\link{TargetExperiment-class}}
}
\examples{
## loading TargetExperiment object
data(ampliPanel, package="TarSeqQC")
## Defining bam file, bed file and fasta file names and paths
setBamFile(ampliPanel)<-system.file("extdata", "mybam.bam", 
    package="TarSeqQC", mustWork=TRUE)
setFastaFile(ampliPanel)<-system.file("extdata", "myfasta.fa", 
    package="TarSeqQC", mustWork=TRUE)
# Exploring the nucleotide percentages compositions of the read counts for a 
# particular amplicon
g<-plotNtdPercentage(ampliPanel,featureID="AMPL20")
if(interactive()){
g
}
}
\author{
Gabriela A. Merino \email{gmerino@bdmg.com.ar}, Cristobal Fresno
\email{cfresno@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
\email{efernandez@bdmg.com.ar}
}
\seealso{
\code{\link{plotFeature}}
}

