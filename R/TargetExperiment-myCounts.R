#'A pileup matrix example for use the TarSeqQC R package.
#'
#'The pileup matrix obtained using pileupCounts. It is built on the non-real
#' dataset containing amplicon sequencing results to test the TarSeqQC package.
#'
#'\describe{
#'\item{pos}{genomic positions of the explored features.} 
#'\item{seqnames}{chromosomes of the explored features.} 
#'\item{seq}{reference nucleotide corresponding to the genomic position.}
#'\item{A,C,G,T,N}{number of nucleotide read.}
#'\item{=}{Amount of read nucleotides matching the reference nucleotide.}
#'\item{-}{Amount of read deletions.}
#'\item{which_label}{feature location.}
#'\item{counts}{Total read counts}
#'}
#'
#'@include TargetExperiment-readFrequencies.R
#'@docType data
#'@format A data.frame object
#'@source see \code{\link{TargetExperiment-class}}
#'@name myCounts
#'@family TargetExperiment
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
NULL