#'A set of two amplicon panels example for use the TarSeqQC R package.
#'
#'A non-real dataset containing amplicon sequencing results to test the TarSeqQC
#'package, principally the use of the TargetExperimentList class.
#'
#'\describe{
#'\item{bedFile}{Bed file containing 29 amplicons and 8 genes in 2 PCR pools.}
#'\item{panels}{GRanges obtaining amplicon coverage for two targeted sequencing
#' experiment performed using the same bed file}
#'\item{feature}{Character "amplicon" indicating that the analyzed features 
#'are amplicon sequences} 
#'\item{attribute}{Character "coverage"}
#'}
#'
#'@include TargetExperiment-plotMetaDataExpl.R
#'@docType data
#'@format A TargetExperimentList object
#'@source see \code{\link{TargetExperimentList-class}}
#'@name TEList
#'@family TargetExperimentList
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
NULL