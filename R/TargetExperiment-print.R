#'Print a TargetExperiment/TargetExperimentList object.
#'
#'Generic \code{print} method for TargetExperiment and TargetExperimentList 
#' classes and descendants.
#'
#'@param x TargetExperiment/TargetExperimentList class object.
#'@param ... Included for generic print compatibility.
#'
#'@return console output of the object.
#'
#'@include TargetExperiment-show.R
#'@exportMethod print
#'@docType methods
#'@name print
#'@rdname TargetExperiment-print
#'@aliases print,TargetExperiment-method
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@author Gabriela Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar} and Elmer A. Fernandez 
#'\email{efernandez@@bdmg.com.ar}
#'@examples
#'## Loading the TargetExperiment object
#'data(ampliPanel,package="TarSeqQC")
#'print(ampliPanel)
setMethod(f="print",signature=signature(x="TargetExperiment"),
definition = function(x, ...){
    show(x)
}
)
