#'@include TargetExperimentList-show.R
#'@exportMethod print
#'@docType methods
#'@name print
#'@rdname TargetExperiment-print
#'@aliases print,TargetExperimentList-method
#'@note see full example in \code{\link{TargetExperimentList-class}}
#'@examples
#'## Loading the TargetExperimentList object
#'data(TEList,package="TarSeqQC")
#'print(TEList)
setMethod(f="print",signature=signature(x="TargetExperimentList"),
definition = function(x, ...){
    show(x)
}
)
