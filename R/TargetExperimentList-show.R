#'@include TargetExperimentList-setters.R
#'@exportMethod show
#'@docType methods
#'@name show
#'@rdname TargetExperiment-show
#'@aliases show,TargetExperimentList-method
#'@note see full example in \code{\link{TargetExperimentList-class}}
#'@author Gabriela Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar} and Elmer A. Fernandez 
#'\email{efernandez@@bdmg.com.ar}
#'@examples
#'## Loading the TargetExperimentList object
#'data(TEList, package="TarSeqQC")
#'show(TEList)
setMethod(f="show", signature=signature(object="TargetExperimentList"),
definition=function(object){
    cat("TargetExperimentList \n")
    cat(paste(getFeature(object), "panels: \n", sep=" "))
    if(length(getPanels(object)) > 0 ){
        cat("\t");show(getPanels(object)[1:min(3, length(getPanels(object)))]);
        cat("\n") 
    }else{
        getPanels(object)
    }
    cat("selected attribute: \n")
    cat("\t", getAttribute(object),"\n") 
})
