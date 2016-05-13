#'Show method for the TargetExperiment and TargetExperimentList classes.
#'
#'\code{show} a TargetExperiment/TargetExperimentList object
#'
#'Generic show method for TargetExperiment and TargetExperimentList classes 
#'output visualization.
#'
#'@param object TargetExperiment/TargetExperimentList class object
#'
#'@return console output of the object
#'
#'@include TargetExperiment-setters.R
#'@exportMethod show
#'@docType methods
#'@name show
#'@rdname TargetExperiment-show
#'@aliases show,TargetExperiment-method
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@author Gabriela Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar} and Elmer A. Fernandez 
#'\email{efernandez@@bdmg.com.ar}
#'@examples
#'## Loading the TargetExperiment object
#'data(ampliPanel, package="TarSeqQC")
#'show(ampliPanel)
setMethod(f="show", signature=signature(object="TargetExperiment"),
definition=function(object){
    cat("TargetExperiment \n")
    cat(paste(getFeature(object), "panel: \n", sep=" "))
    if(length(getFeaturePanel(object)) > 0 ){
        cat("\t");show(getFeaturePanel(object)[1:min(3, length(getFeaturePanel(
            object)))]);cat("\n") 
    }else{
        getFeaturePanel(object)
    }
    cat("gene panel: \n")
    if(length(getGenePanel(object)) > 0){
        cat("\t");show(getGenePanel(object)[1:min(3, length(getGenePanel(
            object)))]);cat("\n") 
    }else{
        getGenePanel(object)
    }
    cat("selected attribute: \n")
    cat("\t", getAttribute(object),"\n") 
})
