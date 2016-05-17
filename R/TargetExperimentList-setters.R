#'Setters for the TargetExperimentList slots
#'
#'Set TargetExperimentList slots, according to the given function call.
#'
#'@include TargetExperimentList-getters.R
#'@exportMethod setFeature<-
#'@docType methods
#'@name setFeature<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature
#'@aliases setFeature<-,TargetExperimentList,character-method
#'@note see full example in \code{\link{TargetExperimentList-class}}
#'@examples
#'## loading TargetExperimentList object
#'data(TEList, package="TarSeqQC")
#'## Set feature slot value
#'setFeature(TEList)<-"amplicon"
#'
setReplaceMethod(f="setFeature", signature=signature(object=
    "TargetExperimentList", value="character"),  definition= function(object, 
    value){
        object@feature<-value
        validObject(object)
        return(object) 
})
