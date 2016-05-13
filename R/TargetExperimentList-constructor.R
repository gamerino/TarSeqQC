#'TargetExperimentList constructor
#'
#'\code{TargetExperimentList} creates a TargetExperimentList object containing
#'a set of targeted sequencing experiment results, all those, carried out using
#'the same bed file. Feature parameter specifies what represent each panel 
#'element (bed file row). Attribute parameter indicates which attribute would 
#'be analyzed, 'coverage' or 'medianCounts' and should be specified in order
#'to indicate which coverage or medianCounts should be conserved.
#'
#'@param TEList List containing all the TargetExperiment objects corresponding
#'to the experiments that will be compared.
#'@param feature Character indicating the name of the feature that will
#'be explored (e.g 'amplicon', 'transcript', 'gene').
#'@param attribute Character indicating the name of the attribute that will
#'be explored. Should be 'coverage' or 'medianCounts'.
#'
#'@return TargetExperimentList object.
#'
#'@include TargetExperimentList-initialize.R
#'@export TargetExperimentList
#'@docType methods
#'@name TargetExperimentList
#'@importClassesFrom IRanges IRanges
#'@importFrom IRanges IRanges
#'@rdname TargetExperimentList-constructor
#'@aliases TargetExperimentList-methods
#'@seealso \code{\link{TargetExperimentList-class}}1
#'@note see full example in \code{\link{TargetExperimentList-class}}
#'@family TargetExperimentList
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar} and Elmer A. Fernandez 
#'\email{efernandez@@bdmg.com.ar}
#'@examples
#'# Defining the set of TargetExperiment objects
#'data(ampliPanel, package="TarSeqQC")
#'data(ampliPanel2, package="TarSeqQC")
#'ampliList<-list(ampliPanel, ampliPanel2)
#'# Defining feature parameter
#'feature<-"amplicon"
#'# Defining attribute parameter
#'attribute<-"coverage"
#'##Calling the constructor
#'object<-TargetExperimentList(TEList=ampliList, attribute=attribute,
#'feature=feature)
#'
TargetExperimentList<-function(TEList, feature=NULL, attribute="coverage"){
    .Object<-new("TargetExperimentList")
    if(nargs() >=1 ){
        .Object<-initialize(.Object, TEList, feature, attribute)
    }
    ##Check the object's validity
    validObject(.Object)
    return(.Object)
}
