#'TargetExperimentList S4 class implementation in R
#' 
#'This S4 class represents a collection of Targeted Sequencing Experiments in 
#'R. All these experiments are characterized by a 'bed file' containing the 
#'specification of the explored 'features', as a 'feature panel'. These 
#'features could be amplicons, exons, transcripts, among others. In general
#'each feature is associated to one gene but a gene could be related to many
#'features. This class allows the representation and quality control of a set
#'of Targeted Sequencing Experiment made over the same or different subjects
#'but using always the same bed file'.
#'
#'@slot bedFile GRanges object that models the bed file.
#'@slot panels GRanges object containing  the feature/gene panels.
#'@slot attribute character indicates which attribute, 'coverage' or 
#''medianCounts' will be used to the analysis.
#'@slot feature character indicates the name of the analyzed features. E.g
#''amplicon', 'exon', 'transcript', 'gene'.
#'
#'@section Features: 
#'\enumerate{
#'  \item Model sets of targeted sequencing experiments in R. 
#'  \item Evaluate the performance of the targeted sequencing technique across
#'  several experiments using coverage/read counts information. 
#'  \item Detect in early stage sequencing or library preparation errors.
#'  \item Report quality control results.
#'}
#'
#'@section Functions:
#'TargetExperimentList S4 class includes the following functions:
#'\describe{
#' \item{initialize}{constructor of TargetExperimentList to generate the 
#' feature panel starting from at least two TargetExperiment objects}
#' \item{getBedFile, getPanels, getAttribute, getFeature}{return the 
#' respective TargetExperimentList slots}
#' \item{setFeature}{set the respective TargetExperimentList slot}
#' \item{show}{generic output of the object}
#' \item{print}{generic output of the object}
#' \item{summary}{print statistics summary for the set attribute}
#' \item{plot}{plot a summarized view of the attribute values achieved by each
#' feature in each sample}
#' \item{plotGlobalAttrExpl}{plot the attribute distribution for each feature}
#' \item{plotAttrExpl}{plot the attribute distribution in each panel}
#' \item{plotpoolPerformance}{plot the attribute distribution in each or pool}
#' }
#'
#'@include TargetExperiment.R
#'@import methods GenomicRanges Rsamtools
#'@importMethodsFrom BiocGenerics strand
#'@importClassesFrom S4Vectors Rle
#'@importMethodsFrom S4Vectors Rle
#'@importClassesFrom IRanges IRanges
#'@importFrom IRanges IRanges
#'@name TargetExperimentList-class
#'@rdname TargetExperimentList-class
#'@exportClass TargetExperimentList
#'@seealso Rsamtools
#'@family TargetExperimentList
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
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
#'## load the example dataset
#'data(TEList, package="TarSeqQC")
#'setFeature(TEList)<-"amplicon"
#'## Early exploration
#'# show/print
#'TEList
#'# summary
#'summary(TEList)
#'## Controlling low counts features
#'# Definition of the interval extreme values
#'attributeThres<-c(0,1,50,200,500, Inf)
#'# Do a frequency table for the attribute intervals
#'summaryIntervals(TEList, attributeThres)
#'# getting low counts features at gene level
#'getLowCtsFeatures(TEList, level=gene, threshold=50)
#'# exploring panel performance along several samples
#'g<-plot(TEList, attributeThres=attributeThres, featureLabs =TRUE)
#'if(interactive()){
#'g
#'}
#'g<-plotGlobalAttrExpl(TEList,log=FALSE)
#'# x11(type="cairo")
#'if(interactive()){
#'g
#'}
#'g<-plotPoolPerformance(TEList, log=FALSE)
#'if(interactive()){
#'g
#'}
setClass(Class="TargetExperimentList", slots=list(
    bedFile ="GRanges",
    panels= "GRanges",
    feature = "character",
    attribute = "character"),
validity=function(object){

    ## Check bedFile
    bedFile<-object@bedFile
    if(class(bedFile)[[1]] != "GRanges"){
        stop(paste("The bedFile slot must be a GRanges object"))
    }
    ## Check attribute
    attribute<-object@attribute
    if(length(attribute) > 0 ){
        if(length(strsplit(attribute,"")[[1]]) > 1){
            if(!(attribute %in% c("coverage", "medianCounts"))){
                stop("The 'attribute' parameter should be 'coverage' or 
                    'medianCounts'")
            }
        }
    }
}, prototype=list(
    bedFile =GRanges(),
    panels= GRanges(),
    feature = character(0),
    attribute = character(0))
)
