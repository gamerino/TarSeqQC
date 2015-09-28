#'Function to build a feature panel based on specific genomic regions.
#'
#'\code{buildFeaturePanel} builds panel slots of a TargetExperiment object.
#'Input can be a bam file or a pileup matrix. If the bed file contains a high
#'number of amplicons, the bam file as input is recommended in order to
#'diminish memory requirements. The resulting object is a GRanges instance
#'having panel and counts/coverage information.
#'
#'@param object TargetExperiment class object.
#'@param BPPARAM An optional BiocParallelParam instance defining the parallel
#'back-end to be used during evaluation.
#'
#'@return GRanges object.
#'
#'@include TargetExperiment-pileupCounts.R
#'@exportMethod buildFeaturePanel
#'@docType methods
#'@name buildFeaturePanel
#'@rdname TargetExperiment-buildFeaturePanel
#'@aliases buildFeaturePanel-methods
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar} and Elmer A. Fernandez 
#'\email{efernandez@@bdmg.com.ar}
#'@examples
#'## loading TargetExperiment object
#'data(ampliPanel, package="TarSeqQC")
#'## Defining bam file, bed file and fasta file names and paths
#'setBamFile(ampliPanel)<-system.file("extdata", "mybam.bam", 
#'    package="TarSeqQC", mustWork=TRUE)
#'setFastaFile(ampliPanel)<-system.file("extdata", "myfasta.fa", 
#'    package="TarSeqQC", mustWork=TRUE)
#'
#'myFeaturePanel<-buildFeaturePanel(ampliPanel)
setGeneric(name="buildFeaturePanel", def=function(object, 
BPPARAM=bpparam()){
    standardGeneric("buildFeaturePanel")
})
#'
#'@name buildFeaturePanel
#'@rdname TargetExperiment-buildFeaturePanel
#'@importMethodsFrom BiocParallel bpmapply bplapply bpworkers bpprogressbar<-
#'@import BiocParallel
#'@aliases buildFeaturePanel,TargetExperiment-method
#'@inheritParams buildFeaturePanel
setMethod(f="buildFeaturePanel", signature="TargetExperiment",
definition=function(object,BPPARAM=bpparam()){
#     BPPARAM2<-bpparam()
#     bpprogressbar(BPPARAM2)<-TRUE
    bed_file<-getBedFile(object)
#     feature_summary<-bplapply(as.list(1:length(bed_file)), 
#     function(feature){
#         bed_feature<-bed_file[feature]
#         subObject<-object
#         subObject@bedFile<-bed_feature
#         setScanBamP(subObject)<-ScanBamParam(which=bed_feature)
        
    feature_counts<-pileupCounts(bed=getBedFile(object), 
        bamFile=path(getBamFile(object)), fastaFile=path(getFastaFile(object)),
        scanBamP=getScanBamP(object), pileupP=getPileupP(object), 
        BPPARAM=BPPARAM)
    features<-as.character(unique(feature_counts[,"which_label"]))
#     remove duplicated features cased by overlapped features
    features<-unique(do.call(c,lapply(features, function(x){
        return(strsplit(x, split=";")[[1]])
        })))
    feature_summary <- bplapply(features,
        function(x) {
        counts <- feature_counts[feature_counts[,"which_label"] == x, "counts"] 
        median_counts <- round(median(counts))
        IQR_counts<-round(IQR(counts))
        coverage<-round(mean(counts))
        sd_coverage<-round(sd(counts))
        return(c(medianCounts=median_counts, IQRCounts=IQR_counts,
            coverage=coverage, sdCoverage=sd_coverage))
        },
        BPPARAM=BPPARAM)
    feature_summary<-as.data.frame(do.call(rbind,feature_summary))
    featurePanel<-getBedFile(object)  
    mcols(featurePanel)<-cbind(mcols(bed_file), feature_summary) 
    return(featurePanel)
})
