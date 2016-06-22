#'Function to summarize a featurePanel slot at a gene level.
#'
#'\code{summarizePanel} helps the initialization of a TargetExperiment object.
#'Is useful to summarize the featurePanel slot at a gene level, building the
#'genePanel slot.
#'
#'@param object TargetExperiment class object.
#'@param BPPARAM An optional BiocParallelParam instance defining the parallel
#'back-end to be used during evaluation.
#'
#'@return TargetExperiment object
#'
#'@include TargetExperiment-buildFeaturePanel.R
#'@exportMethod summarizePanel
#'@docType methods
#'@name summarizePanel
#'@rdname TargetExperiment-summarizePanel
#'@aliases summarizePanel-methods
#'@seealso \code{\link{TargetExperiment}},\code{\link{buildFeaturePanel}}
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
#'
#'@examples
#'## Loading the TargetExperiment object
#'data(ampliPanel, package="TarSeqQC")
#'
#'mySummarizedPanel<-summarizePanel(ampliPanel)
setGeneric(name="summarizePanel", def=function(object, BPPARAM=bpparam()){
    standardGeneric("summarizePanel")
})
#'
#'@name summarizePanel
#'@rdname TargetExperiment-summarizePanel
#'@importClassesFrom S4Vectors Rle
#'@importMethodsFrom S4Vectors Rle
#'@importClassesFrom IRanges IRanges
#'@importFrom IRanges IRanges
#'@import BiocParallel
#'@importMethodsFrom BiocParallel bpprogressbar<-
#'@importMethodsFrom GenomeInfoDb seqnames seqnames<- seqlevels
#'@aliases summarizePanel,TargetExperiment-method
#'@inheritParams summarizePanel
#'
setMethod(f="summarizePanel", signature="TargetExperiment", definition=function(
object, BPPARAM=bpparam()){
#     bpprogressbar(BPPARAM)<-TRUE
    colID<-"gene"
    featurePanel<-getFeaturePanel(object)
    col_ids<-as.character(unique(mcols(featurePanel)[,colID]))
    info<-bplapply(1:length(col_ids), function(ID){
        chr<-as.character(runValue(
        seqnames(featurePanel)[mcols(featurePanel)[,colID] == col_ids[[ID]]]))
        start<-min(start(featurePanel)[mcols(featurePanel)[,colID]==col_ids[[
            ID]]])
        end<-max(end(featurePanel)[mcols(featurePanel)[,colID]== col_ids[[ID]]])
        median_counts<-round(median(
        mcols(featurePanel)[mcols(featurePanel)[,colID]== col_ids[[ID]],
            "medianCounts"]))
        IQR_counts<-round(IQR(mcols(featurePanel)[mcols(featurePanel)[,colID]== 
            col_ids[[ID]], "medianCounts"]))
        coverage<-round(mean(mcols(featurePanel)[mcols(featurePanel)[,colID] ==
            col_ids[[ID]],"coverage"]))
        sd_coverage<-round(sd(mcols(featurePanel)[mcols(featurePanel)[,colID] ==
            col_ids[[ID]],"coverage"]))
        if(is.na(sd_coverage)) sd_coverage<-0
        return(c(chr=chr, start=start, end=end, ID=col_ids[[ID]], 
            medianCounts=median_counts, IQRCounts=IQR_counts, coverage=coverage,
            sdCoverage=sd_coverage))
    }, BPPARAM=BPPARAM)
    info<-as.data.frame(do.call(rbind,info), stringsAsFactor=FALSE)
    info[, "chr"]<-factor(info[, "chr"], levels=seqlevels(featurePanel))
    info[,c("medianCounts", "IQRCounts", "coverage", "sdCoverage")]<-apply(
    info[,c("medianCounts", "IQRCounts", "coverage", "sdCoverage")],2,
        as.numeric)
    info<-info[order(info[,"chr"]),]
    newPanel<-GRanges(seqnames=Rle(info$chr), ranges=IRanges(start=as.integer(
        as.character(info$start)), end=as.integer(as.character(info$end)),
        names=info$ID))
    mcols(newPanel)<-info[,c("medianCounts", "IQRCounts", "coverage", 
        "sdCoverage")]
    return(newPanel)
})
