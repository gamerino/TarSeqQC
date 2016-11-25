#'Getters for TargetExperiment and TargetExperimentList objects.
#'
#'Obtain information about TargetExperiment and TargetExperimentList slots, 
#'according to the given function call.
#'
#'@param object TargetExperiment/TargetExperimentList class object.
#'@param level Character indicating 'gene' or 'feature'. Useful to getRegion
#'function
#'@param ID Character indicating the feature name that getRegion should be 
#'found. 
#'@param collapse Logical. Should the region be collapsed?.
#'@param threshold Numeric what should be the minimum attribute value?.
#'
#'@return according to the call one of the following objects can be returned
#' \item{GRanges}{bed file of the experiment}
#' \item{BamFile}{reference to the BAM file}
#' \item{FaFile}{reference to the fasta file}
#' \item{GRanges}{feature panel with statistical information}
#' \item{GRanges}{summarized version of the feature panel at gene level}
#' \item{character}{name of the explored features (e.g 'amplicon', 'exon')}
#' \item{character}{name of the analyzed attribute ('coverage' or 
#''medianCounts')}
#' \item{ScanBamParam}{parameters for the scan of the BAM file}
#' \item{PileupParam}{parameters for the pileup building}
#' \item{data.frame}{regions or low counts features}
#' \item{data.frame}{regions definition for overlapped features}
#'@include TargetExperiment-ampliPanel2.R
#'@exportMethod getBedFile
#'@docType methods
#'@name getBedFile
#'@rdname TargetExperiment-getters
#'@aliases getBedFile-methods
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@seealso \code{\link{TargetExperiment-class}}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
#'@examples
#'## Loading the TargetExperiment object
#'data(ampliPanel,package="TarSeqQC")
#'## Get the bedFile slot
#'getBedFile(ampliPanel)

setGeneric(name="getBedFile", def=function(object){
    standardGeneric("getBedFile")
})
#'
#'@name getBedFile
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getBedFile,TargetExperiment-method
setMethod(f="getBedFile", signature="TargetExperiment", 
definition=function(object){
    return(object@bedFile)
})
#'
#'@exportMethod getBamFile
#'@name getBamFile
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getBamFile-methods
#'@examples
#'## Get the bamFile slot
#'getBamFile(ampliPanel)
setGeneric(name="getBamFile", def=function(object){
    standardGeneric("getBamFile")
})
#'
#'@name getBamFile
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getBamFile,TargetExperiment-method
setMethod(f="getBamFile", signature="TargetExperiment", 
definition=function(object){
    return(object@bamFile)
})
#'
#'@exportMethod getFastaFile
#'@name getFastaFile
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getFastaFile-methods
#'@examples
#'## Get the fastaFile slot
#'getFastaFile(ampliPanel)
setGeneric(name="getFastaFile", def=function(object){
    standardGeneric("getFastaFile")
})
#'
#'@name getFastaFile
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getFastaFile,TargetExperiment-method
setMethod(f="getFastaFile", signature="TargetExperiment", 
definition=function(object){
    return(object@fastaFile)
})
#'
#'@exportMethod getFeaturePanel
#'@name getFeaturePanel
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getFeaturePanel-methods
#'@examples
#'## Get the feateurePanel slot
#'getFeaturePanel(ampliPanel)
setGeneric(name="getFeaturePanel", def=function(object){
    standardGeneric("getFeaturePanel")
})
#'
#'@name getFeaturePanel
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getFeaturePanel,TargetExperiment-method
setMethod(f="getFeaturePanel", signature="TargetExperiment", 
definition=function(object){
    return(object@featurePanel)
})
#'
#'@exportMethod getGenePanel
#'@name getGenePanel
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getGenePanel-methods
#'@examples
#'## Get the genePanel slot
#'getGenePanel(ampliPanel)
setGeneric(name="getGenePanel", def=function(object){
    standardGeneric("getGenePanel")
})
#'
#'@name getGenePanel
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getGenePanel,TargetExperiment-method
setMethod(f="getGenePanel", signature="TargetExperiment", 
definition=function(object){
    return(object@genePanel)
})
#'
#'@exportMethod getFeature
#'@name getFeature
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getFeature-methods
#'@examples
#'## Get the Feature slot
#'getFeature(ampliPanel)
setGeneric(name="getFeature", def=function(object){
    standardGeneric("getFeature")
})
#'
#'@name getFeature
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getFeature,TargetExperiment-method
setMethod(f="getFeature", signature="TargetExperiment", 
definition=function(object){
    return(object@feature)
})
#'
#'@exportMethod getAttribute
#'@name getAttribute
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getAttribute-methods
#'@examples
#'## Get the attribute slot
#'getAttribute(ampliPanel)
setGeneric(name="getAttribute", def=function(object){
    standardGeneric("getAttribute")
})
#'
#'@name getAttribute
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getAttribute,TargetExperiment-method
setMethod(f="getAttribute", signature="TargetExperiment", 
definition=function(object){
    return(object@attribute)
})
#'
#'@exportMethod getScanBamP
#'@name getFeature
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getScanBamP-methods
#'@examples
#'## Get the scanBamP slot
#'getScanBamP(ampliPanel)
setGeneric(name="getScanBamP", def=function(object){
    standardGeneric("getScanBamP")
})
#'
#'@name getScanBamP
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getScanBamP,TargetExperiment-method
setMethod(f="getScanBamP", signature="TargetExperiment", 
definition=function(object){
    return(object@scanBamP)
})
#'
#'@exportMethod getPileupP
#'@name getPileupP
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getPileupP-methods
#'@examples
#'## Get the pileupP slot
#'getPileupP(ampliPanel)
setGeneric(name="getPileupP", def=function(object){
    standardGeneric("getPileupP")
})
#'
#'@name getRegion
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getPileupP,TargetExperiment-method
setMethod(f="getPileupP", signature="TargetExperiment", 
definition=function(object){
    return(object@pileupP)
})
#'@exportMethod getRegion
#'@name getRegion
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getRegion-methods
#'@examples
#'## Get the region related to a feature or a gene
#'getRegion(ampliPanel, level="gene", ID="gene7", collapse=FALSE)
setGeneric(name="getRegion", def=function(object, level, ID, collapse=TRUE){
    standardGeneric("getRegion")
})
#'
#'@name getRegion
#'@rdname TargetExperiment-getters
#'@inheritParams getRegion
#'@aliases getRegion,TargetExperiment-method
setMethod(f="getRegion", signature="TargetExperiment", 
definition=function(object, level, ID, collapse=TRUE){
    if(!(level %in% c("feature", "gene")) | all(ID != c(names(
        object@featurePanel), names(object@genePanel)))){
        stop("'level' can be 'feature' or 'gene' and 'ID' should be a feature or
            gene name contained in the BED file")
    }
    if(level == "feature"){
        if(all(names(object@featurePanel) !=  ID)){
            stop('The provided ID is not a valid name for features')
        }
        aux<-object@featurePanel[names( object@featurePanel) ==  ID]
        region<-data.frame(names=ID, seqname=as.character(seqnames(aux)),
            start = start(aux),  end=end(aux), gene=mcols(aux)[,"gene"])
    }else{
        aux<-object@featurePanel[mcols(object@featurePanel)[,"gene"] ==  ID]
        if(collapse){
            features<-paste(names(aux), sep= "", collapse=", ")
            startReg<-min(start(aux))
            endReg<-max(end(aux))
        }else{
            features<-names(aux)
            startReg<-(start(aux))
            endReg<-(end(aux))
        }
        region<-data.frame(names=features, seqname=unique(as.character(
        seqnames(aux))),start = startReg, end=endReg, gene=ID)
    }
    return(region)
})
#'@exportMethod getLowCtsFeatures
#'@name getLowCtsFeatures
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getLowCtsFeatures-methods
#'@examples
#'## Get the low counts features
#'getLowCtsFeatures(ampliPanel, level="feature")
setGeneric(name="getLowCtsFeatures", def=function(object, level, threshold=50){
    standardGeneric("getLowCtsFeatures")
})
#'
#'@name getLowCtsFeatures
#'@rdname TargetExperiment-getters
#'@inheritParams getLowCtsFeatures
#'@aliases getLowCtsFeatures,TargetExperiment-method
setMethod(f="getLowCtsFeatures", signature="TargetExperiment", 
definition=function(object, level, threshold=50){
    if(!(level %in% c("feature", "gene")) | threshold < 0){
        stop("'level' can be 'feature' or 'gene' and 'threshold' can't be a 
            negative number")
    }
    if(level == "feature"){
        aux<-object@featurePanel[mcols(object@featurePanel)[,
            object@attribute] < threshold]
    }else{
        aux<-object@genePanel[mcols(object@genePanel)[,
        object@attribute] < threshold]
    }
    lowCts<-data.frame(names=names(aux), seqname=as.character(seqnames(aux)),
        start = start(aux), end=end(aux), mcols(aux))
    return(lowCts)
})
#'@exportMethod getOverlappedRegions
#'@name getOverlappedRegions
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getOverlappedRegions-methods
#'@examples
#'## Get the regions defined by overlapped features
#'getOverlappedRegions(ampliPanel)
setGeneric(name="getOverlappedRegions", def=function(object, collapse=FALSE){
    standardGeneric("getOverlappedRegions")
})
#'
#'@name getOverlappedRegions
#'@rdname TargetExperiment-getters
#'@aliases getOverlappedRegions,TargetExperiment-method
#'@inheritParams getOverlappedRegions
#'
setMethod(f="getOverlappedRegions", signature=signature(
    object= "TargetExperiment"), definition=function(object, collapse=FALSE){
        
        featurePanel<-getFeaturePanel(object)
        df_panel<-as.data.frame(featurePanel)
        attribute<-getAttribute(object)
        feature<-getFeature(object)
        df_panel<-cbind(feature_id=row.names(df_panel), df_panel)
        names(df_panel)[names(df_panel) == attribute]<-"attr"
        overlapPanel<-reduce(featurePanel)
        dfover<-lapply(1:length(overlapPanel), function(x){
            start<-which(df_panel[,"start"] == start(overlapPanel)[x])
            end<-which(df_panel[,"end"] == end(overlapPanel)[x])
            width<-max(end)-min(start)+1
            chr<-as.character(df_panel[start[1], "seqnames"])
            gene<-df_panel[start[1], "gene"]
            starReg<-df_panel[start,"start"]
            endReg<-df_panel[end,"end"]
            attr<-mean(df_panel[c(start,end), "attr"])
            region_id<-paste("Region",x,sep="")
            return(cbind(chr=chr, start=starReg, end=endReg, gene=gene, 
                region_id=rep(region_id, width),attr_region=rep(attr, width)))
        })
        dfover<-as.data.frame(do.call(rbind,dfover))
        dfover[,"region_id"]<-factor(as.character(dfover[, "region_id"]),
            levels=paste("Region", 1:length(overlapPanel),sep=""))
        dfover[,"attr_region"]<-as.numeric(as.character(dfover[, 
            "attr_region"]))
        names(dfover)[names(dfover) == "attr_region"]<-paste("region_",
            attribute, sep="")
        df_panel<-cbind(df_panel, dfover[,c("region_id", paste("region_",
            attribute, sep=""))])
        if(collapse){
            df_panel<-unique(dfover)
        }
        return(df_panel)
}   
)
    
        
