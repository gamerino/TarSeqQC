#'TargetExperimentList object constructor.
#'
#'\code{initialize} creates the TargetExperimentList object containing the 
#'experiment results of several targeted sequencing experiments carried out
#'using a unique bed file.
#'
#'@param .Object TargetExperimentList class.
#'@param TEList List containing all the TargetExperiment objects corresponding
#'to the experiments that will be compared.
#'@param feature Character indicating the name of the feature that will
#'be explored (e.g 'amplicon', 'exon', 'gene').
#'@param attribute Character indicating the name of the attribute that will
#'be explored. Should be 'coverage' or 'medianCounts'.
#'@return TargetExperimentList object.
#'
#'@include TargetExperimentList-TEList.R
#'@exportMethod initialize
#'@docType methods
#'@rdname TargetExperimentList-initialize
#'@import Rsamtools 
#'@import GenomicRanges
#'@importFrom stats IQR
#'@importFrom stats aggregate
#'@importFrom stats median
#'@importFrom stats sd
#'@importMethodsFrom BiocGenerics strand
#'@importClassesFrom S4Vectors Rle
#'@importMethodsFrom S4Vectors Rle
#'@importClassesFrom IRanges IRanges
#'@importFrom IRanges IRanges
#'@aliases TargetExperimentList-method
#'@seealso \code{\link{TargetExperimentList}}
#'@note see full example in \code{\link{TargetExperimentList-class}}
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
#'
setMethod(f="initialize", signature=signature(.Object="TargetExperimentList"),
definition=function(.Object, TEList,feature=NULL, attribute="coverage"){
    ##Set the different slots
    if(nargs() >=2){
        if(is.null(names(TEList))){
            names(TEList)<-paste("subject", 1:length(TEList), sep="_")
        }
        if (length(TEList) <2){
            stop("The TEList should contain at least two TargetExperiment 
                objects")
        }
    # bedFile slot
        bed<-getBedFile(TEList[[1]])
        bed_df<-as.data.frame(bed)
        for (i in 2:length(TEList)){
            bed_i<-getBedFile(TEList[[i]])
            if(any(bed != bed_i)){
                stop("The TargetExperiments should have the same bed file")
            }
        
        }
        
        if(!(all(c("seqnames", "start","end","gene") %in% names(bed_df)))){
            stop("Error! Bed file should contain at least five columns: 
                'seqnames', 'start', 'end' and 'gene'")
        }
        if(any(duplicated(rownames(bed_df)))){
            stop("Error! Bed file should one row per feature. The provided file 
                has duplicated 'name' IDs.")
        }
        if (any(duplicated(bed_df[,"start"]) & duplicated(bed_df[,"end"]))){
        warning("Your bed file have duplicated features. Removing 
            duplicated ...")
        index<-which(duplicated(bed_df[,"start"])& duplicated(bed_df[,"end"]))
        bed_df<-bed_df[-index,]
        }
        .Object@bedFile<-bed
        # panels slot
        GRpanel<-as.data.frame(getFeaturePanel(TEList[[1]]))
        panel<-GRpanel[,!(colnames(GRpanel) %in% c("coverage", "sdCoverage", 
            "medianCounts", "IQRCounts"))]
        panelNames<-colnames(panel)
        attribute<-getAttribute(TEList[[1]]) 
        if(!(attribute %in% c("coverage", "medianCounts"))){
            stop("Attribute slot should be defined in order to call the
                function")
        }
        
        for (i in 1:length(TEList)){
            panel<-cbind(panel, mcols(getFeaturePanel(TEList[[i]]))[,
                attribute])
        }
#         colnames(panel)<-c(panelNames, paste(attribute,"subject", 
#             1:length(TEList), sep="_"))
        colnames(panel)<-c(panelNames, paste(attribute,names(TEList), sep="_"))
        finalPanel<-GRanges(seqnames=panel[,"seqnames"], ranges=IRanges(start=
            panel[,"start"], end=panel[,"end"], names=rownames(panel)), 
            strand=panel[,"strand"])
        mcols(finalPanel)<-panel[,!(names(panel) %in% c("seqnames", "start", 
            "end", "strand", "width")), drop=FALSE]
        .Object@panels <-finalPanel
        # feature slot
        .Object@feature<-getFeature(TEList[[1]])
        # attribute slot
        .Object@attribute<-attribute
        
    }else{

        .Object@bedFile<-GRanges()
        .Object@panels<-GRanges()
        .Object@feature<-character(0)
        .Object@attribute<-character(0)
    }
    ##Check the object's validity
    validObject(.Object)
    return(.Object)
})
