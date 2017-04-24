#'@return
#' \item{GRanges}{feature panels with statistical information}
#'
#'@include TargetExperimentList-initialize.R
#'@exportMethod getBedFile
#'@rdname TargetExperiment-getters
#'@aliases getBedFile,TargetExperimentList-method
#'@examples
#'## Loading the TargetExperimentList object
#'data(TEList,package="TarSeqQC")
#'## Get the bedFile slot
#'getBedFile(TEList)

setMethod(f="getBedFile", signature="TargetExperimentList", 
definition=function(object){
    return(object@bedFile)
})
#'
#'@exportMethod getPanels
#'@name getPanels
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@examples
#'## Get the panels slot
#'getPanels(TEList)
setGeneric(name="getPanels", def=function(object){
    standardGeneric("getPanels")
})
#'
#'@name getPanels
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getPanels,TargetExperimentList-method
setMethod(f="getPanels", signature="TargetExperimentList", 
definition=function(object){
    return(object@panels)
})
#'
#'@exportMethod getFeature
#'@name getFeature
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getFeature,TargetExperimentList-method
#'@examples
#'## Get the Feature slot
#'getFeature(TEList)
setMethod(f="getFeature", signature="TargetExperimentList", 
definition=function(object){
    return(object@feature)
})
#'
#'@exportMethod getAttribute
#'@name getAttribute
#'@rdname TargetExperiment-getters
#'@aliases getAttribute,TargetExperimentList-method
#'@inheritParams getBedFile
#'@examples
#'## Get the attribute slot
#'getAttribute(TEList)
setMethod(f="getAttribute", signature="TargetExperimentList", 
definition=function(object){
    return(object@attribute)
})
#'@exportMethod getRegion
#'@name getRegion
#'@rdname TargetExperiment-getters
#'@inheritParams getBedFile
#'@aliases getRegion,TargetExperimentList-method
#'@examples
#'## Get the region related to a feature or a gene
#'getRegion(TEList, level="gene", ID="gene7", collapse=FALSE)
#'
setMethod(f="getRegion", signature="TargetExperimentList", 
definition=function(object,level, ID, collapse=TRUE){
    if(!(level %in% c("feature", "gene")) | all(ID != c(names(
        object@panels), mcols(object@panels)[,"gene"]))){
        stop("'level' can be 'feature' or 'gene' and 'ID' should be a feature
            or gene name contained in the BED file")
    }
    if(level == "feature"){
        aux<-object@panels[names(object@panels) ==  ID]
        region<-data.frame(names=ID, seqname=as.character(seqnames(aux)),
            start = start(aux),  end=end(aux), gene=mcols(aux)[,"gene"])
    }else{
        aux<-object@panels[mcols(object@panels)[,"gene"] ==  ID]
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
#'@aliases getLowCtsFeatures,TargetExperimentList-method
#'@examples
#'## Get the low counts features
#'getLowCtsFeatures(TEList, level="feature")
setMethod(f="getLowCtsFeatures", signature="TargetExperimentList", 
definition=function(object, level, threshold=50){
    if(!(level %in% c("feature", "gene")) | threshold < 0){
        stop("'level' can be 'feature' or 'gene' and 'threshold' can't be a 
            negative number")
    }
    panel<-object@panels
    index<-which(lapply(colnames(mcols(panel)), function(i){
        return(strsplit(i, split="_")[[1]][1])}) == object@attribute)
    if(level == "feature"){
        rowIndex<-unique(do.call(c,lapply(1:ncol(mcols(panel)[,index]), 
            function(i){return(which(mcols(panel)[,index][,i] < threshold))})))
        lowCts<-as.data.frame(panel[rowIndex])
    }else{
        genes<-unique(mcols(panel)[,"gene"])
        lowCts<-as.data.frame(do.call(rbind,lapply(1:length(genes), 
        function(j){
            aux2<-mcols(panel)[mcols(panel)[,"gene"] == genes[j],index]
            attrG<-do.call(c,lapply(1:ncol(aux2), function(i){
                if (object@attribute=="coverage"){
                    geneAtt<-round(mean(aux2[,i]))
                }else{
                    geneAtt<-round(median(aux2[,i]))
                }

            }))
            if (any(attrG < threshold)){
                resultLCF<-c(genes[j], attrG)
                names(resultLCF)<-c("gene", colnames(mcols(panel)[index]))
                return(resultLCF)
            } 
        })))
        }

    return(lowCts)
})
