#'TargetExperiment and TargetExperimentList summary methods.
#'
#'Explore the TargetExperiment and TargetExperimentList's attribute values 
#'at feature and/or gene level.
#'
#'@param object TargetExperiment/TargetExperimentList class object.
#'@param attributeThres numeric indicating the intervals extreme values
#'required by summaryIntervals.
#'@param pool logical indicating if the summary should be performed for each
#'pool separately
#'@param ... required by summary.
#'
#'@return according to the call one of the following objects can be returned
#' \item{data.frame}{statistics of the analyzed attribute}
#' \item{data.frame}{Frequency table of the feature occurrence in the selected
#' intervals}
#'
#'@include TargetExperiment-constructor.R
#'@exportMethod summaryFeatureLev
#'@docType methods
#'@name summaryFeatureLev
#'@rdname TargetExperiment-statistics
#'@aliases summaryFeatureLev-methods
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
#'@examples
#'## Loading the TargetExperiment object
#'data(ampliPanel, package="TarSeqQC")
#'
#'# Summary at feature level
#'summaryFeatureLev(ampliPanel)
setGeneric(name="summaryFeatureLev", def=function(object){
    standardGeneric("summaryFeatureLev")
})
#'@name summaryFeatureLev
#'@import plyr
#'@importFrom stats IQR
#'@importFrom stats aggregate
#'@importFrom stats median
#'@importFrom stats sd
#'@rdname TargetExperiment-statistics
#'@aliases summaryFeatureLev,TargetExperiment-method
#'@inheritParams summaryFeatureLev
setMethod(f="summaryFeatureLev", signature=signature(object="TargetExperiment"),
definition=function(object){
    df_panel<-as.data.frame(getFeaturePanel(object))
    if(!(getAttribute(object) %in% c("coverage", "medianCounts"))){
        stop("Attribute slot should be defined in order to call the
            function")
    }
    attrSumm<-round(summary(df_panel[,getAttribute(object)]))
    if( "pool" %in% names(df_panel)){
        pool_summary<-(ddply(df_panel, "pool", function(x){
            round(summary(x[,getAttribute(object)]))
        }))
        pool_names<-paste("pool", pool_summary[,1], sep=" ")
        df<-data.frame(cbind(attrSumm),t(pool_summary[,-1]))
        names(df)<-c(getFeature(object), pool_names)
    }else{
        df<-data.frame(as.matrix(attrSumm))
        names(df)<-c(getFeature(object))
    }
    return(t(df))
})
#'@exportMethod summaryGeneLev
#'@name summaryGeneLev
#'@inheritParams summaryFeatureLev
#'@rdname TargetExperiment-statistics
#'@aliases summaryGeneLev-methods
#'@examples
#'# Summary at gene level
#'summaryGeneLev(ampliPanel)
setGeneric(name="summaryGeneLev", def=function(object){
    standardGeneric("summaryGeneLev")
})
#'@name summaryGeneLev
#'@rdname TargetExperiment-statistics
#'@inheritParams summaryFeatureLev
#'@aliases summaryGeneLev,TargetExperiment-method
setMethod(f="summaryGeneLev", signature=signature(object="TargetExperiment"),
definition=function(object){
    object@featurePanel<-getGenePanel(object)
    setFeature(object)<-"gene"
    return(summaryFeatureLev(object))
})
#'@name summary
#'@rdname TargetExperiment-statistics
#'@inheritParams summary
#'@exportMethods summary
#'@aliases summary,TargetExperiment-method
#'# Object summary
#'summary(ampliPanel)
setMethod(f="summary", signature=signature(object="TargetExperiment"),
definition=function(object,...){
    summaryDF<-rbind(summaryGeneLev(object),
        summaryFeatureLev(object))
    return(summaryDF)
})
#'@exportMethod summaryIntervals
#'@name summaryIntervals
#'@inheritParams summaryFeatureLev
#'@rdname TargetExperiment-statistics
#'@aliases summaryIntervals-methods
#'@examples
#'# Defining the attribute interval extreme values
#'attributeThres<-c(0,1,50,200,500, Inf)
#'# Doing a frequency table for the attribute intervals
#'summaryIntervals(ampliPanel, attributeThres=attributeThres)

setGeneric(name="summaryIntervals", def=function(object,attributeThres=c(0, 
1, 50, 200, 500, Inf), pool=FALSE){
    standardGeneric("summaryIntervals")
})
#'@name summaryIntervals
#'@inheritParams summaryFeatureLev
#'@rdname TargetExperiment-statistics
#'@aliases summaryIntervals,TargetExperiment-method
setMethod(f="summaryIntervals",signature=signature(object="TargetExperiment"),
definition=function(object,attributeThres= c(0, 1, 50, 200, 500, Inf), 
pool=FALSE){

    if(pool & !("pool" %in% names(mcols(getFeaturePanel(object))))){
        stop("'pool' was set as TRUE but the bed file doesn't contain a 'pool'
            column")
    }
    if(pool){
    df_panel<-as.data.frame(getFeaturePanel(object))[,c("pool", 
        getAttribute(object)), drop=FALSE]
    }else{
    df_panel<-as.data.frame(getFeaturePanel(object))[,getAttribute(object), 
        drop=FALSE]
    
    }
    if(attributeThres[length(attributeThres)] < Inf){
        attributeThres<-c(attributeThres, Inf)
    }
    # creating a new variable 'score' that groups the features according to 
    # their attribute value and defined intervals
    df_panel[,"score"]<-cut(df_panel[,getAttribute(object)], 
        breaks=attributeThres, include.lowest=TRUE, right=FALSE, dig.lab = 6)
    if(pool){
        att_table<-as.data.frame(table(df_panel[,"score"], df_panel[,"pool"]))

    }else{
        att_table<-as.data.frame(table(df_panel[,"score"]))
    }
    interval_names<-sapply(1:length(attributeThres[attributeThres != "Inf"]),
    function(x){
        if(x < length(attributeThres[attributeThres != "Inf"])) {
            return((paste(attributeThres[x], " <= ", getAttribute(object)," < ",
            attributeThres[x+1])))
        }else{
            paste(  getAttribute(object), " >= ", attributeThres[x])
        }
    })
    att_table[,"Var1"]<-interval_names
    
    if(pool){
        poolLevels<-levels(as.factor(df_panel[,"pool"]))
        tabla<-lapply(1:length(poolLevels), function(i){
            att_tableP<-att_table[att_table[,"Var2"] == poolLevels[i],c("Var1",
                "Freq")]
            tabla <- cbind(att_tableP,cumsum(att_tableP[,"Freq"]), 
                round(100*att_tableP[,"Freq"]/sum(att_tableP[,"Freq"]),1))  

            colnames(tabla) <- c(paste(getFeature(object), "_", 
                getAttribute(object), "_intervals", sep=""),"abs","cum_abs",
                "rel")
            tabla[,"cum_rel"]<-cumsum(tabla[, "rel"])
            if(tabla[nrow(tabla),"cum_rel"] != 100 ){
                tabla[tabla[,"cum_rel"]==tabla[nrow(tabla),"cum_rel"], 
                    "cum_rel"] <-100
            }
            return(tabla)
        })
    names(tabla)<-poolLevels
    }else{
    tabla <- cbind(att_table,cumsum(att_table[,"Freq"]),round(
                    100*att_table[,"Freq"]/sum(att_table[,"Freq"]),1))  
    colnames(tabla) <- c(paste(getFeature(object), "_", getAttribute(object), 
                        "_intervals", sep=""),"abs","cum_abs","rel")
    tabla[,"cum_rel"]<-cumsum(tabla[, "rel"])
    if(tabla[nrow(tabla),"cum_rel"] != 100 ){
        tabla[tabla[,"cum_rel"]==tabla[nrow(tabla),"cum_rel"], "cum_rel"] <-100
    }
    
    
    }
    return(tabla)
})
