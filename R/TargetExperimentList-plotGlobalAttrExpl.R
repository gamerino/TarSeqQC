#'Plot attribute exploration of a TargetExperimentList object.
#'
#'\code{plotGlobalAttrExpl} displays box-plot of the analyzed achieved 
#'attribute values along all samples and at a feature level. This graphic
#'could include density plot together the corresponding  box-plot using the 
#'ggplot2 geom_violin method.
#'
#'@param object TargetExperimentList class object.
#'@param attributeThres Numeric indicating the attribute interval extreme 
#'values.
#'@param dens Logical indicating if boxplot and density function should be 
#'plotted together using the ggplot2 geom_violin method or only the boxplot 
#'(dens=FALSE) should be displayed.
#'@param log Logical indicating if the attribute should be considered in 
#'log10 scale.
#'@param pool Logical indicating if the plots should be performed for each
#'pool separately
#'@param featureLabs logical indicating if feature names should be plotted
#'@param medianMarg numeric indicating the percentage of the median attribute 
#'value to be plotted as lines. If it is NULL no line will be displayed
#'@return ggplot2 graphics.
#'
#'@include TargetExperimentList-plot.R
#'@exportMethod plotGlobalAttrExpl
#'@docType methods
#'@name plotGlobalAttrExpl
#'@rdname TargetExperimentList-plotGlobalAttrExpl
#'@aliases plotGlobalAttrExpl-methods
#'@seealso \code{\link{plot}}
#'@note see full example in \code{\link{TargetExperimentList-class}}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
#'@examples
#'## Loading the TargetExperimentList object
#'data(TEList, package="TarSeqQC")
#'
#'# Attribute boxplot and density plot exploration
#'g<-plotGlobalAttrExpl(TEList,log=FALSE)
#'# x11(type="cairo")
#'if(interactive()){
#'g
#'}
setGeneric(name="plotGlobalAttrExpl", def=function(object, attributeThres=
c(0, 1, 50, 200, 500, Inf), dens=FALSE, log=FALSE, pool=FALSE, featureLabs=
FALSE,medianMarg=NULL){
    standardGeneric("plotGlobalAttrExpl")
})
#'
#'@name plotGlobalAttrExpl
#'@rdname TargetExperimentList-plotGlobalAttrExpl
#'@importFrom cowplot plot_grid
#'@aliases plotGlobalAttrExpl,TargetExperimentList-method
#'@inheritParams plotGlobalAttrExpl
setMethod(f="plotGlobalAttrExpl", signature=signature(object=
"TargetExperimentList"), definition=function(object,attributeThres=c(0, 1,
50, 200, 500, Inf),dens=FALSE, log=FALSE, pool=FALSE, featureLabs=FALSE, 
medianMarg=NULL){
    if(attributeThres[1] !=0){
        attributeThres<-c(0,attributeThres)
    }
    if(attributeThres[length(attributeThres)] !=Inf){
        attributeThres<-c(attributeThres, Inf)
    }
    df<-as.data.frame(getPanels(object))
    if(pool & !("pool" %in% colnames(df))){
        stop("'pool' is TRUE but the panels don't contain pool information")
    }
    if(length(medianMarg)!=2 & !is.null(medianMarg)){
        stop("only  2 median margins can be added")
    }
    attribute<-getAttribute(object)
    feature<-getFeature(object)
    y_lab<-capitalize(attribute)
    index<-do.call(c, lapply(1:ncol(df), function(i){
        if(strsplit(colnames(df)[i], split="_")[[1]][[1]] == 
            attribute){
            return(i)}
    }))
    df2<-cbind(df[, index, drop=FALSE], names=rownames(df))
    dfMelt<-melt(df2, id.vars="names")
    dfMelt[,"names"]<-factor(dfMelt[,"names"], levels=unique(df2[,"names"]))
    dfMelt<-dfMelt[order(dfMelt[, "names"]),]
    if(pool){
        dfInfo<-do.call(rbind, lapply(rownames(df), function(feature){
            i<-which(dfMelt[,"names"] == feature)
            return(cbind(median=rep(median(dfMelt[i, "value"]), 
                times=length(i)), pool=rep(df[rownames(df)==feature, "pool"],
                times=length(i))))
        
        }))
        dfInfo[,"pool"]<-factor(dfInfo[,"pool"], levels=sort(unique(df[,"pool"]
            )))
        dfMelt<-cbind(dfMelt, dfInfo)
    
    }else{
        dfInfo<-do.call(c, lapply(rownames(df), function(feature){
            i<-which(dfMelt[,"names"] == feature)
            return(median=c(rep(median(dfMelt[i, "value"]), 
                times=length(i))))
    
        }))
        dfMelt<-cbind(dfMelt, median=dfInfo)

    }
    scores<-cut(dfMelt[,"median"], breaks=attributeThres, include.lowest=TRUE,
        right=FALSE, dig.lab = 6)
    interval_names<-sapply(1:length(attributeThres[attributeThres != "Inf"]),
    function(x){
        if(x < length(attributeThres[attributeThres != "Inf"])) {
            return((paste(attributeThres[x], " <= ", getAttribute(object),
                " < ", attributeThres[x+1])))
        }else{
            paste(  getAttribute(object), " >= ", attributeThres[x])
        }
    })
    levels(scores)<-interval_names      
    dfMelt<-cbind(dfMelt, scores=scores)
    
    if(log){
        dfMelt["value"]<-log10(dfMelt["value"]+1)
        y_lab<-paste("log10(", attribute,"+1)", sep="")
    }
    value<-NULL
    g<-ggplot(dfMelt)
    x_lab<-capitalize(feature)

    #if dens=TRUE the boxplot and density plot are drawing together as a 
    #violin plot else only the bolxplot is displayed
    if(dens){
        g<-g+geom_violin(aes(x=names, y=value,fill=scores), alpha=0.5, 
            draw_quantiles = c(0.25, 0.5,0.75), trim=FALSE)+labs(
            title="", x=x_lab,y=y_lab)+theme(
            plot.title =element_text(size=rel(1.5), colour ="black",hjust=0.5), 
            title=element_text(size=22), axis.title=element_text(size=22),
            legend.text = element_text(size = 18))
        g<-g+guides(fill= guide_legend(title=paste(capitalize(attribute),
            "_intervals",  sep="")))
    }else{
        g<-g+geom_boxplot(aes(x=names, y=value,fill=scores))+labs(title="", 
            x=x_lab,y=y_lab)+theme(axis.title=element_text(size=14), 
            legend.text = element_text(size = 11), axis.text.x= element_text( 
            angle=90),plot.title =element_text(size=rel(1.5), colour ="black",
            hjust=0.5))+guides(fill=FALSE)
        g<-g+guides(fill=guide_legend(title=paste(capitalize(attribute), 
            "_intervals",  sep="")))
    }
    if(pool){
        g<-g+facet_grid(~ pool, scales="free_x", space="free")
        if( !is.null(medianMarg)){
            if(log){
            lowLine<-as.data.frame(do.call(rbind, lapply(unique(dfMelt[,
                "pool"]), function(i){
                    return(c(pool=i, median=log10(medianMarg[1]/100)+median(
                        dfMelt[dfMelt[,"pool"]==i, "value"])))
            })))
                        
            upLine<-as.data.frame(do.call(rbind, lapply(unique(dfMelt[,
                "pool"]), function(i){
                    return(c(pool=i, median=log10(medianMarg[2]/100)+median(
                        dfMelt[dfMelt[,"pool"]==i, "value"])))
            })))
            }else{
            lowLine<-as.data.frame(do.call(rbind, lapply(unique(dfMelt[,
                "pool"]), function(i){
                    return(c(pool=i, median=medianMarg[1]/100*median(dfMelt[
                        dfMelt[,"pool"]==i, "value"])))
            })))
            upLine<-as.data.frame(do.call(rbind, lapply(unique(dfMelt[,
                "pool"]), function(i){
                    return(c(pool=i, median=medianMarg[2]/100*median(dfMelt[
                        dfMelt[,"pool"]==i, "value"])))
            })))
            
            
            }
            g<-g+geom_hline(data = lowLine, aes(yintercept = median, color=
                as.factor(1)))+ geom_hline(data = upLine, aes(yintercept = 
                median, color=as.factor(2)))+ 
                scale_color_manual(name="median %", breaks=c(1, 2), labels=c(
                medianMarg[1], medianMarg[2]), values=c("red", "blue"))
        }
    }else{
        if( !is.null(medianMarg)){
            if(!log){
                lowLine<-medianMarg[1]/100*median(dfMelt[, "value"])
                upLine<-medianMarg[2]/100*median(dfMelt[, "value"])
            }else{
                lowLine<-log10(medianMarg[1]/100)+median(dfMelt[, "value"])
                upLine<-log10(medianMarg[2]/100)+median(dfMelt[, "value"])
            }
            g<-g+geom_hline(aes(yintercept = lowLine, color= as.factor(1)))+ 
                geom_hline(aes(yintercept = upLine, color=as.factor(2)))+ 
                scale_color_manual(name="median %", breaks=c(1, 2), labels=c(
                medianMarg[1], medianMarg[2]), values=c("red", "blue"))
        }
    }
    if (!featureLabs){
        g<-g+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
    }
    colors<-colorRampPalette(c("red", "green"))(length(interval_names))
    names(colors)<-interval_names
    g<-g+scale_fill_manual(name=paste(attribute, "interval", sep=" "),
        breaks=interval_names, values=colors)+guides(fill=
        guide_legend(title=paste(capitalize(attribute), "_intervals", sep="")))
    
    return(g)
})
