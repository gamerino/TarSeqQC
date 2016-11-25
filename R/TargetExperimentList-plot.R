#'Plot TargetExperimentList object overview.
#'
#'For TargetExperimentList objects, \code{plot} allows a fast and simple 
#'representation of several feature panels using a heatmap plot. Along the
#'x-axis, the features are represented and patients/samples along the y-axis. 
#'Finally, each cell is colored according to the attribute interval.
#'
#'@param pool Logical indicating if the plots should be performed for each
#'pool separately
#'@param sampleLabs Logical. Sample names must be plotted?.
#'@param featureLabs Logical. Feature names must be plotted?.
#'@param ... not used but necessary for redefining the generic function.
#'@include TargetExperimentList-statistics.R
#'@docType methods
#'@name plot
#'@rdname TargetExperiment-plot
#'@import ggplot2
#'@importFrom reshape2 melt
#'@aliases plot,TargetExperimentList,plot.TargetExperimentList
#'@note see full example in \code{\link{TargetExperimentList-class}}
#'@examples
#'## Loading the TargetExperimentList object
#'data(TEList, package="TarSeqQC")
#'# Definition of the interval extreme values
#'attributeThres<-c(0,1,50,200,500, Inf)
#'
#'## Plot panel overview
#'g<-plot(TEList, attributeThres=attributeThres, featureLabs =TRUE)
#'if(interactive()){
#'g
#'}
#'@export plot.TargetExperimentList
plot.TargetExperimentList<-function(x, y, attributeThres=
c(0, 1, 50, 200, 500, Inf),pool=FALSE, sampleLabs=TRUE, 
featureLabs=FALSE, ...){
    if(attributeThres[1] !=0){
        attributeThres<-c(0,attributeThres)
    }
    if(attributeThres[length(attributeThres)] !=Inf){
        attributeThres<-c(attributeThres, Inf)
    }
    df_panel<-as.data.frame(getPanels(x))
    df_panel[,"names"]<-rownames(df_panel)
    attribute<-getAttribute(x)
    if(pool & !(any(names(df_panel) == "pool"))){
        stop("pool was set as TRUE but the object doesn't contain pool 
            information")
    }
    index<-do.call(c, lapply(1:ncol(df_panel), function(i){
        if(strsplit(colnames(df_panel)[i], split="_")[[1]][[1]] == 
            attribute){
            return(i)}
    }))
    listNames<-do.call(c,lapply(colnames(df_panel[,index]), function(j){
        return(strsplit(j, split=paste(attribute, "_", sep="")
            )[[1]][2])
    }))
    interval_names<-sapply(1:length(attributeThres[attributeThres != "Inf"]),
    function(x){
        if(x < length(attributeThres[attributeThres != "Inf"])) {
            return((paste(attributeThres[x], " <= ", attribute,
                " < ", attributeThres[x+1])))
        }else{
            paste(  attribute, " >= ", attributeThres[x])
        }
    })
    interval_names<-cbind(interval_names, nmb=1:length(interval_names))

    # creating a 'score' variable to group features according to the attribute
    #'intervals
    scores<-as.data.frame(do.call(cbind,lapply(1:length(index), function(i){
        return(cut(df_panel[,index[i]], 
            breaks=attributeThres, include.lowest=TRUE, right=FALSE, 
            dig.lab = 6))
    
    })))
    
    scores<-as.data.frame(do.call(cbind, lapply(1:ncol(scores), function(i){
        return(interval_names[match(scores[,i], interval_names[,"nmb"]),
            "interval_names"])
    
    })))

    colnames(scores)<-paste(listNames)
    score_levels<-interval_names[,"interval_names"]
    df_panel<-cbind(df_panel, scores)
    
    for (i in which(names(df_panel) %in% listNames)){
        df_panel[,i]<-factor(df_panel[,i], levels=score_levels)
    }
    dfMelt<-melt(df_panel, id.vars=which(names(df_panel) %in% c(
        "names", "pool")), measure.vars =which(names(df_panel) %in% names(
        scores)))
        
    dfMelt[,"value"]<-factor(dfMelt[,"value"], levels=interval_names[,
        "interval_names"])
    dfMelt[,"names"]<-factor(dfMelt[,"names"], levels=rownames(df_panel))
    names<-variable<-value<-NULL
    g<-ggplot(dfMelt, aes(x=names, y=variable))+geom_tile(aes(fill=value))

    colors<-colorRampPalette(c("red", "green"))(length(interval_names[,
        "interval_names"]))
    names(colors)<-interval_names[,"interval_names"]
    g<-g+scale_fill_manual(name=paste(attribute, "interval", sep=" "), breaks=
        interval_names[,"interval_names"], values=colors)+guides(fill=
        guide_legend(title=paste( capitalize(attribute), "_intervals", sep="")))
    if(!featureLabs){
        g<-g+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
    }else{
        g<-g+theme(axis.text.x=element_text(angle=90))
    }
    if(!sampleLabs){
        g<-g+theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
    }
    g<-g+labs(x=capitalize(getFeature(x)), y= "Samples")
    if (pool){
    g<-g+facet_grid(~ pool, scales="free")
    }
    g
}
#'@S3method plot TargetExperimentList
## S4 method dispatches to S3
setMethod("plot", "TargetExperimentList", plot.TargetExperimentList)
