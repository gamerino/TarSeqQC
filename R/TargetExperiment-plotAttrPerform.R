#'Plot feature performance of a TargetExperiment object.
#'
#'\code{plotAttrPerform} plots the achieved performance for the selected 
#'attribute. The resulting graph shows one bar per each attribute interval and
#'its height is defined according to the amount of features achieving attribute
#'values within that interval.
#'
#'@param object TargetExperiment class object.
#'@param attributeThres Numeric indicating the intervals extreme values.
#'
#'@return ggplot2 graphics
#'
#'@include TargetExperiment-buildReport.R
#'@exportMethod plotAttrPerform
#'@docType methods
#'@name plotAttrPerform
#'@rdname TargetExperiment-plotAttrPerform
#'@aliases plotAttrPerform-methods
#'@seealso \code{\link{plot}}
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar} and Elmer A. Fernandez 
#'\email{efernandez@@bdmg.com.ar}
#'@examples
#'## Loading the TargetExperiment object
#'data(ampliPanel, package="TarSeqQC")
#'
#'# Definition of the interval extreme values
#'attributeThres<-c(0,1,50,200,500, Inf)
#'
#'# Plot panel overview in a feature performance plot
#'g<-plotAttrPerform(ampliPanel, attributeThres=attributeThres)
#'if(interactive()){
#'g
#'}
setGeneric(name="plotAttrPerform", def=function(object, 
attributeThres=c(0,1,50,200,500, Inf)){
    standardGeneric("plotAttrPerform")
})
#'
#'@name plotAttrPerform
#'@rdname TargetExperiment-plotAttrPerform
#'@aliases plotAttrPerform,TargetExperiment-method
#'@inheritParams plotAttrPerform

setMethod(f="plotAttrPerform", signature=signature(object="TargetExperiment"),
definition=function(object, attributeThres=c(0,1,50,200,500, Inf)){
    # definition of the data
    pool<-"pool" %in% names(mcols(getFeaturePanel(object)))
    df_panel<-summaryIntervals(object,attributeThres, pool)
    intervals<-cum_rel<-NULL
    if(is.data.frame(df_panel)){
        intervalName<-names(df_panel)[1]
        names(df_panel)[1]<-"intervals"
        df_panel[,"intervals"]<-factor(df_panel[,"intervals"], 
            levels=df_panel[,"intervals"])
        colors<-colorRampPalette(c("red", "green"))(length(levels(df_panel[,
            "intervals"])))
        names(colors)<-levels(df_panel[,"intervals"])    
        g<-ggplot(df_panel)+geom_bar(stat="identity", aes(x=intervals, 
            y=rel, fill=intervals))+geom_point(aes(x=intervals, y=cum_rel))+
            geom_line(aes(x=intervals, y=cum_rel, group=1, color=as.factor(
            "cumulative frequency")))+xlab(paste(intervalName))+ylab(
            "frequency (%)")+scale_fill_manual(values=colors, name=
            intervalName)+scale_color_manual(name="", values="black")
    
    
    }else{
        stopifnot(is.list(df_panel))
        poolNames<-names(df_panel)
        poolValues<-rep(poolNames[1], nrow(df_panel[[1]]))
        for(i in 2:length(poolNames)){
            poolValues<-c(poolValues, rep(poolNames[i],  nrow(df_panel[[i]])))
        }
        df_panel<-cbind(do.call(rbind, df_panel), poolValues)
        intervalName<-names(df_panel)[1]
        names(df_panel)[1]<-"intervals"
        df_panel[,"intervals"]<-factor(df_panel[,"intervals"], 
            levels=unique(df_panel[,"intervals"]))
        colors<-ggplotColours(object, n=(length(poolNames)))
        g<-ggplot(df_panel)+geom_point(aes(x=intervals, y=cum_rel, color=
            poolValues))+geom_line(aes(x=intervals, y=cum_rel, 
            group=poolValues, color=poolValues))+xlab(paste(intervalName))+
            ylab("frequency (%)")+scale_color_manual(name="pool", values=colors)
        }
    
    return(g)    
})
