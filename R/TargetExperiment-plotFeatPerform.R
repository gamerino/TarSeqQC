#'Plot feature performance of a TargetExperiment object.
#'
#'\code{plotFeatPerform} plots the achieved performance for each feature/gene.
#'The resulting graph shows one bar per each feature/gene with height 
#'according to its attribute value. If complete is set as TRUE, two bar plots
#'(feature and gene level) will be stored in the resulting ggplot object.
#'
#'@param object TargetExperiment class object.
#'@param attributeThres Numeric indicating the intervals extreme values.
#'@param complete Logical indicating if the gene and feature level exploration
#' should be plotted.
#'@param log Logical indicating if the attribute should be considered in 
#'log10 scale.
#'@param featureLabs Logical indicating if feature labels should be plotted.
#'@param sepChr Logical indicating if the plot should show chromosome 
#'divisions.
#'@param legend Logical indicating if legend should be plotted. 
#'
#'@return ggplot2 graphics
#'
#'@include TargetExperiment-plotAttrExpl.R
#'@exportMethod plotFeatPerform
#'@docType methods
#'@name plotFeatPerform
#'@rdname TargetExperiment-plotFeatPerform
#'@aliases plotFeatPerform-methods
#'@seealso \code{\link{plot}}
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
#'@examples
#'## Loading the TargetExperiment object
#'data(ampliPanel, package="TarSeqQC")
#'
#'# Definition of the interval extreme values
#'attributeThres<-c(0,1,50,200,500, Inf)
#'
#'# Plot panel overview in a feature performance plot
#'g<-plotFeatPerform(ampliPanel, attributeThres=attributeThres, log=FALSE, 
#'featureLabs=TRUE, sepChr=TRUE, legend=TRUE)
#'if(interactive()){
#'g
#'}
setGeneric(name="plotFeatPerform", def=function(object, 
attributeThres=c(0,1,50,200,500, Inf), complete=TRUE, log=TRUE, 
featureLabs=FALSE, sepChr=FALSE, legend=TRUE){
    standardGeneric("plotFeatPerform")
})
#'
#'@name plotFeatPerform
#'@rdname TargetExperiment-plotFeatPerform
#'@aliases plotFeatPerform,TargetExperiment-method
#'@inheritParams plotFeatPerform
#'@importFrom cowplot plot_grid
setMethod(f="plotFeatPerform", signature=signature(object="TargetExperiment"),
definition=function(object, attributeThres=c(0,1,50,200,500, Inf),
complete=TRUE, log=TRUE, featureLabs=FALSE, sepChr=FALSE, legend=TRUE){
    # definition of the data 
    if(attributeThres[1] !=0){
        attributeThres<-c(0,attributeThres)
    }
    if(attributeThres[length(attributeThres)] !=Inf){
        attributeThres<-c(attributeThres, Inf)
    }
    if(!(getAttribute(object) %in% c("coverage", "medianCounts"))){
        stop("Attribute slot should be defined in order to call the
            function")
    }
    df_panel<-cbind(as.data.frame(getFeaturePanel(object)), feature=names(
        getFeaturePanel(object)))
    names(df_panel)[names(df_panel) == getAttribute(object)] <- "attribute"
    df_panel[,"score"]<-cut(df_panel[,"attribute"], breaks=attributeThres, 
        include.lowest=TRUE, right=FALSE,dig.lab = 6)
    score_levels<-levels(df_panel[,"score"])
    pool<-"pool" %in% names(df_panel)
    attribute<-getAttribute(object)
    feature<-getFeature(object)
#     colors<-ggplotColours(object, n=(length(attributeThres)-1))
    
    # if complete then feature and gene panel should be plotted
    if (complete) {
        data2<-ddply(df_panel, c("seqnames","gene"), summarize, 
            attribute=mean(attribute))
        data2[,"score"]<-cut(data2[,"attribute"], breaks=attributeThres, 
            include.lowest=TRUE, right=FALSE,dig.lab = 6)

    }
    # if log the attribute scale will be in log10
    if(log){
        df_panel[, "attribute"]<-log10(df_panel[, "attribute"]+1)
        y_lab1<-paste("log10(", feature, "_",attribute,"+1)", sep="")
        y_lab2<-paste("log10(gene_",attribute,"+1)", sep="")
        if(complete){
            data2[, "attribute"]<-log10(data2[, "attribute"]+1)
        }
        
    }else{
        y_lab1<-paste(capitalize(feature),attribute,sep=" ")
        y_lab2<-paste("Gene",attribute,sep=" ")
    }  
    data1<-df_panel[,c("seqnames","feature", "attribute", "score")]
    data1<-arrange(data1,seqnames,feature)
    # plotting gene level
    if(complete){
        x_lab2<-"Gene"
        gene<-score<-NULL
        g2<-ggplot(data2,aes(x=as.factor(gene), y=attribute, 
            fill=score)) + geom_bar(stat="identity") + 
            facet_grid(~ seqnames, scales="free_x", space="free") + theme( 
            axis.text.x = element_text(angle = 90,size=10))
            colors<-colorRampPalette(c("red", "green"))(length(score_levels))
            names(colors)<-score_levels
            g2<-g2+scale_fill_manual(
            values=colors, breaks=score_levels)+ 
            geom_hline(yintercept=mean(data2[, "attribute"]), colour="red") + 
            labs(title="", x=x_lab2,y=y_lab2)
        if(legend) g2<-g2+guides(fill=guide_legend(title=paste(capitalize(
            attribute), "_intervals", sep="")))
    }
    # plot feature level
    # if pool the features ar grouped by pool
    if(pool) {
        data1<-cbind(data1, pool=df_panel[,"pool"])
        x_lab1<-"Pool"
        g1<-ggplot(data1,aes(x=as.factor(feature), y= attribute, 
            fill=score))+geom_bar(stat="identity")+facet_grid(~ pool,
            scales="free_x", space="free")+ theme(axis.text.x = element_blank(),
            axis.title.x=element_blank())+scale_y_continuous(limits = c(0,
            max(data1[,"attribute"])))+labs(title="", x=x_lab1,y=y_lab1) + 
            geom_hline(yintercept=mean(data1[,"attribute"]), colour="red", 
            show.legend=TRUE)
    }else {
        g1<-ggplot( data1,aes(x=as.factor(feature), y= attribute, 
            fill=score))+geom_bar(stat="identity")+ 
            scale_y_continuous(limits = c(0, max(data1[,"attribute"])))+labs(
            title="", y=y_lab1)+geom_hline(yintercept=mean(data1[,"attribute"]),
            colour="red", show.legend=TRUE)
    }
    if(featureLabs){
        g1<-g1+theme( axis.title.x=element_blank(),axis.text.x = element_text(
            angle = 90,size=8))
    }else{
        g1<-g1+theme(axis.text.x = element_blank(), axis.title.x=element_blank(
        ))
    }
    colors<-colorRampPalette(c("red", "green"))(length(score_levels))
    names(colors)<-score_levels
    g1<-g1+scale_fill_manual(values=colors, breaks=score_levels)
    if(legend){
        g1<-g1+guides(fill=guide_legend(title=paste(capitalize(attribute),
            "_intervals", sep="")))
    }
    # arrange the plots together, with appropriate height and width for each 
    # row and column
    if(complete){
        g<-cowplot::plot_grid(g1,g2, ncol=1, nrow=2)
    }else {
        if(sepChr){
            g<-g1+facet_grid(~ seqnames, scales="free_x", space="free")
        }else g<-g1
    }
    return(g)    
})
