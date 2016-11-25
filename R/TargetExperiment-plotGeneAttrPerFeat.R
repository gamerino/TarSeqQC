#'Plot the attribute value for all the features of a selected gene.
#'
#'\code{plotGeneAttrPerFeat} plots the achieved performance for each feature
#'for a particular gene. The resulting graph shows one bar per each gene 
#'feature with heights according to its attribute value.
#'
#'@param object TargetExperiment object.
#'@param geneID Character indicating the ID of the selected gene.
#'@param overlap Logical indicating if the amplicons should be collapsed in 
#'overlapped regions. 
#'@param level Character indicating the level of the plot. Can be 'feature', to
#'plot the features' attribute; 'region', to plot overlapped regions'
#'attribute; 'both'to generate the two previos plots. 
#'@return ggplot2 graphics.
#'
#'@include TargetExperiment-plotFeature.R
#'@exportMethod plotGeneAttrPerFeat
#'@docType methods
#'@name plotGeneAttrPerFeat
#'@rdname TargetExperiment-plotGeneAttrPerFeat
#'@aliases plotGeneAttrPerFeat-methods
#'@seealso \code{\link{plotAttrExpl}}
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
#'@examples
#'## Loading the TargetExperiment object
#'data(ampliPanel, package="TarSeqQC")
#'
#'## Exploring amplicon attribute values for a particular gene
#'# Ignoring amplicon overlapping 
#'g<-plotGeneAttrPerFeat(ampliPanel, geneID="gene4")
#'# Adjust text size
#'g<-g+theme(title=element_text(size=16), axis.title=element_text(size=16),
#'legend.text=element_text(size=14))
#'if(interactive()){
#'g
#'}
#'# Considering amplicon overlapping 
#'g<-plotGeneAttrPerFeat(ampliPanel, geneID="gene4", overlap=TRUE, level="both")
#'# Adjust text size
#'g<-g+theme(title=element_text(size=16), axis.title=element_text(size=16),
#'legend.text=element_text(size=14))
#'if(interactive()){
#'g
#'}


setGeneric(name="plotGeneAttrPerFeat", def=function(object, geneID,
    overlap=FALSE, level="feature"){
        standardGeneric("plotGeneAttrPerFeat")
})
#'
#'@name plotGeneAttrPerFeat
#'@rdname TargetExperiment-plotGeneAttrPerFeat
#'@aliases plotGeneAttrPerFeat,TargetExperiment-method
#'@inheritParams plotGeneAttrPerFeat
#'
setMethod(f="plotGeneAttrPerFeat", signature=signature(
object= "TargetExperiment"), definition=function(object, geneID,
overlap=FALSE, level="feature"){
level<-match.arg(level, c("feature", "region","both"))
if(!overlap & level !="feature"){
    warning("When 'overlap' is set to FALSE, 'level' should be 'feature'")
    level<-"feature"
}

featurePanel<-getFeaturePanel(object)
if(!(geneID %in% mcols(featurePanel)[,"gene"])){
    stop(paste("The selected gene: ", geneID, ", does not exist in the bed 
        file",sep=""))
}
attribute<-getAttribute(object) 
if(!(attribute %in% c("coverage", "medianCounts"))){
    stop("Attribute slot should be defined in order to call the
        function")
}
df_panel<-as.data.frame(featurePanel[mcols(featurePanel)[,"gene"] ==geneID])
feature<-getFeature(object)
df_panel<-cbind(feature_id=row.names(df_panel), df_panel)
names(df_panel)[names(df_panel) == attribute]<-"attr"
pool<-"pool" %in% names(df_panel)

if(overlap){
    df_over<-getOverlappedRegions(object,collapse=FALSE)
    df_panel<-df_over[df_over[,"gene"]==geneID,]
    names(df_panel)[names(df_panel) == paste("region", attribute, sep="_")]<-
        "attr_region"
    dfover<-getOverlappedRegions(object,collapse=TRUE)
    dfover<-dfover[dfover[,"gene"]==geneID,]
    names(dfover)[names(dfover) == paste("region", attribute, sep="_")]<-
        "attr_region"
    lev<-do.call(rbind,strsplit(levels(df_panel[,"region_id"]), split=
                                    "Region"))[,2]
    df_panel[,"region_id"]<-as.character(df_panel[,"region_id"])
    df_panel[,"region_id"]<-do.call(rbind, strsplit(df_panel[,"region_id"], 
                                                    split="Region"))[,2]
    lev<-lev[lev %in% unique(df_panel[,"region_id"])]
    df_panel[,"region_id"]<-factor(df_panel[,"region_id"], levels=lev)
    
}   
feature_id<-NULL
if(pool & level == c("feature")){
    df_panel[,"feature_id"]<-factor(df_panel[,"feature_id"],levels=
        as.character(df_panel[,"feature_id"]))
    g<-ggplot(df_panel, aes(x=as.factor(feature_id), y=attr, fill=as.factor(
        pool)))+geom_bar(stat = "identity")+theme(
        axis.text.x = element_text(size = 18, colour = "black", angle = 90)
        )+scale_fill_hue(name = "Pool")+labs(x=capitalize(feature), 
        y=paste(capitalize(feature), attribute, sep =" "), 
        title=paste(capitalize(feature), attribute, "of ", geneID, "gene",
        sep=" "))+theme(plot.title=element_text(size=22,hjust=0.5), axis.title=
        element_text(size=22), legend.text = element_text(size = 18))
    if(overlap){
        g<-g+facet_grid(~region_id,scales= "free_x", space="free")
    }
}else{
    g<-ggplot(df_panel) 
    if(!overlap){
        g<-g+ geom_bar(aes(x=as.factor(feature_id), y=attr, fill= 
            as.factor(feature_id)),stat = "identity") + scale_fill_hue() +labs(
            x=capitalize(feature), y=paste(capitalize(feature), 
            attribute, sep =" "), title=paste(capitalize(feature), 
            attribute, "of ",geneID, "gene", sep=" ")) + theme(
            axis.text.x = element_text(size=18, colour= "black", angle =90
            ), plot.title=element_text(size=22,hjust=0.5), axis.title=
            element_text(size=22),legend.position="none")
    }else{
        region_id<-attr_region<-NULL
        g<-ggplot(df_panel) 
        g1<-g+ geom_bar(aes(x=as.factor(feature_id), y=attr, fill=
            region_id),stat = "identity")+facet_grid(~region_id,scales=
            "free_x", space="free") + scale_fill_hue() + labs(x=capitalize(
            feature), y=paste(capitalize(feature), attribute, sep =" "), title=
            paste(capitalize( feature), attribute, "of ", geneID, "gene",
            sep=" ")) + theme(axis.text.x = element_text(size=12, colour=
            "black", angle =90 ), plot.title=element_text(size=16,hjust=0.5), 
            axis.title= element_text(size=14,hjust=0.5),legend.position="none")
        
        g2<-ggplot(unique(dfover),aes(x=as.factor(region_id), y= 
            attr_region, fill=
            as.factor(region_id))) + geom_bar(stat = "identity") + theme(
            axis.text.x = element_text(size=16, colour= "black", angle = 
            90)) +  scale_fill_hue() +labs(x="Region", y=paste(
            "Region", attribute, sep =" "), title="") + theme(
            axis.title=element_text(size=16),legend.position="none")
            if(level=="both"){
                g<-cowplot::plot_grid(g1, g2, ncol=1, nrow=2)
            }else{
                if(level=="feature") g<-g1
                    else g<-g2
            }
    }
}
return(g)
})

