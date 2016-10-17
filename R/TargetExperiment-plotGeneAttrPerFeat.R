#'Plot the attribute value for all the features of a selected gene.
#'
#'\code{plotGeneAttrPerFeat} plots the achieved performance for each feature
#'for a particular gene. The resulting graph shows one bar per each gene feature
#'with heights according to its attribute value.
#'
#'@param object TargetExperiment object.
#'@param geneID Character indicating the ID of the selected gene.
#'
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
#'# Exploring amplicon attribute values for a particular gene
#'g<-plotGeneAttrPerFeat(ampliPanel, geneID="gene4")
#'# Adjust text size
#'g<-g+theme(title=element_text(size=16), axis.title=element_text(size=16),
#'legend.text=element_text(size=14))
#'if(interactive()){
#'g
#'}

setGeneric(name="plotGeneAttrPerFeat", def=function(object,geneID){
    standardGeneric("plotGeneAttrPerFeat")
})
#'
#'@name plotGeneAttrPerFeat
#'@rdname TargetExperiment-plotGeneAttrPerFeat
#'@aliases plotGeneAttrPerFeat,TargetExperiment-method
#'@inheritParams plotGeneAttrPerFeat
#'
setMethod(f="plotGeneAttrPerFeat", signature=signature(
object= "TargetExperiment"), definition=function(object, geneID){
    
    featurePanel<-getFeaturePanel(object)
    if(!(geneID %in% mcols(featurePanel)[,"gene"])){
        stop(paste("The selected gene: ", geneID, ", does not exist in the bed 
            file",sep=""))
    }
    df_panel<-as.data.frame(featurePanel[mcols(featurePanel)[,"gene"] ==geneID])
    
    attribute<-getAttribute(object)
    feature<-getFeature(object)
    df_panel<-cbind(feature_id=row.names(df_panel), df_panel)
    names(df_panel)[names(df_panel) == attribute]<-"attr"
    pool<-"pool" %in% names(df_panel)
    feature_id<-NULL
    if(pool){
        g<-ggplot(df_panel, aes(x=as.factor(feature_id), y=attr, fill=as.factor(
            pool)))+geom_bar(stat = "identity")+theme(
            axis.text.x = element_text(size = 18, colour = "black", angle = 60)
            )+scale_fill_hue(name = "Pool")+labs(x=capitalize(feature), 
            y=paste(capitalize(feature), attribute, sep =" "), 
            title=paste(capitalize(feature), attribute, "of ", geneID, "gene",
            sep=" "))+ scale_x_discrete(limits=as.character(df_panel[,
            "feature_id"]))+theme(title=element_text(size=22), axis.title=
            element_text(size=22), legend.text = element_text(size = 18))
    }else{
        g<-ggplot(df_panel,aes(x=as.factor(feature_id), y=attr, fill=as.factor(
            feature_id))) + geom_bar(stat = "identity") + theme(
            axis.text.x = element_text(size=18, colour= "black", angle = 60)) +
            scale_fill_hue() +labs(x=capitalize(feature), y=paste(capitalize(
            feature), attribute, sep =" "), title=paste(capitalize( feature),
            attribute, "of ", geneID, "gene", sep=" ")) +
            scale_x_discrete(limits=as.character(df_panel[,"feature_id"]))+
            theme(title=element_text(size=22), axis.title=element_text(size=22),
            legend.position="none")
    }
    return(g)
})

