#'Plot attribute density and boxplot for each bias source quartile or 
#'category.
#'
#'\code{biasExploration} plots density and box-plot of the analyzed attribute 
#'for each bias source' quartiles per categories. It helps the identification 
#'of some bias due to high source values, for example, high gc content. This 
#'graphics could plot together using the ggplot2 geom_violin method.
#'
#'@param object TargetExperiment class object.
#'@param source Character 'gc','length', or 'pool' indicating the source bias. 
#'In the case of 'gc' and 'length', it will be categorized in four groups 
#'according to its quartiles. In the case of 'pool', its groups will be
#'conserved. 
#'@param dens Logical indicating if density plot should be added using the 
#'geom_violin ggplot2 method.
#'
#'@return ggplot2 graphics.
#'
#'@include TargetExperiment-plotInOutFeatures.R
#'@exportMethod biasExploration
#'@docType methods
#'@name biasExploration
#'@rdname TargetExperiment-biasExploration
#'@aliases biasExploration-methods
#'@seealso \code{\link{plot}}, \code{\link{plotFeatPerform}}
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
#'@examples
#'## Loading the TargetExperiment object
#'data(ampliPanel, package="TarSeqQC")
#'
#'# Attribute boxplot and density plot exploration
#'g<-biasExploration(ampliPanel,source="gc", dens=TRUE)
#'# x11(type="cairo")
#'if(interactive()){
#'g
#'}
setGeneric(name="biasExploration", def=function(object, source=c("length", 
    "gc","pool"), dens=FALSE){
    standardGeneric("biasExploration")
})
#'
#'@name biasExploration
#'@rdname TargetExperiment-biasExploration
#'@importFrom cowplot plot_grid
#'@aliases biasExploration,TargetExperiment-method
#'@inheritParams biasExploration
setMethod(f="biasExploration", signature=signature(object="TargetExperiment"),
definition=function(object,source=c("length", "gc", "pool"), dens=FALSE){
    myDF<-as.data.frame(getFeaturePanel(object))
    if (source %in% c("length", "gc", "pool")){
        source<-match.arg(source)
    }else{
        source<-source
        if(!(source %in% colnames(myDF))){
            stop("'source' should be one of 'length', 'gc', 'pool' or a valid 
                meta data name")
        }
    }
    attribute<-getAttribute(object)
    if(!(attribute %in% c("coverage", "medianCounts"))){
        stop("Attribute slot should be defined in order to call biasExploration
            function")
    }
    colnames(myDF)[colnames(myDF) =="width"]<-"length"
    if(source == "pool"){
        myDF[,"pool"]<-as.factor(myDF[,"pool"])
    }
    # select only one bias source
    myDF<-myDF[,colnames(myDF)[colnames(myDF) %in% c(attribute, source)]]
    colnames(myDF)[colnames(myDF) == attribute]<-"attribute"
    colnames(myDF)[colnames(myDF) == source]<-"source"
    if(is.character(myDF[,"source"])){
        myDF[,"source"]<-factor(myDF[,"source"], levels=unique(myDF[,"source"]
            ))
    }
    # built source intervals according to the four quartiles
    if(is.factor(myDF[,"source"])){
        colnames(myDF)[colnames(myDF) == "source"] <-"intervals"
    }else{
        myDF[,"intervals"]<-cut(myDF[,"source"],breaks=summary(myDF[,
            "source"])[c(1:3, 5:6)],include.lowest=TRUE, dig.lab = 6)
        
    }
    intervals<-NULL
    g<-ggplot(myDF,aes(x=intervals, y=attribute, fill=intervals))
    if(dens){
        g<-g+geom_violin(alpha=0.5,draw_quantiles = c(0.25, 0.5,0.75),
            trim=FALSE) 
    }else{
        g<-g+geom_boxplot()
    }
    labName<-source
    if(labName=="gc") labName<-"GC"
    g<-g+labs(x = "", y = capitalize(attribute))+ theme(axis.title=
        element_text(size=22), legend.text =  element_text(size=16),
        legend.title=element_text(size=18)) +scale_fill_hue(name=paste(
        capitalize(labName), "groups", sep=" "))
    return(g)    
    
})
