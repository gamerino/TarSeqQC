#'Plot attribute exploration of a TargetExperiment/TargetExperimentList object.
#'
#'\code{plotAttrExpl} plots density and/or box-plot of the analyzed attribute 
#'at a feature level. These graphics could be displayed together using the 
#'ggplot2 geom_violin method. If panel's pools are present, one facet for each
#' pool  will be showed.
#'
#'@param object TargetExperiment/TargetExperimentList class object.
#'@param level Character 'feature' or 'gene' indicating at which level should 
#'be analyzed the attribute.
#'@param join Logical indicating if boxplot and density function should be 
#'plotted together using the ggplot2 geom_violin method.
#'@param log Logical indicating if the attribute should be considered in 
#'log10 scale.
#'@param color A character indicating a valid name color.
#'@param ... necessary arguments
#'
#'@return ggplot2 graphics.
#'
#'@include TargetExperiment-plotAttrPerform.R
#'@exportMethod plotAttrExpl
#'@docType methods
#'@name plotAttrExpl
#'@rdname TargetExperiment-plotAttrExpl
#'@aliases plotAttrExpl-methods
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
#'g<-plotAttrExpl(ampliPanel,level="feature",join=TRUE, log=FALSE, color="blue")
#'# x11(type="cairo")
#'if(interactive()){
#'g
#'}
setGeneric(name="plotAttrExpl", def=function(object,dens=FALSE,join=FALSE, 
log=TRUE, pool=FALSE,...){
    standardGeneric("plotAttrExpl")
})
#'
#'@name plotAttrExpl
#'@rdname TargetExperiment-plotAttrExpl
#'@importFrom cowplot plot_grid
#'@aliases plotAttrExpl,TargetExperiment-method
#'@inheritParams plotAttrExpl
setMethod(f="plotAttrExpl", signature=signature(object="TargetExperiment"),
definition=function(object,level="feature",join=TRUE, log=TRUE, color="blue"){
    if(!(level %in% c("feature", "gene"))){
        stop("'level' can be 'feature' or 'gene'")
    }
    #selecting the panel
    if(level=="feature"){
        df<-as.data.frame(getFeaturePanel(object))
        level<-getFeature(object)
    }else{
        df<-as.data.frame(getGenePanel(object))
    }
    pool<-"pool" %in% names(df)
    if(!(getAttribute(object) %in% c("coverage", "medianCounts"))){
        stop("Attribute slot should be defined in order to call the
            function")
    }
    y_lab<-capitalize(getAttribute(object))
    attribute<-getAttribute(object)
    # computting statistics
    mean_attr<-round(mean(df[, attribute]), digits=1)
    sd_attr<-round(sd(df[, attribute]), digits=1)
    median_attr<-round(median(df[, attribute]), digits=1)
    IQR_attr<-round(IQR(df[, attribute]), digits=1)
    #if log the plot will be in log10 scale
    if(log){
        df[, attribute]<-log10(df[, attribute]+1)
        y_lab<-paste("log10(", level,"_", attribute,"+1)", sep="")
    }else{
        y_lab<-paste(capitalize(level), attribute,sep=" ")
    }
    #if pool will appear one graph per pool
    if(pool) {
        df<-df[,c(attribute, "pool"), drop=FALSE]
        names(df)<-c("attribute", "pool")
        g<-ggplot( df,aes(x=as.factor(pool), y=attribute,fill=as.factor(pool)))
        dens.plot<-ggplot(df,aes(attribute,fill=as.factor(pool)))
        x_lab<-"Pool"
    }else {
        df<-df[, attribute, drop=FALSE]
        names(df)<-"attribute"
        g<-ggplot( df,aes(as.factor(1), attribute, fill=as.factor(1)))
        dens.plot<-ggplot(df,aes( attribute, fill=as.factor(1)))
        x_lab<-paste("Mean = ", mean_attr,",  sd = ", sd_attr,",  median = ", 
            median_attr, ",  IQR = ", IQR_attr, sep="")
    }
    #if join the boxplot and density plot are drawing together as a violin plot
    if(join){
        g<-g+geom_violin(alpha=0.5,draw_quantiles = c(0.25,0.5,0.75), 
            trim=FALSE)+labs(title=paste(capitalize(level), attribute, 
            sep=" "), x=x_lab, y=y_lab)+theme(plot.title =element_text(size=
            rel(1.5), colour="black", hjust=0.5),title=element_text(size=22),
            axis.title=element_text(size=22),legend.text = element_text(size =
            18))
        if(pool){
        g<-g+scale_fill_hue(name="Pool")
        }else{
        g<-g+guides(fill=FALSE)+scale_fill_manual(values=c(color))
        }
        return(g)
    }else{
        box.plot<-g+geom_boxplot()+labs(title=paste(capitalize(level), 
            attribute,sep=" "), x=x_lab,y=y_lab)+theme(plot.title = 
            element_text(size = rel(1.5), colour = "black", hjust=0.5), title=
            element_text(size=22), axis.title=element_text(size=22), 
            legend.text = element_text(size = 18))+guides(fill=FALSE)
        #marginal density of y - plot on the right
        dens.plot<- dens.plot+geom_density(alpha=0.5)+coord_flip()+ 
            labs(title="Density",x=y_lab, y="")+theme(plot.title = element_text(
            size = rel(1.5), colour = "black", hjust=0.5), title=element_text(
            size=22), axis.title=  element_text(size=22),legend.text = 
            element_text( size = 18))
        if(pool){
            dens.plot<-dens.plot+scale_fill_hue(name="Pool")
        }else{
            dens.plot<-dens.plot+guides(fill=FALSE)
            box.plot<-box.plot+scale_fill_manual(values=color)
            dens.plot<-dens.plot+scale_fill_manual(values=color)
        }
        #arrange the plots together, with appropriate height and width for each
        #row and column
        g<-cowplot::plot_grid(box.plot, dens.plot, ncol=2, nrow=1)
        return(g)    
    }
    
})
