#'Graphical exploration of a specific metadata column.
#'
#'\code{plotMetaDataExpl} plots density and box-plot of a specific metadata 
#'column. If the characteristic is nonnumerical, then  a frequency plot is 
#'built.
#'@param object TargetExperiment class object.
#'@param name a character indicating the metadata column name that should
#'be analyzed.
#'@param join Logical only for numerical variables. It indicates if boxplot and
#'density function should be plotted together using the ggplot2 geom_violin 
#'method.
#'@param log Logical indicating if the numerical metadata column should be 
#'considered in log10 scale.
#'@param absolute Logical indicating if the frequencies of the selected 
#'categorical metadata column should be in absolute scale. If absolute is FALSE
#'the frequencies are in relative percentages.
#'@param color A character indicating a valid name color.
#'
#'@return ggplot2 graphics.
#'
#'@include TargetExperiment-biasExploration.R
#'@exportMethod plotMetaDataExpl
#'@docType methods
#'@name plotMetaDataExpl
#'@rdname TargetExperiment-plotMetaDataExpl
#'@aliases plotMetaDataExpl-methods
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
#'g<-plotMetaDataExpl(ampliPanel, name="length")
#'if(interactive())
#'{
#'# x11(type="cairo")
#'g
#'}
#'# Explore amount of amplicons per gene
#'g<-plotMetaDataExpl(ampliPanel, name="gene", absolute=TRUE)
#'if(interactive())
#'{
#'# x11(type="cairo")
#'g
#'}
setGeneric(name="plotMetaDataExpl", def=function(object,name=c("length", "gc",
"pool"), log=FALSE,join=TRUE, absolute=FALSE, color="blue"){
    standardGeneric("plotMetaDataExpl")
})
#'
#'@name plotMetaDataExpl
#'@rdname TargetExperiment-plotMetaDataExpl
#'@importFrom cowplot plot_grid
#'@aliases plotMetaDataExpl,TargetExperiment-method
#'@inheritParams plotMetaDataExpl
setMethod(f="plotMetaDataExpl", signature=signature(object="TargetExperiment"),
definition=function(object,name=c("length", "gc", "pool"),log=FALSE, join=TRUE,
absolute=FALSE, color="blue"){
    myDF<-as.data.frame(getFeaturePanel(object))
    if(any(names(myDF) == "pool")){
        myDF[,"pool"]<-as.factor(myDF[,"pool"])
    }
    if (name %in% c("length", "gc", "pool")){
        source<-match.arg(name)
    }else{
        if(name %in% colnames(myDF)){
            source<-name 
        }else {
            stop("name should be one of 'length', 'gc', 'pool' or a valid meta 
            data name")
        }
    }
    x_lab<-paste(capitalize(source))
    
    if(source == "length") source<-"width"
    colnames(myDF)[colnames(myDF) == source]<-"source"
    if(is.factor(myDF[,"source"]) & log) {
        warning("'log' specification is valid only for numerical meta data 
            information")
    }
    if(!is.factor(myDF[,"source"]) & !is.character(myDF[, "source"])){
        mean_attr<-round(mean(myDF[, "source"]), digits=1)
        sd_attr<-round(sd(myDF[, "source"]), digits=1)
        median_attr<-round(median(myDF[, "source"]), digits=1)
        IQR_attr<-round(IQR(myDF[, "source"]), digits=1)
        stat_lab<-paste("Mean = ", mean_attr,",  sd = ",sd_attr,",  median = ",
            median_attr, ",  IQR = ", IQR_attr, sep="")
        if(log){
            myDF[, "source"]<-log10(myDF[, "source"]+1)
#             y_lab<-paste("log10(",x_lab,"+1)", sep="")
        }else{
            y_lab<-paste(x_lab)
        }
        g<-ggplot(myDF)
        if(join){

            p<-g+geom_violin(aes(x=as.factor(1), y=source, fill=as.factor(1)), 
                alpha=0.5, draw_quantiles = c(0.25, 0.5,0.75),trim=FALSE)+ 
                labs(title=paste(x_lab, "distribution", sep=" "), x=stat_lab, 
                y=y_lab)+theme(plot.title = element_text(size=rel(1.5), colour=
                "black", hjust = 0.5), title= element_text(size=14),axis.title= 
                element_text(size=13), legend.text = element_text(size = 13))+
                guides(fill=FALSE) +  scale_fill_manual(values= color)+ 
                scale_x_discrete(breaks=c(1),  labels=c(""))
        }else{
            box.plot<-g+geom_boxplot(aes(x=1, y=source, fill=as.factor(1)))+
                labs(title= "", y=y_lab,x= stat_lab)+ theme(plot.title = 
                element_text(size = rel(1.5), colour = "black",hjust = 0.5),
                title=element_text(size=14),axis.title= element_text(size=13), 
                legend.position="none")+guides(fill=FALSE) + scale_fill_manual(
                values=color) + scale_x_continuous(breaks=1,labels="")
            dens.plot<- g+geom_density(aes(x=source, fill=as.factor(1)), 
                alpha=0.5)+ coord_flip()+ labs(title="", x=y_lab, y="Density")+
                theme(plot.title = element_text(size = rel(1.5), colour = 
                "black",hjust=0.5), title= element_text(size=14), axis.title= 
                element_text(size=13), legend.position = "none")+
                scale_fill_manual(values=color)
            p<-cowplot::plot_grid(box.plot, dens.plot, ncol=2, nrow=1)
        }
    }else{
        DF<-as.data.frame(table(myDF[,"source"]))
        if (absolute){
            names(DF)[names(DF) == "Freq"]<-"PorcFreq"
            y_lab<-"Absoulte frequency"
        }else{
            DF<-cbind(DF, PorcFreq=round(DF[,"Freq"]/sum(DF[,"Freq"])*100,2) )
            y_lab<-"Relative frequency (%)"
        }
        Var1<-PorcFreq<-NULL
        p<-ggplot(DF, aes(x=Var1, y=PorcFreq, fill=as.factor(1))) +geom_bar(
        stat="identity")+ labs(title=paste(x_lab, "frequencies", sep=" "),x=
        x_lab, y=y_lab)+ theme(plot.title =  element_text(size = rel(1.5), 
        colour = "black",hjust=0.5), title=element_text(size=14), axis.title=  
        element_text(size=13), legend.position = "none")+scale_fill_manual(
        values=color)    
    }
    
    return(p)    
    
})
