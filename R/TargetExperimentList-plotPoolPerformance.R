#'Plot pool performance of a TargetExperimentList object.
#'
#'\code{plotPoolPerformance} plots density and/or box-plot of the analyzed
#'attribute achieved in each PCR pool. These graphics could be displayed 
#'together using the ggplot2 geom_violin method. 
#'
#'@param object TargetExperimentList class object.
#'@param dens Logical indicating if density plot should be included
#'@param join Logical indicating if boxplot and density function should be 
#'plotted together using the ggplot2 geom_violin method. For it uses, dens 
#'should be TRUE.
#'@param log Logical indicating if the attribute should be considered in 
#'log10 scale.
#'@param attributeThres Numeric indicating the attribute interval extreme 
#'values. It is not a mandatory parameter but if it is specified,then the 
#'plots will be colored according to the interval in which falls the attribute
#'median values.
#'
#'@return ggplot2 graphics.
#'
#'@include TargetExperimentList-plotAttrExpl.R
#'@exportMethod plotPoolPerformance
#'@docType methods
#'@name plotPoolPerformance
#'@rdname TargetExperimentList-plotPoolPerformance
#'@aliases plotPoolPerformance-methods
#'@examples
#'## Loading the TargetExperimentList object
#'data(TEList, package="TarSeqQC")
#'
#'# Attribute boxplot and density plot exploration
#'g<-plotPoolPerformance(TEList, log=FALSE)
#'# x11(type="cairo")
#'if(interactive()){
#'g
#'}
setGeneric(name="plotPoolPerformance", def=function(object, dens=FALSE,
join=FALSE,log=TRUE, attributeThres=NULL){
    standardGeneric("plotPoolPerformance")
})
#'
#'@name plotPoolPerformance
#'@rdname TargetExperimentList-plotPoolPerformance
#'@importFrom cowplot plot_grid
#'@aliases plotPoolPerformance,TargetExperimentList-method
#'@inheritParams plotPoolPerformance
setMethod(f="plotPoolPerformance", signature=signature(object=
"TargetExperimentList"),definition=function(object,dens=FALSE,join=FALSE,
log=TRUE, attributeThres=NULL){
    if( (!dens) & join ){
        stop("'dens' is FALSE and join is TRUE. Please if you want to plot
            density and box plots together set both dens and join as TRUE")
    }
    
    #selecting the panel
    df<-as.data.frame(getPanels(object))
    stopifnot(any(names(df) == "pool"))
    y_lab<-getAttribute(object)
    attribute<-getAttribute(object)
    #if log the plot will be in log10 scale
    index<-do.call(c, lapply(1:ncol(df), function(i){
        if(strsplit(colnames(df)[i], split="_")[[1]][[1]] == 
            attribute){
            return(i)}
    }))


    x_lab<-""
    #if pool will appear one graph per pool
    
    if(!is.null(attributeThres)){
        if(attributeThres[1] !=0){
            attributeThres<-c(0,attributeThres)
        }
        if(attributeThres[length(attributeThres)] !=Inf){
            attributeThres<-c(attributeThres, Inf)
        }
    interval_names<-sapply(1:length(attributeThres[attributeThres != 
            "Inf"]), function(x){
            if(x < length(attributeThres[attributeThres != "Inf"])) {
                return((paste(attributeThres[x], " <= ", attribute,
                    " < ", attributeThres[x+1])))
            }else{
                paste(  attribute, " >= ", attributeThres[x])
            }
        })

        # creating a 'score' variable to group features according to 
        #the attributeintervals
        medians<-as.data.frame(do.call(rbind,lapply(levels(as.factor(df[,
            "pool"])), function(x){
                data<-as.numeric(as.matrix(df[which(df[,
                    "pool"] == x),index]))
                return(c(pool=x, median=median(data)))
            
        })))
        medians[,"median"]<-as.numeric(as.character(medians[,"median"]))
        scores<-cut(medians[,"median"], breaks=attributeThres, 
            include.lowest=TRUE,   right=FALSE, dig.lab = 6)
        levels(scores)<-interval_names
        medians[,"scores"]<-scores
        df$scores<-df$pool
        df$scores<-medians[match(df$pool, medians[,"pool"]), 
            "scores"]

    }
    if(log){
        df[, index]<-log10(df[, index]+1)
        y_lab<-paste("log10(", attribute,"+1)", sep="")
    }else{
        y_lab<-paste(attribute,sep=" ")
    }
    df<-df[,c(index, which("pool" == names(df)), which("scores" == names(df))),
        drop=FALSE]
    if(is.null(attributeThres)){
        dfMelt<-melt(df, id.vars="pool")
    }else{
        dfMelt<-melt(df, id.vars=c("pool", "scores"))
    }
            
    dfMelt[,"pool"]<-as.factor(dfMelt[,"pool"])
    pool<-variable<-value<-NULL
    g<-ggplot(dfMelt,aes(x=pool, y=value,fill=pool))
    if(!is.null(attributeThres)){
        g<-ggplot(dfMelt,aes(x=pool, y=value,fill=scores))
    }

    if(dens){
        dens.plot<-ggplot(dfMelt,aes(value,fill=pool))
        if(!is.null(attributeThres)){
            dens.plot<-ggplot(dfMelt,aes(value,fill=scores, group=pool))
        }

        #if join the boxplot and density plot are drawing together as a violin
        #plot
        if(join){
            g<-g+geom_violin( alpha=0.5)+geom_boxplot( width=0.2)+labs(
                x=x_lab,y=y_lab)+theme(legend.text = element_text(size = 18),
                legend.title = element_text(size = 18))
            if(!is.null(attributeThres)){
                colors<-colorRampPalette(c("red", "green"))(length(
                    interval_names))
                names(colors)<-interval_names
                g<-g+scale_fill_manual(name=paste(attribute, "interval", 
                    sep=" "), breaks=interval_names, values=colors)+theme(
                    axis.text=element_text(size=16),axis.ticks=element_line()
                    )+labs(x="pool")
            }else{
                g<-g+scale_fill_discrete("pool")+ theme(axis.text=
                    element_text(size=16), axis.text.x=element_blank(), 
                    axis.ticks.x=element_blank(), axis.title=element_text(
                    size=22))
            }    

        }else{
            box.plot<-g+geom_boxplot()+labs(x=x_lab,y=y_lab)+theme( 
                axis.title=element_text(size=22), legend.text = element_text(
                size = 18),axis.text=element_text(size=16))+guides(fill=FALSE)
            #marginal density of y - plot on the right
            dens.plot<- dens.plot+geom_density(alpha=0.5)+coord_flip()+ 
                labs(x=y_lab, y="")+theme(axis.title=  element_text(size=22),
                legend.text = element_text(size = 18), legend.title=
                element_text(size=18),axis.text=element_text(size=16))
            if(!is.null(attributeThres)){
                colors<-colorRampPalette(c("red", "green"))(length(
                    interval_names))
                names(colors)<-interval_names  
                dens.plot<-dens.plot+scale_fill_manual(name=paste(attribute, 
                    "interval", sep=" "), breaks=interval_names, values=colors)
                box.plot<-box.plot+labs(x="pool")+scale_fill_manual(values=
                    colors)

            }else{
                dens.plot<-dens.plot+scale_fill_discrete("pool")
            }    
            #arrange the plots together, with appropriate height and width for
            # each row and column
            g<-cowplot::plot_grid(box.plot, dens.plot, ncol=2, nrow=1)
        }
    }else{
        g<-g+geom_boxplot()+labs(x=x_lab,y=y_lab)+theme(
            axis.title=element_text(size=22), legend.text = element_text(
            size = 18),axis.text=element_blank(), axis.ticks=element_blank(),
            legend.title=element_text(size=18))
        if(!is.null(attributeThres)){
            colors<-colorRampPalette(c("red", "green"))(length(interval_names))
            names(colors)<-interval_names
            g<-g+scale_fill_manual(name=paste(attribute, "interval", sep=" "),
                breaks=interval_names, values=colors)+theme(axis.text= 
                element_text(size=16), axis.ticks=element_line(),axis.title=
                element_text(size=22))+labs(x="pool")
        }else{
            g<-g+scale_fill_discrete("pool")
        }    

        
    }
    return(g)
    
})
