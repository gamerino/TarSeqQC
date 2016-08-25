#'@param dens Logical indicating if density plot should be included
#'@param pool Logical indicating if plots should be displayed for each pool
#'separately
#'@param attributeThres Numeric indicating the attribute interval extreme 
#'values. It is not a mandatory parameter but if it is specified,then the 
#'plots will be colored according to the interval in which falls the attribute
#'median values.
#'
#'@include TargetExperimentList-plotGlobalAttrExpl.R
#'@docType methods
#'@name plotAttrExpl
#'@importFrom cowplot plot_grid
#'@aliases plotAttrExpl,TargetExperimentList-method
#'@inheritParams plotAttrExpl
#'@rdname TargetExperiment-plotAttrExpl
#'@examples
#'## Loading the TargetExperimentList object
#'data(TEList, package="TarSeqQC")
#'
#'# Attribute boxplot and density plot exploration
#'g<-plotAttrExpl(TEList, log=FALSE, pool=FALSE)
#'# x11(type="cairo")
#'if(interactive()){
#'g
#'}
#'
setMethod(f="plotAttrExpl", signature=signature(object="TargetExperimentList"),
definition=function(object,dens=FALSE,join=FALSE, log=TRUE, pool=FALSE,
attributeThres=NULL){
    if( dens==FALSE & join ){
        stop("'dens' is FALSE and join is TRUE. Please if you want to plot
            density and box plots together set both dens and join as TRUE")
    }
    #selecting the panel
    df<-as.data.frame(getPanels(object))
    if (pool & !(any(names(df) == "pool"))){
        stop("'pool' is TRUE but panels do not contain pool information")
    } 
    y_lab<-getAttribute(object)
    attribute<-getAttribute(object)
    #if log the plot will be in log10 scale
    index<-do.call(c, lapply(1:ncol(df), function(i){
        if(strsplit(colnames(df)[i], split="_")[[1]][[1]] == 
            attribute){
            return(i)}
    }))
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
        if(pool){
            poolIndex<-which(colnames(df)=="pool")
            medians<-as.data.frame(do.call(rbind,lapply(index, function(x){
                Sample<-strsplit(colnames(df[,x,drop=FALSE]), split=paste(
                    attribute, "_",sep=""))[[1]][2]
                res<-do.call(rbind,lapply(unique(df[,poolIndex]), function(j){
                    return(c(pool=j, median=median(df[df[,"pool"]==j, x])))}))
                return(cbind(Sample=Sample, res=res))
            })))
            medians[,"median"]<-as.numeric(as.character(medians[,"median"]))

        }else{
            medians<-as.data.frame(do.call(rbind,lapply(index, function(x){
                Sample<-strsplit(colnames(df[,x,drop=FALSE]), split=paste(
                    attribute, "_",sep=""))[[1]][2]
                return(c(Sample=Sample, median=median(df[,x])))
            
            })))
            medians[,"median"]<-as.numeric(as.character(medians[,"median"]))
            
        }
        scores<-cut(medians[,"median"], breaks=attributeThres, 
                include.lowest=TRUE,   right=FALSE, dig.lab = 6)
        levels(scores)<-interval_names
        medians[,"scores"]<-scores
        colnames(df[,index])<-unique(medians[,"Sample"])
    }
    if(log){
        df[, index]<-log10(df[, index]+1)
        y_lab<-paste("log10(", attribute,"+1)", sep="")
    }else{
        y_lab<-capitalize(attribute)
    }
    x_lab<-"Sample"
    #if pool will appear one graph per pool
    if(pool) {
        df<-df[,c(index, which("pool" == names(df))), drop=FALSE]
        dfMelt<-melt(df, id.vars="pool")
        dfMelt[,"pool"]<-as.factor(dfMelt[,"pool"])
    }else {
        df<-df[, index, drop=FALSE]
        dfMelt<-melt(df)
    }
    levels(dfMelt[,"variable"])<-do.call(c, lapply(levels(dfMelt[,"variable"]), 
        function(panel){
            return(strsplit(panel, split=paste(attribute, "_",sep=""))[[1]][2])
    }))
    if(!(is.null(attributeThres))){
        if(pool){
            m1<-paste(dfMelt[,"pool"],dfMelt[,"variable"],sep="")
            m2<-paste(medians[,"pool"], medians[,"Sample"],sep="")
        }else{
            m1<-dfMelt[,"variable"]
            m2<-medians[,"Sample"]
        }
        dfMelt[,"scores"]<-medians[match(m1,m2),"scores"]

    }
    variable<-value<-NULL
    g<-ggplot(dfMelt,aes(x=variable, y=value,fill=variable))
    dens.plot<-ggplot(dfMelt,aes(value,fill=variable))
    if(!is.null(attributeThres)){
        g<-ggplot(dfMelt,aes(x=variable, y=value,fill=scores))
        dens.plot<-ggplot(dfMelt,aes(value,fill=scores, group=variable))
    }
    
    if(dens){
        if(join){
            g<-g+geom_violin( alpha=0.5,draw_quantiles = c(0.25, 0.5,0.75),
                trim=FALSE)+labs(x=x_lab,y=y_lab)+theme(axis.text= 
                element_text(size=12), axis.ticks=element_line(), 
                axis.title=element_text(size=16), legend.text=
                element_text(size=16), legend.title=element_text(size=16))
            if(pool){
                g<-g+facet_grid(~pool)+theme(strip.text.x = element_text(
                    size= 10))
            }
            if(!is.null(attributeThres)){
                colors<-colorRampPalette(c("red", "green"))(length(
                    interval_names))
                names(colors)<-interval_names
                g<-g+scale_fill_manual(name=paste(capitalize(attribute), 
                    "intervals", sep=" "), breaks=interval_names, values=
                    colors)+theme(axis.text=element_text(size=12), axis.ticks=
                    element_line())+labs(x="Sample")
            }else{
                g<-g+scale_fill_discrete("Sample")+ theme(axis.text=
                    element_text(size=12), axis.text.x=element_blank(), 
                    axis.ticks.x=element_blank(), axis.title=element_text(
                    size=16))
            }    

        }else{
            box.plot<-g+geom_boxplot()+labs(x=x_lab,y=y_lab)+theme(axis.text=
                element_text(size=12), axis.ticks= element_line(), axis.title=
                element_text(size=16), legend.text=element_text(size=16), 
                legend.title=element_text(size=16))+guides(fill=FALSE)
            #marginal density of y - plot on the right
            dens.plot<- dens.plot+
                geom_density(alpha=0.5)+coord_flip()+ labs(x=y_lab, y="")+
                theme(axis.text=element_text(size=12), axis.ticks= 
                element_line(), axis.title=element_text(size=16), legend.text=
                element_text(size=16), legend.title=element_text(size=16))

            if(pool){
                box.plot<-box.plot+facet_grid(~pool)+theme(strip.text.x = 
                    element_text(size= 10))
                dens.plot<-dens.plot+facet_grid(~pool)+theme(strip.text.x = 
                    element_text(size= 10))
            }
            if(!is.null(attributeThres)){
                colors<-colorRampPalette(c("red", "green"))(length(
                    interval_names))
                names(colors)<-interval_names
                dens.plot<-dens.plot+scale_fill_manual(name=paste(capitalize(
                    attribute), "intervals", sep=" "), breaks=interval_names, 
                    values=colors)
                box.plot<-box.plot+labs(x="Sample")+ scale_fill_manual(name=
                    paste(capitalize(attribute), "intervals", sep=" "), breaks=
                    interval_names, values=colors)
            }else{
                dens.plot<-dens.plot+scale_fill_discrete("Sample")+ theme( 
                    axis.text.x=element_blank(), axis.ticks.x=element_blank())
                box.plot<-box.plot+theme(axis.text.x=element_blank(), 
                    axis.ticks.x=element_blank())
                
            }
            
        #arrange the plots together, with appropriate height and width for each
        #row and column
            g<-cowplot::plot_grid(box.plot, dens.plot, ncol=2, nrow=1)
        }
    }else{
        g<-g+geom_boxplot()+labs(x=x_lab,y=y_lab)+theme( axis.title=
            element_text(size=16), legend.text = element_text(size = 16),
            legend.title=element_text(size=16),axis.text=element_text(size=12))

        if(pool){
            g<-g+facet_grid(~pool)+theme(strip.text.x = element_text(size= 10))
        }
        if(!is.null(attributeThres)){
            colors<-colorRampPalette(c("red", "green"))(length(interval_names))
            names(colors)<-interval_names
            g<-g+scale_fill_manual(name=paste(capitalize(attribute), 
                "intervals", sep=" "), breaks=interval_names, values=colors)+
                theme( axis.text= element_text(size=12), axis.ticks= 
                element_line(), axis.title= element_text(size=16), legend.text=
                element_text( size=16), legend.title= element_text(size=16))+
                labs(x="Sample")
        }else{
            g<-g+scale_fill_discrete("")+theme(axis.text.x=element_blank(), 
                axis.ticks.x=element_blank())
        }

    }
    return(g)    
    
})
