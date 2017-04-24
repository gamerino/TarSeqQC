#'Plot TargetExperiment object overview.
#'
#'\code{plot} allows a fast and simple representation of one feature panel
#'using a polar histogram plot. Histogram bar reflects the percentage of 
#'features that have shown the analyzed attribute in a user set interval.
#'The resulting graph can be busy and might be better off saved.
#'
#'@param x TargetExperiment/TargetExperimentList class object.
#'@param y not used but necessary for redefining the generic function.
#'@param attributeThres Numeric indicating the interval extreme values.
#'@param binSize Numeric indicating bin width. Should probably be left
#'as 1, as other parameters are relative to it.
#'@param spaceGene Numeric. Space between bins.
#'@param spaceChr Numeric. Space between chromosomes.
#'@param innerRadius Numeric. Radius of the inner circle.
#'@param outerRadius Numeric. Radius of the outer circle. 
#'@param guides A vector with percentages to use for the white guide lines.
#'@param alphaStart Numeric offset from 12 o'clock in radians.
#'@param circleProportion Numeric proportion of the circle to cover.
#'@param direction Character indicating if the increasing count goes from or
#'to the center.
#'@param chrLabels Logical. Chromosome names must be plotted?.
#'
#'@return a ggplot2 graph.
#'
#'@include TargetExperiment-statistics.R
#'@docType methods
#'@name plot
#'@rdname TargetExperiment-plot
#'@import ggplot2
#'@importFrom grDevices colorRampPalette
#'@importFrom grDevices hcl
#'@importFrom graphics plot
#'@importFrom Hmisc capitalize
#'@aliases plot,TargetExperiment,plot.TargetExperiment
#'@seealso \code{\link{plotFeatPerform}}
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
#'@references
#'\url{http://www.r-bloggers.com/polar-histogram-pretty-and-useful/}
#'@examples
#'## Loading the TargetExperiment object
#'data(ampliPanel, package="TarSeqQC")
#'# Definition of the interval extreme values
#'attributeThres<-c(0,1,50,200,500, Inf)
#'
#'## Plot panel overview
#'g<-plot(ampliPanel, attributeThres, chrLabels =TRUE)
#'if(interactive()){
#'g
#'}
#'@export plot.TargetExperiment
plot.TargetExperiment <- function(x, y, attributeThres=c(0, 1, 50, 200, 500, 
Inf),binSize=1, spaceGene=0.2,  spaceChr=1.2,  innerRadius=0.3,
outerRadius=1, guides=c(20,40,60,80),  alphaStart=-0.3,  
circleProportion=0.95, direction="inwards",  chrLabels=FALSE,...){
    if(attributeThres[1] !=0){
        attributeThres<-c(0,attributeThres)
    }
    if(attributeThres[length(attributeThres)] !=Inf){
        attributeThres<-c(attributeThres, Inf)
    }
    df_panel<-as.data.frame(getFeaturePanel(x))
    df_panel[,"names"]<-rownames(df_panel)
    attribute<-getAttribute(x)
    if(!(attribute %in% c("coverage", "medianCounts"))){
        stop("Attribute slot should be defined in order to call biasExploration
            function")
    }
    # creating a 'score' variable to group features according to the attribute
    #'intervals
    df_panel[,"score"]<-cut(df_panel[,attribute], breaks=attributeThres, 
        include.lowest=TRUE, right=FALSE, dig.lab = 6)
    score_levels<-levels(df_panel[,"score"])

    df_panel<-df_panel[order(df_panel[,"seqnames"], df_panel[,"gene"], 
        df_panel[,"score"]),]
    geneFeat<-sapply(unique(df_panel[,"gene"]), function(gene){
        paste(gene, "(",length(which(df_panel[,"gene"] == gene)), ")", sep="")
    })
    gene_idx<-match(df_panel[, "gene"],unique(df_panel[,"gene"]))
    gene<-geneFeat[gene_idx]
#     df_panel<-cbind(df_panel, gene)
    df_panel[,"gene"]<-gene

    aux<-sapply(unique(df_panel[,"seqnames"]), function(chr){
        info<-t(sapply(unique(df_panel[df_panel[,"seqnames"] == chr,"gene"]), 
        function(gene){
            return(as.matrix(table(df_panel[which(df_panel[,"gene"] == gene &
            df_panel[,"seqnames"]== chr),"score"])))
        }))
        return(list(info))
    })
    aux<-do.call(rbind,aux)
    rownames(aux)<-unique(df_panel[, "gene"])

    gene_names<-as.list(unique(as.character(df_panel[, "gene"])))
    df_panel<-do.call(rbind,lapply(unique(df_panel[, "gene"]), function(gene){
        seqnames<-rep(as.character(unique(df_panel[df_panel[,"gene"] == gene, 
            "seqnames"])), times=ncol(aux))
        score<-levels(df_panel[, "score"])
        values<-aux[rownames(aux) == gene,]
        gene<-rep(gene, times=ncol(aux))
        return(data.frame(seqnames=seqnames, gene=gene, score=score, 
            values=values))
    }))
    levels(df_panel[,"score"])<-score_levels
    df_panel<-ddply(df_panel,c("seqnames","gene"),transform,values= cumsum(
        values/(sum(values))))
    if(any(is.na(df_panel[,"values"]))){
        df_panel[is.na(df_panel[, "values"]), "values"]<-0
    }

    df_panel<-ddply(df_panel, c("seqnames", "gene"), transform, previous= c(0, 
        head(values, length(values)-1)))
    df_panel<-df_panel[order(df_panel[,"seqnames"], df_panel[,"gene"]),]
    df_panel$indexGene<-factor(df_panel[, "gene"], levels= as.character(
        unique(df_panel[, "gene"])))
    levels(df_panel$indexGene)<-1:length(levels(df_panel[, "indexGene"]))

    df_panel$indexChr<-factor(df_panel[,"seqnames"], levels=as.character(
        unique(df_panel[,"seqnames"])))
    levels(df_panel$indexChr)<-1:length(levels(df_panel[,"indexChr"]))
    df_panel[,which(names(df_panel) %in% c("indexGene", "indexChr"))]<-apply(
        df_panel[,which(names(df_panel) %in% c("indexGene", "indexChr"))],
            2,as.numeric)

    affine<-switch(direction,
        'inwards'= function(y) (outerRadius-innerRadius)*y+innerRadius,
        'outwards'=function(y) (outerRadius-innerRadius)*(1-y)+innerRadius,
        stop(paste("Unknown direction value")))
    xmin<-(df_panel[,"indexGene"]-1)*binSize + 
        (df_panel[,"indexGene"]-1)*spaceGene +
        (df_panel[,"indexGene"]-1)*(spaceChr-spaceGene)
    
    xmax<-xmin+binSize
    ymin<-affine(1-df_panel[,"previous"])
    ymax<-affine(1-df_panel[,"values"])
    df_panel<-cbind(df_panel, xmin, xmax, ymin, ymax)
    guidesDF<-data.frame(
        xmin=rep(xmin, times=length(guides)), y=rep(1-guides/100, 
        times=1, each=nrow(df_panel)))
    xend<-guidesDF[,"xmin"]+binSize
    guidesDF<-cbind(guidesDF, xend)
    guidesDF$y<-affine(guidesDF[,"y"])
    # Building the ggplot object
    totalLength<-tail(df_panel[, "xmin"]+binSize+spaceChr,1)/circleProportion-0
    p<-ggplot(df_panel)+geom_rect(aes( xmin=xmin, xmax=xmax, ymin=ymin, 
        ymax=ymax,fill=score))
    colors<-colorRampPalette(c("red", "green"))(length(score_levels))
    names(colors)<-score_levels
    p<-p+scale_fill_manual(name=paste(attribute, "interval", sep=" "),
                breaks=score_levels, values=colors)
    # names labels
    readableAngle<-function(x){
        angle<-x*(-360/totalLength)-alphaStart*180/pi+90
        angle+ifelse(sign(cos(angle*pi/180))+sign(sin(angle*pi/180))==-2,180,0)
    }
    readableJustification<-function(x){
        angle<-x*(-360/totalLength)-alphaStart*180/pi+90
        ifelse(sign(cos(angle*pi/180))+sign(sin(angle*pi/180))==-2,1,0)
    }

    df_panelItemLabels<-ddply(df_panel, "gene", summarize, xmin=xmin[1])
    df_panelItemLabels$x<-df_panelItemLabels[,"xmin"]+binSize/2
    angle<-readableAngle(df_panelItemLabels[,"xmin"] + binSize/2)
    hjust<-readableJustification(df_panelItemLabels[,"xmin"] + binSize/2)
    df_panelItemLabels<-cbind(df_panelItemLabels, angle, hjust)
    p<-p+geom_text(aes(x=x, label=gene, angle=angle, hjust=hjust), y=1.02, 
        size=3, vjust=0.5, data=df_panelItemLabels)
    # guides  
    p<-p+geom_segment(aes( x=xmin, xend=xend, y=y, yend=y), colour="white",
        data=guidesDF)

    # label for guides
    label<-paste(guides, "% ", sep='')
    guideLabels<-data.frame( x=0, y=affine(1-guides/100), label=label)

    p<-p+geom_text( aes(x=x, y=y, label=label), data=guideLabels, 
        angle=-alphaStart*180/pi, hjust=1, size=4)

    # gene labels
    if(chrLabels){
        chrLabelsDF<-aggregate(formula=xmin~seqnames,data=df_panel,
            FUN=function(s){
                mean(s+binSize)
            })
#         chrLabelsDF<-within(chrLabelsDF,{
#             x<-xmin
#             angle<-xmin*(-360/totalLength)-alphaStart*180/pi
#         })
        chrLabelsDF$x<-chrLabelsDF[,"xmin"]
        chrLabelsDF$angle<-chrLabelsDF[,"xmin"]*(-360/totalLength) - 
            alphaStart*180/pi
        p<-p+geom_text( aes( x=x, label=seqnames, angle=angle), 
            data=chrLabelsDF, y=1.4)
    }  
    p<-p+theme(panel.background=element_blank(), axis.title.x=element_blank(),
        axis.title.y=element_blank(), panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank(), axis.ticks=element_blank() )
    p<-p+xlim(0,tail(df_panel[,"xmin"]+binSize+spaceChr,1)/circleProportion)
    p<-p+ylim(0,outerRadius+0.2)
    p<-p+guides(fill=guide_legend(title=paste(capitalize(attribute), 
        "_intervals", sep="")))
    p<-p+coord_polar(start=alphaStart)
    p

}
#'@S3method plot TargetExperiment
## S4 method dispatches to S3
setMethod("plot", "TargetExperiment", plot.TargetExperiment)
