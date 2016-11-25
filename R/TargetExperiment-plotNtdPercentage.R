#'Plot nucleotide read percentages for a particular feature. 
#'
#'\code{plotNtdPercentage} plots the percentages of the ocurrence of each
#'nucleotide in each position for a selected feature. 
#'
#'@param object a TargetExperiment object.
#'@param featureID a character indicating the feature ID.
#'@param BPPARAM An optional BiocParallelParam instance defining the parallel
#'back-end to be used during evaluation.
#'
#'returned by the function.
#'
#'@return ggplot2 graphics
#'
#'@include TargetExperiment-plotGeneAttrPerFeat.R
#'@exportMethod plotNtdPercentage
#'@docType methods
#'@name plotNtdPercentage
#'@rdname TargetExperiment-plotNtdPercentage
#'@aliases plotNtdPercentage-methods
#'@seealso \code{\link{plotFeature}}
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
#'@examples
#'## loading TargetExperiment object
#'data(ampliPanel, package="TarSeqQC")
#'## Defining bam file, bed file and fasta file names and paths
#'setBamFile(ampliPanel)<-system.file("extdata", "mybam.bam", 
#'    package="TarSeqQC", mustWork=TRUE)
#'setFastaFile(ampliPanel)<-system.file("extdata", "myfasta.fa", 
#'    package="TarSeqQC", mustWork=TRUE)
#'# Exploring the nucleotide percentages compositions of the read counts for a 
#'# particular amplicon
#'g<-plotNtdPercentage(ampliPanel,featureID="AMPL20")
#'if(interactive()){
#'g
#'}

setGeneric(name="plotNtdPercentage", def=function(object, featureID, 
BPPARAM=bpparam()){
    standardGeneric("plotNtdPercentage")
})
#'
#'@name plotNtdPercentage
#'@rdname TargetExperiment-plotNtdPercentage
#'@importFrom reshape2 melt
#'@import BiocParallel
#'@importMethodsFrom BiocParallel bpprogressbar<-
#'@aliases plotNtdPercentage,TargetExperiment-method
#'@inheritParams plotNtdPercentage
setMethod(f="plotNtdPercentage", 
signature=signature(object="TargetExperiment"), definition=function(object, 
featureID, BPPARAM=bpparam()){
    bpprogressbar(BPPARAM)<-TRUE
    bedFile<-getBedFile(object)
    if(!(featureID %in% names(bedFile))){
        stop(paste("The feature",featureID,"is not part of the bed file", 
            sep =" "))
    }
    bed_feature<-bedFile[names(bedFile) == featureID, drop=FALSE]
    #if(is.null(getScanBamP(object))) scanBamP<-ScanBamParam()
    scanBamP<-getScanBamP(object)
    bamWhich(scanBamP)<-bed_feature
    cts<-pileupCounts(bed=bed_feature, bamFile=path(getBamFile(object)),
        fastaFile=path(getFastaFile(object)), scanBamP=scanBamP, 
        pileupP=getPileupP(object), BPPARAM=BPPARAM)
    if(!(any(c("A", "C", "T", "G") %in% colnames(cts)))) {
        stop("In order to count the nucleotide percentages for each position
            'distinguish_nucleotide' should be set as TRUE in the pileupP 
            parameter")
    }
    cts[,c("A","C","T","G")]<-100*cts[,c("A", "C", "T", "G") ]/cts[,c("counts")]
    cts[cts[,"counts"] == 0,c("A", "C", "T", "G")]<-0

    melt_perc<-melt(cts,id.vars="pos", measure.vars=names(cts)[names(cts) %in% 
        c("A", "C", "G","T")])
    title<-paste("Ntd percentages of ", featureID,sep="")
    pos<-value<-variable<-NULL
    if("seq" %in% colnames(cts)){
        refSeq<-(cts$seq)
        melt_perc[,"seq"]<-as.character(rep(refSeq,4))
        melt_perc[,"value"][melt_perc[,"variable"] == melt_perc[,"seq"]]<-0
    }
    p<-ggplot(data=melt_perc,aes(x=as.factor(pos), y=value, fill=variable))+
        geom_bar(stat="identity")
    color<-c(A="green", C="blue", T="red", G="brown")
    p<-p+scale_fill_manual(name="Nucleotides", values=color, breaks=names(
        color)) + labs(title=title, x="", y="Ntd %")+theme(plot.title=
        element_text(size=22, hjust=0.5), axis.title=element_text(size=22),
        legend.text=element_text(size=18))
    if("seq" %in% colnames(cts)){
        refSeq<-(cts$seq)
        melt_perc$seq<-as.character(rep(refSeq,4))
        p<-p+scale_x_discrete(labels=melt_perc[,"seq"])
    }
    p<-p+ylim(c(0,100))
    return(p)
})
