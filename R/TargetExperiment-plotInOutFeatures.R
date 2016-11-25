#'Function to explore read percentages in targeted regions and out 
#'targeted regions.
#'
#'\code{plotInOutFeatures} allows the graphical exploration of the data 
#'frame obtained using \code{readFrequencies}. This data frame contains
#'information about the amount of reads mapped to the targeted regions and out
#'of them. This information is presented in rows, one for each chromosome and
#'in absolute and relative amounts. After it invocation, a bar plot built as a 
#'ggplot object is returned
#'
#'@param object a data frame or a TargetExperiment.
#'@param absolute logical indicating if absolute frequency should be used.
#'@param ... additional parameters according to the function call
#'
#'@return ggplot object.
#'
#'@include TargetExperiment-myCounts.R
#'@exportMethod plotInOutFeatures
#'@docType methods
#'@name plotInOutFeatures
#'@rdname plotInOutFeatures
#'@aliases plotInOutFeatures
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
#'
#'g<-plotInOutFeatures(ampliPanel)
#'
setGeneric(name="plotInOutFeatures", def=function(object,...){
    standardGeneric("plotInOutFeatures")
})
#'@rdname plotInOutFeatures
#'@aliases plotInOutFeatures
#'@inheritParams plotInOutFeatures
setMethod(f="plotInOutFeatures", signature="data.frame",
definition=function(object, absolute=FALSE){
    if(!is.data.frame(object)) {
        stop("object should be a data frame")
    }
    if( any(!c("chr", "In", "Out", "InPerc", "OutPerc") %in% colnames(object
        ))){ stop ("object should contain 'chr', 'In', 'Out', 'InPerc' and
            'OutPerc' columns")
    }
    if (absolute){
        object<-object[,colnames(object) %in% c("chr","In", "Out")]
        colnames(object)<-c("Chromosome", "In targets", "Out targets")
    }else {
        object<-object[,colnames(object) %in% c("chr", "InPerc", 
            "OutPerc")]
        colnames(object)<-c("Chromosome", "In targets (%)", 
            "Out targets (%)")
    }
    reshapeInfo<-melt(object, id.vars="Chromosome")
    Chromosome<-value<-variable<-NULL
    g<-ggplot(reshapeInfo, aes(x=Chromosome, y=value, fill=variable))+geom_bar(
        stat="identity", position="dodge")+theme(legend.title=element_blank())+
        scale_fill_manual(values=c("blue", "red"))+labs(y="Value")
    return(g)
})
#'@rdname plotInOutFeatures
#'@aliases plotInOutFeatures
#'@inheritParams plotInOutFeatures
#'@param BPPARAM An optional BiocParallelParam instance defining the parallel
#'back-end to be used during evaluation.
setMethod(f="plotInOutFeatures", signature="TargetExperiment",
definition=function(object, absolute=FALSE, BPPARAM=bpparam()){
    readsInfo<-readFrequencies(object, BPPARAM)
    g<-plotInOutFeatures(object=readsInfo, absolute)
    return(g)
})
