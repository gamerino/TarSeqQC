#'Plot read profiles for a particular feature.
#'
#'\code{plotFeature} plots the read profiles for a selected feature The minAAF 
#'parameter set the minimum proportion value
#'to call an SNP and the minRD the minimum read depth. They are combined to 
#'obtain a minimum read count value at each position used to distinguish 
#'between possible SNPs and background noise. If SNPs is
#'set as 'TRUE', colored bars will appear indicating the occurrence of possible 
#'SNPs surpassing the minAAF and minRD, at each genomic position.
#'
#'@param object TargetExperiment object.
#'@param featureID Character indicating the ID of the feature.
#'@param SNPs Logical flag indicating if SNPs should be plotted.
#'@param minAAF Numeric indicating the minimum alternative allele proportion 
#'necessary to call a SNP.
#'@param minRD Numeric indicating the minimum read depth of alternative 
#'alleles necessary to call a SNP.
#'@param xlab Character containing the axis x label.
#'@param title Character containing the plot title.
#'@param size Numeric indicating the size of line plots.
#'@param BPPARAM An optional BiocParallelParam instance defining the parallel
#'back-end to be used during evaluation.
#'
#'@return ggplot2 graphics.
#'
#'@include TargetExperiment-plotRegion.R
#'@exportMethod plotFeature
#'@docType methods
#'@name plotFeature
#'@rdname TargetExperiment-plotFeature
#'@aliases plotFeature-methods
#'@seealso \code{\link{plotRegion}}
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
#'# Exploring the read count profile for a particular amplicon
#'g<-plotFeature(ampliPanel, featureID="AMPL20")
#'if(interactive()){
#'g
#'}
setGeneric(name="plotFeature", def=function(object,featureID, SNPs=TRUE, 
minAAF=0.05, minRD=10, xlab="", title="",size=0.5, BPPARAM=bpparam()){
    standardGeneric("plotFeature")
})
#'
#'@name plotFeature
#'@rdname TargetExperiment-plotFeature
#'@import BiocParallel
#'@importMethodsFrom BiocParallel bpprogressbar<-
#'@aliases plotFeature,TargetExperiment-method
#'@inheritParams plotFeature
setMethod(f="plotFeature", signature=signature(object="TargetExperiment"),
definition=function(object, featureID, SNPs=TRUE,  minAAF=0.05, minRD=10, 
xlab="", title=featureID, size=0.5, BPPARAM=bpparam()){
    bpprogressbar(BPPARAM)<-TRUE
    bed_file<-getBedFile(object)
    
    if(!(featureID %in% names(bed_file))) {
        stop(paste("The feature", featureID, "is not contained in the bed file",
        sep=" "))
    }
    
    bed_feature<-bed_file[names(bed_file) == featureID]
    if(length(bed_feature) >1){
        stop(paste("The bed file contains two regions with the same name: ",
        featureID, sep =""))
    }
    p<-plotRegion(object,region=c(start(bed_feature), end(bed_feature)), 
        seqname=runValue(seqnames(bed_feature)), SNPs,  minAAF, minRD, xlab,
        title, size, 
        BPPARAM)

    return(p)

})

