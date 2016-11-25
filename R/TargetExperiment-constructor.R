#'TargetExperiment constructor
#'
#'\code{TargetExperiment} creates a TargetExperiment object with the
#'architecture specified by the bed and alignment BAM files. If 'scanBamP' 
#'and/or 'pileupP' parameters are not specified, default values of their 
#'constructors will be used. attribute and feature parameters can be set 
#'after constructor calling.
#'
#'@param bedFile Character indicating the bed file full path.
#'@param bamFile Character indicating the alignment and index bam
#'files full path.
#'@param fastaFile Character indicating the full path to the genome reference
#'and index files.
#'@param scanBamP ScanBamParam indicating the parameters the BAM file read.
#'@param pileupP PileupParam indicating the parameters for the pileup build.
#'@param feature Character indicating the name of the feature that will
#'be explored (e.g 'amplicon', 'exon').
#'@param attribute Character indicating the name of the attribute that will
#'be explored. Should be 'coverage' or 'medianCounts'.
#'@param BPPARAM An optional BiocParallelParam instance defining the parallel
#'back-end to be used during evaluation.
#'
#'@return TargetExperiment object.
#'
#'@include TargetExperiment-initialize.R
#'@export TargetExperiment
#'@docType methods
#'@name TargetExperiment
#'@importMethodsFrom BiocGenerics strand
#'@importClassesFrom S4Vectors Rle
#'@importMethodsFrom S4Vectors Rle
#'@importClassesFrom IRanges IRanges
#'@importFrom IRanges IRanges
#'@import BiocParallel
#'@importMethodsFrom BiocParallel bpprogressbar<-
#'@rdname TargetExperiment-constructor
#'@aliases TargetExperiment-methods
#'@seealso \code{\link{TargetExperiment-class}}1
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@family TargetExperiment
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
#'@examples
#'## Defining bam file, bed file and fasta file names and paths
#'bamFile<-system.file("extdata", "mybam.bam", package="TarSeqQC",
#'mustWork=TRUE)
#'bedFile<-system.file("extdata", "mybed.bed", package="TarSeqQC",
#'mustWork=TRUE)
#'fastaFile<-system.file("extdata", "myfasta.fa", package="TarSeqQC", 
#'mustWork=TRUE)
#'
#'## Creating a TargetExperiment object
#'
#'# Defining feature parameter
#'feature<-"amplicon"
#'# Defining attribute parameter
#'attribute<-"coverage"
#'##Calling the constructor
#'object<-TargetExperiment(bedFile, bamFile, fastaFile, attribute=attribute,
#'feature=feature)
#'
TargetExperiment<-function(bedFile, bamFile, fastaFile, scanBamP=NULL, 
pileupP=NULL, feature=NULL, attribute=NULL, BPPARAM=bpparam()){
#     bpprogressbar(BPPARAM)<-TRUE
    .Object<-new("TargetExperiment")
    if(nargs() >=1 ){
        .Object<-initialize(.Object, bedFile, bamFile, fastaFile, scanBamP, 
            pileupP, feature, attribute, BPPARAM)
    }
    ##Check the object's validity
    validObject(.Object)
    return(.Object)
}