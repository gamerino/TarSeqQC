#'TargetExperiment object constructor.
#'
#'\code{initialize} creates the TargetExperiment object architecture for the
#'specified bed and alingment BAM files. If 'scanBamP' and/or 'pileupP' 
#'parameters are not specified, default values of their constructors will be
#'used.
#'
#'@param .Object TargetExperiment class.
#'@param bedFile Character indicating the bed file full path.
#'@param bamFile Character indicating the alignment and index bam
#'files full paths.
#'@param fastaFile Character indicating the full path to the genome reference
#'and index files.
#'@param scanBamP ScanBamParam indicating the parameters for read the BAM 
#'file.
#'@param pileupP PileupParam indicating the parameters for pileup building.
#'@param feature Character indicating the name of the feature that will
#'be explored (e.g 'amplicon', 'exon').
#'@param attribute Character indicating the name of the attribute that will
#'be explored. Should be 'coverage' or 'medianCounts'.
#'@param BPPARAM An optional BiocParallelParam instance defining the parallel
#'back-end to be used during evaluation.
#'
#'@return TargetExperiment object.
#'
#'@include TargetExperiment-summarizePanel.R
#'@exportMethod initialize
#'@docType methods
#'@name initialize
#'@rdname TargetExperiment-initialize
#'@import Rsamtools 
#'@import GenomicRanges
#'@importMethodsFrom BiocGenerics strand
#'@importClassesFrom S4Vectors Rle
#'@importMethodsFrom S4Vectors Rle
#'@importClassesFrom IRanges IRanges
#'@importFrom IRanges IRanges
#'@import BiocParallel
#'@importMethodsFrom BiocParallel bpprogressbar<-
#'@aliases initialize,TargetExperiment-method
#'@seealso \code{\link{TargetExperiment}}, \code{\link{buildFeaturePanel}}
#'\code{\link{summarizePanel}}
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@family TargetExperiment
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar} and Elmer A. Fernandez 
#'\email{efernandez@@bdmg.com.ar}
#'@examples
#'## Defining bam file, bed file and fasta file names and paths
#'if(interactive()){
#'  bamFile<-system.file("extdata", "mybam.bam", package="TarSeqQC",
#'      mustWork=TRUE)
#'  bedFile<-system.file("extdata", "mybed.bed", package="TarSeqQC", 
#'      mustWork=TRUE)
#'  fastaFile<-system.file("extdata", "myfasta.fa", package="TarSeqQC", 
#'      mustWork=TRUE)
#'
#'  ## Creating a TargetExperiment object
#'
#'  ## Defining feature parameter
#'  feature<-"amplicon"
#'  ## Defining attribute parameter
#'  attribute<-"coverage"
#'  ##Calling the constructor
#'  ampliPanel<-TargetExperiment(bedFile, bamFile, fastaFile, 
#'      attribute=attribute, feature=feature)
#'}
setMethod(f="initialize", signature=signature(.Object="TargetExperiment"),
definition=function(.Object, bedFile, bamFile, fastaFile, scanBamP=NULL, 
pileupP=NULL, feature=NULL, attribute=NULL, BPPARAM=bpparam()){
#     bpprogressbar(BPPARAM)<-TRUE
    ##Set the different slots
    if(nargs() >=2){
    # bedFile slot
        bed_df<-read.delim(bedFile, header=TRUE, stringsAsFactor=FALSE)
        if(!(all(c("chr", "start","end", "name", "gene") %in% names(bed_df)))){
            stop("Error! Bed file should containt at least five columns: 'chr', 
                'start', 'end', 'name' and 'gene'")
        }
        if(any(duplicated(bed_df[,"name"]))){
            stop("Error! Bed file should one row per feature. The provided file 
                has duplicated 'name' IDs.")
        }
        if (any(duplicated(bed_df[,"start"]) & duplicated(bed_df[,"end"]))){
        warning("Your bed file have duplicated features. Removing 
            duplicated ...")
        index<-which(duplicated(bed_df[,"start"])& duplicated(bed_df[,"end"]))
        bed_df<-bed_df[-index,]
        }
        strand<-"strand" %in% names(bed_df)
        bed<-GRanges(seqnames=Rle(bed_df[,"chr"]), ranges=IRanges(start=
        bed_df[,"start"], end=bed_df[,"end"], names=bed_df[,"name"]))
        if(strand) strand(bed)<-bed_df[,"strand"]
        mcols(bed)<-bed_df[,!(names(bed_df) %in% c("chr", "start", "end", 
            "name", "strand")), drop=FALSE]

        .Object@bedFile<-bed
        # bamFile slot
        bam<-BamFile(file=bamFile,index=paste(bamFile, ".bai", sep=""))
        if (!file.exists(index(bam))){
            print("The index of your BAM file doesn't exist")
            print("Building BAM file index")
            dummy<-indexBam(bam)
        }
        .Object@bamFile<-bam
        # fastaFile slot
        fasta<-FaFile(file=fastaFile,index=paste(fastaFile, ".fai", sep="") )
        if(!file.exists(index(fasta))){
            print("The index of your fasta file doesn't exist")
            print("Building fasta file index")
            dummy<-indexFa(fasta)
        }
        .Object@fastaFile<-fasta
        # scanBamP slot
        if(is.null(scanBamP)) scanBamP<-ScanBamParam(which=bed)
        if(length(bamWhich(scanBamP)) == 0) {
            warning("The 'scanBamP' should have the regions to scan in its 
                'which' parameter \n", immediate. =TRUE)
            warning("Specifying scanBamP 'which' parameter based on the bed file
                \n", immediate. =TRUE)
            bamWhich(scanBamP)<-bed
        }
        setScanBamP(.Object)<-scanBamP
        # pileupP slot
        if(is.null(pileupP)) pileupP<-PileupParam()
        setPileupP(.Object)<-pileupP
        # featurePanel slot
        .Object@featurePanel <- buildFeaturePanel(.Object, BPPARAM)
        # genePanel slot
        .Object@genePanel <- summarizePanel(.Object,  BPPARAM)
        # feature slot
        if (!is.null(feature)) setFeature(.Object)<-feature
        # attribute slot
        if(is.null(attribute)){
            setAttribute(.Object)<-""
            warning("The 'attribute' slot is empty", immediate. =TRUE)
        }else{
            if(!(attribute %in% c("coverage", "medianCounts"))){
                stop("'attribute' should be 'coverage' or 'medianCounts'")}
            setAttribute(.Object)<-attribute
        }
        if(is.null(feature)){
        setFeature(.Object)<-""
        warning("The 'feature' slot is empty", immediate. =TRUE)
        }
    }else{

        .Object@bedFile<-GRanges()
        .Object@bamFile<-BamFile(file=character(1))
        .Object@fastaFile<-FaFile(file=character(1))
        setScanBamP(.Object)<-ScanBamParam()
        setPileupP(.Object)<-PileupParam()
        .Object@featurePanel<-GRanges()
        .Object@genePanel<-GRanges()
        setFeature(.Object)<-character(0)
        setAttribute(.Object)<-character(0)
    }
    ##Check the object's validity
    validObject(.Object)
    return(.Object)
})
