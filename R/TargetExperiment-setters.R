#'Setters for the TargetExperiment slots
#'
#'Set TargetExperiment slots, according to the given function call.
#'
#'@param object TargetExperiment class object.
#'@param value value to set the slot.
#'
#'@return a TargetExperiment object
#'
#'@include TargetExperiment-getters.R
#'@exportMethod setFeature<-
#'@docType methods
#'@name setFeature<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature
#'@aliases setFeature<--methods
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar} Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar} and Elmer A. Fernandez 
#'\email{efernandez@@bdmg.com.ar}
#'@examples
#'## loading TargetExperiment object
#'if (interactive()){
#'  data(ampliPanel, package="TarSeqQC")
#'  ## Defining bam file, bed file and fasta file names and paths
#'  setBamFile(ampliPanel)<-system.file("extdata", "mybam.bam", 
#'      package="TarSeqQC", mustWork=TRUE)
#'  setBedFile(ampliPanel)<-system.file("extdata", "mybed.bed", 
#'      package="TarSeqQC", mustWork=TRUE)
#'  setFastaFile(ampliPanel)<-system.file("extdata", "myfasta.fa", 
#'      package="TarSeqQC", mustWork=TRUE)
#'
#'  ## Set feature slot value
#'  setFeature(ampliPanel)<-"amplicon"
#'}
setGeneric(name="setFeature<-", def=function(object, value){
    standardGeneric("setFeature<-")
})
#'
#'@name setFeature<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature<-
#'@aliases setFeature<-,TargetExperiment,character-method
setReplaceMethod(f="setFeature", signature=signature(object="TargetExperiment",
value="character"),  definition= function(object, value){
    object@feature<-value
    validObject(object)
    return(object) 
})
#'
#'@exportMethod setAttribute<-
#'@name setAttribute<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature
#'@aliases setAttribute<--methods
#'@examples
#'## Set attribute slot value
#'setAttribute(ampliPanel)<-"coverage"
setGeneric(name="setAttribute<-", def=function(object, value){
    standardGeneric("setAttribute<-")
})
#'
#'@name setAttribute<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature<-
#'@aliases setAttribute<-,TargetExperiment,character-method
setReplaceMethod(f="setAttribute", signature=signature(
object="TargetExperiment",value="character"), definition=function(
object, value){
    object@attribute<-value
    validObject(object)
    return(object) 
})
#'
#'@exportMethod setScanBamP<-
#'@name setScanBamP<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature
#'@aliases setScanBamP<--methods
#'@examples
#'## Set scanBamP slot value
#'setScanBamP(ampliPanel)<-ScanBamParam()
setGeneric(name="setScanBamP<-", def=function(object, value){
    standardGeneric("setScanBamP<-")
})
#'
#'@name setScanBamP<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature<-
#'@aliases setScanBamP<-,TargetExperiment,ScanBamParam-method
setReplaceMethod(f="setScanBamP", signature=signature(
object="TargetExperiment",value="ScanBamParam"),definition=function(
object, value){
    object@scanBamP<-value
    validObject(object)
    return(object) 
})
#'
#'@exportMethod setPileupP<-
#'@name setPileupP<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature<-
#'@aliases setPileupP<--methods
#'@examples
#'## Set pileupP slot value
#'setPileupP(ampliPanel)<-PileupParam()
setGeneric(name="setPileupP<-", def=function(object, value){
    standardGeneric("setPileupP<-")
})
#'
#'@name setPileupP<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature<-
#'@aliases setPileupP<-,TargetExperiment,PileupParam-method
setReplaceMethod(f="setPileupP",signature=signature(object="TargetExperiment",
value= "PileupParam"),  definition=function(object, value){
    object@pileupP<-value
    validObject(object)
    return(object) 
})
#'
#'@exportMethod setFeaturePanel<-
#'@name setFeaturePanel<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature
#'@aliases setFeaturePanel<--methods
#'@examples
#'## Set featurePanel slot value
#'setFeaturePanel(ampliPanel)<-buildFeaturePanel(ampliPanel)
setGeneric(name="setFeaturePanel<-", def=function(object, value){
    standardGeneric("setFeaturePanel<-")
})
#'
#'@name setFeaturePanel<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeaturePanel<-
#'@aliases setFeaturePanel<-,TargetExperiment,GRanges-method
setReplaceMethod(f="setFeaturePanel",signature=signature(
object="TargetExperiment", value= "GRanges"),  definition=function(object,
value){
    object@featurePanel<-value
    validObject(object)
    return(object) 
})
#'
#'@exportMethod setGenePanel<-
#'@name setGenePanel<-
#'@rdname TargetExperiment-setters
#'@inheritParams setGenePanel<-
#'@aliases setGenePanel<--methods
#'@examples
#'## Set genePanel slot value
#'setGenePanel(ampliPanel)<-summarizePanel(ampliPanel)
setGeneric(name="setGenePanel<-", def=function(object, value){
    standardGeneric("setGenePanel<-")
})
#'
#'@name setGenePanel<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature<-
#'@aliases setGenePanel<-,TargetExperiment,GRanges-method
setReplaceMethod(f="setGenePanel",signature=signature(object="TargetExperiment",
value= "GRanges"),  definition=function(object, value){
    if(length(value) == 0){
        stop("value should have at least one element")
    }
    if(!(all(names(value) %in% unique(mcols(object@featurePanel)[,"gene"])))){
        stop("value should contain all the same gene names that are contained
            in the 'gene' metadata column of the featurePanel slot")
    }
    object@genePanel<-value
    validObject(object)
    return(object) 
})
#'@exportMethod setBedFile<-
#'@name setBedFile<-
#'@rdname TargetExperiment-setters
#'@inheritParams setBedFile<-
#'@aliases setBedFile<--methods
#'@examples
#'## Set bedFile slot value
#'setBedFile(ampliPanel)<-system.file("extdata", "mybed.bed", 
#'package="TarSeqQC", mustWork=TRUE)
setGeneric(name="setBedFile<-", def=function(object, value){
    standardGeneric("setBedFile<-")
})
#'
#'@name setBedFile<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature<-
#'@aliases setBedFile<-,TargetExperiment,character-method
setReplaceMethod(f="setBedFile",signature=signature(object="TargetExperiment",
value= "character"),  definition=function(object, value){
    if(!file.exists(value)){
        stop("The file ", value, " doesn't exist", sep="")
    }
    bed_df<-read.delim(value, header=TRUE, stringsAsFactor=FALSE)
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
    bedFile<-GRanges(seqnames=Rle(bed_df[,"chr"]), ranges=IRanges(start=
    bed_df[,"start"], end=bed_df[,"end"], names=bed_df[,"name"]))
    if(strand) strand(bedFile)<-bed_df[,"strand"]
    mcols(bedFile)<-bed_df[,!(names(bed_df) %in% c("chr", "start", "end", 
        "name", "strand")), drop=FALSE]
    if (any(duplicated(end(bedFile)) & duplicated(start(bedFile)))){
    warning("Your bed file have duplicated features. Removing duplicated ...")
    bedFile<-bedFile[-which(duplicated(end(bedFile)) & duplicated(start(
        bedFile)))]
    }

    object@bedFile<-bedFile
    validObject(object)
    return(object) 
})
#'@exportMethod setBamFile<-
#'@name setsetBamFile<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature
#'@aliases setsetBamFile<--methods
#'@examples
#'## Set bamFile slot value
#'setBamFile(ampliPanel)<-system.file("extdata", "mybam.bam", 
#'package="TarSeqQC", mustWork=TRUE)
setGeneric(name="setBamFile<-", def=function(object, value){
    standardGeneric("setBamFile<-")
})
#'
#'@name setBamFile<-
#'@rdname TargetExperiment-setters
#'@inheritParams setBamFile<-
#'@aliases setBamFile<-,TargetExperiment,character-method
setReplaceMethod(f="setBamFile", signature=signature(
object="TargetExperiment",value="character"), definition=function(
object, value){
    if(!file.exists(value)){
        stop("The file ", value, " doesn't exist.", sep="")
    }else{
        bam<-BamFile(file=value,index=paste(value, ".bai", sep=""))
        if (!file.exists(index(bam))){
            print("The index of your BAM file doesn't exist")
            print("Building BAM file index")
            dummy<-indexBam(bam)
        }
        object@bamFile<-bam
        validObject(object)
        return(object) 
    }   
})
#'@exportMethod setFastaFile<-
#'@name setFastaFile<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature
#'@aliases setFastaFile<--methods
#'@examples
#'## Set fastaFile slot value
#'setFastaFile(ampliPanel)<-system.file("extdata", "myfasta.fa", 
#'package="TarSeqQC", mustWork=TRUE)
setGeneric(name="setFastaFile<-", def=function(object, value){
    standardGeneric("setFastaFile<-")
})
#'
#'@name setFastaFile<-
#'@rdname TargetExperiment-setters
#'@inheritParams setFeature<-
#'@aliases setFastaFile<-,TargetExperiment,character-method
setReplaceMethod(f="setFastaFile", signature=signature(
object="TargetExperiment",value="character"), definition=function(
object, value){
    if(!file.exists(value)){
        stop("The file ", value, " doesn't exist.", sep="")
    }else{
        fasta<-FaFile(file=value,index=paste(value, ".fai", sep="") )
        if(!file.exists(index(fasta))){
            print("The index of your fasta file doesn't exist")
            print("Building fasta file index")
            dummy<-indexFa(fasta)
        }
        object@fastaFile<-fasta
        validObject(object)
        return(object)
    }
})
