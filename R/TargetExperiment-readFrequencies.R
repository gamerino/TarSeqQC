#'Function to explore read frequencies in targeted regions and out 
#'targeted regions.
#'
#'\code{readFrequencies} builds a data frame containing the read frequencies 
#'falling in targeted regions and out of those, separated by chromosome. 
#'
#'@param object TargetExperiment class object.
#'@param BPPARAM An optional BiocParallelParam instance defining the parallel
#'back-end to be used during evaluation.
#'
#'@return data.frame object.
#'
#'@include TargetExperiment-plotNtdPercentage.R
#'@exportMethod readFrequencies
#'@docType methods
#'@name readFrequencies
#'@rdname TargetExperiment-readFrequencies
#'@aliases readFrequencies-methods
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
#'myReadPercentages<-readFrequencies(ampliPanel)
setGeneric(name="readFrequencies", def=function(object, 
BPPARAM=bpparam()){
    standardGeneric("readFrequencies")
})
#'
#'@name readFrequencies
#'@rdname TargetExperiment-readFrequencies
#'@importMethodsFrom BiocParallel bpmapply bplapply bpworkers bpprogressbar<-
#'@importMethodsFrom GenomeInfoDb seqlengths
#'@importMethodsFrom GenomicAlignments readGAlignments
#'@importFrom Rsamtools path
#'@import BiocParallel
#'@aliases readFrequencies,TargetExperiment-method
#'@inheritParams readFrequencies
setMethod(f="readFrequencies", signature="TargetExperiment",
definition=function(object,BPPARAM=bpparam()){
    bed_file<-getBedFile(object)
    param <- getScanBamP(object)
    bamFile<-getBamFile(object)
    #extract GC content from Fasta File
    chrLengths<-seqlengths(getFastaFile(object))
    chrs <-as.character(seqnames(bed_file))

    chrGaps<-gaps(bed_file)
    # add last gaps
    aux<-(do.call(rbind, bplapply(unique(chrs), 
        function(chrL){
            if (max(end(ranges(chrGaps)[seqnames(chrGaps) == chrL])) < 
                chrLengths[chrL])
                return(cbind(chr= chrL, start=max(end(ranges(bed_file)[
                    seqnames(bed_file) == chrL]))+1, end=chrLengths[chrL])) 
            BPPARAM=BPPARAM})))
    #check if any chromosome doesn't have features
    if(any(!(names(chrLengths)%in% as.character(unique(seqnames(bed_file)))))){
        notFeature<-names(chrLengths)[!(names(chrLengths) %in% as.character(
            unique(seqnames(bed_file))))] 
        aux<-rbind(aux,cbind(notFeature,1,chrLengths[!(names(chrLengths) %in% 
            as.character(unique(seqnames(bed_file))))]))
    }
        
    regionsOut<-GRanges(Rle(c(as.character(seqnames(chrGaps)), as.character(
        aux[,"chr"]))), IRanges( start=c(start(ranges(chrGaps)), as.numeric(
        aux[,"start"])),end=c(end(ranges(chrGaps)), as.numeric(aux[,"end"]))))
        
    # ensure that only those reads NON overlapping targeted regions are counted
    paramOut<-param
    bamWhat(paramOut) <-c("qname", "pos", "qwidth", "rname")
    bamWhich(paramOut)<-regionsOut
    # Count reads overlapping in targeted regions
    bamWhat(param) <-c("qname", "pos", "qwidth", "rname")
    aln<-readGAlignments(path(bamFile))
    outCounts<-countOverlaps(aln,GRanges(bamWhich(paramOut)),type="within")
    inCounts<-countOverlaps(aln,GRanges(bamWhich(param)))

    
    #compute statistics
    alnChr<-as.character(seqnames(aln))
    info<-as.data.frame(do.call(rbind,lapply(unique(alnChr),function(i){
        return(c(chr=i, In=sum(inCounts[alnChr==i]), Out=sum(
            outCounts[alnChr==i])))

    })))
    info[,"In"]<-as.numeric(as.character(info[,"In"]))
    info[,"Out"]<-as.numeric(as.character(info[,"Out"]))
    totalReads<-sum(info[, "In"], info[, "Out"])
    info[,"InPerc"]<-round(info[,"In"]/totalReads*100,2)
    info[,"OutPerc"]<-round(info[,"Out"]/totalReads*100,2)
    info[,"chr"]<- factor(info[,"chr"], levels=c(unique(as.character(info[,
        "chr"]))))
    
    return(info)
    
})

