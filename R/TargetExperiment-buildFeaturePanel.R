#'Function to build a feature panel based on specific genomic regions.
#'
#'\code{buildFeaturePanel} builds panel slots of a TargetExperiment object.
#'Input can be a bam file or a pileup matrix. If the bed file contains a high
#'number of amplicons, the bam file as input is recommended in order to
#'diminish memory requirements. The resulting object is a GRanges instance
#'having panel and counts/coverage information.
#'
#'@param object TargetExperiment class object.
#'@param BPPARAM An optional BiocParallelParam instance defining the parallel
#'back-end to be used during evaluation.
#'
#'@return GRanges object.
#'
#'@include TargetExperiment-pileupCounts.R
#'@exportMethod buildFeaturePanel
#'@docType methods
#'@name buildFeaturePanel
#'@rdname TargetExperiment-buildFeaturePanel
#'@aliases buildFeaturePanel-methods
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar} and Elmer A. Fernandez 
#'\email{efernandez@@bdmg.com.ar}
#'@examples
#'## loading TargetExperiment object
#'data(ampliPanel, package="TarSeqQC")
#'## Defining bam file, bed file and fasta file names and paths
#'setBamFile(ampliPanel)<-system.file("extdata", "mybam.bam", 
#'    package="TarSeqQC", mustWork=TRUE)
#'setFastaFile(ampliPanel)<-system.file("extdata", "myfasta.fa", 
#'    package="TarSeqQC", mustWork=TRUE)
#'
#'myFeaturePanel<-buildFeaturePanel(ampliPanel)
setGeneric(name="buildFeaturePanel", def=function(object, 
BPPARAM=bpparam()){
    standardGeneric("buildFeaturePanel")
})
#'
#'@name buildFeaturePanel
#'@rdname TargetExperiment-buildFeaturePanel
#'@importMethodsFrom BiocParallel bpmapply bplapply bpworkers bpprogressbar<-
#'@import BiocParallel
#'@importFrom Biostrings letterFrequency
#'@aliases buildFeaturePanel,TargetExperiment-method
#'@inheritParams buildFeaturePanel
setMethod(f="buildFeaturePanel", signature="TargetExperiment",
definition=function(object,BPPARAM=bpparam()){
#     BPPARAM2<-bpparam()
#     bpprogressbar(BPPARAM2)<-TRUE
    bed_file<-getBedFile(object)
    param <- getScanBamP(object)
    #extract GC content from Fasta File
    aux<-scanFa(getFastaFile(object), param=bed_file)
    gc <- round(rowSums(letterFrequency(aux, letters= c("G","C")))/width(aux),3)
    mcols(bed_file)<-cbind(mcols(bed_file),GC=as.data.frame(gc))
    rm(aux)
    # ensure that only those reads overlapping targeted regions are counted
    bamWhat(param) <-c("qname", "pos", "qwidth", "rname")
    bamFile<-getBamFile(object)
    aln <- scanBam(bamFile, param=param)
    # create RangedData object
    reads<- do.call(rbind, bplapply(1:length(aln), function(x){
    with(aln[[x]], data.frame(pos=pos, width=qwidth, ID=qname,
        seqnames=rname))}, BPPARAM=BPPARAM))
    reads<-GRanges(IRanges(reads[,"pos"], width=reads[,"width"]), 
        ID=reads[,"ID"], seqnames=reads[,"seqnames"])
   
#     reads<- do.call(c, bplapply(1:length(aln), function(x){
#     with(aln[[x]], GRanges(IRanges(pos, width=qwidth), ID=qname,
#         seqnames=rname))}, BPPARAM=BPPARAM))
#     
    chrs <-as.character(unique(seqnames(bed_file)))
    
    # compute counts of all reads
    counts <- coverage(reads)
    #compute statistics
    info<-do.call(rbind,lapply(chrs, function(chr){  
        index<-which(as.character(seqnames(bed_file)) == chr)
    # select cromosome counts and features
        chrCounts <- counts[[chr]]
        featRange <- ranges(bed_file[index])
        aux <- bplapply(featRange, function(x){
                if(all(x <= length(chrCounts))){
                    return(chrCounts[x])
                }else{ return(c(0,0))}
            }) 
    # compute average median SD and IQR for counts per feature
            cov <- floor(sapply(aux, mean))
            sdcov <- floor(sapply(aux, sd))
            medcount<-floor(sapply(aux, median))
            iqrcount<-floor(sapply(aux, IQR))
            return(cbind(coverage=cov, sdCoverage=sdcov,
                medianCounts= medcount,  IQRCounts=iqrcount))
    }))
    
    m<-match(names(bed_file), rownames(info))
    mcols(bed_file)<-cbind(mcols(bed_file),as.data.frame(info[m,]))
    return(bed_file)
})

