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
    chrs <-as.character(unique(seqnames(bed_file)))
    alName<-do.call(rbind,strsplit(names(aln), ":"))
    info<-do.call(rbind,lapply(1:length(chrs), function(x){
        chr<-chrs[x]
        index<-which(alName[,1] == chr)
        dfInfo<-do.call(rbind,lapply(1:length(index), function(i){
            df<-as.data.frame(aln[[index[i]]])
            reads<-GRanges(IRanges(df[,"pos"], width=df[,"qwidth"]), 
                ID=df[,"qname"], seqnames=df[,"rname"])
        # compute counts of all reads
            count<-coverage(reads)

        #compute statistics

            index2<-which(as.character(seqnames(bed_file)) == chr)
        # select cromosome counts and features
            chrCounts <- count[[chr]]
            featRange <- ranges(bed_file[index2])[i]
            aux1 <- lapply(1:length(featRange), function(j){
                aux<-start(featRange)[j]:end(featRange)[j]
                if(all(aux <= length(
                    chrCounts))){
                        return(chrCounts[aux])
                }else{ return(c(0,0))}
            }) 
    # compute average median SD and IQR for counts per feature
            cov <- floor(as.numeric(lapply(aux1, mean)))
            sdcov <- floor(as.numeric(lapply(aux1, sd)))
            medcount<-floor(as.numeric(lapply(aux1, median)))
            iqrcount<-floor(as.numeric(lapply(aux1, IQR)))
            return(cbind(names=names(featRange),coverage=cov, sdCoverage=sdcov,
                medianCounts= medcount,  IQRCounts=iqrcount))}))
        return(as.data.frame(dfInfo))
    })
    )
    rownames(info)<-info[,1]
    info<-info[,-1]
    info[,"coverage"]<-as.numeric(as.character(info[,"coverage"]))
    info[,"sdCoverage"]<-as.numeric(as.character(info[,"sdCoverage"]))
    info[,"medianCounts"]<-as.numeric(as.character(info[,"medianCounts"]))
    info[,"IQRCounts"]<-as.numeric(as.character(info[,"IQRCounts"]))
    m<-match(names(bed_file), rownames(info))
    mcols(bed_file)<-cbind(mcols(bed_file),info[m,])
    return(bed_file)
})

