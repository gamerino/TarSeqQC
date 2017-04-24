#'Function to obtain the pileup counts for a bam file.
#'
#'\code{pileupCounts} waits for a TargetExperiment object containing the bed
#' file information in order to obtain pileup counts only for the specified 
#'genomic regions. The resulting object is a data.frame instance, in which each
#'row represents one position of the specified features across the bed file. 
#'The first three columns called 'pos', 'seqnames' and 'which_label,' represent
#'the position in the seqnames (e.g. pos=10183795 and seqnames=chr3) and the 
#'associated feature. According to the 'pileupP' parameters set before, the 
#'number of next columns could change. If 'distinguish_nucleotide' was set to
#'TRUE, then one column per ntd will appear containing the counts obtained for
#'each of them. Same will occur when 'distinguish_strands' is set to TRUE. The
#'last column, called 'counts', contains the total counts obtained for the
#'corresponding position.
#'
#'@param bed a Granges object containing the bed file information.
#'@param bamFile Character indicating the alignment and index bam
#'files full path.
#'@param fastaFile Character indicating the full path to the genome reference
#'and index files.
#'@param scanBamP ScanBamParam indicating the parameters the BAM file read.
#'@param pileupP PileupParam indicating the parameters for the pileup build.
#'@param BPPARAM An optional BiocParallelParam instance defining the parallel
#'back-end to be used during evaluation.
#'
#'@return data.frame object.
#'
#'@include TargetExperiment-print.R
#'@export pileupCounts
#'@docType methods
#'@name pileupCounts
#'@rdname pileupCounts
#'@aliases pileupCounts-methods
#'@seealso \code{Rsamtools-pileup}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
#'@references 
#'\enumerate{
#'\item Morgan M, Pages H, Obenchain V and Hayden N. Rsamtools: 
#'Binary alignment (BAM), FASTA, variant call (BCF), and tabix file import. 
#'R package version 1.20.1}
#'@examples
#'##Obtain the pileup matrix for the first amplicon
#'data(ampliPanel, package="TarSeqQC")
#'bed<-getBedFile(ampliPanel)[1]
#'## Defining bam file and fasta file names and paths
#'bamFile<-system.file("extdata", "mybam.bam", 
#'    package="TarSeqQC", mustWork=TRUE)
#'fastaFile<-system.file("extdata", "myfasta.fa", 
#'    package="TarSeqQC", mustWork=TRUE)
#'## extracting the pileup matrix
#'myCounts<-pileupCounts(bed, bamFile, fastaFile)
#'head(myCounts)
#'
#'@name pileupCounts
#'@importMethodsFrom GenomeInfoDb seqnames seqnames<- seqlevels
#'@importMethodsFrom S4Vectors runValue
#'@import BiocParallel
#'@importMethodsFrom BiocParallel bpprogressbar<- bplapply bpmapply
pileupCounts<-function(bed, bamFile, fastaFile, scanBamP=NULL, pileupP=NULL,
    BPPARAM=bpparam()){
    #checking for duplicated features
    if (any(duplicated(end(bed)) & duplicated(start(bed)))){
    ("Your bed file have duplicated features. Removing duplicated ...")
    bed<-bed[-which(duplicated(end(bed)) & duplicated(start(bed)))]
    }
    #setting values
    bam<-BamFile(file=bamFile,index=paste(bamFile, ".bai", sep=""))
    if (!file.exists(index(bam))){
        print("The index of your BAM file doesn't exist")
        print("Building BAM file index")
        dummy<-indexBam(bam)
    }
    fasta<-FaFile(file=fastaFile,index=paste(fastaFile, ".fai", sep="") )
    if(!file.exists(index(fasta))){
        print("The index of your fasta file doesn't exist")
        print("Building fasta file index")
        dummy<-indexFa(fasta)
    }
    if(is.null(scanBamP)) scanBamP<-ScanBamParam(which=bed)
    if(length(bamWhich(scanBamP)) == 0) {
        warning("The 'scanBamP' should have the regions to scan in its 
            'which' parameter \n", immediate. =TRUE)
        warning("Specifying scanBamP 'which' parameter based on the bed file
            \n", immediate. =TRUE)
        bamWhich(scanBamP)<-bed
    }
    
    if(is.null(pileupP)) pileupP<-PileupParam()
    
    # building pileup using RSamtools pileup function
    res <- pileup(bam,scanBamParam=scanBamP,  pileupParam=pileupP)
    strand<-"strand" %in% names(res)
    if(!strand) res[,"strand"]<-"*"
    # adding reference sequence
    fasta_seq<-getSeq(fasta, param=bed)
    sequence<-unlist(strsplit(paste(sapply(1:length(fasta_seq), function(x){
        return(as.character(fasta_seq[[x]]))
    }), collapse=''), split=""))
#     positions<-bplapply(as.list(1:length(bed_file)), function(x){
#         chr<-as.character(runValue(seqnames(bed_file)[x]))
#         start<-start(bed_file)[x]
#         end<-end(bed_file)[x]
#         which_label<-paste(chr,":", start, "-", end, sep="")
#         return(cbind(seqnames=rep(chr, times=length(start:end)),pos=start:end,
#         which_label=which_label))
#     }, BPPARAM=BPPARAM)
    chr <- as.character(seqnames(bed))
    start_bed <- start(bed)
    end_bed <- end(bed)
    which_label <- paste0(chr, ":", start_bed, "-", end_bed)
    fun <- function(chr, start_bed, end_bed, which_label) {
        cbind(seqnames=chr, pos=start_bed:end_bed, which_label=which_label)
    }
    positions <- mapply(fun, chr, start_bed, end_bed,  which_label, 
                SIMPLIFY=FALSE)

    if(is.matrix(positions)){
        fasta_info<-data.frame(cbind(matrix(positions, ncol=3),seq=sequence))
        names(fasta_info)<-c("seqnames", "pos", "which_label", "seq")
    }else{
        fasta_info<-data.frame(do.call(rbind, positions), seq=sequence)
    }

    if("nucleotide" %in% names(res)  ){
        nucleotides<-as.factor(res$nucleotide)
        ntd_types<-levels(nucleotides)
        ntd_matrix<-as.data.frame(matrix(
        rep(0,(length(ntd_types)*length(nucleotides))),nrow=length(nucleotides),
            ncol=length(ntd_types)) )
        names(ntd_matrix)<-ntd_types
        count_matrix<-sapply(ntd_types, function(ntd){
            if( length(nucleotides) > 0){
                    return(as.numeric(as.character(nucleotides) == ntd)) 
            }else{
                return(0)
            }
        })
    }else{
    count_matrix<-1
    }
    final_counts<-cbind(res[,!(names(res) %in% 
        c("nucleotide", "count", "which_label"))], count_matrix*res$count, 
            which_label=res$which_label)
    simples<-unique(fasta_info[,1:2])        
    final_counts<-as.data.frame(do.call(rbind, bplapply(1:nrow(simples), 
        function(x){
        index<-(final_counts[,"pos"] == simples[,"pos"][x] & 
#         as.character(final_counts[,"which_label"]) == as.character(
#         fasta_info[,"which_label"][x]))
        as.character(final_counts[,"seqnames"]) == as.character(
        simples[,"seqnames"][x]))
        if(length(which(index))==0){
            seqnames<-as.character(simples[x, "seqnames"])
            pos<-as.character(simples[x, "pos"])
            which_label<-paste0(fasta_info[fasta_info[,"seqnames"] == 
            simples[x,1] & fasta_info[,"pos"] == simples[x,2], "which_label"], 
            collapse=";")
            if(length(count_matrix) ==1){
                counts<-0
            }else{
                counts<-rep(0,7)
                names(counts)<-c("A", "C", "G", "T", "N" , "=","-")
            }
        }else{
            seqnames<-unique(as.character(final_counts[index, "seqnames"]))
            pos<-as.character(unique(final_counts[index,"pos"]))
            which_label<-paste0(fasta_info[fasta_info[,"seqnames"] == 
            simples[x,1] & fasta_info[,"pos"] == simples[x,2], "which_label"], 
            collapse=";")
            if (length(count_matrix) ==1){
                counts<-sum(unique(final_counts[index, 
                    "count_matrix * res$count"]))
                names(counts)<-"counts"
            }else{
                final_counts[index, c("A", "C", "G", "T", "N" , "=", 
                    "-")]<-apply(final_counts[index, c("A", "C", "G", "T", "N",
                    "=", "-")], 2, as.numeric)
                counts<-colSums(unique(final_counts[index, c("strand", "A", 
                    "C", "G", "T", "N" , "=", "-")])[,-1])

            }
        }
        return(c(seqnames=seqnames, pos=pos,counts,which_label= which_label))
    }, BPPARAM=BPPARAM)))
    if(!("counts" %in% names(final_counts))){
    final_counts[, c("A", "C", "G", "T", "N" , "=", "-")]<-apply(
    final_counts[, c("A", "C", "G", "T", "N" , "=", "-")], 2, as.numeric)
    final_counts$counts<-rowSums(
        final_counts[,c("A", "C","G", "T", "N", "=", "-")])

    }

    final_pileup<-merge(x=unique(fasta_info[,-3]), y = final_counts, by=sort(
        intersect(names(final_counts), names(fasta_info[,-3]))), all.x=TRUE,
        all.y=FALSE)
    if(all(c("=", "seq") %in% names(final_pileup))){
        final_pileup[,"="]<-bpmapply(function(x){final_pileup[x,
            names(final_pileup) %in% as.character(final_pileup$seq[x])]},
            1:nrow(final_pileup), BPPARAM=BPPARAM)
    }
    final_pileup[,"seqnames"]<-factor(final_pileup[,"seqnames"], 
        levels=seqlevels(bed))
    
    final_pileup[,"pos"]<-as.numeric(as.character(final_pileup[,"pos"]))
    final_pileup<-final_pileup[order(final_pileup[,"seqnames"],
    final_pileup[,"pos"]),]
    return(final_pileup)
}
