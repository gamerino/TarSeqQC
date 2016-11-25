#'Plot read profiles for a particular genomic region.
#'
#'\code{plotRegion} plots the read profiles for a selected region. The minAAF 
#'parameter set the minimum proportion value
#'to call an SNP and the minRD the minimum read depth. They are combined to 
#'obtain a minimum read count value at each position used to distinguish 
#'between possible SNPs and background noise. If SNPs is
#'set as 'TRUE', colored bars will appear indicating the occurrence of possible 
#'SNPs surpassing the minAAF and minRD, at each genomic position.
#'
#'@param object TargetExperiment object.
#'@param region Numeric of length two indicating the selected genomic region.
#'@param seqname Character indicating the chromosome of the genomic region.
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
#'include TargetExperiment-FeatPerform.R
#'@exportMethod plotRegion
#'@docType methods
#'@name plotRegion
#'@rdname TargetExperiment-plotRegion
#'@aliases plotRegion-methods
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
#'
#'# getting and exploring a sequenced region of a particular gene
#'getRegion(ampliPanel, level="gene", ID="gene7", collapse=FALSE)
#'# plot a particular genomic region
#'g<-plotRegion(ampliPanel,region=c(4500,6800), seqname="chr10", SNPs=TRUE,  
#'xlab="", title="gene7 amplicons",size=0.5)
#'# x11(type="cairo")
#'if(interactive()){
#'g
#'}

setGeneric(name="plotRegion", def=function(object, region, seqname, SNPs=TRUE, 
minAAF=0.05, minRD=10,xlab="", title="", size=0.5, BPPARAM=bpparam()){
    standardGeneric("plotRegion")
})
#'
#'@name plotRegion
#'@rdname TargetExperiment-plotRegion
#'@import BiocParallel
#'@importMethodsFrom BiocParallel bpprogressbar<-
#'@aliases plotRegion,TargetExperiment-method
#'@inheritParams plotRegion
setMethod(f="plotRegion", signature=signature(object="TargetExperiment"), 
definition=function(object,region, seqname, SNPs=TRUE,  minAAF=0.05, minRD=10,
xlab="", title="", size=0.5, BPPARAM=bpparam()){
#   bpprogressbar(BPPARAM)<-TRUE
    bedFile<-getBedFile(object)
    scanBamP<-getScanBamP(object)
    pileupP<-getPileupP(object)
    # keep only the region of interest
    bamWhich<-GRanges(seqnames=seqname, ranges=IRanges(start=region[1], 
        end=region[2]), strand="*")
    if(is.null(scanBamP)) scanBamP<-ScanBamParam()
    bed_region<-bamWhich  
    bamWhich(scanBamP)<-bed_region
    # obtaing the pileup for the selected region
    cts<-pileupCounts(bed=bed_region, bamFile=path(getBamFile(object)),
        fastaFile=path(getFastaFile(object)), scanBamP=scanBamP, 
        pileupP=pileupP, BPPARAM=BPPARAM)
    if(!all(c("pos", "seqnames") %in% colnames(cts))){
        stop(paste("pileupDF should contain at least 'seqnames' and 'pos' 
            columns", "Please consider to use the pileupCounts() function to 
            build it", sep="\n"))
    }
    cts[,"pos"]<-as.numeric(as.character(cts[,"pos"]))
    cts[,"seqnames"]<-as.character(cts[,"seqnames"])
    cts<-cts[order(cts[,"seqnames"], cts[,"pos"]),]
    if(SNPs & ! all(c("A", "C", "G", "T", "=") %in% names(cts))){
        stop("'SNPs was set as TRUE but 'distinguish_nucleotide' was set as 
            FALSE")
    }
    if(any("-" == names(cts) )){
        names(cts)["-" == names(cts)]<-"deletion"
    }
    if(any("=" == names(cts) )){
        names(cts)["=" == names(cts)]<-"ref_consensus"
    }
    pos<-NULL
    counts<-NULL
    ref_consensus<-NULL
    A<-G<-C<-Te<-NULL
    #building the plot
    p<-ggplot(data=cts)+geom_area(data=cts,aes(x=pos, y=counts, fill="counts"),
        alpha=0.5)
    p<-p+geom_line(data=cts,aes(x=pos, y=ref_consensus, color="ref_consensus"),
        size=size)
    color<-c(ref_consensus="darkmagenta")
    # adding SNP information
    if(SNPs){
        nosnps<-cts[,c("A", "C", "T", "G")] == cts[,"ref_consensus"]
        cts[nosnps[,1], "A"]<-0
        cts[nosnps[,2], "C"]<-0
        cts[nosnps[,3], "T"]<-0
        cts[nosnps[,4], "G"]<-0
        if(minAAF > 0 | minRD >0){
        idx<-(rowSums(cts[,c("A","C","G","T")])>0)
        if(any(idx) ){
            idx<-which(idx)
            possibleSNP<-do.call(rbind, bplapply(idx, function(x){
                aux<-cts[x,]
                thres<-max(c(minRD, round(aux[,"counts"]*minAAF)))
                ref<-as.character(aux[,"seq"])
                alts<-c("A","C","G","T")[c("A","C","G","T") != ref]
                alts<-alts[aux[,alts] !=0]
                for(i in alts){
                    if(aux[,i] < thres) aux[,i]<-0
                }
                return(aux)
            }, BPPARAM=BPPARAM))
            cts[idx,]<-possibleSNP
        }
        }
        names(cts)[names(cts) == "T"]<-"Te"
        #only for BiocCheck T detection
        cts[nosnps[,4], "G"]<-0
        A<-G<-C<-Te<-NULL
        p<-p+geom_bar(data=cts,aes(x=pos, y=A, fill="A"), size=size, 
            stat="identity")
        p<-p+geom_bar(data=cts,aes(x=pos, y=C, fill="C"), size=size,
            stat="identity")
        p<-p+geom_bar(data=cts,aes(x=pos, y=Te, fill="T"), size=size,
            stat="identity")
        p<-p+geom_bar(data=cts,aes(x=pos, y=G, fill="G"), size=size, 
            stat="identity")
        fill<-c( A="green", C="blue", T="red", G="brown")
        p<-p+scale_color_manual(name="Profiles", values=color, breaks=names(
            color))+scale_fill_manual(name="", values=c(counts="grey29", fill),
            breaks=names(c(counts="grey29", fill)))
    }else{
    p<-p+scale_color_manual(name="Profiles", values=color, breaks=names(color))+
        scale_fill_manual(name="", values=c(counts="grey29"), breaks=
        names(c(counts="grey29")))
    }
    bed<-as.data.frame(bedFile[seqnames(bedFile) == seqname])
    aux<-bed[(bed$end >=region[1]  & bed$end <region[2]) | (bed$start >=
        region[1] & bed$start <region[2]),c("start", "end")]
    aux[aux[,"start"]< region[1] & aux[,"end"]<=region[2] ,"start"]<-region[1]
    aux[aux[,"end"] > region[2] & aux[,"start"]>= region[1] ,"end"]<-region[2]
#     x_start<-which 
#     x_end<-which(cts[,"pos"] %in% end(bedFile[seqnames(bedFile) == seqname]))
#     # take into account several cases
#     # i) length(x_start) = length(x_end) and x_start > x_end
#     if (length(x_start) == length(x_end) & any(x_start > x_end)){
#         x_start<-c(1,x_start)
#         x_end<-c(x_end, nrow(cts))    
#     }
#     # ii) length(x_start) > length(x_end) 
#     if (length(x_start) != length(x_end)){
#         if(length(x_start) > length(x_end)) {
#             x_end<-c(x_end, nrow(cts))
#         }else x_start<-c(1, x_start) # iii) length(x_start) > length(x_end)
#     }
#    
# adding feature information
    p<-p+geom_segment(data=aux, aes(x=start,  y=-10, xend=end, yend=-10),
        color="darkcyan", alpha=rep(0.7,nrow(aux)), size=2)
    p<-p+labs(title = title, x = xlab, y="Counts")
    p<-p+theme(panel.background=element_rect(fill="white", color="black"),
        legend.key=element_rect( fill="white", color="white"), 
        panel.grid.major=element_line(color="lightgoldenrod3"), 
        panel.grid.minor=element_line(color="tomato", linetype = "dashed"),
        legend.text = element_text(size = 18),legend.title = element_text(
        size = 16, face = "bold"), axis.text=element_text(size=12), 
        plot.title=element_text(size=22, hjust=0.5), axis.title=element_text(
        size=22)) + guides(color = guide_legend(order = 1), fill = 
        guide_legend(order = 2))

    return(p)
})
