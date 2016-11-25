#'Function to control Bed and FASTA files compatibility.
#'
#'\code{checkBedFasta} checks the compatibility of a Bed file and a Fasta file.
#' The functions first will control the consistency of the Bed file in terms of
#'duplicated positions or feature's IDs and correct definition of start-end
#'values. Then, the method will control the consistency between the specified
#'features and the reference file. During its execution, several testing 
#'messages will be printed.
#'
#'@param bedFile Character indicating the bed file full path.
#'@param fastaFile Character indicating the full path to the genome reference 
#'file.
#'
#'@return NULL
#'
#'@include TargetExperiment-print.R
#'@export checkBedFasta
#'@docType methods
#'@name checkBedFasta
#'@rdname checkBedFasta
#'@aliases checkBedFasta-methods
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar}, Yanina Murua \email{ymurua@leloir.org.ar},
#'Andrea S. Llera \email{allera@leloir.org.ar} and Elmer A. Fernandez 
#'\email{efernandez@bdmg.com.ar}
#'@examples
#'##Define the bed and fasta file full paths
#'bedFile<-system.file("extdata", "mybed.bed", package="TarSeqQC",
#'    mustWork=TRUE)
#'fastaFile<-system.file("extdata", "myfasta.fa", 
#'    package="TarSeqQC", mustWork=TRUE)
#'##Checking the bed-fasta consistency
#'checkBedFasta(bedFile, fastaFile)
#'
#'@name checkBedFasta
#'@import Rsamtools
checkBedFasta<-function(bedFile, fastaFile){
    bed_df<-read.delim(bedFile, header=TRUE, stringsAsFactor=FALSE)
    cat("- Checking bed file, mandatory columns:\n")
    if(!(all(c("chr", "start","end", "name", "gene") %in% names(bed_df)))){
        stop(paste("Bed file should contain at least five columns: 'chr', ", 
            "'start', 'end', 'name' and 'gene'", sep=""))
    }else{cat(" OK")}
    cat("\n- Checking bed file, duplicated feature IDs:\n")
    if(any(duplicated(bed_df[,"name"]))){
        stop(paste( "Bed file should have one row per feature. The provided ", 
            "file has duplicated 'name' IDs.", sep=""))
    }else{cat(" OK")}
    cat("\n- Checking bed file, duplicated feature locations:\n")
    if (any(duplicated(bed_df[,"start"]) & duplicated(bed_df[,"end"]))){
        cat(paste("Warning: Your bed file have duplicated features sharing ", 
            "chromosome, start and end positions in rows ", which(duplicated(
            bed_df[,"start"])& duplicated(bed_df[,"end"])), ". If duplicates ", 
            "persist, they will be removed when TargetExperiment constructor ",
            "is called", sep=""))
    }else{cat(" OK")}
    cat("\n- Checking bed file, start and end values:\n")
    if(any(bed_df[,"start"] < 1 | bed_df[,"end"] <1) ){
        stop(paste(" The start and end could not be lower than 1"))
    }else{cat(" OK")}
    strand<-"strand" %in% names(bed_df)
    bed<-GRanges(seqnames=Rle(bed_df[,"chr"]), ranges=IRanges(start=
        bed_df[,"start"], end=bed_df[,"end"], names=bed_df[,"name"]))
    if(strand) strand(bed)<-bed_df[,"strand"]
    mcols(bed)<-bed_df[,!(names(bed_df) %in% c("chr", "start", "end", 
        "name", "strand")), drop=FALSE]
    
    #get seqinfo Fasta file
    fasta<-FaFile(fastaFile)
    fastaInfo<-as.data.frame(seqinfo(fasta))
    
    #check chromosome consistency
    cat("\n- Checking bed-fasta chromosome consistency:\n")
    chrs<-levels(seqnames(bed))
    if(!all(chrs %in% rownames(fastaInfo))){
        notPres<-chrs[!chrs %in% rownames(fastaInfo)]
        cat(paste("Warning: Your bed file has the chromosomes: '",notPres, 
            "', not found in the reference FASTA file",sep=""))
    }else{cat(" OK")}
    #check genomic coordinates consistency
    cat("\n- Checking bed-fasta genomic coordinates consistency:\n")
    chrsSE<-as.matrix(do.call(rbind, lapply(1:length(chrs), function(x){
        ids<-which(as.character(seqnames(bed)) ==chrs[x])
        rangeChr<-ranges(bed[ids])
        minR<-min(c(start(rangeChr), end(rangeChr)))
        maxR<-max(c(start(rangeChr), end(rangeChr)))
        return(c(chr=chrs[x], min=minR,max=maxR))
    })))
    chrsSE<-cbind(chrsSE, length=fastaInfo[chrsSE[,"chr"],"seqlengths"])
    if(any(as.numeric(chrsSE[,"min"]) > as.numeric(chrsSE[,"length"])) |any( 
        as.numeric(chrsSE[,"max"]) > as.numeric(chrsSE[,"length"] ))){
            cat(paste("Warning: Your bed file has features out of the ", 
                "chromosomes' ranges from the Fasta file", sep=""))
        cat("\n")
        idx<-which(as.numeric(chrsSE[,"min"]) > as.numeric(chrsSE[,"length"])| 
            as.numeric(chrsSE[,"max"]) > as.numeric(chrsSE[,"length"]))
            warning(paste("Please, check feature's definitions of chromosome ",
                paste(chrsSE[idx,"chr"],collapse = ", "),".",sep="" ))
    }else{cat(" OK\n")}
    return(invisible(NULL))
}    
