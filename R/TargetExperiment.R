#'TargetExperiment S4 class implementation in R
#' 
#'This S4 class represents a Targeted Sequencing Experiment in R. 
#'Targeted Sequencing Experiments are characterized by a 'bed file' that
#'contains the specification of the explored 'features' as a 'panel'. This  
#'features could be amplicons, exons, transcripts, among others. In general 
#'each feature is associated to one gene. A gene could be related to many 
#'features. This class allows the representation and quality control of a 
#'Targeted Sequencing Experiment.
#'
#'@slot scanBamP ScanBamParam containing 
#'the information to scan the BAM file.
#'@slot pileupP PileupParam containing 
#'the information to build the pileup.
#'@slot bedFile GRanges object that models the bed file.
#'@slot bamFile BamFile object that is a reference to the BAM file.
#'@slot fastaFile FaFile object that is a 
#'reference to the reference sequence.
#'@slot featurePanel GRanges object that models 
#'the feature panel and related statistics.
#'@slot genePanel GRanges object that models
#'the analyzed panel and related statistics at a gene level.
#'@slot attribute character indicates which attribute 'coverage' or 
#''medianCounts' will be used to the analysis.
#'@slot feature character indicates the name of the analyzed features. E.g
#''amplicon', 'exon', 'transcript'.
#'
#'@section Features: 
#'\enumerate{
#'  \item Model Targeted Sequencing Experiments in R. 
#'  \item Obtain coverage and read counts per sequenced feature.
#'  \item Evaluate the performance of a targeted sequencing experiment using
#'   coverage/read counts information. 
#'  \item Detect in early stage sequencing or library preparation errors.
#'  \item Explore read profiles for particular features or genomic regions.
#'  \item Explore any kind of experiment in which 'feature' definition is
#'   possible for several genes. E.g RNA-seq experiments in which transcripts
#'   could be the 'features'.
#'  \item Report quality control results.
#'}
#'
#'@section Functions:
#'TargetExperiment S4 class includes the following functions:
#'\describe{
#' \item{pileupCounts}{calculate pileup statistics for the BAM file}
#' \item{buildFeaturePanel}{build and model a feature panel as a GRanges
#'  object and compute read statistics}
#' \item{summarizePanel}{summarize the feature panel to a gene panel and
#'  compute read statistics}
#' \item{initialize}{constructor of TargetExperiment to generate the feature
#' and gene panels starting from an alignment BAM file and the bed file} 
#' \item{getBedFile, getBamFile, getFeaturePanel, getGenePanel,  getAttribute,
#' getFeature, getScanBamP, getPileupP}{return the respective TargetExperiment
#' slot}
#' \item{setAttribute,setFeature, setScanBamP, setPileupP}{set the 
#' respective TargetExperiment slots}
#' \item{show}{generic output of the object}
#' \item{print}{generic output of the object}
#' \item{summary}{print statistics summary for the set attribute}
#' \item{freqTable}{build a frequency table of the attribute occurrence in
#'  user configured intervals}
#' \item{plot}{plot a summarized view of the feature panel performance}
#' \item{plotAttrExpl}{plot the density and distribution of the attribute}
#' \item{plotFeatPerform}{plot the sequencing performance for each 
#'  feature and/or gene}
#' \item{plotFeature}{plot the reads profile for a particular feature}
#' \item{plotGeneAttrPerFeat}{plot the explored attribute for each feature of
#'  a particular gene}
#' \item{plotNtdPercentages}{plot nucleotide percentages for each position of
#'  a particular feature}
#' \item{plotRegion}{plot the reads profile for a particular genomic region}
#' \item{readFrequencies}{calculate frequencies of reads fall in and 
#'  out of targeted regions}
#' \item{plotInOutFeatures}{plot frequencies of reads fall in and out of 
#'  targeted regions}
#' \item{biasExploration}{plot attribute distributions along groups of bias
#'  sources}
#' \item{plotMetaDataExpl}{plot density and box plots or frequency bar plot of
#'  metadata columns}
#' \item{addStatSummSheet}{internal function to add the first sheet of xlsx 
#'  reports}
#' \item{buildReport}{build the experiment report as an xlsx file}
#' }
#'
#'@include TarSeqQC-package.R
#'@import methods GenomicRanges Rsamtools
#'@importMethodsFrom BiocGenerics strand
#'@importClassesFrom S4Vectors Rle
#'@importMethodsFrom S4Vectors Rle
#'@importClassesFrom IRanges IRanges
#'@importFrom IRanges IRanges
#'@name TargetExperiment-class
#'@rdname TargetExperiment-class
#'@exportClass TargetExperiment
#'@seealso Rsamtools
#'@family TargetExperiment
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar} and Elmer A. Fernandez 
#'\email{efernandez@@bdmg.com.ar}
#'@examples
#'
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
#'ampliPanel<-TargetExperiment(bedFile, bamFile, fastaFile, attribute=attribute,
#'feature=feature)
#'
#'## Alternative object creation
#'# Creating the TargetExperiment object
#'ampliPanel<-TargetExperiment(bedFile, bamFile, fastaFile)
#'# Set feature slot value
#'setFeature(ampliPanel)<-"amplicon"
#'# Set attribute slot value
#'setAttribute(ampliPanel)<-"coverage"
#'# Set pileupP slot value in order to set the maximum depth at 1000
#'setPileupP(ampliPanel)<-PileupParam(max_depth=1000)
#'# Set the featurePanel slot but now using the new pileupP definition
#'setFeaturePanel(ampliPanel)<-buildFeaturePanel(ampliPanel)
#'## Early exploration
#'# show/print
#'ampliPanel
#'# summary
#'summary(ampliPanel)
#'# summary at feature level
#'summaryFeatureLev(ampliPanel)
#'# summary at gene level
#'summaryGeneLev(ampliPanel)
#'# attribute boxplot and density plot exploration
#'g<-plotAttrExpl(ampliPanel,level="feature",join=TRUE, log=FALSE, color="blue")
#'if(interactive()){
#'g
#'}
#'# explore amplicon length distribution
#'g<-plotMetaDataExpl(ampliPanel, "length", log=FALSE, join=FALSE, color=
#'"blueviolet")
#'if(interactive()){
#'g
#'}
#'# explore gene's relative frequencies
#'g<-plotMetaDataExpl(ampliPanel, "gene", abs=FALSE)
#'if(interactive()){
#'g
#'}
#'## Deep exploration and Quality Control
#'myfrequencies<-readFrequencies(ampliPanel)
#'g<-plotInOutFeatures(readFrequencies(ampliPanel))
#'if(interactive()){
#'g
#'}
#'# definition of the interval extreme values
#'attributeThres<-c(0,1,50,200,500, Inf)
#'# plot panel overview
#'g<-plot(ampliPanel, attributeThres, chrLabels =TRUE)
#'if(interactive()){
#'g
#'}
#'# plot panel overview in a feature performance plot
#'g<-plotFeatPerform(ampliPanel, attributeThres, complete=TRUE, log=FALSE, 
#'featureLabs=TRUE, sepChr=TRUE, legend=TRUE)
#'if(interactive()){
#'g
#'}
#'# explore possible attribute bias
#'g<-biasExploration(ampliPanel, source="gc", dens=TRUE)
#'if(interactive()){
#'x11(type="cairo");g
#'}
#'## Controlling low counts features
#'# Do a frequency table for the attribute intervals
#'summaryIntervals(ampliPanel, attributeThres)
#'#plotting attribute intervals
#'g<-plotAttrPerform(ampliPanel)
#'if(interactive()){
#'x11(type="cairo");g
#'}
#'# getting low counts features at gene level
#'getLowCtsFeatures(ampliPanel, level="gene", threshold=50)
#'# getting low counts features at feature level
#'getLowCtsFeatures(ampliPanel, level="feature", threshold=50)
#'# exploring amplicon attribute values for a particular gene
#'g<-plotGeneAttrPerFeat(ampliPanel, geneID="gene4")
#'# adjust text size
#'g<-g+theme(title=element_text(size=16), axis.title=element_text(size=16),
#'legend.text=element_text(size=14))
#'if(interactive()){
#'x11(type="cairo");g
#'}
#'##Obtain the pileup matrix for the first amplicon
#'bed<-getBedFile(ampliPanel)[1]
#'## extracting the pileup matrix
#'myCounts<-pileupCounts(bed, bamFile, fastaFile)
#'head(myCounts)
#'
#'# getting and exploring a sequenced region of a particular gene
#'getRegion(ampliPanel, level="gene", ID="gene7", collapse=FALSE)
#'# plot a particular genomic region
#'g<-plotRegion(ampliPanel,region=c(4500,6800), seqname="chr10", SNPs=TRUE,  
#'xlab="", title="gene7 amplicons",size=0.5)
#'if(interactive()){
#'x11(type="cairo");g
#'}
#'# exploring the read count profile for a particular amplicon
#'g<-plotFeature(ampliPanel, featureID="AMPL20")
#'# x11(type="cairo")
#'if(interactive()){
#'g
#'}
#'# exploring the nucleotide percentages compositions of the read counts for a 
#'# particular amplicon
#'g<-plotNtdPercentage(ampliPanel,featureID="AMPL20")
#'if(interactive()){
#'g
#'}
#'## Building the XLSX report
#'imageFile<-system.file("extdata", "plot.pdf", package="TarSeqQC",
#'mustWork=TRUE)
#'buildReport(ampliPanel, attributeThres, imageFile ,file="Results.xlsx")
setClass(Class="TargetExperiment", slots=list(
    bedFile ="GRanges",
    bamFile="BamFile", 
    fastaFile="FaFile",
    scanBamP= "ScanBamParam",
    pileupP= "PileupParam",
    featurePanel= "GRanges",
    genePanel = "GRanges",
    feature = "character",
    attribute = "character"),
validity=function(object){

    ## Check bedFile
    bedFile<-getBedFile(object)
    if(class(bedFile)[[1]] != "GRanges"){
        stop(paste("The bedFile slot must be a GRanges object"))
    }
    ## Check bam file
    bamFile<-getBamFile(object)
    if(!( file.exists(path(bamFile))) & length(path(bamFile)
        ) > 1 ){
        stop(paste("The BAM file ", path(bamFile), " does not exist", 
        sep=""))
    }
    if(length(index(bamFile)) & length(path(bamFile)) > 1 ){
        warning("The index of the BAM file does not exist and will be build", 
            immediate. =TRUE)
        indexBam(path(bamFile))
    }
    ## Check fasta file
    fastaFile<-getFastaFile(object)
    if(!( file.exists(path(fastaFile))) & length(path(fastaFile)
        ) > 1 ){
        stop(paste("The Fasta file ", path(fastaFile), " does not exist",
            sep=""))
    }
    if(length(index(fastaFile)) & length(path(fastaFile)
        ) > 1 ){
        warning("The index of the Fasta file does not exist and will be build", 
            immediate. =TRUE)
        indexFa(path(fastaFile))
    }
    ## Check attribute
    attribute<-getAttribute(object)
    if(length(attribute) > 0 ){
        if(length(strsplit(attribute,"")[[1]]) > 1){
            if(!(attribute %in% c("coverage", "medianCounts"))){
                stop("The 'attribute' parameter should be 'coverage' or 
                    'medianCounts'")
            }
        }
    }
}, prototype=list(
    bedFile =GRanges(),
    bamFile=BamFile(file=character(1)), 
    fastaFile=FaFile(file=character(1)),
    scanBamP= ScanBamParam(),
    pileupP= PileupParam(),
    featurePanel= GRanges(),
    genePanel = GRanges(),
    feature = character(0),
    attribute = character(0))
)
