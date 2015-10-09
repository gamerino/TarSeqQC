# rm(list=ls())
# library("RUnit")
# library("TarSeqQC")
# library("Rsamtools")
library(BiocParallel)

##-----------------------------------------------------------------------------
##TargetExperiment-class Tests 
##-----------------------------------------------------------------------------
##Empty object test: Does TargetExperiment work without any parameter?
test_TargetExperiment<-function(){
    checkTrue(validObject(TargetExperiment()), 
    msg="TargetExperiment works without any parameter: OK.")  
}

##Build TargetExperiment object with user data
test_TargetExperimentWithData<-function(){
    # Defining bam file, bed file and fasta file names and paths
    bamFile<-system.file("extdata", "mybam.bam", package="TarSeqQC", 
        mustWork=TRUE)
    bedFile<-system.file("extdata", "mybed.bed", package="TarSeqQC",
        mustWork=TRUE)
    fastaFile<-system.file("extdata", "myfasta.fa", package="TarSeqQC", 
        mustWork=TRUE)
    attribute<-"coverage"    
    feature<-"amplicon"
    object<-TargetExperiment(bedFile=bedFile, bamFile=bamFile, 
        fastaFile=fastaFile, feature=feature, attribute=attribute)
    # Checking slots
        #bedFile
        checkEquals(class(getBedFile(object))[1],"GRanges", 
            msg="TargetExperiment bedFile slot type: OK.")
        #bamFile
        checkTrue(file.exists(bamFile),  
            msg="TargetExperiment bamFile existence: OK.")
        checkTrue(class(getBamFile(object))=="BamFile",  
            msg="TargetExperiment bamFile slot type: OK.")
        #fastaFile
        checkTrue(file.exists(fastaFile),  
            msg="TargetExperiment fastaFile existence: OK.")
        checkEquals(class(getFastaFile(object))[1],"FaFile",  
            msg="TargetExperiment fastaFile slot type: OK.")
        #scanBamP
        checkEquals(class(getScanBamP(object))[1],"ScanBamParam",
            msg="TargetExperiment scanBamP slot type= OK.")
        aux<-lapply(levels(seqnames(getBedFile(object))), function(x){
            ranges(getBedFile(object)[seqnames(getBedFile(object)) == x])
        })
        names(aux)<-levels(seqnames(getBedFile(object)))
                
        checkEquals(as.list(bamWhich(getScanBamP(object))),aux,
            msg="Which parameter setting= OK.")
        #pileupP
        checkEquals(class(getPileupP(object))[1],"PileupParam",
            msg="TargetExperiment pileupP slot type= OK.")
        #featurePanel
        checkEquals(class(getFeaturePanel(object))[1],"GRanges",
            msg="TargetExperiment featurePanel slot type= OK.")
        checkEquals(names(getFeaturePanel(object)),names(getBedFile(object)),
            msg="TargetExperiment featurePanel names= OK.")
        checkTrue(all(c("medianCounts", "IQRCounts", "coverage", 
            "sdCoverage") %in% names(mcols(object@featurePanel))),
            msg="TargetExperiment featurePanel metadata= OK.")
        checkEquals(length(getFeaturePanel(object)),length(getBedFile(object)),
            msg="TargetExperiment featurePanel dimension= OK.")
        #genePanel
        checkEquals(class(getGenePanel(object))[1],"GRanges",
            msg="TargetExperiment genePanel slot type= OK.")
        checkTrue(all(names(getGenePanel(object)) %in% unique(mcols(getBedFile(
            object))[,"gene"])), msg="TargetExperiment genePanel names= OK.")
        checkTrue(all(c("medianCounts", "IQRCounts", "coverage", 
            "sdCoverage") %in% names(mcols(object@genePanel))),
            msg="TargetExperiment featurePanel metadata= OK.")
        #feature
        checkTrue(is.character(getFeature(object)),
            msg="TargetExperiment feature slot type= OK.")
        #attribute
        checkTrue(is.character(getAttribute(object)),
            msg="TargetExperiment attribute slot type= OK.")
        checkTrue(any(getAttribute(object) %in% c("medianCounts", "coverage")),
            msg="TargetExperiment attribute slot value= OK.")
                  
    
}
##Getters/setters: check getBedFile
test_getBedFile<-function(){
    data(ampliPanel, package="TarSeqQC")
    checkEquals(class(getBedFile(ampliPanel))[1], "GRanges",
    msg="bedFile getter: OK.")
}
##Getters/setters: check setBedFile
test_setBedFile<-function(){
    data(ampliPanel, package="TarSeqQC")
    setBamFile(ampliPanel)<-system.file("extdata", "mybam.bam", 
        package="TarSeqQC", mustWork=TRUE)
    setBedFile(ampliPanel)<-system.file("extdata", "mybed.bed", 
        package="TarSeqQC", mustWork=TRUE)
    setFastaFile(ampliPanel)<-system.file("extdata", "myfasta.fa", 
        package="TarSeqQC", mustWork=TRUE)
    bedFile<-system.file("extdata", "mybed.bed", 
        package="TarSeqQC", mustWork=TRUE)
    #bedFile must be setted starting from a character object containing a valid
    #file path and name
    checkException(setBedFile(ampliPanel)<-"", 
        msg="bedFile setter: OK.", silent=TRUE)
    #bedFile must be setted starting from a character object
    setBedFile(ampliPanel)<-bedFile
    checkTrue(validObject(ampliPanel),msg="bedFile setter: OK.")
}

##Getters/setters: check getBamFile
test_getBamFile<-function(){
    data(ampliPanel, package="TarSeqQC")
    checkEquals(class(getBamFile(ampliPanel))[1], "BamFile",
        msg="bamFile getter: OK.")
}
##Getters/setters: check setBamFile
test_setBamFile<-function(){
    data(ampliPanel, package="TarSeqQC")
    bamFile<-system.file("extdata", "mybam.bam", 
        package="TarSeqQC", mustWork=TRUE)                
    #bamFile cannot be defined as a non-existing file
    checkException(setBamFile(ampliPanel)<-"", 
        msg="bamFile setter: OK.", silent=TRUE)
    #bamFile must be setted starting from a character object containing a valid
    #file path and name
    setBamFile(ampliPanel)<-bamFile
    checkTrue(validObject(ampliPanel),msg="bamFile setter: OK.")
}
##Getters/setters: check getFastaFile
test_getFastaFile<-function(){
    data(ampliPanel, package="TarSeqQC")
    checkEquals(class(getFastaFile(ampliPanel))[1], "FaFile", 
        msg="fastaFile getter: OK.")
}
##Getters/setters: check setFastaFile
test_setFastaFile<-function(){
    data(ampliPanel, package="TarSeqQC")
    #fastaFile cannot be defined as a non-existing file
    checkException(setFastaFile(ampliPanel)<-"", 
        msg="fastaFile setter: OK.", silent=TRUE)
    
    #file path and name
    fastaFile<-system.file("extdata", "myfasta.fa", 
        package="TarSeqQC", mustWork=TRUE)
    setFastaFile(ampliPanel)<-fastaFile
    checkTrue(validObject(ampliPanel),msg="fastaFile setter: OK.")
}
##Getters/setters: check getScanBamP
test_getScanBamP<-function(){
    data(ampliPanel, package="TarSeqQC")
    checkEquals(class(getScanBamP(ampliPanel))[1], "ScanBamParam", 
        msg="scanBamP getter: OK.")
}
##Getters/setters: check setScanBamP
test_setScanBamP<-function(){
    data(ampliPanel, package="TarSeqQC")
    #scanBamP should be ScanBamParam class
    checkException(setScanBamP(ampliPanel)<-"", 
        msg="scanBamP setter: OK.", silent=TRUE)
    #scanBamP must be setted using a ScanBamParam object
    setScanBamP(ampliPanel)<-ScanBamParam()
    checkTrue(validObject(ampliPanel),msg="scanBamP setter: OK.")
}
##Getters/setters: check getPileupP
test_getPileupP<-function(){
    data(ampliPanel, package="TarSeqQC")
    checkEquals(class(getPileupP(ampliPanel))[1], "PileupParam", 
        msg="pileupP getter: OK.")
}
##Getters/setters: check setPileupP
test_setPileupP<-function(){
    data(ampliPanel, package="TarSeqQC")
    #pileupP must be a PileupParam object
    checkException(setPileupP(ampliPanel)<-"", 
        msg="pileupP setter: OK.", silent=TRUE)
    #pileupP must be setted using a PileupParam object
    setPileupP(ampliPanel)<-PileupParam()
    checkTrue(validObject(ampliPanel),msg="pileupP setter: OK.")
}
##Getters/setters: check getFeaturePanel
test_getFeaturePanel<-function(){
    data(ampliPanel, package="TarSeqQC")
    checkEquals(class(getFeaturePanel(ampliPanel))[1], "GRanges",
        msg="featurePanel getter: OK.")
}
##Getters/setters: check setFeaturePanel
test_setFeaturePanel<-function(){
    data(ampliPanel, package="TarSeqQC")
    #featurePanel cannot be setted as a data.frame
    checkException(ampliPanel@featurePanel<-data.frame(), 
        msg="featurePanel setter: OK.", silent=TRUE)
    #featurePanel must be setted starting from a GRanges object
    setFeaturePanel(ampliPanel)<-GRanges()
    checkTrue(validObject(ampliPanel),msg="featurePanel setter: OK.")
    checkEquals(length(getFeaturePanel(ampliPanel)), 0,
        msg="featurePanel setter: OK.")
}
##Getters/setters: check getGenePanel
test_getGenePanel<-function(){
    data(ampliPanel, package="TarSeqQC")
    checkEquals(class(getGenePanel(ampliPanel))[1],"GRanges",
        msg="genePanel getter: OK.")
}
##Getters/setters: check setGenePanel
test_setGenePanel<-function(){
    data(ampliPanel, package="TarSeqQC")
    setBamFile(ampliPanel)<-system.file("extdata", "mybam.bam", 
        package="TarSeqQC", mustWork=TRUE)
    setBedFile(ampliPanel)<-system.file("extdata", "mybed.bed", 
        package="TarSeqQC", mustWork=TRUE)
    setFastaFile(ampliPanel)<-system.file("extdata", "myfasta.fa", 
        package="TarSeqQC", mustWork=TRUE)
    #genePanel cannot be setted as a data.frame
    checkException(ampliPanel@genePanel<-data.frame(), 
        msg="genePanel setter: OK.", silent=TRUE)
    #genePanel cannot be setted starting from an empty GRanges object
    checkException(setGenePanel(ampliPanel)<-GRanges(), 
        msg="genePanel setter: OK.", silent=TRUE)
    setGenePanel(ampliPanel)<-summarizePanel(ampliPanel)
    checkTrue(validObject(ampliPanel), msg="genePanel setter: OK.")
}
##Getters/setters: check getFeature
test_getFeature<-function(){
    data(ampliPanel, package="TarSeqQC")
    feature<-"amplicon"
    checkEquals(getFeature(ampliPanel),feature, msg="feature getter: OK.")
}
##Getters/setters: check setFeature
test_setFeature<-function(){
    data(ampliPanel, package="TarSeqQC")
    feature<-getFeature(ampliPanel)
    checkEquals(getFeature(ampliPanel), feature, msg="feature setter: OK.")
    feature2<-factor(feature)    
    checkException(setFeature(ampliPanel)<-feature2, msg="feature setter: OK.", 
        silent=TRUE)
}
##Getters/setters: check getAttribute
test_getAttribute<-function(){
    data(ampliPanel, package="TarSeqQC")
    attribute<-"coverage"
    checkEquals(getAttribute(ampliPanel),attribute, msg="feature getter: OK.")
}
##Getters/setters: check setFeature
test_setAttribute<-function(){
    data(ampliPanel, package="TarSeqQC")
    attribute<-"medianCounts"
    setAttribute(ampliPanel)<-attribute
    checkEquals(getAttribute(ampliPanel), attribute, msg="feature setter: OK.")
    attribute2<-"mean"    
    checkException(setAttribute(ampliPanel)<-attribute2, 
        msg="feature setter: OK.", silent=TRUE)
      
}

##Test statistics:summaryFeatureLev
##
test_summaryFeatureLev<-function(){
    data(ampliPanel, package="TarSeqQC")
    checkTrue(is.matrix(summaryFeatureLev(ampliPanel)), 
        msg="summaryFeatureLev returned type: OK.")
    checkTrue(all(colnames(summaryFeatureLev(ampliPanel)) %in% c("Min.", 
        "1st Qu.", "Median","Mean", "3rd Qu.", "Max.")), 
        msg="summaryFeatureLev returned matrix colnames: OK.")
        
}
##Test statistics:summaryGeneLev
test_summaryGeneLev<-function(){
    data(ampliPanel, package="TarSeqQC")
    checkTrue(is.matrix(summaryGeneLev(ampliPanel)), 
        msg="summaryGeneLev returned type: OK.")
    checkTrue(all(colnames(summaryGeneLev(ampliPanel)) %in% c("Min.", 
        "1st Qu.", "Median","Mean", "3rd Qu.", "Max.")), 
        msg="summaryGeneLev returned matrix colnames: OK.")
}
##Test statistics:summaryIntervals
test_summaryIntervals<-function(){
    data(ampliPanel, package="TarSeqQC")
    checkTrue(is.data.frame(summaryIntervals(ampliPanel)), 
        msg="summaryIntervals returned type: OK.")
    checkTrue(all(colnames(summaryIntervals(ampliPanel)) %in% c(paste(
        getFeature(ampliPanel), getAttribute(ampliPanel), "intervals", sep="_"),
        "abs", "cum_abs","rel", "cum_rel")), 
        msg="summaryIntervals returned data.frame colnames: OK.")
    checkEquals(dim(summaryIntervals(ampliPanel)), c(5,5), 
        msg="summaryIntervals returned data.frame dimension: OK.")
    checkEquals(summaryIntervals(ampliPanel)[nrow(summaryIntervals(ampliPanel)),
        "cum_abs"], length(getFeaturePanel(ampliPanel)), 
        msg="summaryIntervals returned data.frame consistency: OK.")
}

##Test pileupCounts
test_pileupCounts<-function(){
    data(ampliPanel, package="TarSeqQC")
    bamFile<-system.file("extdata", "mybam.bam", 
        package="TarSeqQC", mustWork=TRUE)
    bed<-getBedFile(ampliPanel)
    fastaFile<-system.file("extdata", "myfasta.fa", 
        package="TarSeqQC", mustWork=TRUE)
    myCounts<-pileupCounts(bed=bed, bamFile=bamFile, fastaFile=fastaFile)
    checkEquals(class(myCounts), "data.frame", 
        msg="pileupCounts returned object type: OK.") 
    checkTrue(all(c("pos", "seqnames", "which_label", "counts") %in% 
        colnames(myCounts)), msg="pileupCounts returned data.frame colnames: 
        OK.")
    checkEquals(dim(myCounts), c(2165,12), 
        msg="pileupCounts returned data.frame dimension: OK.")
 
    auxData<-unique(as.character(myCounts[,"which_label"]))
    strsplit(auxData, ":")[1]
    checkTrue(all(unique(sapply(1:length(auxData),function(x){ strsplit(
        auxData[x], ":")[[1]][1]})) %in% c("chr1", "chr3", "chr7", "chr10")),
        msg="pileupCounts returned data.frame seqnames: OK.")
}
##Test buildFeaturePanel
test_buildFeaturePanel<-function(){
    if (.Platform$OS.type != "windows")
        BPPARAM <- MulticoreParam(2)
    else
        BPPARAM <- SerialParam()
 
    data(ampliPanel, package="TarSeqQC")
    setBamFile(ampliPanel)<-system.file("extdata", "mybam.bam", 
        package="TarSeqQC", mustWork=TRUE)
    setBedFile(ampliPanel)<-system.file("extdata", "mybed.bed", 
        package="TarSeqQC", mustWork=TRUE)
    setFastaFile(ampliPanel)<-system.file("extdata", "myfasta.fa", 
        package="TarSeqQC", mustWork=TRUE)

    myFeaturePanel<-buildFeaturePanel(ampliPanel, BPPARAM=BPPARAM)
    checkEquals(class(myFeaturePanel)[1], "GRanges", 
        msg="buildFeaturePanel returned object type: OK.") 
    checkTrue(all(colnames(mcols(myFeaturePanel)) %in% c("gene", 
        "medianCounts", "IQRCounts", "coverage", "sdCoverage")), 
        msg="buildFeaturePanel returned metadata colnames: OK.")
    checkEquals(length(myFeaturePanel), 29, 
        msg="buildFeaturePanels returned GRanges dimension: OK.")
    checkTrue(all(names(myFeaturePanel) %in% names(getFeaturePanel(ampliPanel
        ))), msg="buildFeaturePanel returned GRanges names: OK.")
}
##Test summarizePanel
test_summarizePanel<-function(){
    data(ampliPanel, package="TarSeqQC")
    myGenePanel<-summarizePanel(ampliPanel)
    checkEquals(class(myGenePanel)[1], "GRanges", 
        msg="summarizePanel returned object type: OK.") 
    checkTrue(all(colnames(mcols(myGenePanel)) %in% c("medianCounts", 
        "IQRCounts", "coverage", "sdCoverage")), 
        msg="summarizePanel returned metadata colnames: OK.")
    checkEquals(length(myGenePanel), 8, 
        msg="summarizePanel returned GRanges dimension: OK.")
    checkTrue(all(names(myGenePanel) %in% names(getGenePanel(ampliPanel))), 
        msg="summarizePanel returned GRanges names: OK.")
}
#### Test plot
##test_plot<-function(){
##    data(ampliPanel, package="TarSeqQC")
##    g<-plot(ampliPanel)
##    checkTrue(is.ggplot(g), msg="returned plot type: OK.")
##    checkEquals(dim(g$data), c(40,11), msg="returned plot data dimension: OK.")
##}
#### Test plotAttrExpl
##test_plotAttrExpl<-function(){
##    data(ampliPanel, package="TarSeqQC")
##    g<-plotAttrExpl(ampliPanel)
##    checkTrue(is.ggplot(g), msg="returned plot type: OK.")
##    checkEquals(dim(g$data), c(29,1), msg="returned plot data dimension: OK.")
##}
#### Test plotFeatPerform
##test_plotFeatPerform<-function(){
##    data(ampliPanel, package="TarSeqQC")
##    g<-plotFeatPerform(ampliPanel)
##    checkTrue(class(g)[1] == "gg" , msg="returned plot type: OK.")
##}
#### Test plotFeature
##test_plotFeature<-function(){
##    data(ampliPanel, package="TarSeqQC")
##    setBamFile(ampliPanel)<-system.file("extdata", "mybam.bam", 
##        package="TarSeqQC", mustWork=TRUE)
##    setBedFile(ampliPanel)<-system.file("extdata", "mybed.bed", 
##        package="TarSeqQC", mustWork=TRUE)
##    setFastaFile(ampliPanel)<-system.file("extdata", "myfasta.fa", 
##        package="TarSeqQC", mustWork=TRUE)
##    checkException(plotFeature(ampliPanel, featureID="nopresent"), 
##        msg="Testing feature ID: OK.", silent=TRUE)
##    featureID<-"AMPL20"
##    checkTrue(featureID %in% names(getFeaturePanel(ampliPanel)),
##        msg="featureID present in the featurePanel: OK.")
##    g<-plotFeature(ampliPanel, featureID)
##    checkTrue(is.ggplot(g), msg="returned plot type: OK.")
##    checkEquals(dim(g$data), c(63,12), msg="returned plot data dimension: OK.")
##}
#### Test plotGeneAttrPerFeat
##test_plotGeneAttrPerFeat<-function(){
##    data(ampliPanel, package="TarSeqQC")
##    checkException(plotGeneAttrPerFeat(ampliPanel, geneID="nopresent"), 
##        msg="Testing gene ID: OK.", silent=TRUE)
##    geneID<-"gene1"
##    checkTrue(geneID %in% names(getGenePanel(ampliPanel)),
##        msg="geneID present in the featurePanel: OK.")
##    g<-plotGeneAttrPerFeat(ampliPanel, geneID)
##    checkTrue(is.ggplot(g), msg="returned plot type: OK.")
##    checkEquals(dim(g$data), c(1,11), msg="returned plot data dimension: OK.")
##}
#### Test plotNtdPercentage
##test_plotNtdPercentage<-function(){
##    data(ampliPanel, package="TarSeqQC")
##    setBamFile(ampliPanel)<-system.file("extdata", "mybam.bam", 
##        package="TarSeqQC", mustWork=TRUE)
##    setBedFile(ampliPanel)<-system.file("extdata", "mybed.bed", 
##        package="TarSeqQC", mustWork=TRUE)
##    setFastaFile(ampliPanel)<-system.file("extdata", "myfasta.fa", 
##        package="TarSeqQC", mustWork=TRUE)
##    checkException(plotNtdPercentage(ampliPanel, featureID="nopresent"), 
##        msg="Testing feature ID: OK.", silent=TRUE)
##    featureID<-"AMPL20"
##    checkTrue(featureID %in% names(getFeaturePanel(ampliPanel)),
##        msg="featureID present in the featurePanel: OK.")
##    g<-plotNtdPercentage(ampliPanel, featureID)
##    checkTrue(is.ggplot(g), msg="returned plot type: OK.")
##    checkEquals(dim(g$data), c(252,3), msg="returned plot data dimension: OK.")
##}
#### Test plotRegion
##test_plotRegion<-function(){
##    data(ampliPanel, package="TarSeqQC")
##    setBamFile(ampliPanel)<-system.file("extdata", "mybam.bam", 
##        package="TarSeqQC", mustWork=TRUE)
##    setBedFile(ampliPanel)<-system.file("extdata", "mybed.bed", 
##        package="TarSeqQC", mustWork=TRUE)
##    setFastaFile(ampliPanel)<-system.file("extdata", "myfasta.fa", 
##        package="TarSeqQC", mustWork=TRUE)
##    region<-c(4500,6000)
##    seqname<-"chr10"
##    checkException(plotRegion(ampliPanel), 
##        msg="Testing function calling: OK.", silent=TRUE)
##    checkException(plotRegion(ampliPanel, seqname="chr0"), 
##        msg="Testing function calling: OK.", silent=TRUE)
##    checkException(plotRegion(ampliPanel, seqname="chr10"), 
##        msg="Testing function calling: OK.", silent=TRUE)
##    checkTrue(seqname %in% levels(seqnames(getBedFile(ampliPanel))),
##        msg="seqname present in the featurePanel: OK.")
##    g<-plotRegion(ampliPanel, region, seqname)
##    checkTrue(is.ggplot(g), msg="returned plot type: OK.")
##    checkEquals(dim(g$data), c(1501,12), 
##        msg="returned plot data dimension: OK.")
##}
#### Test xlsx report creation
##test_reportCreation<-function(){
##    data(ampliPanel, package="TarSeqQC")
##    imageFile<-system.file("extdata", "plot.pdf", package="TarSeqQC",
##        mustWork=TRUE)
##    buildReport(ampliPanel, imageFile=imageFile, file="test.xlsx")        
##    checkTrue(file.exists("test.xlsx"), msg="Xlsx file creation: OK.")
##}

##-----------------------------------------------------------------------------
##Test functions
##-----------------------------------------------------------------------------
##TargetExperiment class test
# test_TargetExperiment()
# test_TargetExperimentWithData()
# test_getBamFile()
# test_setBamFile()
# test_getBedFile()
# test_setBedFile()
# test_getFastaFile()
# test_setFastaFile()
# test_getScanBamP()
# test_setScanBamP()
# test_getPileupP()
# test_setPileupP()
# test_getFeature()
# test_setFeature()
# test_getAttribute()
# test_setAttribute()
# test_getFeaturePanel()
# test_setFeaturePanel()
# test_getGenePanel()
# test_setGenePanel()
# test_summaryFeatureLev()
# test_summaryGeneLev()
# test_summaryIntervals()
# test_pileupCounts()
# test_buildFeaturePanel()
# test_summarizePanel()
# test_plot()
# test_plotAttrExpl()
# test_plotFeatPerform()
# test_plotFeature()
# test_plotGeneAttrPerFeat()
# test_plotNtdPercentage()
# test_plotRegion()
# test_reportCreation()
