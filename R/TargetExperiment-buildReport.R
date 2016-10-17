#'Build an excel report of the Target Experiment.
#'
#'\code{buildReport} builds an excel file containing some statistical results.
#' These are computed to the selected attribute (e.g. 'coverage') along features
#' (e.g. 'amplicon') and genes. If 'imageFile' is null, the graph generated
#' calling the generic plot function will be used.
#'
#'@param file Character indicating the name of the report.
#'
#'@return NULL.
#'
#'@include TargetExperiment-addStatSummSheet.R 
#'@exportMethod buildReport
#'@docType methods
#'@name buildReport
#'@inheritParams addStatSummSheet
#'@rdname TargetExperiment-buildReport
#'@aliases buildReport-methods
#'@examples
#'## Loading the TargetExperiment object
#'data(ampliPanel,package="TarSeqQC")
#'# definition of the interval extreme values
#'attributeThres<-c(0,1,50,200,500, Inf)
#'
#'## Building the XLSX report
#'imageFile<-system.file("extdata", "plot.pdf", package="TarSeqQC",
#'mustWork=TRUE)
#'buildReport(ampliPanel, attributeThres=attributeThres, imageFile=imageFile,
#'file="results.xlsx")
#'

setGeneric(name="buildReport", def=function(object,attributeThres=c(0, 1, 50,
200, 500, Inf),imageFile=NULL, file="Results.xlsx"){
    standardGeneric("buildReport")
})
#'
#'@name buildReport
#'@rdname TargetExperiment-buildReport
#'@import openxlsx 
#'@aliases buildReport,TargetExperiment-method
#'@inheritParams buildReport
setMethod(f="buildReport", signature="TargetExperiment", definition=function(
object, attributeThres=c(0, 1, 50, 200, 500, Inf),
imageFile=NULL, file="Results.xlsx"){
    colors.fill<-colorRampPalette(c("red", "green"))(length(attributeThres)-1)
    featurePanel<-getFeaturePanel(object)
    genePanel<-getGenePanel(object)
    attribute<-getAttribute(object)
    feature<-getFeature(object)
    pool<-"pool" %in% names(mcols(featurePanel))
    df_panel<-as.data.frame(cbind(chr=as.data.frame(featurePanel)[,1],
        names(featurePanel), mcols(featurePanel)[,1, drop=FALSE], 
        as.data.frame(featurePanel)[,-c(1, which(names(as.data.frame(
        featurePanel)) == names(mcols(featurePanel)[,1,  drop=FALSE]))) ]))
    names(df_panel)[2]<-feature
    if(attribute== "coverage") {
        rm_cols<-which(names(df_panel) %in% c("medianCounts", "IQRCounts"))
    }else{
        rm_cols<-which(names(df_panel) %in% c("coverage", "sdCoverage"))
    }
    df_panel<-df_panel[,-rm_cols]
    gene_panel<-cbind(chr=as.data.frame(genePanel)[,1], gene=names(genePanel),
        as.data.frame(genePanel)[,-1])
    if(attribute== "coverage") {
        rm_cols<-which(names(gene_panel) %in% c("medianCounts", "IQRCounts"))
    }else{
        rm_cols<-which(names(gene_panel) %in% c("coverage", "sdCoverage"))
    }
    gene_panel<-gene_panel[,-rm_cols]
    if(is.null(imageFile)){
        g<-plot(object, attributeThres=attributeThres)
        imageFile<-"plotTargetExp.png"
        ggplot2::ggsave(g,filename=imageFile)
        
    }
    # creating the xlsx file
    wb <- createWorkbook()
    # create sheet1
    wb<-addStatSummSheet(object, wb,attributeThres, imageFile=imageFile )

    # creating sheet2
    sheet2<- addWorksheet(wb, "genes")
    headerStyles<-createStyle(fontSize=12, textDecoration="bold")

    writeDataTable(wb,x=gene_panel, sheet=sheet2, colNames=TRUE, rowNames=FALSE,
        startRow=1, startCol=1,withFilter =FALSE, tableStyle = "None", 
        headerStyle=headerStyles)

    # coloring cells
    gene_panel[,"score"]<-cut(gene_panel[,attribute], breaks=attributeThres,
        include.lowest=TRUE, right=FALSE, dig.lab = 6)

    cov_intervals<-levels(gene_panel[,"score"])
    sapply(1:length(cov_intervals), function(no_interval){
        rows<-which(gene_panel[,"score"] == cov_intervals[no_interval])+1
        cols<-rep(which(names(gene_panel) == attribute), length(rows))
        style<-createStyle(fgFill=colors.fill[no_interval], )  
        addStyle(wb,sheet=sheet2, style=style, rows=rows, cols= cols) 

    })  
    setColWidths(wb, sheet2, cols=1:ncol(gene_panel), widths="auto")

    sheet3<- addWorksheet(wb, feature)

    writeDataTable(wb,x=df_panel, sheet=sheet3, colNames=TRUE, rowNames=FALSE,
        startRow=1, startCol=1,withFilter =FALSE, tableStyle = "None", 
        headerStyle=headerStyles)
    df_panel[,"groups"]<-cut(df_panel[,attribute], breaks=attributeThres,
        include.lowest=TRUE, right=FALSE,dig.lab = 6)

    # coloring cells

    intervals<-levels(df_panel[,"groups"])

    sapply(1:length(intervals), function(no_interval){
        rows<-which(df_panel[,"groups"] == intervals[no_interval])+1
        cols<-rep(which(names(df_panel) == attribute), length(rows))
        style<-createStyle(fgFill=colors.fill[no_interval], )  
        addStyle(wb,sheet=sheet3, style=style, rows=rows, cols= cols) 
    })  

    setColWidths(wb, sheet3, cols=1:ncol(df_panel), widths="auto")

        #determining last feature of each gene
        gene_names<-unique(df_panel[,"gene"])
        gene_index<-sapply(gene_names, function(x){
        index<-which(df_panel[,"gene"] ==x)
        return(index[length(index)])
    })
# changing the format of amplicon coverage cells
# adding a border line to the last amplicon
    sapply(1:length(gene_names), function(x){

        index<-which(df_panel[,"gene"] == gene_names[x])
        row.genes <-rep(1+index[length(index)], ncol(df_panel)-1)
        col.genes <- 1:(ncol(df_panel)-1)
        color.fill<-colors.fill[which(df_panel[index[length(index)],
            "groups"] == levels(df_panel[,"groups"]))]
        geneStyle<-createStyle(border="Bottom",borderStyle="medium")
        covStyle<-createStyle(border="Bottom",fgFill=color.fill, 
            borderStyle="medium")
        col.cov<-which(names(df_panel) == attribute)
        addStyle(wb,sheet=sheet3, style=geneStyle, rows=row.genes ,
            cols= col.genes)
        addStyle(wb,sheet=sheet3, style=covStyle, rows=unique(row.genes) , 
            cols= col.cov)
    })
    #save the excel file
    saveWorkbook(wb, file=file, overwrite = TRUE)
})
