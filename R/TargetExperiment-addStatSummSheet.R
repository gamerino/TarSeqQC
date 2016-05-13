#'Build excel report of the Target Experiment.
#'
#'\code{addStatSummSheet} adds the statistics summary sheet to the workbook
#'that contains the Target Experiment Report.
#'
#'@param object TargetExperiment class object.
#'@param wb A workbook object that will contain the report.
#'@param attributeThres Numeric indicating the intervals extreme values.
#'@param imageFile Character indicating the name of the file that contains
#'the plot that could be insert in the report. 
#'
#'@return Workbook object.
#'
#'@include TargetExperiment-ggplotColours.R
#'@exportMethod addStatSummSheet
#'@docType methods
#'@name addStatSummSheet
#'@rdname TargetExperiment-buildReport
#'@aliases addStatSummSheet-methods
#'@seealso \code{\link{TargetExperiment-class}}
#'@note see full example in \code{\link{TargetExperiment-class}}
#'@author Gabriela A. Merino \email{gmerino@@bdmg.com.ar}, Cristobal Fresno
#'\email{cfresno@@bdmg.com.ar} and Elmer A. Fernandez 
#'\email{efernandez@@bdmg.com.ar}

setGeneric(name="addStatSummSheet", def=function(object, wb,
attributeThres=c(0,1,50,200,500, Inf), imageFile){
    standardGeneric("addStatSummSheet")
})
#'
#'@name addStatSummSheet
#'@rdname TargetExperiment-buildReport 
#'@import openxlsx 
#'@inheritParams addStatSummSheet
#'@aliases addStatSummSheet,TargetExperiment-method
setMethod(f="addStatSummSheet", signature=signature(object="TargetExperiment"),
definition=function(object, wb, attributeThres=c(0, 1, 50, 200, 500, Inf),
imageFile){
featurePanel<-getFeaturePanel(object)
    genePanel<-getGenePanel(object)
    pool<-any("pool" == names(mcols(featurePanel)))
#     colors.fill<-ggplotColours(object,n=(length(attributeThres)-1))
    colors.fill<-colorRampPalette(c("red", "green"))(length(attributeThres)-1)
    statistics_summary<-as.data.frame(rbind(summaryGeneLev(object), 
        summaryFeatureLev(object)))
    #   Computing frequency tables for genes and feature
    ftable.ampli <- summaryIntervals(object=object,
        attributeThres=attributeThres)
    object2<-object
    object2@featurePanel<-getGenePanel(object)
    setFeature(object2)<-"gene"
    ftable.genes <- summaryIntervals(object=object2,
        attributeThres=attributeThres)
    sheet_name="Summary"
    sheet1<-addWorksheet(wb, sheet_name)
    df_panel<-as.data.frame(featurePanel)
    attribute<-getAttribute(object)
    feature<-getFeature(object)
    df_panel[,"score"]<-cut(df_panel[,attribute], breaks=attributeThres, 
        include.lowest=TRUE, right=FALSE,dig.lab = 6)
    first_feat_ids<-data.frame(sapply(levels(df_panel[,"score"]), function(x){
        aux<-row.names(df_panel[df_panel[,"score"]== x,])[1:5]
        out<-paste(aux[!is.na(aux)], sep="", collapse=", ")
        if(length(row.names(df_panel[df_panel[,"score"]== x,])) > 5){
            out<-paste(out, "...", sep=", ")
        }
        return(out)
    }))
    names(first_feat_ids)[1]<-paste(feature, "_id", sep="")
    gene_panel<-as.data.frame(genePanel)
    gene_panel[,"score"]<-cut(gene_panel[,attribute], breaks=attributeThres, 
        include.lowest=TRUE, right=FALSE,dig.lab = 6)
    first_gene_ids<-data.frame(gene_id=sapply(levels(gene_panel[,"score"]),
    function(x){
        aux<-row.names(gene_panel[gene_panel[,"score"]== x,])[1:5]
        out<-paste(aux[!is.na(aux)], sep="", collapse=", ")
        if(length(row.names(gene_panel[gene_panel[,"score"]== x,])) > 5){
            out<-paste(out, "...", sep=", ")
        }
        return(out)
    }))
    # adding table titles
    titleStyle<-createStyle(fontSize=14, textDecoration = c("bold", "italic"),
        fgFill="#FFFF99")
    
    mergeCells(wb, sheet=sheet_name, cols=1:(ncol(statistics_summary)+1), 
        rows=1)
    headerStyles<-createStyle(fontSize=12, textDecoration="bold",
        halign="center", valign="center")
    writeData(wb,x=paste("statistics_summary for attribute ", attribute, 
        sep=""), sheet=sheet_name, colNames=TRUE, rowNames=FALSE,startRow=1, 
        startCol=1,  withFilter =FALSE)

    mergeCells(wb, sheet=sheet_name, cols=1:ncol(ftable.genes), 
        rows= 5+nrow(statistics_summary))
    writeData(wb,x="gene_level", sheet=sheet_name, colNames=TRUE, 
        rowNames=FALSE, startRow=5+nrow(statistics_summary), startCol=1,
        withFilter =FALSE)
    mergeCells(wb, sheet=sheet_name, cols=1:ncol(ftable.ampli), rows=(10+
        nrow(statistics_summary)+nrow(ftable.genes)))
    writeData(wb,x="amplicon_level", sheet=sheet_name, colNames=TRUE, 
        rowNames=FALSE, startRow=(10+nrow(statistics_summary)+nrow(
        ftable.genes)), startCol=1,withFilter =FALSE)
    addStyle(wb,sheet=sheet_name, style=titleStyle, rows=c(1,5+nrow(
        statistics_summary),(10+nrow(statistics_summary)+nrow(ftable.genes))), 
        cols=rep(1, 3))

    #adding statistics_summary table
    writeData(wb,x=rownames(statistics_summary), sheet=sheet_name,
        colNames=FALSE, rowNames=FALSE,startRow=3, startCol=1,withFilter =FALSE)
    writeDataTable(wb,x=statistics_summary, sheet=sheet_name, colNames=TRUE,
        rowNames=FALSE,startRow=2, startCol=2,withFilter=FALSE, 
        tableStyle="None", headerStyle=headerStyles)
    
    # adding frequency tables
    mergeCells(wb, sheet=sheet_name, cols=2:ncol(ftable.genes), rows=6+nrow(
        statistics_summary))
    writeData(wb,x="freq", sheet=sheet_name, colNames=FALSE, rowNames=FALSE,
        startRow=6+nrow(statistics_summary), startCol=2,withFilter =FALSE)
    addStyle(wb,sheet=sheet_name, style=headerStyles, rows=6+
        nrow(statistics_summary), cols= 2)
    writeDataTable(wb,x= ftable.genes, sheet=sheet_name, colNames=TRUE, 
        rowNames=FALSE, startRow=7+nrow(statistics_summary), startCol=1,
        withFilter=FALSE, tableStyle = "None", headerStyle=headerStyles)

    mergeCells(wb, sheet=sheet_name, cols=2:ncol(ftable.genes), rows=11+nrow(
        statistics_summary)+nrow(ftable.genes))
    writeData(wb,x="freq", sheet=sheet_name, colNames=FALSE, rowNames=FALSE,
        startRow=11+nrow(statistics_summary)+nrow(ftable.genes), startCol=2,
        withFilter =FALSE)
    addStyle(wb,sheet=sheet_name,style=headerStyles, rows=11+nrow(
        statistics_summary)+nrow(ftable.genes), cols= 2)
    writeDataTable(wb,x= ftable.ampli, sheet=sheet_name, colNames=TRUE, 
        rowNames=FALSE, startRow=(12+nrow(statistics_summary)+nrow(
        ftable.genes)), startCol=1, withFilter =FALSE, tableStyle = "None", 
        headerStyle=headerStyles)
    setColWidths(wb, sheet=sheet_name, cols=c(1:(ncol(statistics_summary)-1), 
        ncol(statistics_summary)+1), widths="auto")
    setColWidths(wb, sheet=sheet_name, cols=6, widths=10)

    mergeCells(wb, sheet=sheet_name, cols=1:3+ncol(ftable.genes), rows=8+nrow(
        statistics_summary))
    writeData(wb, x=first_gene_ids, sheet=sheet_name, colNames=TRUE,
        startRow=7+nrow(statistics_summary), startCol=1+ncol(ftable.genes), 
        withFilter=FALSE, headerStyle=headerStyles)
    
    mergeCells(wb, sheet=sheet_name, cols=1:3+ncol(ftable.ampli), rows=13+nrow(
        statistics_summary+nrow(ftable.genes)))
    writeData(wb, x=first_feat_ids, sheet=sheet_name,  colNames=TRUE,
        startRow=12+nrow(statistics_summary)+nrow(ftable.genes),startCol=1+ncol(
        ftable.ampli), withFilter=FALSE, headerStyle=headerStyles)
    insertImage(wb, sheet=sheet_name,file=imageFile,width=10, height=10,
        startRow=(15+nrow(statistics_summary)+nrow(ftable.genes)+nrow(
        ftable.ampli)), startCol=1)

    sapply(1:length(colors.fill), function(x){
        row_freq_tab<-rep(c((7+nrow(statistics_summary)+x),
        (12+nrow(statistics_summary)+nrow(ftable.genes)+x)),
        length(1:ncol(ftable.genes)))
        col_tables<-c(sapply(1:ncol(ftable.genes), function(x){
            rep(x, length(unique(row_freq_tab)))
        }))
        style<-createStyle(fgFill=colors.fill[x])  
        addStyle(wb,sheet=sheet1, style=style, rows=row_freq_tab, 
            cols= col_tables)
    })
    return(wb)
})
