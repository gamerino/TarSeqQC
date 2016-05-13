#'TargetExperiment auxiliar function.
#'
#'\code{ggplotColours} is a function to know what color is used when  
#'ggplot is called.
#'
#'@param n amount of colors.
#'
#'@return colours
#'
#'@include TargetExperiment-plot.R
#'@exportMethod ggplotColours
#'@docType methods
#'@name ggplotColours
#'@rdname TargetExperiment-buildReport
#'@aliases ggplotColours-methods
setGeneric(name="ggplotColours", def=function(object,n){
    standardGeneric("ggplotColours")
})
#'
#'@name ggplotColours
#'@rdname TargetExperiment-buildReport
#'@inheritParams ggplotColours
#'@aliases ggplotColours,TargetExperiment-method
#'@examples
#'## Loading the TargetExperimentList object
#'data(ampliPanel, package="TarSeqQC")
#'colors<-ggplotColours(ampliPanel, n=5)
setMethod(f="ggplotColours", signature="TargetExperiment",
definition=function(object,n){
    
    h<-c(0, 360) +15
    if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
    return(hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65))
})
#'@name ggplotColours
#'@rdname TargetExperiment-buildReport
#'@inheritParams ggplotColours
#'@aliases ggplotColours,TargetExperimentList-method
#'@examples
#'## Loading the TargetExperimentList object
#'data(TEList, package="TarSeqQC")
#'colors<-ggplotColours(object, n=5)
setMethod(f="ggplotColours", signature="TargetExperimentList",
definition=function(object,n){
    
    h<-c(0, 360) +15
    if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
    return(hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65))
})
