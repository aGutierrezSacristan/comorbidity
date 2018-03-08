#' Plot the comorbidity analysis results in a heatmap.
#'
#' Given an object of class \code{cAnalysis} or \code{molecularcAnalysis} obtained from a comorbidity analysis, 
#' a heatmap is obtained. 
#'
#' @param input A \code{cAnalysis} or \code{molecularcAnalysis} object, obtained 
#' by applying the \code{comorbidityAnalysis} or \code{comorbidityAnalysisMolecular} 
#' function
#' @param selectValue By default \code{"score"} variable will be selected. Change
#' it to any of the other possible variables (\code{'correctedPvalue'},\code{'odds ratio'}, 
#' \code{'phi'}, \code{'rr'}).  
#' @param cutOff By default \code{'0.05'}. The value of the argument can be changed 
#' to any other numeric variable, according to the range of the selected value.
#' @param npairs by default \code{'0'}.  The value of the argument can be changed
#' to any other numeric variable to show in the network only those comorbidities 
#' suffered by at least \code{npairs} of patients.
#' @param interactive Determines if the output heatmap is interactive or not. 
#' By default the \code{interactive} argument is set up as \code{FALSE}. The value 
#' of the argument can be changed to \code{TRUE}, as a result an interactive 
#' heatmap will be obtained.
#' @param lowColor Determines the heatmap color for the lowest value.
#' By default it is set to \code{"#0000FF"}. 
#' @param highColor Determines the heatmap color for the highest value.
#' By default it is set to \code{"yellow"}. 
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return A heatmap
#' @examples
#' load(system.file("extdata", "comorMale.RData", package="comoRbidity"))
#' htmp <- heatmapPlot( input = comorMale, 
#'               selectValue  = "score", 
#'               cutOff       = 0.5, 
#'               npairs       = 2
#'               )
#' htmp
#' @export heatmapPlot

heatmapPlot <- function( input , selectValue = "score", cutOff = 0.05, npairs = 0, interactive = FALSE, lowColor = "#0000FF", highColor = "yellow", verbose = FALSE ) {

    message("Checking the input object")
    checkClass <- class(input)[1]
    
    if(checkClass != "cAnalysis" & checkClass != "molecularcAnalysis"){
        message("Check the input object. Remember that this
                    object must be obtained after applying the query
                    function to your input file. The input object class must
                    be:\"cAnalysis\" or \"molecularcAnalysis\"")
        stop()
    }
    
    if(class(input)[1]== "cAnalysis"){
        obj <- input@result
        obj <- obj[as.numeric(obj$AB) >= npairs, ]
        column <- which(colnames(obj )==selectValue)
        obj$value <- obj[,column]
        
        if( selectValue == "correctedPvalue"  | selectValue == "fisher"){
            obj  <- obj [obj$value <= cutOff,]
        }else{
            obj  <- obj [obj$value >= cutOff,]
        }
        
       
        
        casted_rr <- reshape::cast(obj , disAcode ~ disBcode, value= selectValue  )
        l <- dim(casted_rr)[2]        
        Rowm  <- rowMeans(casted_rr[2:l], na.rm = T)       
        casted_rr<- cbind(casted_rr, Rowm)
        ordered<-  casted_rr[order(casted_rr["Rowm"]), "disAcode" ]          
        obj$disAcode  <- factor(obj $disAcode  , levels= as.factor(ordered))
        m <- max(obj [,column])
        n <- min(obj [,column])
        
        if( interactive == FALSE ){
            p <- ggplot2::ggplot(obj , ggplot2::aes ( disAcode, disBcode ) ) +
                ggplot2::geom_tile(ggplot2::aes(fill = value), colour = "black") +
                ggplot2::scale_fill_gradient(limits = c(n,m), low = lowColor,   high = highColor, na.value = "black") +
                ggplot2::theme_grey(base_size = 9) +
                ggplot2::labs ( title = "Comorbidity between diseases", x = "diagnosis code under study", y = "disease comorbidities") +
                ggplot2::scale_x_discrete(expand = c(0, 0)) +
                ggplot2::theme( plot.margin = grid::unit ( x = c ( 5, 15, 5, 15 ), units = "mm" ),
                                axis.text.x = ggplot2::element_text ( angle = 45, size = 11, hjust = 1 ), panel.background = ggplot2::element_blank() )
            p
        }else if(interactive == TRUE){
            tt <- reshape2::acast(obj, disBcode~disAcode, value.var=selectValue)
            d3heatmap::d3heatmap(tt, dendrogram = "none", scale = "none",
                                 colors =c("mediumorchid4","black"), 
                                 yaxis_font_size = 10,
                                 xaxis_font_size = 10
            )
        }
    }
    else if(class(input)[1]=="molecularcAnalysis"){
        obj <- input@result
        column <- which(colnames(obj )==selectValue)
        obj$value <- obj[,column]
        obj  <- obj [as.numeric(obj$value) > cutOff,]
        
        casted_rr <- reshape::cast(obj , V1 ~ V2, value= selectValue  )
        l <- dim(casted_rr)[2]        
        Rowm  <- rowMeans(casted_rr[2:l], na.rm = T)       
        casted_rr<- cbind(casted_rr, Rowm)
        ordered<-  casted_rr[order(casted_rr["Rowm"]), "V1" ]          
        obj$V1  <- factor(obj $V1  , levels= as.factor(ordered))
        m <- max(obj [,column])
        n <- min(obj [,column])
        
        if( interactive == FALSE ){
            p <- ggplot2::ggplot(obj , ggplot2::aes ( V1, V2 ) ) +
                ggplot2::geom_tile(ggplot2::aes(fill = value), colour = "black") +
                ggplot2::scale_fill_gradient(limits = c(n,m), low = lowColor,   high = highColor, na.value = "black") +
                ggplot2::theme_grey(base_size = 11) +
                ggplot2::labs ( title = "Molecular comorbidity between diseases", x = "diseases", y = "disease comorbidities") +
                ggplot2::scale_x_discrete(expand = c(0, 0)) +
                ggplot2::theme( plot.margin = grid::unit ( x = c ( 5, 15, 5, 15 ), units = "mm" ),
                                axis.text.x = ggplot2::element_text ( angle = 45, size = 11, hjust = 1 ), panel.background = ggplot2::element_blank() )
            p
        }else if(interactive == TRUE){
            tt <- reshape2::acast(obj, V2~V1, value.var=selectValue)
            d3heatmap::d3heatmap(tt, dendrogram = "none", scale = "none",
                                 colors =c("mediumorchid4","black"), 
                                 yaxis_font_size = 10,
                                 xaxis_font_size = 10
            )
        }
    }

}





