#' Plot the directionality results in a heatmap.
#'
#' Given an \code{dataframe} with the directionality results obtained after apply
#' the \code{directionality} function, a heatmap plot of the results is obtained.  
#'
#' @param input A \code{data.frame} obtained after applying the \code{directionality}
#' function
#' @param interactive Determines if the output heatmap is interactive or not. 
#' By default the \code{interactive} argument is set up as \code{FALSE}. The value 
#' of the argument can be changed to \code{TRUE}, as a result an interactive 
#' heatmap will be obtained.
#' @param fromAtoBColor Determines the heatmap color for the direction A to B.
#' By default it is set to \code{"olivegreen"}. 
#' @param fromBtoAColor Determines the heatmap color for the direction B to A.
#' By default it is set to \code{"orange"}. 
#' @param noDirectionColor Determines the heatmap color when there is no preferred direction.
#' By default it is set to \code{"grey"}. 
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return A heatmap
#' @examples
#' load(system.file("extdata", "directionalityExample.RData", package="comoRbidity"))
#' htmpDirection <- heatmapDirection( input            = directionalityEx, 
#'                                    fromAtoBColor    = "darkgreen", 
#'                                    fromBtoAColor    = "orange", 
#'                                    noDirectionColor = "grey")
#' htmpDirection
#' @export heatmapDirection


heatmapDirection <- function( input , interactive = FALSE, fromAtoBColor = "darkgreen", fromBtoAColor = "orange", noDirectionColor = "grey", verbose = FALSE ) { 
    
    message("Checking the input object")
    checkClass <- class(input)[1]
    
    if(checkClass != "data.frame" ){
        message("Check the input object. Remember that this
                    object must be obtained after applying the directionality
                    function. The input object class must
                    be:\"data.frame\"")
        stop()
    }
    
    checkColnames <- colnames(input)

    if(checkColnames[1] != "disAcode" | checkColnames[2] != "disBcode"|
       checkColnames[3] != "AtoB" | checkColnames[4] != "BtoA" |
       checkColnames[5] != "test" | checkColnames[6] != "correctedPvalue" |
       checkColnames[7] != "result" ){
        message("Check the input object. Remember that this
                    object must be obtained after applying the directionality
                    function. The input data.frame must contain the next 6 columns:
                \"disAcode\", \"disBcode\", \"AtoB\", \"BtoA\", \"test\" and \"result\" ")
        stop()
    }
    
    if( interactive == FALSE ){
    p <- ggplot2::ggplot(input , ggplot2::aes ( disAcode, disBcode ), colour = factor(result) ) +
        ggplot2::geom_tile(ggplot2::aes(fill = as.factor(result)), colour = "white") +
        ggplot2::scale_fill_manual(values = c(fromAtoBColor, fromBtoAColor, noDirectionColor), name="Direction", labels=c(paste0(fromAtoBColor,"=From A to B"), paste0(fromBtoAColor,"= From B to A"), paste0(noDirectionColor,"= No directionality")))+
        ggplot2::theme_grey(base_size = 9) +
        ggplot2::labs ( title = "Comorbidity directionality", x = "diagnosis code under study (A)", y = "disease comorbidities (B)") +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::theme( plot.margin = grid::unit ( x = c ( 5, 15, 5, 15 ), units = "mm" ),
                        axis.line = ggplot2::element_line ( size = 0.7, color = "black" ), text = ggplot2::element_text ( size = 11 ) ,
                        axis.text.x = ggplot2::element_text ( angle = 45, size = 11, hjust = 1 ), panel.background = ggplot2::element_blank() )
    p
    }else if(interactive == TRUE){
        tt <- reshape2::acast(input, disBcode~disAcode, value.var="test")
        d3heatmap::d3heatmap(tt, dendrogram = "none", scale = "none",
                             colors =c(fromAtoBColor,fromBtoAColor, noDirectionColor), 
                             yaxis_font_size = 10,
                             xaxis_font_size = 10
        )
    }
    
}
