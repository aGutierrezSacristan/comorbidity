#' A graphical summary of the code uses. 
#'
#' Given an object of class \code{comorbidity}, and the file which contains the
#' codes.
#'
#' @param codesPth Determines the path where the file with the index diseases is 
#' located (indexDiseaseCode) 
#' @param input A comorbidity object, obtained applying query function.
#' @param cutOff By default \code{0}. Change it to show only those codes with a 
#' percentage of use equal or higher than the cutOff value. 
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param barColor Determines the bar color. By default \code{darkblue}. 
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return A barplot with the diagnosis code uses. 
#' @examples
#' load(system.file("extdata", "comorbidity.RData", package="comoRbidity"))
#' diagnosticUse( codesPth    = system.file("extdata", package="comoRbidity"), 
#'                input       = comor_obj, 
#'                cutOff      = 0.5, 
#'                barColor    = "darkcyan"
#'            )
#' @export diagnosticUse

diagnosticUse <- function( codesPth, input, barColor = "darkblue",  cutOff=0, interactive = FALSE ){
    
    message("Checking the input object")
    checkClass <- class(input)[1]
    
    if(checkClass != "comorbidity"){
        message("Check the input object. Remember that this
                    object must be obtained after applying the query
                    function to your input file. The input object class must
                    be:\"comorbidity\"")
        stop()
    }

    
    codes <- read.delim(paste0(codesPth, "/indexDiseaseCode.txt"), 
                        header = TRUE, 
                        sep = "\t", 
                        colClasses = "character" )
    inputObj <- input@qresult
    
    if(input@aggregated == FALSE){
        codesUnderStudy <- as.character(unique(codes$Code))
        codesNumber <- length(unique(codesUnderStudy))
            
    }else if(input@aggregated == TRUE){
            codesUnderStudy <- as.character(unique(codes$Agg))
            codesNumber <- length(unique(codesUnderStudy))
    }
    
    if( as.numeric(codesNumber) < 2 ) {
        stop( "This analysis is only performed when more than 1 code for 
                      disease is available" )
    }
    
    tt <- inputObj[ inputObj$diagnosis_code %in% codesUnderStudy, ]
    rr <- as.data.frame( table( tt$diagnosis_code) )
    rr$perc <- rr$Freq/nrow( tt )*100
    rr <- rr[ rr$perc >= cutOff, ]
        
    ord <- rr[order( -rr [ "perc" ] ), "Var1" ]
    rr$Var1 <- factor( rr$Var1, levels= as.factor( ord ) )
        
    p <- ggplot2::ggplot ( rr, ggplot2::aes ( x = Var1, y = perc ), order = ord )  + 
        ggplot2::geom_bar ( stat = "identity", fill = barColor ) +
        ggplot2::labs ( title = "Codes' usage" , x = "diagnosis codes", y = "percentage of times used")
    
    p <- p + ggplot2::theme_classic( ) + 
        ggplot2::theme( plot.margin = ggplot2::unit ( x = c ( 5, 15, 5, 15 ), units = "mm" ), 
                       axis.line = ggplot2::element_line ( size = 1, color = "black" ), text = ggplot2::element_text ( size = 12 ) ,
                       axis.text.x = ggplot2::element_text ( angle = 45, size = 12, hjust = 1 ))
    
    if( interactive == FALSE ){
        return ( p )
    }else if ( interactive == TRUE ){
        plotly::ggplotly(p)
    }           
    
}

