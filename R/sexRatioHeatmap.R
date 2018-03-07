#' Sex Ratio Heatmap 
#'
#' Given the data frame obtained after apply the \code{sexRatio} function, a 
#' heatmap representing the sex ratio results is obtained. 
#'
#' @param input Data frame obtained after applying the \code{sexRatio} function. 
#' @param maleColor Determines the heatmap color for those comorbidities that
#' are more likely in men than women. By default \code{"blue"}. 
#' @param femaleColor Determines the heatmap color for those comorbidities that
#' are more likely in women than men. By default \code{"red"}. 
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return A dataframe with the sex ratio estimation is obtained
#' @examples
#' load(system.file("extdata", "srEstimation.RData", package="comoRbidity"))
#' heatmapSexRatio( input = srEstimation )
#' @export heatmapSexRatio

heatmapSexRatio <- function( input, maleColor = "blue", femaleColor = "red", verbose = FALSE, warnings = TRUE, interactive = FALSE ){
    
    input$SR <- as.numeric(input$SR)
    casted_rr <- reshape::cast(input, disA ~ disB, value= "SR"  )
    l <- dim(casted_rr)[2]        
    Rowm  <- rowMeans(casted_rr[2:l], na.rm = T)       
    casted_rr<- cbind(casted_rr, Rowm)
    ordered<-  casted_rr[order(casted_rr["Rowm"]), "disA" ]          
    input$disA  <- factor(input$disA  , levels= as.factor(ordered))
    
    if( interactive == FALSE ){
    p <- ggplot2::ggplot(input, ggplot2::aes(disA, disB)) + 
         ggplot2::geom_tile(ggplot2::aes(fill = SR), colour = "black") + 
        ggplot2::scale_fill_gradient2(midpoint=0, low= maleColor, high= femaleColor)
    
    p <- p + ggplot2::theme_grey(base_size = 9) + 
        ggplot2::labs(title = "Sex Ratio Comorbidity", x = "diagnosis under study",  y = "disease comorbidities") + 
        ggplot2::scale_x_discrete(expand = c(0, 0))
    
    p <- p + ggplot2::theme( plot.margin = ggplot2::unit(x=c(5,15,5,15),units="mm") , 
                             axis.line = ggplot2::element_line(size = 0.5, color = "black"), 
                             text = ggplot2::element_text(size = 9), 
                             axis.text.x = ggplot2::element_text(angle=45, size = 9, hjust = 1), 
                             panel.background = ggplot2::element_blank())
    p
    }else if(interactive == TRUE){
        tt <- reshape2::acast(input, as.character(disB) ~ as.character(disA), value.var="SR")
        d3heatmap::d3heatmap(tt, dendrogram = "none", scale = "none",
                             colors =c(maleColor, femaleColor), 
                             yaxis_font_size = 10,
                             xaxis_font_size = 10
        )
    }
    
}