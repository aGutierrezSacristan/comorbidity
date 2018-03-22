#' A graphical summary of the database content. 
#'
#' Given an object of class \code{comorbidity}, and the characters used to 
#' specify the sex, a figure containing 3 plots, one of the age distribution, 
#' another one of the age distribution and the third one with the relation 
#' between age and sex distribution. 
#'
#' @param input Object of \code{comorbidity} class. 
#' @param maleCode Characters(s) used to determine the male condition of a patient. 
#' Depending on the database it can be determined, for example, as \code{Male}, .
#' \code{MALE}, \code{M}, with digits as \code{0} or \code{1}. 
#' @param femaleCode Characters(s) used to determine the female condition of a patient. 
#' Depending on the database it can be determined, for example, as \code{Female}, .
#' \code{FEMALE}, \code{F}, with digits as \code{0} or \code{1}. 
#' @param barColor Determines the bar color. By default \code{darkblue}. 
#' @param femaleColor Determines the color for female values. By default \code{gold}. 
#' @param maleColor Determines the color for male values. By default \code{lightblue}. 
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return A multiple graph containing a barplot with age distribution,  
#' a boxplot representing age distribution by sex and a pie chart representing 
#' sex distribution.
#' @examples
#' load(system.file("extdata", "comorbidity.RData", package="comoRbidity"))
#' summaryDB( input      = comor_obj, 
#'            maleCode   = "Male", 
#'            femaleCode = "Female"
#'            )
#' @export summaryDB



summaryDB <- function( input, maleCode = "Male", femaleCode ="Female", barColor = "darkblue", femaleColor = "gold", maleColor = "lightblue", verbose = FALSE, warnings = TRUE) {
    
    message("Checking the input object")
    checkClass <- class(input)[1]
    
    if(checkClass != "comorbidity"){
        message("Check the input object. Remember that this
                object must be obtained after applying the query
                function to your input file. The input object class must
                be:\"comorbidity\"")
        stop()
    }
    
    if( verbose == TRUE ){
        message( "Creating a summary table with the main characteristics of population" )
    }
    
    tt <- input@qresult
    tt <- tt[ !duplicated( tt$patient_id ), ]
    
    male<- tt[ tt$patient_sex == maleCode, ]
    female<- tt[ tt$patient_sex == femaleCode, ]
    tt$age    <- as.numeric( tt$age )
    
    ### age distribution
    p <-ggplot2::ggplot( tt ) +
        ggplot2::geom_bar( ggplot2::aes(age) , fill = barColor) +
        ggplot2::labs ( title = "Age distribution" , x = "age", y = "# of patients") +
        ggplot2::scale_x_continuous(breaks = seq(min(tt$age), max(tt$age), by= 3),labels = seq(min(tt$age), max(tt$age), by=3))
    
    
    p <- p + ggplot2::theme_classic( ) + ggplot2::theme( plot.margin = ggplot2::unit ( x = c ( 5, 15, 5, 15 ), units = "mm" ), 
                                                         axis.line = ggplot2::element_line ( size = 0.7, color = "black" ), text = ggplot2::element_text ( size = 14 ) ,
                                                         axis.text.x = ggplot2::element_text ( angle = 45, size = 9, hjust = 1 ))
    
    #Boxplot
    male1<-male[,c(1,6,4)]
    female1<-female[,c(1,6,4)]
    
    tot<-rbind(male1, female1)
    
    bp <- ggplot2::ggplot(tot, ggplot2::aes(x=patient_sex, y=age)) + 
        ggplot2::geom_boxplot( fill= c( femaleColor, maleColor)) +
        ggplot2::stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
        ggplot2::labs ( title = "Age distribution by Sex" , x = "sex", y = "age")
    
    bp <- bp + ggplot2::theme_classic( ) +
        ggplot2::theme( plot.margin = ggplot2::unit ( x = c ( 5, 15, 5, 15 ), units = "mm" ), 
                        axis.line = ggplot2::element_line ( size = 0.7, color = "black" ), 
                        text = ggplot2::element_text ( size = 11 ),
                        axis.text.x = ggplot2::element_text ( angle = 45, size = 11, hjust = 1 ))
    
    
    ### sex distribution
    slices<-c(length ( unique( female$patient_id ) ), length ( unique( male$patient_id ) ) )
    lbl1 <- paste0(round((length (unique( female$patient_id ) ) ) /length ( unique( tt$patient_id ) ) *100, 2 ), "%")
    lbl2 <- paste0(round((length (unique( male$patient_id ) ) ) /length ( unique( tt$patient_id ) ) *100, 2 ), "%")
    
    grp1 <- paste("Female (",lbl1, ")" )    
    grp2 <- paste("Male (",lbl2, ")" )
    
    
    df <- data.frame(
        group = c( grp1, grp2 ),
        value = c(slices[1], slices[2])
    )
    
    bp2<- ggplot2::ggplot(df, ggplot2::aes(x="", y=value, fill=group))+
        ggplot2::geom_bar(width = 1, stat = "identity")
    
    
    blank_theme <- ggplot2::theme_minimal()+
        ggplot2::theme(
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            panel.border = ggplot2::element_blank(),
            panel.grid   = ggplot2::element_blank(),
            axis.ticks   = ggplot2::element_blank(),
            plot.title   = ggplot2::element_text(size=14, face="bold")
        )
    
    pie <- bp2 + ggplot2::coord_polar("y", start=0)
    pie <- pie + ggplot2::scale_fill_manual(values=c(femaleColor, maleColor)) +  
        blank_theme + 
        ggplot2::theme(axis.text.x = ggplot2::element_blank())
    
    
    #plot them together
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
    vplayout <- function(x, y) grid::viewport(layout.pos.row = x, layout.pos.col = y)
    print(p, vp = vplayout(1, 1:2))
    print(bp, vp = vplayout(2, 1))
    print(pie, vp = vplayout(2, 2))
}
