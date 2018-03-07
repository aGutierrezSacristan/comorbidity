#' A graphical summary of the disease prevalence. 
#'
#' Given an object of class \code{comorbidity}, and the characters used to 
#' specify the gender, a barplot showing the disease prevalence in general 
#' and according to gender is obtained. 
#'
#' @param input Object of \code{comorbidity} class. 
#' @param maleCode Characters(s) used to determine the male condition of a patient. 
#' Depending on the database it can be determined, for example, as \code{Male}, .
#' \code{MALE}, \code{M}, with digits as \code{0} or \code{1}. 
#' @param femaleCode Characters(s) used to determine the female condition of a patient. 
#' Depending on the database it can be determined, for example, as \code{Female}, .
#' \code{FEMALE}, \code{F}, with digits as \code{0} or \code{1}. 
#' @param databasePth Determines the path where the three required input files 
#' (patientData, diagnosisData, admissionData) are located. 
#' @param barColor Determines the bar color. By default \code{lightblue}. 
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return A barplot with disease prevalence in all the population and according to 
#' the gender.
#' @examples
#' load(system.file("extdata", "comorbidity.RData", package="comoRbidity"))
#' diseasePrevalence( input      = comor_obj, 
#'            maleCode   = "Male", 
#'            femaleCode = "Female",
#'            databasePth      = system.file("extdata", package="comoRbidity"), 
#'            barColor    = "darkcyan"
#'            )
#' @export diseasePrevalence



diseasePrevalence <- function( input, maleCode = "Male", femaleCode ="Female", databasePth, barColor = "lightblue", verbose = FALSE, warnings = TRUE) {
    
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
        message( "Creating the disease prevalence barplot" )
    }
    
    tt <- input@qresult
    
    maleDisease<- tt[ tt$patient_gender == maleCode, ]
    femaleDisease<- tt[ tt$patient_gender == femaleCode, ]

    #load comorbidity object with the total data for the estimations
    all <- load(paste0(databasePth, "/allData.RData"))
    allPopulation <- allData@qresult

    maleAll<- allPopulation[ allPopulation$patient_gender == maleCode, ]
    femaleAll<- allPopulation[ allPopulation$patient_gender == femaleCode, ]
    
    disPrev <- as.data.frame(matrix(ncol=4, nrow=3))
    colnames(disPrev) <- c("Group", "All population", "Disease population", "Prevalence")
    disPrev[,1] <- c("All population", "Male", "Female")
    disPrev[,2] <- c(length(unique(allPopulation$patient_id)), length(unique(maleAll$patient_id)), length(unique(femaleAll$patient_id)))
    disPrev[,3] <- c(length(unique(tt$patient_id)), length(unique(maleDisease$patient_id)), length(unique(femaleDisease$patient_id)))
    disPrev[,4] <- (disPrev[,3]/disPrev[,2])*100

    ### disease prevalence
    p <- ggplot2::ggplot(disPrev, ggplot2::aes(x = Group, y = Prevalence)) +
    ggplot2::geom_bar(stat="identity", position="dodge", color="black", fill= barColor) +
    ggplot2::geom_text(ggplot2::aes(x=Group, y=Prevalence, label=paste0(round(Prevalence,2), "% "),
                                    hjust=ifelse(sign(Prevalence)>0, 1, 0)),
                       position = ggplot2::position_dodge(width=1)) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs ( title = "Prevalence of the index disease" , y = "prevalence (%)", x = "population group")

    return(p)
   
}

 