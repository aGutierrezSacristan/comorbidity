#' Age distribution comparison: patients subset vs all population. 
#'
#' Given an object of class \code{comoRbidity} the age distribution of your patients
#' of interest (those that have at least one code of the ones contained in indexDiseaseCodes.txt 
#' file) is represent comparing with the age distribution of the your whole patients 
#' database.
#'
#' @param input A comorbidity object, obtained with the query function. 
#' @param codesPth Determines the path where the file with the index diseases is 
#' located (indexDiseaseCode) 
#' @param databasePth Determines the path where the three required input files 
#' (patientData, diagnosisData, admissionData) are located.
#' @param type by Befault \code{'together'}. This argument allows the user to 
#' select the output barplot type. The other possible option is \code{'separate'}
#' @param allColor Determines the bar color representing the whole population. 
#' By default \code{darkblue}. 
#' @param disorderColor Determines the bar color representing the disorder patients. 
#' By default \code{gold}.
#' @param interactive Determines if the output barplot is interactive or not. 
#' By default the \code{interactive} argument is set up as \code{FALSE}. The value 
#' of the argument can be changed to \code{TRUE}, as a result an interactive 
#' barplot will be obtained. 
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return A barplot comparing the ages distributions
#' @examples
#' load(system.file("extdata", "comorbidity.RData", package="comoRbidity"))
#' popAge <- populationAge( input       = comor_obj, 
#'                          codesPth    = system.file("extdata", package="comoRbidity"), 
#'                          databasePth = system.file("extdata", package="comoRbidity"),
#'                          type        = "separate"
#'               )
#' @export populationAge



populationAge <-function( input, codesPth, databasePth, type="together", disorderColor = "gold", allColor = "darkblue",interactive = FALSE, verbose = FALSE, warnings = TRUE ){
    
    
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
    
    if(input@aggregated == TRUE){
        codesUnderStudy    <- as.character(unique(codes$Agg))
    }else{
        codesUnderStudy    <- as.character(unique(codes$Code))
    }
    
    
    input <- input@qresult
    input <- input[ input$diagnosis_code %in% codesUnderStudy, ]
    input <- input[ with ( input, order(admissionStartDate, patient_id)), ]
    input <- input[!duplicated(input$patient_id), ]
    
    all <- load(paste0(databasePth, "/allData.RData"))
    all <- allData@qresult
    all <- all[ with ( all, order(admissionStartDate, patient_id)), ]
    all <- all[!duplicated(all$patient_id), ]
    
    
    countAll<-as.data.frame(table(all$age))
    colnames(countAll) <- c("Ages", "Total")
    countAll$perc <- ((countAll$Total)/sum(countAll$Total))*100
    countAll$class <- "all"
    
    countDisorder<-as.data.frame(table(input$age))
    colnames(countDisorder) <- c("Ages", "Total")
    countDisorder$perc <- ((countDisorder$Total)/sum(countDisorder$Total))*100
    countDisorder$class <- "disorder"
    
    
    final <- rbind(countAll, countDisorder)
    final$Ages <- as.numeric(final$Ages)
    
    
    if(type=="together"){
        if(interactive == FALSE){
            p <- ggplot2::ggplot(data = final, 
                                 ggplot2::aes(x = as.numeric(Ages), y = perc)) +
                ggplot2::scale_fill_manual( values = c( allColor, disorderColor )) +
                ggplot2::geom_bar(stat = "identity", position = "dodge",
                                  ggplot2::aes(fill = class)) 
            return( p )
        }
        if(interactive == TRUE){
            p <- ggplot2::ggplot(data = final, 
                                 ggplot2::aes(x = as.numeric(Ages), y = perc)) +
                ggplot2::scale_fill_manual( values = c( allColor, disorderColor )) +
                ggplot2::geom_bar(stat = "identity", position = "dodge",  
                                  ggplot2::aes(fill = class)) 
            
            plotly::ggplotly(p)
        }
        
        
        
    }
    else if(type=="separate"){
        if(interactive == FALSE){
            p <- ggplot2::ggplot(final, 
                                 ggplot2::aes(x=as.numeric(Ages), y=perc, fill= class)) +
                ggplot2::geom_bar(stat = "identity", position = "dodge") +
                ggplot2::scale_fill_manual( values = c( allColor, disorderColor )) +
                ggplot2::facet_grid(class ~ .)
            return( p )
        }
        if(interactive == TRUE){
            p <- ggplot2::ggplot(final, 
                                 ggplot2::aes(x=as.numeric(Ages), y=perc, fill= class)) +
                ggplot2::geom_bar(stat = "identity", position = "dodge") +
                ggplot2::scale_fill_manual( values = c( allColor, disorderColor )) +
                ggplot2::facet_grid(class ~ .)
            plotly::ggplotly(p)
        }
    }
}
