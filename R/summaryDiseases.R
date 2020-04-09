#' A graphical summary of the diseases of interest. 
#'
#' Given an object of class \code{molecularComorbidity}, a graphical summary of 
#' the main characteristic of the disorders is obtained. 
#'
#' @param input Object of \code{molecularComorbidity} class.
#' @param database The default value is set to \code{'CURATED'}. User can select 
#' any of the databases available in DisGeNET.
#' @param type Determines if the analysis is performed for genes or diseases. By 
#' default the argument value is set to \code{dis_barplot}, allowing to analyze 
#' the number of genes associated to each one of the index diseases as well as
#' the number of disorders that share some genes with them. It can be changed to 
#' \code{gene_barplot} to analyze and characterize the genes associated to the 
#' index diseases.
#' @param interactive Determines if the output barplot is interactive or not. 
#' By default the \code{interactive} argument is set up as \code{FALSE}. The value 
#' of the argument can be changed to \code{TRUE}, as a result an interactive 
#' barplot generated with Shiny will be obtained. This argument is only available 
#' when \code{gene_barplot} type is selected.
#' @param assocGeneColor Determines the bar color that represents the number of associated genes. 
#' By default it is set to "#E69F00".  
#' @param assocDiseaseColor Determines the bar color that represents the number of associated diseases. 
#' By default it is set to "#136593".
#' @param dsiColor Determines the dot color that represents the Disease Specificity Index (DSI). 
#' By default it is set to "red".
#' @param dpiColor Determines the dot color that represents the Disease Pleiotropy Index (DPI). 
#' By default it is set to "black".
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return A barplot characterizing the genes or the disease in terms of gene-disease
#' association. 
#' @examples
#' load(system.file("extdata", "mc.RData", package="comoRbidity"))
#' summaryDiseases( input    = mc, 
#'                  type     = "dis_barplot",
#'                  database = "CURATED" 
#'            )
#' @export summaryDiseases

summaryDiseases <- function( input, database = "CURATED", type = "dis_barplot", assocGeneColor = "#E69F00", assocDiseaseColor = "#136593", dsiColor = "red", dpiColor = "black", interactive = FALSE, verbose = FALSE, warnings = TRUE) {
    

    #check if the input object is of class molecularComorbidity
    message("Checking the input object")
    checkClass <- class(input)[1]
    
    if(checkClass != "molecularComorbidity"){
        message("Check the input object. Remember that this
                    object must be obtained after applying the queryMolecular
                    function to your input file. The input object class must
                    be:\"molecularComorbidity\"")
        stop()
    }else{
        if( input@userInput == TRUE){
            message("Sorry, this function is only available for object 
                    obtained after applying the queryMolecular
                    function.")
            stop()  
        }
        
    }
    
    input <- input@qresult
    diseases <- unique( as.character( input$diseaseId ) )
    
    if( type == "dis_barplot"){
        
        disGdata <- disgenetDiseaseData( db = database )
        disGdata <- disGdata[disGdata$diseaseid %in% diseases, ]
        disGdata <- disGdata[,c( "diseaseid", "num_genes", "num_diseases" ) ]
        colnames( disGdata ) <- c( "disease", "assocGenes", "assocDiseases" )
        disGdata$assocDiseases <- gsub( "null", 0, disGdata$assocDiseases)
        disGdata$disease <- gsub( "umls:", "", disGdata$disease)      
        
        disGdata.m <- reshape2::melt(disGdata, id.vars='disease')
        disGdata.m$value <- as.numeric( as.character(disGdata.m$value ))
        
        p <- ggplot2::ggplot(disGdata.m, ggplot2::aes(disease, value)) +   
            ggplot2::geom_bar(ggplot2::aes(fill = variable), 
                              position = "dodge", 
                              stat="identity", 
                              colour = "black")
        
        p <- p + ggplot2::scale_fill_manual(values=c(assocGeneColor, assocDiseaseColor))
        p <- p + ggplot2::theme_classic( ) + ggplot2::theme( plot.margin = ggplot2::unit ( x = c ( 5, 15, 5, 15 ), units = "mm" ), 
                                                             axis.line = ggplot2::element_line ( size = 0.7, color = "black" ), text = ggplot2::element_text ( size = 14 ) ,
                                                             axis.text.x = ggplot2::element_text ( angle = 45, size = 10, hjust = 1 ))
        return( p )  
        
    }
    
    if( type == "gene_barplot"){
        input$pair <- paste( input$geneId, input$diseaseId, sep="-" )
        gdaPairs   <- unique( as.character( input$pair ) )
        
        disGDAdata <- disgenetgDaData( db = database )
        disGDAdata$pair <- paste( disGDAdata$geneid, disGDAdata$diseaseid, sep="-" )
        disGDAdata <- disGDAdata[ disGDAdata$pair %in% gdaPairs, ]
        
        result <- data.frame(matrix(0, ncol = 4, nrow = length(unique(disGDAdata$gene_symbol))))
        colnames(result) <- c("genes", "diseases","dsi", "dpi")
        result$genes <- as.character( unique(disGDAdata$gene_symbol) )
        
        for( i in 1:nrow(result)){
            selection <- disGDAdata[ disGDAdata$gene_symbol == result$genes[i], ]
            result[i,2] <- length( unique( selection$diseaseid ))
            result[i,3] <- selection$gene_dpi[1]
            result[i,4] <- selection$gene_dsi[1]
        }
        if( interactive == FALSE ){
            
            if(nrow(result) > 100){         
                data <- result[with(result, order(-result[,2])), ]
                data <- data[1:10, ]
                message("The top 10 genes with more diseases associated are shown")
            }else{
                data <- result
            }
            p <- ggplot2::ggplot(data, ggplot2::aes(x = genes,  y = diseases))+
                ggplot2::geom_bar(stat = "identity", position = "dodge", fill="#136593",
                                  ggplot2::aes(fill = diseases)) 
            
            p <- p + ggplot2::theme_classic( ) + ggplot2::theme( plot.margin = ggplot2::unit ( x = c ( 5, 15, 5, 15 ), units = "mm" ), 
                                                                 axis.line = ggplot2::element_line ( size = 0.7, color = "black" ), text = ggplot2::element_text ( size = 14 ) ,
                                                                 axis.text.x = ggplot2::element_text ( angle = 45, size = 10, hjust = 1 ))
            p <- p+ ggplot2::geom_point(ggplot2::aes(y=as.numeric(dsi)),
                               stat="identity",
                               alpha=.8,
                               size=3, 
                               color = dsiColor)
            p <- p+ ggplot2::geom_point(ggplot2::aes(y=as.numeric(dpi)),
                               stat="identity",
                               alpha=.8,
                               size=3, 
                               color = dpiColor)

            return( p )
            
            
        }
        
        if (interactive == TRUE) {
            
            # Define UI
            ui <- shiny::fluidPage(
                
                # Application title
                shiny::titlePanel("Gene Attributes"),
                
                # Show a plot of the generated distribution
                shiny::mainPanel(
                    shiny::plotOutput("distPlot")
                ), 
                
                shiny::sidebarLayout(
                    
                    # Sidebar with a slider input
                    shiny::sidebarPanel(
                        shiny::sliderInput("obs",
                                    "Number of observations:",
                                    min = 1,
                                    max = 500,
                                    value = 10)
                    ), 
                    
                    shiny::withTags(div(class='row-fluid',
                                 div(class='span3', checkboxInput(inputId = "dsi", label = "Disease Similarity Index (DSI)",value=FALSE)),
                                 div(class='span3', checkboxInput(inputId = "dpi", label = "Disease Pleiotropy Index (DPI)",value=FALSE))
                    ))
                )
            )
            
            # Server logic
            server <- function(input, output) {
                output$distPlot <- renderPlot({
                    
                    data <- result[1:input$obs,]
                    data <- data[with(data, order(-data[,2])), ]
                    colnames(data) <- c("genes", "diseases", "dsi", "dpi")
                    
                    p <- ggplot2::ggplot(data, ggplot2::aes(x = genes,  y = diseases))+
                        ggplot2::geom_bar(stat = "identity", position = "dodge", fill="#136593",
                                          ggplot2::aes(fill = diseases)) 
                    
                    p <- p + ggplot2::theme_classic( ) + ggplot2::theme( plot.margin = ggplot2::unit ( x = c ( 5, 15, 5, 15 ), units = "mm" ), 
                                                                         axis.line = ggplot2::element_line ( size = 0.7, color = "black" ), text = ggplot2::element_text ( size = 14 ) ,
                                                                         axis.text.x = ggplot2::element_text ( angle = 45, size = 10, hjust = 1 ))
                    
                    if (input$dsi) {
                        
                        p <- p+ ggplot2::geom_point(ggplot2::aes(y=dsi),
                                           stat="identity",
                                           #position_dodgen="dodge",
                                           alpha=.8,
                                           size=3, 
                                           color = "red")
                        p
                    }
                    if(input$dpi) {
                        
                        p <- p+ ggplot2::geom_point(ggplot2::aes(y=dpi),
                                           stat="identity",
                                           #position="dodge",
                                           alpha=.8,
                                           size=3, 
                                           color = "black")
                        p
                    }else{
                        p
                    }
                })
            }
            
            # Complete app with UI and server components
            shiny::shinyApp(ui, server)
        }
        
    }
 
}

disgenetDiseaseData <- function(db){
    
    url <- paste0( "https://www.disgenet.org/api/disease/source/", db, "?format=tsv" )
    dataTsv <- RCurl::getURLContent( url )
    myTextConnection <- textConnection( dataTsv )
    data <- read.csv( myTextConnection, header = TRUE, sep = "\t", colClasses=c("character"))
    close(myTextConnection  )
    return( data )
}

disgenetgDaData <- function( db ){
    
    url <- paste0( "https://www.disgenet.org/api/gda/source/", db, "?format=tsv" )
    dataTsv <- RCurl::getURLContent( url )
    myTextConnection <- textConnection( dataTsv )
    data <- read.csv( myTextConnection, header = TRUE, sep = "\t", colClasses=c("character"))
    close(myTextConnection  )
    return( data )
    
}
