#' From a file containing a table with gene and disease information generates a \code{molecularComorbidity}
#' object for comorbidity analysis
#'
#' Given a file with gene and diseases, the \code{dataframe2molecularcomorbidity}
#' generates a \code{molecularComorbidity} object.
#'
#' @param filePth The file name and with the complete path where the file with
#' the genes and disorders of interest is located
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return An object of class \code{molecularComorbidity}
#' @examples
#' ex1 <- dataframe2molecularcomorbidity(
#'               filePth = system.file("extdata", "genediseaseTable.txt", package="comoRbidity"),
#'               verbose  = FALSE
#'               )
#' @export dataframe2molecularcomorbidity

dataframe2molecularcomorbidity <- function ( filePth, verbose = FALSE, warnings = TRUE ){
    
    #load the input file with the gene and diseases
    genedisease <- read.delim ( filePth,
                          header=TRUE,
                          sep="\t",
                          colClasses = "character")
    message("Checking the gene-disease input file structure \n")
    
    colnamesRequired <- c("gene", "diseaseName")
    
    checkA <- colnamesRequired[colnamesRequired %in% colnames(genedisease)]
    
    if(length(checkA) != length(colnamesRequired) ){
        message("Check the input file structure. Remember that this
                file must contain at least two columns with the column
                    names as follows:\n -> gene \n -> diseaseName\n")
        stop()
    }
    

    message("Creating the molecularComorbidity object \n")
    
    
    result <- new( "molecularComorbidity", 
                   search       = ifelse( length(genedisease$diseaseName) > 1, "list", "single" ), 
                   aggregated   = FALSE, 
                   indexDis     = length(unique( genedisease$diseaseName )),                   
                   nfDisease    = length(unique( genedisease$diseaseName )), 
                   nGenes       = length( unique ( genedisease$gene) ), 
                   indexDisList = unique( genedisease$diseaseName ),
                   qresult      = genedisease  
    )
    
    return(result)
    
}
