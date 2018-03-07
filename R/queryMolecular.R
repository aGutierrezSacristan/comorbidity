#' Query DisGeNET for given disease(s) and generates an \code{molecularComorbidity}
#' object for comorbidity analysis
#'
#' Given a file with diseases, a comorbidity analysis is perform according to
#' the genes shared between the diseases. The comorbidity analysis can be performed
#' at cui level or in a higher category specified in the input file. It
#' generates a \code{molecularComorbidity} object.
#'
#' @param filePth The file name and with the complete path where the file with
#' disorders of interest is located
#' @param unify By default it is set to \code{FALSE}. If the argument is set to 
#' \code{TRUE}, the \code{name} colum from the cui disease file will be selected 
#' for doing the comorbidity analysis. 
#' @param database Name of the database that will be queried. It can take the values \code{'CTD_human'} to use Comparative
#' Toxicogenomics Database, human data; \code{'UNIPROT'} to use Universal
#' Protein Resource;\code{'CLINVAR'} to use ClinVar, a public archive of relationships
#' among sequence variation and human phenotype; \code{'GWASCAT'} to use
#' the NHGRI-EBI GWAS Catalog; \code{'ORPHANET'}, to use
#' Orphanet, the portal for rare diseases and orphan drugs;
#' \code{'CURATED'} to use expert curated, human databases;
#' \code{'RGD'}, to use Rat Genome Database; \code{'MGD'}, to use the Mouse Genome Database;
#' \code{'CTD_rat'} to use Comparative Toxicogenomics Database, rat data;
#' \code{'CTD_mouse'} to use Comparative Toxicogenomics Database, mouse data;
#' \code{'PREDICTED'} to use the expert curated, animal models data;
#' \code{'ALL'} to use all these databases. Default \code{'CURATED'}.
#' @param score A vector with two elements: 1) character with greather
#' \code{'>'} or with lower \code{'<'} meaing greather or equal and lower or
#' equal; 2) the score to be compared. By default it is set to \code{c('>', 0)}. 
#' It means that all the data available in DisGeNET will be used for the comorbidity analysis. 
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return An object of class \code{molecularComorbidity}
#' @examples
#' ex1 <- querymolecular(
#'               filePth = system.file("extdata", "cuiDiseaseList.txt", package="comoRbidity"),
#'               unify    = TRUE,
#'               verbose  = FALSE
#'               )
#' @export querymolecular

querymolecular <- function ( filePth,  unify = FALSE, database = "CURATED", score = c(">", 0), verbose = FALSE, warnings = TRUE ){

  #load the input file with the disorders
  #filePth <- "/home/alba/Desktop/exampleInput.txt"
  codes <- read.delim ( filePth,
                        header=TRUE,
                        sep="\t",
                        colClasses = "character")
  message("Checking the input file structure \n")

  colnamesRequired <- c("identifier", "name")

  checkA <- colnamesRequired[colnamesRequired %in% colnames(codes)]
  if(length(checkA) != length(colnamesRequired)){
    message("Check the input file structure. Remember that this
                    file must contain at least two columns with the column
                    names as follows:\n -> identifier \n -> name\n")
    stop()
  }
  
 
  
  if("disgenet2r" %in% rownames(installed.packages()) == FALSE) {
      message("disgenet2r package needed for molecular comorbidity analysis")
      message("For installing the last version type: 
                library(devtools)
                install_bitbucket(\"ibi_group/disgenet2r\", force = TRUE)
                library(disgenet2r)")
      stop()
  }else{
      library(disgenet2r)    
      }

  message("Extracting the disorders of interes for the comorbidity analysis\n")

  if( unify == TRUE){
    codesUnderStudy    <- as.character(unique(codes$name))
    message("unify = TRUE \nComorbidity analysis will be performed at the name level\n")

    for( i in 1:length(codesUnderStudy)){
      message(paste0("Extracting genes associated to ", codesUnderStudy[i], "\n"))
      diseases <- unique(as.character(codes[codes$name == codesUnderStudy[i],1]))
      results  <- disgenet2r::disgenetDisease( disease  = diseases,
                                   database = database,
                                   score    = score)
      if( i == 1){
        comorb <- results@qresult
        comorb$diseaseCategory <- codesUnderStudy[i]
        comorb <- comorb[,c("c2.geneId", "c2.symbol", "c0.diseaseId", "c1.name")]
        comorb <- comorb[!duplicated( comorb ), ]
      }else if( nrow(results@qresult)!=0){
            comorbResult <- results@qresult
            comorbResult$diseaseCategory <- codesUnderStudy[i]
            comorbResult <- comorbResult[,c("c2.geneId", "c2.symbol", "c0.diseaseId", "c1.name")]
            comorbResult <- comorbResult[!duplicated( comorbResult ), ]
            
            comorb <- rbind( comorb, comorbResult )
        }

      }

    }else if( unify == FALSE){
    codesUnderStudy    <- as.character(unique(codes$identifier))
    message("unify = FALSE \nComorbidity analysis will be performed using the disease list\n")
    diseases <- unique(as.character(codes$identifier))
    
    results  <- disgenet2r::disgenetDisease( disease  = diseases,
                                 database = database,
                                 score    = score, 
                                 verbose  = FALSE, 
                                 warnings  = FALSE)
    comorb <- results@qresult
    comorb <- comorb[,c("c2.geneId", "c2.symbol", "c0.diseaseId", "c1.name")]
    
    
   
  }

  colnames(comorb) <- c("geneId", "geneSymbol", "diseaseId", "diseaseName")
  
  result <- new( "molecularComorbidity", 
                 search       = ifelse( length(codesUnderStudy) > 1, "list", "single" ), 
                 aggregated   = unify, 
                 indexDis     = length(unique( codesUnderStudy )),                   
                 nfDisease    = length(unique(comorb$diseaseName)), 
                 nGenes       = length( unique ( comorb$geneId ) ), 
                 indexDisList = unique( codesUnderStudy ),
                 qresult      = comorb  
  )
  
  return(result)

}
