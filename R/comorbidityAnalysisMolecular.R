#' Comorbidity Analysis \code{comorbidityAnalysisMolecular}
#'
#' Given an object of type \code{molecularComorbidity}, a comorbidity analysis is perform, 
#' taking into account the number of genes shared between the disorders.
#' It generates a \code{molecularcAnalysis} object.
#'
#' @param input   A \code{molecularComorbidity} object, obtained with the 
#' \code{queryMolecular} function.
#' @param pValue Determines if the p-value is estimated or not. By default it is 
#' set to \code{'FALSE'}. The \code{pValue} argument can be set to \code{'TRUE'} 
#' in order to estimate the P-value associated to each Jaccard Index. 
#' @param nboot Determines the number of random times that the Jaccard Index is 
#' computed using random sets. By default it is set to \code{100}. The value 
#' of the argument can be changed to any other numeric variable. 
#' @param ncores By default \code{1}. To run parallel computations on machines 
#' with multiple cores or CPUs, the \code{ncores} argument can be changed.
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return An object of class \code{molecularcAnalysis}
#' @examples
#' load(system.file("extdata", "mc.RData", package="comoRbidity"))
#' ex1 <- comorbidityAnalysisMolecular( 
#'               input      = mc, 
#'               pValue     = FALSE,
#'               verbose    = FALSE
#'               )
#' @export comorbidityAnalysisMolecular


comorbidityAnalysisMolecular <- function ( input, pValue = FALSE, nboot = 100, ncores = 1, verbose = FALSE, warnings = TRUE ){

    #check if the input object is of class molecularComorbidity
    message("Checking the input object")
    checkClass <- class(input)[1]
    
    if(checkClass != "molecularComorbidity"){
        message("Check the input object. Remember that this
                    object must be obtained after applying the queryMolecular
                    function to your input file. The input object class must
                    be:\"molecularComorbidity\"")
        stop()
    }
    
    message("Estimating the overlap and jaccard")
    comorb <- input@qresult
    codes.f    <- unique(as.character(comorb$diseaseName))
    codes.pair <- as.data.frame(gtools::combinations(n = length(codes.f),
                                           r = 2,
                                           v = codes.f,
                                           repeats.allowed = FALSE))
    
    codes.pair$V1 <- as.character(codes.pair$V1)
    codes.pair$V2 <- as.character(codes.pair$V2)
    codes.pair$geneV1  <- NA
    codes.pair$geneV2  <- NA
    codes.pair$overlap <- NA
    codes.pair$jaccard <- NA

  for( cc in 1:nrow(codes.pair)){
    selection1 <- unique(comorb[comorb$diseaseName == as.character(codes.pair$V1[cc]), 1])
    selection2 <- unique(comorb[comorb$diseaseName == as.character(codes.pair$V2[cc]), 1])

    overlapV12 <- selection1[ selection1 %in% selection2]
    unionV1V2  <- unique(c(selection1, selection2))

    codes.pair$geneV1[cc]  <- length(selection1)
    codes.pair$geneV2[cc]  <- length(selection2)
    codes.pair$overlap[cc] <- length(overlapV12)
    codes.pair$jaccard[cc] <- round(length(overlapV12)/length(unionV1V2), 3)
  }
  
  codes.pair <- codes.pair[codes.pair$overlap != 0, ]
  
  if( pValue == TRUE){
      message("Estimating the p-value")
      codes.pair <- pValueEstimation(codes.pair, nboot=nboot)
      
      codes.pair$pval <- round( as.numeric( codes.pair$pval ), 3)
      
      
  }else{
      message("p-value will not be estimated. Set the pVal argument to \"TRUE\" for p-value estimation")
  }
  
  codes.pair$jaccard <- round( as.numeric( codes.pair$jaccard ), 3)
  
  resultsGC <- new( "molecularcAnalysis", 
                    ovlpMin   = min(codes.pair$overlap), 
                    ovlpMax    = max(codes.pair$overlap), 
                    jaccardMin= min(codes.pair$jaccard), 
                    jaccardMax= max(codes.pair$jaccard),
                    pValue    = pValue, 
                    tdiseases = length(unique(c(codes.pair$V1, codes.pair$V2))),
                    dispairs  = nrow( codes.pair ),
                    result    = codes.pair 
  )
  return( resultsGC )

}