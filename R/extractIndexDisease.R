#' Obtain the index diseases from a \code{comoRbidity} object
#'
#' Given an object of type \code{comoRbidity}, extract the index diseases 
#' present in the data.
#'
#' @param input An object of class \code{comoRbidity}
#' @return A list with the index diseases present in the input data
#' @examples
#' indexDis <- extractIndexDisease( obj )
#' @export extractIndexDisease
#' 
extractIndexDisease <- function( input ){
    if(class(input)== "comorbidity"){
        diseaseList <- input@indexDisList
        return(diseaseList)
    }
}
    
