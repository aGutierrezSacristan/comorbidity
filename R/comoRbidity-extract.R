#' Obtain the data from a \code{comorbidity}, \code{molecularComorbidity}, 
#' \code{cAnalysis} or \code{molecularcAnalysis} object.
#'
#' @name extract
#' @rdname extract-methods
#' @param object Object of class \code{comorbidity}, \code{molecularComorbidity}, 
#' \code{cAnalysis} or \code{molecularcAnalysis} object.
#' @return A \code{data.frame} containing the raw result
#' @examples
#' \dontrun{
#' #Being x an comoRbidity
#' qr <- extract(x) 
#' }
#' @export
#' 
setMethod( "extract",
   signature = "comorbidity",
   definition = function( object ) {
     return( object@qresult )
   }
)

setMethod( "extract",
           signature = "cAnalysis",
           definition = function( object ) {
               return( object@result )
           }
)

setMethod( "extract",
           signature = "molecularcAnalysis",
           definition = function( object ) {
               return( object@result )
           }
)

setMethod( "extract",
           signature = "molecularComorbidity",
           definition = function( object ) {
               return( object@qresult )
           }
)
