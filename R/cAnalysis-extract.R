#' Obtain the raw query from a \code{cAnalysis} object.
#'
#' @name extract
#' @rdname extract-methods
#' @param object of class \code{cAnalysis}, \code{molecularcAnalysis}, 
#' \code{molecularComorbidity} or \code{comorbidity}
#' @return A \code{data.frame} containing the raw result from the initial data 
#' (comorbidity)or comorbidity analysis (cAnalysis)
#' @examples
#' \dontrun{
#' #Being x an cAnalysis
#' qr <- extract(x) 
#' }
#' @export
#' 
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
           signature = "comorbidity",
           definition = function( object ) {
               return( object@qresult )
           }
)

setMethod( "extract",
           signature = "molecularComorbidity",
           definition = function( object ) {
               return( object@qresult )
           }
)
