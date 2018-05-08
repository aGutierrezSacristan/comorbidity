#' Sex Ratio Analysis. 
#'
#' Given two objects of class \code{comoRbidity}, one for each gender, the sex 
#' ratio estimation is performed.    
#'
#' @param female A \code{comoRbidity} object, containing the comorbidity data for
#' female patients. 
#' @param male A \code{comoRbidity} object, containing the comorbidity data for
#' male patients. 
#' @param fisherTest By default the fisher test is not performed. The fisherTest 
#' default argument is set to \code{0}. It can be changed to 1 in order to perform 
#' the fisher test to each common comorbidity present in male and female.  
#' @param fisherCutOff by default \code{0.05}. The value of the argument can be changed 
#' to any other numeric variable, according to the range of the fisher value.
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return A dataframe with the sex ratio estimation is obtained
#' @examples
#' load(system.file("extdata", "comorMale.RData", package="comoRbidity"))
#' load(system.file("extdata", "comorFemale.RData", package="comoRbidity"))
#' srEstimation <- sexRatioAdapted( female       = comorFemale, 
#'                                  male         = comorMale, 
#'                                 fisherTest    = 0
#'                                 )
#' @export sexRatioAdapted

sexRatioAdapted <- function (female, male, femalePatients, malePatients, fisherTest = 0, fisherCutOff=0.05, verbose = FALSE, warnings = TRUE){

    #check if the input object is of class molecularComorbidity
    message("Checking the input objects")
    checkClass1 <- class(female)[1]
    checkClass2 <- class(male)[1]
    
    if(checkClass1 != "cAnalysis" | checkClass2 != "cAnalysis" ){
        message("Check the input objects. Remember that these
                    objects must be obtained after applying the comorbidityAnalysis
                    function to your input file. The input object class must
                    be:\"cAnalysis\"")
        stop()
    }
    
    fdep   <- female@tpatients
    mdep   <- male@tpatients
    
    female <- female@result
    #female$pairs <- NA
    
    #for(i in 1:nrow(female)){
    #    pairDisF <- sort(c(female[i,"disAcode"],female[i,"disBcode"]))
    #    female$pairs[i] <- paste(pairDisF[1], pairDisF[2], sep="-")
    #}
    female$pairs <- paste(female$disAcode, female$disBcode, sep="-")
    female <- female[,c("disAcode", "disBcode", "pairs", "AB")]
    
    male   <- male@result
    #male$pairs <- NA
    
    #for(i in 1:nrow(male)){
    #    pairDisM <- sort(c(male[i,"disAcode"], male[i,"disBcode"]))
    #    male$pairs[i] <- paste(pairDisM[1], pairDisM[2], sep="-")
    #}
    male$pairs <- paste(male$disAcode, male$disBcode, sep="-")
    male <- male[,c("disAcode", "disBcode", "pairs", "AB")]
    

    maleFilter <- male[ male$pairs %in% female$pairs, ]
    femaleFilter <- female[ female$pairs %in% male$pairs, ]
    
    if(nrow(maleFilter)==0){
        stop("There are no pairs in common")
    }
    
    if(fisherTest==1){
    fisherResults <- as.data.frame(matrix(ncol=3))
    
    for(i in 1:nrow(femaleFilter)){
                ff <- as.numeric(femaleFilter[i,4])
                mm <- as.numeric(maleFilter[maleFilter$pairs==femaleFilter[i,3],4])
                
                group1AB <- ff
                group2AB <- mm
                group1AnotB <- fdep - ff
                group2AnotB <- mdep - mm
                
                mtrx <- matrix( c( group1AB, group1AnotB, group2AB, group2AnotB), nrow = 2 )
    
                tryCatch( {ff <- fisher.test( mtrx )}, error=function(msg) {
                    message(msg)
                    message("code1:", femaleFilter[i,1], " - code2:", femaleFilter[i,2])
                })
                
                fisherResults[i,] <- rbind( c(femaleFilter[i,1], femaleFilter[i,2], ff$p.value ))
        }
        
        fisherResults <- fisherResults[as.numeric(fisherResults$V3) < fisherCutOff,]
        fisherResults$pairs <- paste(fisherResults[,1],fisherResults[,2], sep="-")
            
        maleFilter <- maleFilter[ maleFilter$pairs %in% fisherResults$pairs, ]
        femaleFilter <- femaleFilter[ femaleFilter$pairs %in% fisherResults$pairs, ]  

        
    }
    
    sexRatio <- as.data.frame(matrix(ncol=3))
    
    for(i in 1:nrow(femaleFilter)){
        pairf <- as.numeric(femaleFilter[i,4])
        pairm <- as.numeric(maleFilter[maleFilter$pairs==femaleFilter[i,3],4])
        
        num <- 1+(mdep/pairm)
        den <- 1+(fdep/pairf)
        
        
        sr <- log(num/den)
        sexRatio[i,] <- rbind( c(femaleFilter[i,1], femaleFilter[i,2], round(sr,3)))
        
    }
    colnames(sexRatio) <- c("disA", "disB", "SR")
    
    sexRatio$pairs <- paste( sexRatio$disA, sexRatio$disB, sep="-")
    
    colnames(fisherResults)[3] <- "FisherTest"
    sexRatioFinal <- merge( sexRatio, fisherResults, by="pairs")
    sexRatioFinal <- sexRatioFinal[, c("disA", "disB", "SR", "FisherTest")]
    sexRatioFinal$FisherTest <- round( as.numeric(sexRatioFinal$FisherTest), 3)
    sexRatioFinal <- sexRatioFinal[with(sexRatioFinal, order(SR)),]

    return(sexRatioFinal)
}

