#' Comorbidity Analysis \code{cAnalysis}
#'
#' Given an object of type \code{comoRbidity}, a comorbidity analysis is perform, 
#' for the subset of population under specific conditions of age and gender. It 
#' generates a \code{cAnalysis} object.
#'
#' @param input  A comorbidity object, obtained with the query function. 
#' @param databasePth Determines the path where the three required input files 
#' (patientData, diagnosisData, admissionData) are located. 
#' @param codesPth Determines the path where the file with the index diseases is 
#' located (indexDiseaseCode)
#' @param ageRange Determines what is the age range of interest for
#' performing the comorbidity analysis. By default it is set from 0 to 100 
#' years old. 
#' @param gender Determine what is the gender of interest for 
#' performing the comorbidity analysis. By default \code{ALL}. Change it to the 
#' gender of interest for your comorbidity analysis.
#' @param score The comorbidity score is a measure based on  the observed comorbidities
#' and the expected ones, based on the occurrence of each disease.
#' @param fdr A Fisher exact test for each pair of diseases is performed to assess 
#' the null hypothesis of independence between the two diseases. The Benjamini-Hochberg 
#' false discovery rate method (FDR) is applied to correct for multiple testing.
#' @param oddsRatio The odds ratio represents the increased chance that someone 
#' suffering disease X will have the comorbid disorder Y.
#' @param rr The relative risk refers to the fraction between the number of 
#' patients diagnosed with both diseases and random expectation based on disease 
#' prevalence.
#' @param phi The Pearsons correlation for binary variables (Phi) measures the 
#' robustness of the comorbidity association.
#' @param cores By default \code{1}. To run parallel computations on machines 
#' with multiple cores or CPUs, the cores argument can be changed. 
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return An object of class \code{cAnalysis}
#' @examples
#' load(system.file("extdata", "comorbidity.RData", package="comoRbidity"))
#' ex1 <- comorbidityAnalysis( 
#'               input              = comor_obj,
#'               databasePth      = system.file("extdata", package="comoRbidity"),
#'               codesPth         = system.file("extdata", package="comoRbidity"),
#'               ageRange         = c(0,50),
#'               gender           = "Female", 
#'               score            = 1, 
#'               fdr              = 1
#'               )
#' @export comorbidityAnalysis

comorbidityAnalysis <- function ( input, codesPth, databasePth, ageRange=c(0,100), gender="ALL", score, fdr, oddsRatio, rr, phi, cores = 1, verbose = FALSE, warnings = TRUE ){
    
    message("Checking the input object")
    checkClass <- class(input)[1]
    
    if(checkClass != "comorbidity"){
        message("Check the input object. Remember that this
                    object must be obtained after applying the query
                    function to your input file. The input object class must
                    be:\"comorbidity\"")
        stop()
    }
    
    message( "Staring the comorbidity analysis" )
    message( "Loading the indexDiseaseCode file" )
    #load the codes of interest
    codes <- read.delim ( file.path(codesPth, "indexDiseaseCode.txt"),
                          header=TRUE, 
                          sep="\t", 
                          colClasses="character" ) 
    
    if( input@aggregated == TRUE){
        message( "Aggregating the disorders of interest" )
        codesUnderStudy    <- as.character(unique(codes$Agg))
    }else if( input@aggregated == FALSE){
        codesUnderStudy    <- as.character(unique(codes$Code))
    }
    
    
    #load comorbidity object with the total data for the estimations
    all <- load(paste0(databasePth, "/allData.RData"))
    all <- allData@qresult

    input <- input@qresult
    input$age <- as.numeric(input$age)
    all$age <- as.numeric(all$age)

    
    if ( !missing( ageRange ) ) {
        input <- input[ input$age >= ageRange[ 1 ] & input$age <= ageRange[ 2 ], ]
        all <- all[ all$age >= ageRange[ 1 ] & all$age <= ageRange[ 2 ], ]
        
    }
    
    if ( !missing( gender ) ) {
        if(gender!="ALL"){
            input <- input[ input$patient_gender == gender, ]
            all <- all[ all$patient_gender == gender, ]   
        }
    }
    
    
    ##active patients
    totPatients <- length( unique( all$patient_id ) )
    activePatients <- unique( input[ input$diagnosis_code %in% codesUnderStudy,1] )
    input <- input[ input$patient_id %in% activePatients, ]
    ##
    
    ##
    totPatients
    length( activePatients )
    ##
        
    codePairs <- function ( pt ){
        pp <- input[ input$patient_id == pt, ]
        icd9C <- unique( pp$diagnosis_code )
        codes.f <- codesUnderStudy[codesUnderStudy %in% icd9C]
        icd9C.c <- unique(do.call(c, apply(expand.grid(codes.f, icd9C), 1, combn, m=2, simplify=FALSE)))
        icd9C.f <- icd9C.c[sapply(icd9C.c, function(x) x[1] != x[2])]
    }
    
    message( "Generating the cAnalysis object" )
    
    
    finalCP  <- parallel::mclapply( activePatients, codePairs, mc.preschedule = TRUE, mc.cores = cores )
    
    finalCP <- finalCP[ sapply(finalCP, function(x) { length(x) != 0 }) ]
    
    f <- function( j ){ t( data.frame( j ) ) }
    unnest <-  do.call( f, list( j = finalCP  ) )
    unnest <- unnest[!duplicated(unnest), ]
    unnest <- lapply(1:nrow(unnest), function(ii) unnest[ii, ])
    
    
    
    resultado <- parallel::mclapply( unnest, tableData, mc.preschedule = TRUE, mc.cores = cores, data = all, lenActPa=totPatients)
    resultad2 <- do.call("rbind", resultado )
    resultad2 <- as.data.frame( resultad2, stringsAsFactors=FALSE )
    
    
    
    colnames(resultad2) <- c( "disAcode", "disBcode", "disA", "disB", "AB", "AnotB", "BnotA", "notAnotB", "fisher", "oddsRatio", "relativeRisk", "phi" )
    
    
    resultad2$expect <-  as.numeric( resultad2$disA ) * as.numeric( resultad2$disB ) / totPatients
    resultad2$score  <- log2( ( as.numeric( resultad2$AB ) + 1 ) / ( resultad2$expect + 1) )
    resultad2        <- resultad2[ with( resultad2, order( resultad2$fisher ) ), ]
    resultad2$fdr    <- p.adjust( as.numeric( resultad2$fisher ), method = "fdr", n = nrow( resultad2 ) )
    # 
    # resultad2$pair   <- NA
    # for(cont in 1:nrow(resultad2)){
    #     pairDis <- sort(c(resultad2$disAcode[cont], resultad2$disBcode[cont]))
    #     resultad2$pair[cont] <- paste(pairDis[1], pairDis[2], sep="*")
    # }
    # 
    # resultad2 <- resultad2[!duplicated(resultad2$pair),]
    # resultad2 <- resultad2[,c(1:15)]

    if ( !missing( score ) ) {
        resultad2 <- resultad2[ resultad2$score > score, ]
    }
    if ( !missing( fdr ) ) {
        resultad2 <- resultad2[ resultad2$fdr < fdr, ]
    }
    if ( !missing( oddsRatio ) ) {
        resultad2 <- resultad2[ resultad2$oddsRatio > oddsRatio, ]
    }
    if ( !missing( rr ) ) {
        resultad2 <- resultad2[ resultad2$rr > rr, ]
    }
    if ( !missing( phi ) ) {
        resultad2 <- resultad2[ resultad2$phi > phi, ]        
    }
    
    resultad2$fisher <- round(as.numeric(resultad2$fisher), 3)
    resultad2$oddsRatio <- round(as.numeric(resultad2$oddsRatio), 3)
    resultad2$relativeRisk <- round(as.numeric(resultad2$relativeRisk), 3)
    resultad2$phi <- round(as.numeric(resultad2$phi), 3)
    resultad2$expect <- round(as.numeric(resultad2$expect), 3)
    resultad2$score <- round(as.numeric(resultad2$score), 3)
    resultad2$fdr <- round(as.numeric(resultad2$fdr), 3)
    
    comb_table_rank <- resultad2
    comb_table_rank$scoreRank <- rank(-comb_table_rank$score)
    comb_table_rank$fisherRank <- rank(comb_table_rank$fisher)
    comb_table_rank$fdrRank <- rank(comb_table_rank$fdr)
    comb_table_rank$oddsRatioRank <- rank(-comb_table_rank$oddsRatio)
    comb_table_rank$rrRank <- rank(-comb_table_rank$relativeRisk)
    comb_table_rank$phiRank <- rank(-comb_table_rank$phi)        
    
    ## sum the ranking and create the "sumRank"
    col_num_before <- dim(resultad2)[2]
    col_num_after <- dim(comb_table_rank)[2]
    
    comb_table_rank$sum <- apply(comb_table_rank[, col_num_before:col_num_after], 1, sum)
    comb_table_rank$sumRank <- rank(comb_table_rank$sum)
    
    
    ## get the sort file and return
    comb_table_sort <- comb_table_rank[order(comb_table_rank$sumRank), ]
    
    if( nrow( comb_table_sort ) == 0 ){
        warning("None of the disease pairs has pass the filters") 
    }
    
    
    cAnalysis <- new( "cAnalysis", 
                      ageMin    = ageRange[ 1 ], 
                      ageMax    = ageRange[ 2 ], 
                      gender    = gender, 
                      patients  = totPatients,
                      tpatients = length(activePatients),
                      prevalence= (length(activePatients)/totPatients)*100,
                      rangeOR = paste0("[", round(min(as.numeric(comb_table_sort$oddsRatio)), digits = 3)," , ", round(max(as.numeric(comb_table_sort$oddsRatio)), digits = 3)  ,"]" ),
                      rangeRR = paste0("[", round(min(as.numeric(comb_table_sort$relativeRisk)), digits = 3)," , ", round(max(as.numeric(comb_table_sort$relativeRisk)), digits = 3)  ,"]" ),
                      rangePhi= paste0("[", round(min(as.numeric(comb_table_sort$phi)), digits = 3)," , ", round(max(as.numeric(comb_table_sort$phi)), digits = 3)  ,"]" ),
                      dispairs  = nrow( comb_table_sort ),
                      result    = comb_table_sort[,c(1:15,23)] 
    )
    return( cAnalysis )
}

tableData <- function ( pairCode, data, lenActPa ) {
    
    code1 <- pairCode[[ 1 ]]
    code2 <- pairCode[[ 2 ]]
    
    data$diagnosis_code[is.na(data$diagnosis_code)] <- 0
    
    dis1 <- data[ data$diagnosis_code == code1, ]
    dis2 <- data[ data$diagnosis_code == code2, ]
    
    dis12 <- dis2[ dis2$patient_id %in% dis1$patient_id, ]
    
    disAcode <- code1
    disBcode <- code2
    disA     <- length( unique ( dis1$patient_id ) )
    disB     <- length( unique ( dis2$patient_id ) )
    AB       <- length( unique ( dis12$patient_id ) )
    AnotB    <- disA - AB
    BnotA    <- disB - AB
    notAB    <- as.numeric(lenActPa) - AB - AnotB - BnotA
    
    mm <- matrix( c( AB, AnotB, BnotA, notAB), nrow = 2 )
    
    tryCatch( {ff <- fisher.test( mm )}, error=function(msg) {
      message(msg)
      message("code1:", code1, " - code2:", code2)
    })
    
    relativeRisk <- (as.numeric(AB)*as.numeric(lenActPa))/as.numeric(disA* disB)
    den <- as.numeric(disA*disB)*(as.numeric(lenActPa)-as.numeric(disA))*(as.numeric(lenActPa)-as.numeric(disB))
    num <- (as.numeric(AB)*as.numeric(lenActPa))-as.numeric(disA*disB)
    phi <- ((num)/sqrt(den))
    
    oddsRatio <- (as.numeric(AB)*as.numeric(notAB))/(as.numeric(AnotB)*as.numeric(BnotA))
    
    c( disAcode, disBcode, disA, disB, AB, AnotB, BnotA, notAB, ff$p.value, oddsRatio, relativeRisk, phi )    
    
  }
  








