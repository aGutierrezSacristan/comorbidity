#' Query your data and generates a \code{datafame}
#'
#' Given an object of class \code{cAnalysis} obtained from a comorbidity analysis,
#' the comorbidity direction can be estimated. An object of class \code{dataframe} 
#' is obtained.
#'
#' @param input A \code{cAnalysis} object, containing the comorbidity data for
#' a set of patients. 
#' @param databasePth Determines the path where the intermediate RData objects 
#' have been created. It is the same path where the three required input files 
#' (patientData, diagnosisData, admissionData) are located. 
#' @param minPairs Determines the minimum number of patients that must suffer the
#' comorbidity. By default the \code{minPairs} value is set to \code{1}. The value 
#' of the argument can be changed to any other numeric variable.
#' @param gender Determine what is the gender of interest for 
#' performing the comorbidity analysis. 
#' @param ageRange Determines what is the age range of interest for
#' performing the comorbidity analysis. By default it is set from 0 to 100 
#' years old. 
#' @param days Determines the number of days of difference needed for considering 
#' two diseases as comorbid.  
#' @param dataSep Determines the separator symbol used in the admission date.
#' By default \code{"-"}. 
#' @param correctionMethod A binomial test for each pair of diseases is performed to assess 
#' the null hypothesis of independence between the two diseases. The Bonferroni correction 
#' ("bonferroni") is applied to correct for multiple testing by default. 
#' However user can select the best correction method for the analysis. The adjustment methods 
#' include the Benjamini-Hochberg false discovery rate method ("fdr"),  Holm correction ("holm"), 
#' Hochberg correction ("hochberg"), Hommel ("hommel") and Benjamini & Yekutieli ("BY"). 
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return An object of class \code{dataframe}
#' @examples
#' load(system.file("extdata", "comorMale.RData", package="comoRbidity"))
#' directionEx <- directionality( input       = comorMale, 
#'                                databasePth    = system.file("extdata", package="comoRbidity"), 
#'                                minPairs         = 1, 
#'                                gender           = "Male",
#'                                ageRange         = c(0,80),
#'                                days             = "0",
#'                                dataSep          = "-", 
#'                                correctionMethod = "bonferroni"
#'               )
#' @export directionality

directionality <- function( input, databasePth, minPairs = 1, gender, ageRange=c(0,100), days, dataSep="-", correctedPval = 1, correctionMethod = "bonferroni", verbose ){
    
    #check if the input object is of class molecularComorbidity
    message("Checking the input objects")
    checkClass <- class(input)[1]

    if(checkClass != "cAnalysis" ){
        message("Check the input object. Remember that this
                    object must be obtained after applying the comorbidityAnalysis
                    function to your input file. The input object class must
                    be:\"cAnalysis\"")
        stop()
    }
    
    load(paste0(databasePth, "/direction.RData"))
    
    direction <- direccionality@qresult
    #direction$data <- direction$admissionStartDate
    direction$age <- as.numeric(direction$age)
    direction <- direction[direction$patient_gender == gender, ]
    #direction$admissionStartDate <- do.call(rbind, strsplit(direction$admissionStartDate, dataSep) )[,1]
    #direction$patient_dateBirth <- do.call(rbind, strsplit(direction$patient_dateBirth, dataSep) )[,1]
    #direction$age <- as.numeric(direction$admissionStartDate) - as.numeric(direction$patient_dateBirth)
    direction <- direction[ direction$age > ageRange[ 1 ] & direction$age < ageRange[ 2 ], ]
    
    input <- input@result
    input <- input[as.numeric(input$AB) >=minPairs, ]
    
    if(nrow(input) == 0){
        stop("There are not any disease comorbidities with this minimum pairs")
    }
    
    pairs <- input[,c(1:2)]
    pairs$AtoB <- 0
    pairs$BtoA <- 0
    pairs$test <- 0
    
    for(i in 1:nrow(pairs)){
        
        code1 <- pairs[i,1]
        code2 <- pairs[i,2]
        
        studyCode1 <- direction[direction$diagnosis_code==code1,]
        studyCode1 <- studyCode1[with(studyCode1, order(admissionStartDate, patient_id)), ]
        studyCode1 <- studyCode1[!duplicated(studyCode1$patient_id), ]
        
        studyCode2 <- direction[direction$diagnosis_code==code2,]
        studyCode2 <- studyCode2[with(studyCode2, order(admissionStartDate, patient_id)), ]
        studyCode2 <- studyCode2[!duplicated(studyCode2$patient_id), ]
        
        studyCode1 <- studyCode1[studyCode1$patient_id %in% studyCode2$patient_id,]
        studyCode2 <- studyCode2[studyCode2$patient_id %in% studyCode1$patient_id,]
        
        df <- merge(studyCode1, studyCode2, by="patient_id")
        dfColnames <- c("patient_id", "data.x", "data.y")
        df  <- df[,colnames(df)%in% dfColnames] 
        
        df$diff <- as.Date(as.character(df$data.x), format=paste0("%Y", dataSep, "%m",dataSep ,"%d"))-
            as.Date(as.character(df$data.y), format=paste0("%Y", dataSep, "%m",dataSep ,"%d"))
        
        
        pairs$AtoB[i] <- nrow(df[df$diff<days,])
        pairs$BtoA[i] <- nrow(df[df$diff>days,])
        
        if(pairs$AtoB[i]==0 & pairs$BtoA[i]==0){
            pairs$test[i] <- "NA"
            
        }else{
            Btest <- binom.test( x = pairs$AtoB[i], 
                                 n = pairs$AtoB[i]+pairs$BtoA[i], 
                                 p = 0.5, 
                                 alternative = "two.sided")
            
            pairs$test[i] <- Btest$p.value
        }
        

    }
    pairs <- pairs[pairs$test != "NA",]
    pairs$correctedPvalue <- p.adjust( as.numeric( pairs$test ), method = correctionMethod, n = nrow( pairs ) )
    
    pairs$result <- NA
    
    for(i in 1:nrow(pairs)){
        if(pairs$test[i]>0.05){
            pairs$result[i]<- "No directionality"
        }
        if(pairs$test[i]<=0.05){
            if(pairs$AtoB[i]>pairs$BtoA[i]){
                pairs$result[i] <- "From A to B"
            }
            if(pairs$AtoB[i]<pairs$BtoA[i]){
                pairs$result[i] <- "From B to A"
            }
        }
        
    }
    
    
    return(pairs)
}




