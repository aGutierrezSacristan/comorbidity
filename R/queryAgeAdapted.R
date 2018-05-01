#' Query your data and generates a \code{comoRbidity}
#'
#' Given 4 files (patientData.txt, diagnosisData.txt, admissionData.txt and 
#' indexDiseaseCodes.txt), generates some RData and create an object of type 
#' \code{comoRbidity}.
#'
#' @param databasePth Determines the path where the three required input files 
#' (patientData, diagnosisData, admissionData) are located. 
#' @param codesPth Determines the path where the file with the index diseases is 
#' located (indexDiseaseCode)
#' @param birthDataSep Determines the separator symbol used in the birth date. 
#' @param intraCodes Comorbidities will be estimated only between the index 
#' diseases. By default \code{FALSE} 
#' @param aggregatedDis Data extraction is done using the Agg column from the 
#' index disease file,  that collapse the diseases in a supperior class.
#' By default \code{FALSE} 
#' @param python By default \code{FALSE}. Change it to \code{TRUE} to run the 
#' query quicklier using python. It is necessary to have python installed.
#' @param verbose By default \code{FALSE}. Change it to \code{TRUE} to get a
#' on-time log from the function.
#' @param warnings By default \code{TRUE}. Change it to \code{FALSE} to don't see
#' the warnings.
#' @return An object of class \code{comoRbidity}
#' @examples
#' ex1 <- queryAgeAdapted( databasePth      = system.file("extdata", package="comoRbidity"),
#'               codesPth         = system.file("extdata", package="comoRbidity"),
#'               birthDataSep     = "-",
#'               python           = FALSE)
#' @export queryAgeAdapted


queryAgeAdapted <- function( databasePth, codesPth, birthDataSep, intraCodes = FALSE, aggregatedDis = FALSE, python= FALSE, verbose = FALSE, warnings= TRUE) {
    
    message( "Starting querying the index diseases in the dataset" )
    
    if( python == TRUE ){
        
        message( "Loading the input datasets" )
        
        pythonPth <- system.file("comorbidity.py", package="comoRbidity")
        qq <- paste0("'",pythonPth,"' '", databasePth,"/' '",codesPth,"/' '", admissionDataSep,"' '", birthDataSep, "' '", "FALSE", "' '",intraCodes,"' '",aggregatedDis, "'" )
        
        message( "Starting querying for your index diseases" )
        
        system(qq)
        
        final <- read.delim( file.path( databasePth, "final"  ), 
                             header = FALSE, sep = "\t", 
                             colClasses = "character" )
        
        direct <- read.delim( file.path( databasePth, "direct"  ), 
                              header = FALSE, sep = "\t", 
                              colClasses = "character" )
        
        all <- read.delim( file.path( databasePth, "all"  ), 
                           header = FALSE, sep = "\t", 
                           colClasses = "character" )
        
        codes <- read.delim ( file.path(codesPth, "indexDiseaseCode.txt"), 
                              header=TRUE, sep="\t", 
                              colClasses="character" )
        
        
        if( aggregatedDis == TRUE){
            message( "Aggregating the disorders of interest" )
            codesUnderStudy    <- as.character(unique(codes$Agg))
        }else if( aggregatedDis == FALSE){
            codesUnderStudy    <- as.character(unique(codes$Code))
        }
        
        colnames(final) <- c("patient_id", "admissionStartDate", "diagnosis_code",
                             "patient_sex", "patient_dateBirth", "age")
        
        colnames(all) <- c("patient_id", "admissionStartDate", "diagnosis_code",
                           "patient_sex", "patient_dateBirth", "age")
        
        colnames(direct) <- c("patient_id", "admissionStartDate", "diagnosis_code",
                              "patient_sex", "patient_dateBirth", "age")
        
    }
    else{
        message( "Loading the input datasets" )
        patients <- read.delim( file.path( databasePth, "patientData.txt"  ), 
                                header = TRUE, sep = "\t", 
                                colClasses = "character" )
        
        diagnostic <- read.delim( file.path( databasePth, "diagnosisData.txt"  ), 
                                  header = TRUE, sep = "\t", 
                                  colClasses = "character" )
        
        admission <- read.delim( file.path( databasePth, "admissionData.txt"  ), 
                                 header = TRUE, sep = "\t", 
                                 colClasses = "character" )
        
        codes <- read.delim ( file.path(codesPth, "indexDiseaseCode.txt"),
                              header=TRUE, 
                              sep="\t", 
                              colClasses="character" ) 
        
        message("Checking the patientData file structure")
        colnamesPatients   <- c("patient_id","patient_sex", "patient_dateBirth")   
        check <- colnamesPatients[colnamesPatients %in% colnames(patients)]
        if(length(check) != length(colnamesPatients)){
            message("Check the patientData file structure. Remember that this
                    file must contain at least three columns with the column 
                    names as follows:\n -> patient_id \n -> patient_sex \n -> patient_dateBirth")
            stop()
        }
        
        colnamesDiagnosis   <- c("patient_id", "admission_id", "diagnosis_code") 
        checkD <- colnamesDiagnosis[colnamesDiagnosis %in% colnames(diagnostic)]
        message("Checking the diagnosisData file structure")
        
        if(length(checkD) != length(colnamesDiagnosis)){
            message("Check the diagnosisData file structure. Remember that this
                    file must contain at least three columns with the column 
                    names as follows:\n -> patient_id \n -> admission_id \n -> diagnosis_code")
            stop()
        }
        
        message("Checking the admissionData file structure")
        
        colnamesAdmission   <- c("patient_id", "admission_id", "admissionStartDate")
        checkA <- colnamesAdmission[colnamesAdmission %in% colnames(admission)]
        if(length(checkA) != length(colnamesAdmission)){
            message("Check the admissionData file structure. Remember that this
                    file must contain at least three columns with the column 
                    names as follows:\n -> patient_id \n -> admission_id \n -> admissionStartDate")
            stop()
            
        }
        
        message("Checking the patients")
        
        if( length(unique(patients$patient_id)) != length(unique(diagnostic$patient_id)) |
            length(unique(patients$patient_id)) != length(unique(admission$patient_id))){
            message("The number of patients in the files is not the same")
            stop()
        }
        
        
        message( "Starting querying for your index diseases" )
            
        if( intraCodes == TRUE ){
            diagnostic <- diagnostic[diagnostic$diagnosis_code %in% codes$Code,]
        }
        if( intraCodes == FALSE ){
            diagnostic <- diagnostic
        }
        
        if( aggregatedDis == TRUE){
            
            message( "Aggregating the disorders of interest" )
            
        
            for( i in 1:nrow(diagnostic)){
                if(diagnostic$diagnosis_code[i] %in% codes$Code){
                    aggCode <-  codes[codes$Code == diagnostic$diagnosis_code[i], ]
                    diagnostic$diagnosis_code[i] <- aggCode$Agg
                }else{
                    next
                }
            }

            codesUnderStudy    <- as.character(unique(codes$Agg))
            patientsUnderStudy <- diagnostic[diagnostic$diagnosis_code %in% codesUnderStudy,]
            patientsUnderStudy <- unique(patientsUnderStudy$patient_id)
        }
        
        if( aggregatedDis == FALSE){
            #save data of interest
            codesUnderStudy    <- as.character(unique(codes$Code))
            patientsUnderStudy <- diagnostic[as.character(diagnostic$diagnosis_code) %in% as.character(codesUnderStudy),]
            patientsUnderStudy <- unique(patientsUnderStudy$patient_id)
         }
            
     
        patientsSelection  <- patients[patients$patient_id %in% patientsUnderStudy,] 
        patientsSelection  <- patientsSelection[,colnames(patientsSelection)%in% colnamesPatients] 
        totalPatients      <- patients[,colnames(patients)%in% colnamesPatients] 
        
       
        diagnosisSelection  <- diagnostic[diagnostic$patient_id %in% patientsUnderStudy,]
        diagnosisSelection  <- diagnosisSelection[,colnames(diagnosisSelection)%in% colnamesDiagnosis] 
        diagnosisSelection$pair <- paste(diagnosisSelection$patient_id, diagnosisSelection$admission_id, sep="*")
        totalDiagnosis      <- diagnostic[,colnames(diagnostic)%in% colnamesDiagnosis] 
        totalDiagnosis$pair <- paste(totalDiagnosis$patient_id, totalDiagnosis$admission_id, sep="*")
            
        pairsUnderStudy <- diagnosisSelection$pair
       
        message("Generating the resulting objects")
        
        admissionSelection  <- admission[,colnames(admission)%in% colnamesAdmission] 
        admissionSelection$pair <- paste(admissionSelection$patient_id, admissionSelection$admission_id, sep="*")
        admissionSelection  <- admissionSelection[admissionSelection$pair %in% pairsUnderStudy,] 
        totalAdmission  <- admission[,colnames(admission)%in% colnamesAdmission]
        totalAdmission$pair <- paste(totalAdmission$patient_id, totalAdmission$admission_id, sep="*")
            
        final <- plyr::join(admissionSelection,diagnosisSelection, by="pair")
        final <- final[,c("patient_id","admissionStartDate","diagnosis_code")] 
        all <- plyr::join(totalAdmission,totalDiagnosis, by="pair")
        all <- all[,c("patient_id","admissionStartDate","diagnosis_code")] 
        all <- all[!duplicated(all), ]
        
        final <- plyr::join(final, patientsSelection, by="patient_id" )
        all   <- plyr::join(all, totalPatients, by="patient_id" )
            
        direct <- final
        direct$data <- final$admissionStartDate
            
        #final$admissionStartDate <- do.call(rbind, strsplit(final$admissionStartDate, admissionDataSep) )[,1]
        #final$patient_dateBirth <- do.call(rbind, strsplit(final$patient_dateBirth, birthDataSep) )[,1]
        
        #all$admissionStartDate <- do.call(rbind, strsplit(all$admissionStartDate, admissionDataSep) )[,1]
        #all$patient_dateBirth <- do.call(rbind, strsplit(all$patient_dateBirth, birthDataSep) )[,1]
        
        #direct$admissionStartDate <- do.call(rbind, strsplit(direct$admissionStartDate, admissionDataSep) )[,1]
        #direct$patient_dateBirth <- do.call(rbind, strsplit(direct$patient_dateBirth, birthDataSep) )[,1]
            
        final$age <- as.numeric(as.Date(final$admissionStartDate)-as.Date(final$patient_dateBirth))%/%365
        all$age <- as.numeric(as.Date(all$admissionStartDate)-as.Date(all$patient_dateBirth))%/%365
        direct$age <- as.numeric(as.Date(direct$admissionStartDate)-as.Date(direct$patient_dateBirth))%/%365
        
        }
    
    indexDiseases <- codesUnderStudy[codesUnderStudy %in% final$diagnosis_code]
    
    
    if( verbose ) {
        message( "There are ", length( unique ( final$patient_id)), " patients in our database that suffer at least 1 of the ", length( unique ( codesUnderStudy ) ), " different input diseases of our list.")
    }
    
    #with the data we have, we create a comorbidity object
    result <- new( "comorbidity", 
                   search       = ifelse( length(codesUnderStudy) > 1, "list", "single" ), 
                   intraCode    = intraCodes,
                   aggregated   = aggregatedDis, 
                   tDiseases    = length(unique( codesUnderStudy )),
                   indexDis     = length(unique( indexDiseases )),                   
                   nfDisease    = length(unique(final$diagnosis_code) ), 
                   nPatient     = length( unique ( final$patient_id ) ), 
                   indexDisList = unique( indexDiseases ),
                   qresult      = final  
    )
    
    #we create a second object with all the data of the dataset
    allData <- new( "comorbidity", 
                    search    = ifelse( length(codesUnderStudy) > 1, "list", "single" ), 
                    intraCode = intraCodes,
                    aggregated= aggregatedDis, 
                    tDiseases = length(unique( codesUnderStudy )),
                    indexDis  = length( unique( indexDiseases )) ,                  
                    nfDisease = length(  unique(final$diagnosis_code) ), 
                    nPatient  = length( unique ( all$patient_id ) ), 
                    indexDisList = unique( indexDiseases ),
                    qresult   = all 
    )
    save(allData, file=paste0(databasePth, "allData.RData"))
    
    
    #we create a second object with all the data of the dataset
    direccionality <- new( "comorbidity", 
                           search    = ifelse( length(codesUnderStudy) > 1, "list", "single" ), 
                           intraCode = intraCodes,
                           aggregated= aggregatedDis, 
                           tDiseases = length(unique( codesUnderStudy )),
                           indexDis  = length( unique( indexDiseases )),                   
                           nfDisease = length(  unique(final$diagnosis_code) ),  
                           nPatient  = length( unique ( direct$patient_id ) ), 
                           indexDisList = unique( indexDiseases ),
                           qresult   = direct 
    )
    save(direccionality, file=paste0(databasePth, "direction.RData"))
    
    return( result )

}



