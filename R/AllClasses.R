setClass( "comorbidity",
          representation =
              representation( search       = "character",  # single or list
                              intraCode    = "logical",
                              aggregated   = "logical", 
                              tDiseases    = "numeric",
                              indexDis     = "numeric",
                              nfDisease    = "numeric",    # number of the initial diseases present in our database
                              nPatient     = "numeric",    # number of patients suffering them
                              indexDisList = "character",
                              qresult      = "data.frame"  # result
              ),
          prototype = 
              prototype( search    = "",
                         intraCode = logical(),
                         aggregated= logical(),
                         tDiseases = numeric(),
                         indexDis  = numeric(),
                         nfDisease = numeric(),
                         nPatient  = numeric(),
                         indexDisList =  "",
                         qresult   = data.frame()
              )
)


setClass( "cAnalysis",
          representation =
              representation( ageMin    = "numeric",  # single or list
                              ageMax    = "numeric",    # max age
                              gender    = "character",    # gender
                              patients  = "numeric",    # subsetPatients
                              tpatients = "numeric",    # totalPatients
                              prevalence= "numeric",    # prevalence respet to the total population
                              rangeOR   = "character",    # range value of the OR
                              rangeRR   = "character",    # range value of the relative risk
                              rangePhi  = "character",    # range value of the phi value
                              dispairs  = "numeric",    # number of pairs
                              result    = "data.frame"  # result
              ),
          prototype = 
              prototype( ageMin    = numeric(),
                         ageMax    = numeric(),
                         gender    = character(),
                         patients  = numeric(),
                         tpatients = numeric(),
                         prevalence= numeric(),
                         rangeOR = character(),
                         rangeRR = character(),
                         rangePhi= character(),
                         dispairs  = numeric(),
                         result    = data.frame()
              )
)

setClass( "molecularComorbidity",
          representation =
              representation( search       = "character",  # single or list
                              aggregated   = "logical", 
                              indexDis     = "numeric",
                              nfDisease    = "numeric",    # number of the initial diseases present in our database
                              nGenes       = "numeric",    # number of patients suffering them
                              indexDisList = "character",
                              qresult      = "data.frame"  # result
              ),
          prototype = 
              prototype( search    = "",
                         aggregated= logical(),
                         indexDis  = numeric(),
                         nfDisease = numeric(),
                         nGenes    = numeric(),
                         indexDisList =  "",
                         qresult   = data.frame()
              )
)

setClass( "molecularcAnalysis",
          representation =
              representation( ovlpMin    = "numeric",  # minimum value overlap
                              ovlpMax    = "numeric",  # maximum value overlap
                              jaccardMin = "numeric",  # minimum value jaccard
                              jaccardMax = "numeric",  # maximum value jaccard
                              pValue     = "logical",  # pValue
                              tdiseases  = "numeric",  # total diseases with comorbidities
                              dispairs   = "numeric",  # number of comorbidity pairs
                              result     = "data.frame"# result
              ),
          prototype = 
              prototype( ovlpMin    = numeric(),
                         ovlpMax    = numeric(),
                         jaccardMin = numeric(),
                         jaccardMax = numeric(),
                         pValue     = logical(),
                         tdiseases  = numeric(),
                         dispairs   = numeric(),
                         result     = data.frame()
              )
)
