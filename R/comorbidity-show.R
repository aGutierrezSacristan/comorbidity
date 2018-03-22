setMethod( "show",
    signature = "comorbidity",
    definition = function( object ) {
        cat( "Object of class 'comorbidity'\n" )
        cat( " . Search:                   ", object@search, "\n" )
        cat( " . Only comorbidities between index diseases:", object@intraCode, "\n" )
        cat( " . Aggregate the disease codes:", object@aggregated, "\n" )
        cat( " . Number of Input Index Diseases:  ",   object@tDiseases, "\n" )
        cat( " . Number of Index Diseases Present:", object@indexDis, "\n" )
        cat( " . Number of Concomintant Diseases: ",  object@nfDisease, "\n" )
        cat( " . Number of Patients:              ", length( unique ( object@qresult$patient_id ) ), "\n" )
    }
)

setMethod( "show",
           signature = "cAnalysis",
           definition = function( object ) {
               cat( "Object of class 'cAnalysis'\n" )
               cat( " . Minimum age:", object@ageMin, "\n" )
               cat( " . Maximum age:", object@ageMax, "\n" )
               cat( " . Sex  :", object@sex, "\n" )
               cat( " . Patients in the age and sex interval:", object@patients, "\n" )
               cat( " . Patients suffering the index disease(s):", object@tpatients, "\n" )
               cat( " . Disease Prevalence:", round(as.numeric(object@prevalence), 3), "\n")
               cat( " . Odds ratio range:", object@rangeOR, "\n")
               cat( " . Relative risk range:", object@rangeRR, "\n")
               cat( " . Phi value range:", object@rangePhi, "\n")
               cat( " . Number of comorbidities:", object@dispairs, "\n")
           }
)


setMethod( "show",
           signature = "molecularComorbidity",
           definition = function( object ) {
               cat( "Object of class 'molecularComorbidity'\n" )
               cat( " . Search:                     ", object@search, "\n" )
               cat( " . Aggregate the disease codes:", object@aggregated, "\n" )
               cat( " . Number of Input Diseases:          ", object@indexDis, "\n" )
               cat( " . Number of Index Diseases Present:  ", object@nfDisease, "\n" )
               cat( " . Number of Genes   :                ", object@nGenes, "\n" )
           }
)

setMethod( "show",
           signature = "molecularcAnalysis",
           definition = function( object ) {
               cat( "Object of class 'molecularcAnalysis'\n" )
               cat( " . Minimum number of genes overlaped:", object@ovlpMin, "\n" )
               cat( " . Maximum number of genes overlaped:", object@ovlpMax, "\n" )
               cat( " . Jaccard Minimum value:            ", object@jaccardMin, "\n" )
               cat( " . Jaccard Maximum value:            ", object@jaccardMax, "\n" )
               cat( " . P-value     :                     ", object@pValue, "\n" )
               cat( " . Number of comorbidities:          ", object@dispairs, "\n")
           }
)