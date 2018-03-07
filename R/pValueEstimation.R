pValueEstimation <- function(input, nboot = 100, ncores = 1, verbose = FALSE) {

  #extract universe genes
  oql <- "DEFINE
  c0='/data/gene_disease_summary',
  c1='/data/genes',
  c3='/data/sources'
  ON
  'http://www.disgenet.org/web/DisGeNET'
  SELECT
  c0 (geneId, score, diseaseId, source ),
  c1 (symbol)
  FROM
  c0
  WHERE
  c3 = 'CURATED'
  ORDER BY
  c0.score DESC"
  
  getUrlDis <- function() {
      url <- "http://www.disgenet.org/oql"
      return( url )
  }
  
  dataTsv <- RCurl::getURLContent(
      getUrlDis(),
      readfunction  = charToRaw(oql),
      upload        = TRUE,
      customrequest = "POST"
  )

  res <- read.csv(textConnection(dataTsv), header = TRUE, sep="\t")
  universe <- unique(as.character(res$c0.geneId))

  message(paste0("A total of ", length(universe), 
                 " genes obtained from DisGeNET CURATED database are 
                 being used for the bootstrap process\n")
          )

  
  input$pval <- NA
  message("Pvalue estimation for each comorbidity")
  
  for(i in 1:nrow( input )){

      if(verbose == TRUE){
          message(paste0("\n\t->Comorbidity pair ", i, " of ", nrow(input), " total diseases' pairs."))
      }
      
      
      bb <- ji.internal(input$geneV1[i], input$geneV2[i], universe, nboot, ncores)
      pval <- (sum(bb > input$jaccard[i]) * 1.0) / (nboot+1)
      input$pval[i] <- round(pval, 5)
    }

    return(input)

}


ji.internal <- function(len1, len2, universe, nboot, ncores) {
  if (!requireNamespace("parallel", quietly = TRUE)) {
    pfun <- lapply
  } else {
    pfun <- parallel::mclapply
  }

  unlist(pfun(1:nboot, function(ii) {
    g1 <- sample( universe, len1 )
    g2 <- sample( universe, len2 )
    ja.coefr <- length(intersect(g1, g2)) / length(union(g1, g2))
  }, mc.cores = ncores))
}
