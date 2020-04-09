# comoRbidity

`comoRbidity` is an R package to analyze comorbidities from clinical data


## What is this repository for?

This report is used for `comoRbidity` package distribution while we walk thought publication process. 

## Package' Status

 * __Version__: 1.1.1
 * __Authors__: Alba Gutierrez-Sacristan (GRIB-UPF / HMS)
 * __Maintainer__: Alba Gutierrez-Sacristan (GRIB-UPF / HMS)

## How to start

### Installation

`comoRbidity` can be installed using `devtools` from this repository:

```R
library(devtools)
install_bitbucket("ibi_group/comoRbidity")
```

### comoRbidity:

The following lines show one examples of how `comoRbidity` R package can be used in comorbidity studies:

* __Example Query__

```R
library(comoRbidity)
ex1 <- query( databasePth      = system.file("extdata", package="comoRbidity"),
               codesPth         = system.file("extdata", "indexDiseaseCodes.txt", package="comoRbidity"),
               birthDataSep     = "-",
               admissionDataSep = "-",
               determinedCodes  = FALSE,
               python           = FALSE)
)
```

* __Heatmap Representation__

```R
library(comoRbidity)
load(system.file("extdata", "comorMale.RData", package="comoRbidity"))
htmp <- heatmapPlot( input      = comorMale, 
               selectValue       = "score", 
               cutOff     = 0.5
        )

```
 
```
