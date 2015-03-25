# lavaan.survey.fiml
## Complex survey standard errors for Structural Equation Models in R

This is an experimental package, whose functionality will eventually be merged into lavaan.survey. 
Currently it estimates standard errors for Structural Equation Models (SEM) estimated using lavaan in R, 
accounting for clustering and stratification. It also allows for models with missing data estimated using FIML, 
where lavaan.survey requires multiple imputation to deal with missing data.

It does not (yet) account for weights, so there is no point estimation in this package (see lavaan.survey for that). 
