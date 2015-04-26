# Complex survey standard errors for lavaan fits estimated with FIML
#' @import survey

# Defines a bread method for estfun in sandwich
bread.lavaan <- function(fit) {
  #-lavaan:::computeObservedInformation(lavmodel = fit@Model, lavsamplestats = fit@SampleStats, lavdata = fit@Data, estimator = fit@Options$estimator )
  lavaan:::computeExpectedInformation(lavmodel = fit@Model, lavsamplestats = fit@SampleStats, lavdata = fit@Data, estimator = fit@Options$estimator )
}

#' Calculate complex survey standard errors for lavaan fits estimated with FIML.
#' 
#' @param fit A lavaan Fit object.
#' @param ids Clusters specified as a vector or data frame (not a formula!)
#' @param strata Strata specified as a vector or data frame (not a formula!)
#' @return The estimated variance-covariance matrix of the paramters, accounting for complex sampling design.
#' @export
#' @examples
#' library(lavaan)
#' 
#' ###### A single group example #######
#' 
#' # European Social Survey Denmark data (SRS)
#' data(ess.dk, package = "lavaan.survey")
#' 
#' # A saturated model with reciprocal effects from Saris & Gallhofer
#' dk.model <- "
#'   socialTrust ~ 1 + systemTrust + fearCrime
#'   systemTrust ~ 1 + socialTrust + efficacy
#'   socialTrust ~~ systemTrust
#' "
#' lavaan.fit <- lavaan(dk.model, data=ess.dk, auto.var=TRUE, estimator="MLM")
#' summary(lavaan.fit)
#' # Use vcov_complex to obtain se's instead of lavaan.survey:
#' sqrt(diag(vcov_complex(lavaan.fit, ids = ess.dk$intnum)))
#' 
#' # Create random missings and use FIML to estimate the model
#' set.seed(4897)
#' ess.dk.mis <- ess.dk
#' ess.dk.mis[cbind(sample(1:nrow(ess.dk), size=200), sample(3:6, size=200, replace=TRUE))] <- NA
#' lavaan.fit.mis <- lavaan(dk.model, data=ess.dk.mis, auto.var=TRUE, missing="fiml")
#' summary(lavaan.fit.mis)
#' # Use vcov_complex to obtain se's accounting for clustering
#' sqrt(diag(vcov_complex(lavaan.fit.mis, ids = ess.dk$intnum)))
vcov_complex <- function(fit, ids=~1, strata=NULL) {
  
  # Treat casewisegrad*inv(Hessian) as data
  # This is also how svymle works.
  the_bread <- bread.lavaan(fit)
  if(min(eigen(the_bread, only.values=TRUE)$values) < .Machine$double.eps*10)
    warning("Singular information matrix.")
  db <- as.data.frame(lavaan::estfun.lavaan(fit) %*% MASS::ginv(the_bread))
  names(db) <- paste("V", 1:ncol(db), sep="")
  
  des <- svydesign(ids = ids, probs=~1, strata=strata, data = db)
  
  vcov_fit <- attr(svymean(make.formula(names(db)), design = des), "var")
  colnames(vcov_fit) <- rownames(vcov_fit) <- names(coef(fit)) # Conserve parameter names
  vcov_fit
}

