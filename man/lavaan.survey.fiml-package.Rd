\name{lavaan.survey.fiml-package}
\alias{lavaan.survey.fiml-package}
\alias{lavaan.survey.fiml}
\docType{package}
\title{
Complex survey standard errors (clustering and stratification) for Structural Equation Models
}
\description{
Similarly to lavaan.survey, this package ties together lavaan (for
    SEM) and survey (for complex survey designs). Contrary to lavaan.survey,
    however, this experimental package allows the lavaan model to have been
    fitted by FIML. It does not do estimation, only standard errors. Since the
    point estimates as well as se's will be affected by weights, weights cannot
    be accounted for.
    }
\details{
\tabular{ll}{
Package: \tab lavaan.survey.fiml\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2015-03-25\cr
License: \tab GPL\cr
}


}
\author{
Daniel Oberski
}
\references{

}

\keyword{ package }
\seealso{
}
\examples{
library(lavaan)

###### A single group example #######

# European Social Survey Denmark data (SRS)
data(ess.dk, package = "lavaan.survey")

# A saturated model with reciprocal effects from Saris & Gallhofer
dk.model <- "
  socialTrust ~ 1 + systemTrust + fearCrime
  systemTrust ~ 1 + socialTrust + efficacy
  socialTrust ~~ systemTrust
"
lavaan.fit <- lavaan(dk.model, data=ess.dk, auto.var=TRUE, estimator="MLM")
summary(lavaan.fit)

sqrt(diag(vcov_complex(lavaan.fit, ids = ess.dk$intnum)))
}
