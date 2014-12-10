\name{standrm}
\alias{standrm}
\title{Dose-response models with rstan}
\description{Formula interface for fitting dose-response curves with rstan}
\usage{
standrm(formula, data, fct, curveid = NULL, random = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A two-sided formula with response and dose covariate names}
  \item{data}{A dataframe containing the variables in the model}
  \item{fct}{A drc model function, currelntly only LL.5() is supported.}
  \item{curveid}{A two-sided formula with parameter names (b + c + d) on the left and a factor name on the right hand side, defining several clusters, modelled as fixed effects}
  \item{random}{A two-sided formula with parameter names (b + c + d) on the left and a factor name on the right hand side, defining several clusters, modelled as random effects}
  \item{\dots}{Additional arguments passed to function stan in ackage rstan}
}

\value{
\item{stan_dat}{rstan compatible list containing the data}
\item{model}{A character vector containing the stan model}
\item{stan}{S4 output of rstan containing the samples from the posterior distributions}
\item{pars}{Character vector with fixed effects parameter names}
\item{random}{Character vector with random effects parameter names}
}

\author{
Daniel Gerhard
}

\examples{
\dontrun{
data(spinach)
head(spinach)

m <- standrm(SLOPE ~ DOSE, data=spinach, 
             fct=LL.5(fixed=c(NA, NA, NA, NA, 0)), 
             curveid=b + c + e ~ HERBICIDE, 
             random=c + e ~ CURVE, 
             iter=3000, chains=4)
print(m)
ranef(m)
VarCorr(m)

ED(m, respLev=c(25, 50, 75))

}
}

\keyword{models}
