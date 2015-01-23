\name{standrm_count}
\alias{standrm_count}
\title{Dose-response models with rstan for binary responses}
\description{Formula interface for fitting dose-response curves with rstan for binary responses}
\usage{
standrm_count(formula, data, fct, curveid = NULL, random = NULL,
               priors = standrc_priors(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A two-sided formula with counts as the response (integer >= 0) and dose covariate names}
  \item{data}{A dataframe containing the variables in the model}
  \item{fct}{A drc model function.}
  \item{curveid}{A two-sided formula with parameter names (b + c + d) on the left and a factor name on the right hand side, defining several clusters, modelled as fixed effects}
  \item{random}{A two-sided formula with parameter names (b + c + d) on the left and a factor name on the right hand side, defining several clusters, modelled as random effects}
  \item{priors}{A list with prior definitions, see \code{\link{standrc_priors}}}
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

\details{The counts are assumed to be Poisson distributed. The parameterization of the drc model functions is changed, modelling the slope and the asymmetry parameter on a logarithmic scale and allowing the lower and upper asymptotes to swap parameters. When not fixed, the asymptotes are modelled on the log scale assuming normal priors, also adding random effects on the log scale; samples from the posterior of the asymptotes and predictions are obtained on the original scale of the counts. There might be a problem with three-parameter functions, fixing a specific asymptote by the fixed argument in a more complex model might be a better solution. Random effect parameters can be specified even when the fixed effects counterpart is fixed at a specific value, but might result in convergence problems with small lower asymptotes as asymptote parameters are fixed on the original scale.}


\keyword{models}
