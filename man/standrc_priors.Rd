\name{standrc_priors}
\alias{standrc_priors}
\title{Priors for standrc model parameters}
\description{
  Definition of prior distributions for standrc model parameters
}
\usage{
  standrc_priors(b="normal(pb, 100)",
                 c="normal(pc, 100)",
                 d="normal(pd, 100)",
                 e="normal(pe, 100)",
                 f="normal(pf, 100)",
                 sy="inv_gamma(0.001, 0.001)",
                 sb="inv_gamma(0.001, 0.001)",
                 sc="inv_gamma(0.001, 0.001)",
                 sd="inv_gamma(0.001, 0.001)",
                 se="inv_gamma(0.001, 0.001)",
                 sf="inv_gamma(0.001, 0.001)",
                 pb=NULL,
                 pc=NULL,
                 pd=NULL,
                 pe=NULL,
                 pf=NULL)
}
\arguments{
  \item{b}{prior for the slope parameter (on log scale)}
  \item{c}{prior for the asymptote parameter}
  \item{d}{prior for the 2nd asymptote parameter}
  \item{e}{prior for the effective dose parameter}
  \item{f}{prior for the asymmetry parameter (on log scale)}
  \item{sy}{prior for the residual standard deviation}
  \item{sb}{prior for the slope random effect standard deviation}
  \item{sc}{prior for the asymptote random effect standard deviation}
  \item{sd}{prior for the 2nd asymptote random effect standard deviation}
  \item{se}{prior for the effective doserandom effect standard deviation}
  \item{sf}{prior for the asymmetry random effect standard deviation}
  \item{pb}{prior mean vector for the slope parameter (if NULL, it is drived from the data)}
  \item{pc}{prior mean vector for the asymptote parameter (if NULL, it is drived from the data)}
  \item{pd}{prior mean vector for the 2nd asymptote parameter (if NULL, it is drived from the data)}
  \item{pe}{prior mean vector for the effective dose parameter (if NULL, it is drived from the data)}
  \item{pf}{prior mean vector for the asymmetry parameter (if NULL, it is drived from the data)}
}

\details{
The priors are defined by character strings, containing stan code with the distributional assumptions for each parameter. pb, pc, pd, pe, pf contain the automatically derived prior means and can be replaced by supplying vectors for the arguments pb, ..., pf.
}


\author{
  Daniel Gerhard
}

\seealso{
  \code{\link{standrm}}
}

\keyword{htest}