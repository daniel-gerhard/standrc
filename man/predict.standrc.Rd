\name{predict.standrc}
\alias{predict.standrc}
\alias{mapredict}
\title{Predictions from standrc objects}
\description{
Predict new response values from a standrc object
}
\usage{
\method{predict}{standrc}(object, ..., newdata=NULL)
mapredict(object, ..., newdata=NULL)
}
\arguments{
  \item{object}{A standrc object}
  \item{...}{further arguments}
  \item{newdata}{data.frame with dose levels and curveids as support points, names should be the same as in the original data.frame; if NULL the predictions are obtained at the observed dose levels.}
}


\author{
Daniel Gerhard
}

\seealso{
\code{\link{standrm}}
}

\keyword{methods}