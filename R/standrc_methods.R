print.standrc <- function(x, ...){
  print(x$stan, pars=x$pars) 
}

ranef.standrc <- function(object, ...){
  print(object$stan, pars=paste("r", object$random, sep=""))
}

VarCorr.standrc <- function(x, sigma=1, rdig=3){
  print(x$stan, pars=c(paste("sigma_", x$random, sep=""), "sigma_y"))
}