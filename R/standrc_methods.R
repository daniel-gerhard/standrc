print.standrc <- function(x, ...){
  print(x$stan, pars=x$pars) 
}

fixef.standrc <- function(object, ...){
  print(object$stan, pars=object$fixedpars)
}

ranef.standrc <- function(object, ...){
  print(object$stan, pars=paste("r", object$random, sep=""))
}

VarCorr.standrc <- function(x, sigma=1, rdig=3){
  print(x$stan, pars=c(paste("sigma_", x$random, sep=""), "sigma_y"))
}

predict.standrc <- function(object, ..., newdata=NULL){ 
  if (is.null(newdata)) x <- object$data$x else x <- newdata[,as.character(object$call$formula[3])]
  if (is.null(newdata)){
    if (is.null(object$data$idc)) idc <- rep(1, length(x)) else idc <- object$data$idc
  } else {
    if (is.null(object$data$idc)) idc <- rep(1, nrow(newdata)) else idc <- newdata[,as.character(object$call$curveid)[3]]
  }  
  plist <- extract(object$stan, pars=object$pars)
  fix <- object$fixed
  fid <- logical(length=length(fix))
  fid[is.na(fix)] <- object$curves$pars
  samp <- lapply(1:object$curves$J, function(i){
    sapply(1:length(fix), function(j){
      if (fid[j]) plist[[j]][,i] else plist[[j]]
    })
  })
  pred <- sapply(1:length(x), function(i){
    apply(samp[[as.numeric(idc)[i]]], 1, function(xp){
      xp[2] + (xp[3] - xp[2]) / (1 + exp(-exp(xp[1]) * (log(x[i]/ xp[4]))))^exp(xp[5])
    })
  })
  return(pred)
}
