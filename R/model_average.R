maED <- function(object, ..., respLev=NULL){ 
  if (is.null(respLev)) stop("Please provide response levels!")
  
  modellist <- list(object, ...)
  inc <- sapply(modellist, function(x) waic(x)$waic)
  d <- inc-min(inc)
  wi <- exp(-0.5*d)/sum(exp(-0.5*d))

  isstandrc <- sapply(modellist, function(x) inherits(x, "standrc"))
  mllist <- modellist[isstandrc]  
  Call <- match.call()
  Call$respLev <- NULL
  mnames <- as.character(Call[-1L])[isstandrc]
  names(wi) <- mnames

  nl <- length(respLev)  
  edl <- lapply(1:length(wi), function(i) lapply(ED(mllist[[i]], respLev=respLev), function(x) x*wi[i]))
  EDlist <- edl[[1]]
  for (i in 2:length(mllist)){
    EDlist <- lapply(1:length(EDlist), function(j) apply(EDlist[[j]], 2, sort) + apply(edl[[i]][[j]], 2, sort))
  }  
  names(EDlist) <- respLev
  class(EDlist) <- "EDsamp"
  return(EDlist)  
}


mapredict <- function(object, ..., newdata=NULL){
  modellist <- list(object, ...)
  inc <- sapply(modellist, function(x) waic(x)$waic)
  d <- inc-min(inc)
  wi <- exp(-0.5*d)/sum(exp(-0.5*d))
  
  isstandrc <- sapply(modellist, function(x) inherits(x, "standrc"))
  mllist <- modellist[isstandrc]  
  Call <- match.call()
  Call$newdata <- NULL
  mnames <- as.character(Call[-1L])[isstandrc]
  names(wi) <- mnames
  
  pred <- wi[1]*predict(mllist[[1]], newdata=newdata)
  for (i in 2:length(modellist)){
    pred <- apply(pred, 2, sort) + apply(wi[i]*predict(mllist[[i]], newdata=newdata), 2, sort)
  }
  return(pred)
}

maplot <- function(x, ..., ndose=25, logx=FALSE, lim=NULL){  
  if (is.null(x$data$total)) xcc <- as.character(x$call$formula[[2]]) else xcc <- as.character("p")
  if (x$curves$J > 1){
    dframe <- data.frame(x$data$y, x$data$x, x$curves$names[x$data$idc])
    names(dframe) <- c(xcc,
                       as.character(x$call$formula[[3]]), 
                       as.character(x$call$curveid[[3]]))
  } else {
    dframe <- data.frame(x$data$y, x$data$x)
    names(dframe) <- c(xcc,
                       as.character(x$call$formula[[3]]))
  }
  if (!is.null(x$data$total)) dframe[,1] <- dframe[,1]/x$data$total
  
  if (x$curves$J > 1){
    curvn <- paste(",", as.character(x$call$curveid[[3]]), "=x$curves$names")
    dframe[,3] <- factor(dframe[,3], levels=x$curves$names)
  } else {
    curvn=NULL
  }
  if (is.null(lim)){
    minx <- min(x$data$x)
    maxx <- max(x$data$x)
  } else {
    if (length(lim) != 2) stop("Please provide limits as vector with 2 elements.")
    minx <- lim[1]
    maxx <- lim[2]
  }
  if (logx){    
    newd <- eval(parse(text=paste("expand.grid(", as.character(x$call$formula[[3]]), "=exp(seq(log(", minx, "), log(", maxx, "), length=ndose))", curvn, ")")))
  } else {
    newd <- eval(parse(text=paste("expand.grid(", as.character(x$call$formula[[3]]), "=seq(", minx, ", ", maxx, ", length=ndose)", curvn, ")")))    
  }
  pm <- mapredict(x, ..., newdata=newd)
  newd$p <- apply(pm, 2, function(x) mean(x, na.rm=TRUE))
  newd$pmin <- apply(pm, 2, function(x) quantile(x, 0.025, na.rm=TRUE))
  newd$pmax <- apply(pm, 2, function(x) quantile(x, 0.975, na.rm=TRUE))
  
  if (logx) lxt <- "+ coord_trans(x='log')" else lxt <- NULL
  if (x$curves$J > 1){
    eval(parse(text=paste("ggplot(dframe, aes(x=",  as.character(x$call$formula[[3]]),", y=", xcc, ", colour=", as.character(x$call$curveid[[3]]),")) +  geom_point() + geom_ribbon(data=newd, aes(y=p, ymin=pmin, ymax=pmax, fill=", as.character(x$call$curveid[[3]]),", colour=NULL), alpha=0.2) + geom_line(data=newd, aes(y=p))", lxt)) ) 
  } else {
    eval(parse(text=paste("ggplot(dframe, aes(x=",  as.character(x$call$formula[[3]]),", y=", xcc, ")) +  geom_point() + geom_ribbon(data=newd, aes(y=p, ymin=pmin, ymax=pmax), alpha=0.2) + geom_line(data=newd, aes(y=p))", lxt)) ) 
  }
}
