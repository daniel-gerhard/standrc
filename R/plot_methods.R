plot.standrc <- function(x, ..., ndose=25, logx=FALSE){
  if (x$curves$J > 1){
    dframe <- data.frame(x$data$y, x$data$x, x$curves$names[x$data$idc])
    names(dframe) <- c(as.character(x$call$formula[[2]]),
                       as.character(x$call$formula[[3]]), 
                       as.character(x$call$curveid[[3]]))
  } else {
    dframe <- data.frame(x$data$y, x$data$x)
    names(dframe) <- c(as.character(x$call$formula[[2]]),
                       as.character(x$call$formula[[3]]))
  }
  
  if (x$curves$J > 1) curvn <- paste(",", as.character(x$call$curveid[[3]]), "=x$curves$names") else curvn=NULL
  if (logx){    
    newd <- eval(parse(text=paste("expand.grid(", as.character(x$call$formula[[3]]), "=exp(seq(log(min(x$data$x)), log(max(x$data$x)), length=ndose))", curvn, ")")))
  } else {
    newd <- eval(parse(text=paste("expand.grid(", as.character(x$call$formula[[3]]), "=exp(seq(log(min(x$data$x)), log(max(x$data$x)), length=ndose))", curvn, ")")))    
  }
  pm <- predict(x, newdata=newd)
  newd$p <- apply(pm, 2, function(x) mean(x, na.rm=TRUE))
  newd$pmin <- apply(pm, 2, function(x) quantile(x, 0.025, na.rm=TRUE))
  newd$pmax <- apply(pm, 2, function(x) quantile(x, 0.975, na.rm=TRUE))
  
  if (logx) lxt <- "+ coord_trans(x='log')" else lxt <- NULL
  if (x$curves$J > 1){
    eval(parse(text=paste("ggplot(dframe, aes(x=",  as.character(x$call$formula[[3]]),", y=", as.character(x$call$formula[[2]]), ", colour=", as.character(x$call$curveid[[3]]),")) +  geom_point() + geom_ribbon(data=newd, aes(y=p, ymin=pmin, ymax=pmax, fill=", as.character(x$call$curveid[[3]]),", colour=NULL), alpha=0.2) + geom_line(data=newd, aes(y=p))", lxt)) ) 
  } else {
    eval(parse(text=paste("ggplot(dframe, aes(x=",  as.character(x$call$formula[[3]]),", y=", as.character(x$call$formula[[2]]), ")) +  geom_point() + geom_ribbon(data=newd, aes(y=p, ymin=pmin, ymax=pmax), alpha=0.2) + geom_line(data=newd, aes(y=p))", lxt)) ) 
  }
}