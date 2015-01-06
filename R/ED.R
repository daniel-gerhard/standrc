ED.standrc <- function(object, respLev=NULL){ 
  if (is.null(respLev)) stop("Please provide response levels!")
  plist <- extract(object$stan, pars=object$pars)
  fix <- object$fixed
  fid <- logical(length=length(fix))
  fid[is.na(fix)] <- object$curves$pars
  samp <- lapply(1:object$curves$J, function(i){
    sapply(1:length(fix), function(j){
      if (fid[j]) plist[[j]][,i] else plist[[j]]
    })
  })
  EDlist <- lapply(respLev, function(p){
    smat <- cbind(sapply(1:length(samp), function(i){
      apply(samp[[i]], 1, function(x){
        xt <- x
        xt[1] <- -exp(xt[1])
        xt[5] <- exp(xt[5])
        p <- 100-p
        if (object$fct$name %in% c("LL.5", "LL.4", "LL.3")){
          tempVal <- log((100 - p)/100)
          value <- xt[4] * (exp(-tempVal/xt[5]) - 1)^(1/xt[1])
        }
        if (object$fct$name %in% c("L.5", "L.4", "L.3")){
          tempVal <- 100/p
          value <- xt[4] + (log(tempVal^(1/xt[5]) - 1))/xt[1]
        }
        if (object$fct$name %in% c("W1.4", "W1.3")){
          tempVal <- log(-log((100 - p)/100))
          value <- exp(tempVal/xt[1] + log(xt[4]))
        }
        return(value)
      })
    }))
    colnames(smat) <- object$curves$names
    return(smat)
  })
  names(EDlist) <- respLev
  class(EDlist) <- "EDsamp"
  return(EDlist)  
}

print.EDsamp <- function(x, ...){
  out <- lapply(x, function(xl){
    dat <- as.data.frame(t(apply(cbind(xl), 2, function(xc) quantile(xc, c(0.5, 0.025, 0.975)))))
    colnames(dat) <- c("median", "2.5%", "97.5%")
    return(dat)
  })
  print(out)
}