ED.standrc <- function(object, respLev=NULL){ 
  if (is.null(respLev) & !is.null(object$respLev)){
    print(object$stan, pars="dED") 
  } else {
    if (is.null(respLev)) stop("Please provide response levels!")
    plist <- extract(object$stan, pars=object$pars)
    fix <- object$fixed
    samp <- lapply(1:object$curves$J, function(i){
      sapply(1:length(fix), function(j) if (is.na(fix[j])){
        if (object$curves$pars[sum(is.na(fix[1:j]))] == TRUE) plist[[sum(is.na(fix[1:j]))]][,i] else plist[[sum(is.na(fix[1:j]))]]
      } else {rep(fix[j], dim(plist[[1]])[1])})
    })
    EDlist <- lapply(respLev, function(p){
      smat <- cbind(sapply(1:length(samp), function(i){
        apply(samp[[i]], 1, function(x){
          xt <- x
          xt[1] <- -exp(xt[1])
          xt[5] <- exp(xt[5])
          p <- 100-p
          tempVal <- log((100 - p)/100)
          xt[4] * (exp(-tempVal/xt[5]) - 1)^(1/xt[1])
        })
      }))
      colnames(smat) <- object$curves$names
      return(smat)
    })
    names(EDlist) <- respLev
    class(EDlist) <- "EDsamp"
    return(EDlist)
  }
}

print.EDsamp <- function(x, ...){
  out <- lapply(x, function(xl){
    dat <- as.data.frame(t(apply(cbind(xl), 2, function(xc) quantile(xc, c(0.5, 0.025, 0.975)))))
    colnames(dat) <- c("median", "2.5%", "97.5%")
    return(dat)
  })
  print(out)
}