standrm <- function(formula, data, fct, curveid=NULL, random=NULL, respLev=NULL, ...){
  callDetail <- match.call()
  
  mf <- model.frame(formula, data)
  
  if (!is.null(curveid)){
    cid <- as.character(curveid)[3]
    idc <- as.numeric(data[,cid])
    J <- length(unique(idc))
    stsp <- strsplit(as.character(curveid)[2], "+")[[1]]
    cin <- stsp[!stsp %in% c(" ", "+")]
    pnl <- sapply(fct$names, function(x) x %in% cin)
    curvenames <- levels(as.factor(data[,cid]))
  } else {
    cid <- idc <- NULL
    J <- 1
    pnl <- rep(FALSE, length(fct$names))
    curvenames <- NULL
  }
  
  if (!is.null(random)){
    rid <- as.character(random)[3]
    idr <- as.numeric(data[,rid])
    K <- length(unique(idr))
    stsp <- strsplit(as.character(random)[2], "+")[[1]]
    rin <- stsp[!stsp %in% c(" ", "+")]
    pnlr <- sapply(c("b", "c", "d", "e", "f"), function(x) x %in% rin)
  } else {
    rid <- idr <- NULL
    K <- 1
    pnlr <- rep(FALSE, 5)
  }
 
  fix <- fct$fixed
  isfix <- !is.na(fix)
  
  jv <- rep(J, 5)
  jv[!isfix][!pnl] <- 1
  
  N <- nrow(mf) 
  y <- mf[,1]
  x <- mf[,2] 
  if (fct$name == "LL.5") x[x == 0] <- 0.5*min(x[x > 0])
  
  pb <- rep(0, jv[1])  
  if (mean(y[x == min(x)]) < mean(y[x == max(x)])){
    pc <- rep(min(y), jv[2])
    pd <- rep(max(y), jv[3])
  } else {
    pd <- rep(min(y), jv[3])
    pc <- rep(max(y), jv[2])
  }  
  pe <- rep(median(x), jv[4])
  pf <- rep(0, jv[5])
  
  ped <- respLev/100
  Ned <- length(ped)
  
  stan_dat <- list(N=N, J=J, K=K, idc=idc, idr=idr, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf, Ned=Ned, ped=ped)
  if (is.null(curveid)) stan_dat <- list(N=N, K=K, idr=idr, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf, Ned=Ned, ped=ped)
  if (is.null(random)) stan_dat <- list(N=N, J=J, idc=idc, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf, Ned=Ned, ped=ped)
  if (is.null(curveid) & is.null(random)) stan_dat <- list(N=N, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf, Ned=Ned, ped=ped)
  if (is.null(respLev)) stan_dat <- list(N=N, J=J, K=K, idc=idc, idr=idr, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf)
  if (is.null(curveid) & is.null(respLev)) stan_dat <- list(N=N, K=K, idr=idr, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf)
  if (is.null(random) & is.null(respLev)) stan_dat <- list(N=N, J=J, idc=idc, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf)
  if (is.null(curveid) & is.null(random) & is.null(respLev)) stan_dat <- list(N=N, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf)
  
  # data
  d1 <- "int<lower=0> N; real y[N];  real<lower=0> x[N];"  
  dJ <- "int<lower=0> idc[N]; int<lower=0> J;" 
  dK <- "int<lower=0> idr[N]; int<lower=0> K;" 
  ded <- if (Ned <= 1) "real<lower=0, upper=100> ped;" else "int<lower=0> Ned; real<lower=0, upper=100> ped[Ned];"
  
  # priors
  dJp <- c("real pb[J];", "real pc[J];", "real pd[J];", "real pe[J];", "real pf[J];")[!isfix]
  dp <- c("real pb;", "real pc;", "real pd;", "real pe;", "real pf;")[!isfix]  
  dpc <- paste(sapply(1:length(pnl), function(i) if (pnl[i]) dJp[i] else dp[i]), collapse=" ")

  # define parameters
  paraJ <- c("real slope[J];", "real lasy[J];", "real uasy[J];", "real<lower=0> ed[J];", "real assym[J];")[!isfix]
  para <- c("real slope;", "real lasy;", "real uasy;", "real<lower=0> ed;", "real assym;")[!isfix]
  parac <- paste(sapply(1:length(pnl), function(i) if (pnl[i]) paraJ[i] else para[i]), collapse=" ")
  
  # random effects
  rpara <- c("real rslope[K];", "real rlasy[K];", "real ruasy[K];", "real red[K];", "real rassym[K];")
  rparasig <- c("real<lower=0> sigmasq_slope;", "real<lower=0> sigmasq_lasy;", "real<lower=0> sigmasq_uasy;", "real<lower=0> sigmasq_ed;", "real<lower=0> sigmasq_assym;")
  rparac <- paste(c(rpara[pnlr], rparasig[pnlr]), collapse=" ")
   
  trans1 <- "real<lower=0> sigma_y; real mu[N];" 

  trap <- c("slope", "lasy", "uasy", "ed", "assym")[!isfix]
  traj <- c("slope[idc[i]]", "lasy[idc[i]]", "uasy[idc[i]]", "ed[idc[i]]", "assym[idc[i]]")[!isfix]
  tra <- character(length(fix))
  tra[!isfix] <- sapply(1:length(pnl), function(i) if (pnl[i]) traj[i] else trap[i])
  tra[isfix] <- fix[isfix]
  if (is.null(random)){
    trans <- paste("sigma_y <- sqrt(sigmasq_y); for(i in 1:N){ mu[i] <-", tra[2], " + (", tra[3], "-", tra[2], ") / (1 + exp(-exp(", tra[1], ") * (log(x[i]/", tra[4], "))))^exp(", tra[5], ");}", collapse="")
  } else {
    stra <- paste(c("real<lower=0> sigma_slope;", "real<lower=0> sigma_lasy;", "real<lower=0> sigma_uasy;", "real<lower=0> sigma_ed;", "real<lower=0> sigma_assym;")[pnlr], collapse=" ")
    strasq <- paste(c("sigma_slope <- sqrt(sigmasq_slope);",
                      "sigma_lasy <- sqrt(sigmasq_lasy);",
                      "sigma_uasy <- sqrt(sigmasq_uasy);",
                      "sigma_ed <- sqrt(sigmasq_ed);",
                      "sigma_assym <- sqrt(sigmasq_assym);")[pnlr], collapse=" ")
    trar <- c(paste("(", tra[1], " + rslope[idr[i]] )"), 
              paste("(", tra[2], " + rlasy[idr[i]] )"), 
              paste("(", tra[3], " + ruasy[idr[i]] )"), 
              paste("(", tra[4], " + red[idr[i]] )"), 
              paste("(", tra[5], " + rassym[idr[i]] )"))
    trc <- sapply(1:length(tra), function(i) if (pnlr[i]) trar[i] else tra[i])    
    trans <- paste(stra, strasq, "sigma_y <- sqrt(sigmasq_y); for(i in 1:N){ mu[i] <-", trc[2], " + (", trc[3], "-", trc[2], ") / (1 + exp(-exp(", trc[1], ") * (log(x[i]/", trc[4], "))))^exp(", trc[5], ");}", collapse="")
  }
  
  # derived ED
  rped <- if (Ned <= 1){
    if (is.null(idc)) c("real dED;") else c("real dED[J];")
  } else {
    if (is.null(idc)) c("real dED[Ned];") else c("real dED[Ned, J];")
  }
  edtj <- c("slope[j]", "lasy[j]", "uasy[j]", "ed[j]", "assym[j]")[!isfix]
  edt <- character(length(fix))
  edt[!isfix] <- sapply(1:length(pnl), function(i) if (pnl[i]) edtj[i] else trap[i])
  edt[isfix] <- fix[isfix]
  if (is.null(random)){ 
    if (is.null(idc)){
      if (Ned <= 1){  
        ted <- paste("dED <-", edt[4], " * (exp(-1*log( ped )/ exp(", edt[5], ")) - 1)^(1/(-1*exp(", edt[1], ")));", collapse="")
      } else {
        ted <- paste("for (ne in 1:Ned){ dED[ne] <-", edt[4], " * (exp(-1*log( ped[ne] )/ exp(", edt[5], ")) - 1)^(1/(-1*exp(", edt[1], ")));}", collapse="")
      }
    } else {
      if (Ned <= 1){  
        ted <- paste("for (j in 1:J){ dED[j] <-", edt[4], " * (exp(-1*log( ped )/ exp(", edt[5], ")) - 1)^(1/(-1*exp(", edt[1], ")));}", collapse="")
      } else {
        ted <- paste("for (j in 1:J){ for (ne in 1:Ned){ dED[ne, j] <-", edt[4], " * (exp(-1*log( ped[ne] )/ exp(", edt[5], ")) - 1)^(1/(-1*exp(", edt[1], ")));}}", collapse="")
      }
    }
  } else {  
    # random effects
    rnpara <- c("real rnslope;", "real rnlasy;", "real rnuasy;", "real rned;", "real rnassym;")
    rnparac <- paste(rnpara[pnlr], collapse=" ")
    trnar <- c(paste(edt[1], " + rnslope"), 
              paste("(", edt[2], " + rnlasy )"), 
              paste("(", edt[3], " + rnuasy )"), 
              paste("(", edt[4], " + rned )"), 
              paste(edt[5], " + rnassym"))
    trnc <- sapply(1:length(tra), function(i) if (pnlr[i]) trnar[i] else edt[i]) 
    if (is.null(idc)){
      if (Ned <= 1){  
        ted <- paste("dED <-", trnc[4], " * (exp(-1*log( ped )/ exp(", trnc[5], ")) - 1)^(1/(-1*exp(", trnc[1], ")));", collapse="")
      } else {
        ted <- paste("for (ne in 1:Ned){ dED[ne] <-", trnc[4], " * (exp(-1*log( ped[ne] )/ exp(", trnc[5], ")) - 1)^(1/(-1*exp(", trnc[1], ")));}", collapse="")
      }
    } else {   
      if (Ned <= 1){  
        ted <- paste("for (j in 1:J){ dED[j] <-", trnc[4], " * (exp(-1*log( ped )/ exp(", trnc[5], ")) - 1)^(1/(-1*exp(", trnc[1], ")));}", collapse="")
      } else {
        ted <- paste("for (j in 1:J){ for (ne in 1:Ned){ dED[ne, j] <-", trnc[4], " * (exp(-1*log( ped[ne] )/ exp(", trnc[5], ")) - 1)^(1/(-1*exp(", trnc[1], ")));}}", collapse="")
      }
    }
  }
  
  mod <- paste(c("slope ~ normal(pb, 100);", "lasy ~ normal(pc, 100);", "uasy ~ normal(pd, 100);", "ed ~ normal(pe, 100);", "assym ~ normal(pf, 100);")[!isfix], collapse=" ") 
  
  mody <- "sigmasq_y ~ inv_gamma(0.001, 0.001); y ~ normal(mu, sigma_y);"
  
  modr <- paste(c("rslope ~ normal(0, sigma_slope);", "rlasy ~ normal(0, sigma_lasy);", "ruasy ~ normal(0, sigma_dasy);", "red ~ normal(0, sigma_ed);", "rassym ~ normal(0, sigma_assym);")[pnlr], collapse=" ") 
  modrsig <- paste(c("sigmasq_slope ~ inv_gamma(0.001, 0.001);", "sigmasq_lasy ~ inv_gamma(0.001, 0.001);", "sigmasq_uasy ~ inv_gamma(0.001, 0.001);", "sigmasq_ed ~ inv_gamma(0.001, 0.001);", "sigmasq_assym ~ inv_gamma(0.001, 0.001);")[pnlr], collapse=" ") 
  
  # new random effects for derived ed
  moded <- paste(c("rnslope ~ normal(0, sigma_slope);", "rnlasy ~ normal(0, sigma_lasy);", "rnuasy ~ normal(0, sigma_dasy);", "rned ~ normal(0, sigma_ed);", "rnassym ~ normal(0, sigma_assym);")[pnlr], collapse=" ") 
  
  stancode <- paste("data {", 
                    d1,   
                    if (!is.null(idc)) dJ,
                    if (!is.null(random)) dK,
                    if (!is.null(respLev)) ded,
                    dpc,
                    "} ",
                    "parameters { real<lower=0> sigmasq_y;",
                    parac,
                    if (!is.null(random)) rparac, 
                    if (!is.null(random) & !is.null(respLev)) rnparac,  
                    "} ",
                    "transformed parameters {",
                    trans1,                    
                    trans,                      
                    "} ",
                    "model {",
                    mod,
                    mody,                     
                    if (!is.null(random)) modr,
                    if (!is.null(random)) modrsig,    
                    if (!is.null(random) & !is.null(respLev)) moded,                    
                    "}",
                    if (!is.null(respLev)) "generated quantities {",
                    if (!is.null(respLev)) rped,
                    if (!is.null(respLev)) ted,                                      
                    if (!is.null(respLev)) "}",
                    sep="")
  
  assign("stancode", stancode, envir=.GlobalEnv)
  fit <- stan(model_code = 'stancode', data = stan_dat, ...)
  
  out <- list()
  out$data <- stan_dat
  out$model <- stancode
  out$stan <- fit
  out$pars <- c(trap, "sigma_y")
  if (!is.null(random)) out$random <- c("slope", "lasy", "uasy", "ed", "assym")[pnlr] 
  out$curves <- list(pars=pnl, J=J, names=curvenames)
  out$fixed <- fix
  out$fct <- fct
  out$respLev <- respLev
  class(out) <- "standrc"
  return(out)
}

