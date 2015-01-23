standrm_binary <- function(formula, data, fct, curveid=NULL, random=NULL, priors=standrc_priors(), ...){
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
  if (fct$name %in% c("W1.4", "W1.3", "W2.4", "W2.3", "LN.4", "LN.3")) fix <- c(fix, 0)  
  isfix <- !is.na(fix)
  
  jv <- rep(J, 5)
  jv[!isfix][!pnl] <- 1
  
  N <- nrow(mf) 
  y <- as.integer(mf[,1])
  x <- mf[,2]   
  if (fct$name %in% c("LL.5", "LL.4", "LL.3", "W1.4", "W1.3", "W2.4", "W2.3", "LN.4", "LN.3")) x[x == 0] <- 0.5*min(x[x > 0])
  if (fct$name %in% c("LL.4", "LL.3", "L.4", "L.3")) fix[5] <- 0
  
  if (is.null(priors$pb)) pb <- rep(0, jv[1])  
  if (mean(y[x == min(x)]) < mean(y[x == max(x)])){
    if (is.null(priors$pc)) pc <- rep(-10, jv[2])
    if (is.null(priors$pd)) pd <- rep(10, jv[3])
  } else {
    if (is.null(priors$pd)) pd <- rep(-10, jv[3])
    if (is.null(priors$pc)) pc <- rep(10, jv[2])
  }  
  if (is.null(priors$pe)) pe <- rep(median(x), jv[4])
  if (is.null(priors$pf)) pf <- rep(0, jv[5])
  
  stan_dat <- list(N=N, J=J, K=K, idc=idc, idr=idr, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf)
  if (is.null(curveid)) stan_dat <- list(N=N, K=K, idr=idr, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf)
  if (is.null(random)) stan_dat <- list(N=N, J=J, idc=idc, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf)
  if (is.null(curveid) & is.null(random)) stan_dat <- list(N=N, y=y, x=x, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf)
  
  # data 
  dJ <- "int<lower=0> idc[N]; int<lower=0> J;" 
  dK <- "int<lower=0> idr[N]; int<lower=0> K;" 
  
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
  
  trap <- c("slope", "1/(1 + exp(-1*lasy))", "1/(1 + exp(-1*uasy))", "ed", "assym")[!isfix]
  traj <- c("slope[idc[i]]", "1/(1 + exp(-1*lasy[idc[i]]))", "1/(1 + exp(-1*uasy[idc[i]]))", "ed[idc[i]]", "assym[idc[i]]")[!isfix]
  tra <- character(length(fix))
  tra[!isfix] <- sapply(1:length(pnl), function(i) if (pnl[i]) traj[i] else trap[i])
  tra[isfix] <- fix[isfix]
  if (is.null(random)){
    if (fct$name %in% c("LL.5", "LL.4", "LL.3")) trans <- paste("for(i in 1:N){ mu[i] <-", tra[2], " + (", tra[3], "-", tra[2], ") / (1 + exp(-exp(", tra[1], ") * (log(x[i]/", tra[4], "))))^exp(", tra[5], ");}", collapse="")
    if (fct$name %in% c("L.5", "L.4", "L.3")) trans <- paste("for(i in 1:N){ mu[i] <-", tra[2], " + (", tra[3], "-", tra[2], ") / (1 + exp(-exp(", tra[1], ") * (x[i] - ", tra[4], ")))^exp(", tra[5], ");}", collapse="")
    if (fct$name %in% c("W1.4", "W1.3")) trans <- paste("for(i in 1:N){ mu[i] <-", tra[2], " + (", tra[3], "-", tra[2], ") * exp(-exp(-exp(", tra[1], ") * (log(x[i]) - log(", tra[4], "))));}", collapse="")
    if (fct$name %in% c("W2.4", "W2.3")) trans <- paste("for(i in 1:N){ mu[i] <-", tra[2], " + (", tra[3], "-", tra[2], ") * (1 - exp(-exp(-exp(", tra[1], ") * (log(x[i]) - log(", tra[4], ")))));}", collapse="")
    if (fct$name %in% c("LN.4", "LN.3")) trans <- paste("for(i in 1:N){ mu[i] <-", tra[2], " + (", tra[3], "-", tra[2], ") * normal_cdf(-exp(", tra[1], ") * (log(x[i]) - log(", tra[4], ")), 0, 1);}", collapse="")
  } else {
    stra <- paste(c("real<lower=0> sigma_slope;", "real<lower=0> sigma_lasy;", "real<lower=0> sigma_uasy;", "real<lower=0> sigma_ed;", "real<lower=0> sigma_assym;")[pnlr], collapse=" ")
    strasq <- paste(c("sigma_slope <- sqrt(sigmasq_slope);",
                      "sigma_lasy <- sqrt(sigmasq_lasy);",
                      "sigma_uasy <- sqrt(sigmasq_uasy);",
                      "sigma_ed <- sqrt(sigmasq_ed);",
                      "sigma_assym <- sqrt(sigmasq_assym);")[pnlr], collapse=" ")
    trar <- c(paste("(", tra[1], " + rslope[idr[i]] )"), 
              paste("(1/(1 + exp(-1*(", tra[2], " + rlasy[idr[i]] )))"), 
              paste("(1/(1 + exp(-1*(", tra[3], " + ruasy[idr[i]] )))"), 
              paste("(", tra[4], " + red[idr[i]] )"), 
              paste("(", tra[5], " + rassym[idr[i]] )"))
    trc <- sapply(1:length(tra), function(i) if (pnlr[i]) trar[i] else tra[i])    
    if (fct$name %in% c("LL.5", "LL.4", "LL.3")) trans <- paste(stra, strasq, "for(i in 1:N){ mu[i] <-", trc[2], " + (", trc[3], "-", trc[2], ") / (1 + exp(-exp(", trc[1], ") * (log(x[i]/", trc[4], "))))^exp(", trc[5], ");}", collapse="")
    if (fct$name %in% c("L.5", "L.4", "L.3")) trans <- paste(stra, strasq, "for(i in 1:N){ mu[i] <-", trc[2], " + (", trc[3], "-", trc[2], ") / (1 + exp(-exp(", trc[1], ") * (x[i] -", trc[4], ")))^exp(", trc[5], ");}", collapse="")
    if (fct$name %in% c("W1.4", "W1.3")) trans <- paste(stra, strasq, "for(i in 1:N){ mu[i] <-", trc[2], " + (", trc[3], "-", trc[2], ") * exp(-exp(-exp(", trc[1], ") * (log(x[i]) - log(", trc[4], "))));}", collapse="")
    if (fct$name %in% c("W2.4", "W2.3")) trans <- paste(stra, strasq, "for(i in 1:N){ mu[i] <-", trc[2], " + (", trc[3], "-", trc[2], ") * (1 - exp(-exp(-exp(", trc[1], ") * (log(x[i]) - log(", trc[4], ")))));}", collapse="")
    if (fct$name %in% c("LN.4", "LN.3")) trans <- paste(stra, strasq, "for(i in 1:N){ mu[i] <-", trc[2], " + (", trc[3], "-", trc[2], ") * normal_cdf(-exp(", trc[1], ") * (log(x[i]) - log(", trc[4], ")), 0, 1);}", collapse="")
  }
  
  ### population parameters
  popJ <- c("real pslope[J];", "real plasy[J];", "real puasy[J];", "real ped[J];", "real passym[J];")
  pop <- c("real pslope;", "real plasy;", "real puasy;", "real ped;", "real passym;")
  pop[!isfix] <- sapply(1:length(pnl), function(i) if (pnl[i]) popJ[!isfix][i] else pop[!isfix][i]) 
  popc <- paste(pop, collapse=" ")
  rnpara <- c("real rnslope;", "real rnlasy;", "real rnuasy;", "real rned;", "real rnassym;")
  rnparac <- paste(rnpara[pnlr], collapse=" ")
  ppp <- c("slope", "lasy", "uasy", "ed", "assym")[!isfix]
  ppj <- c("slope[j]", "lasy[j]", "uasy[j]", "ed[j]", "assym[j]")[!isfix]
  pp <- character(length(fix))
  pp[!isfix] <- sapply(1:length(pnl), function(i) if (pnl[i]) ppj[i] else ppp[i])
  pp[isfix] <- fix[isfix]
  po <- c("pslope", "plasy", "puasy", "ped", "passym")
  poj <- paste("for (j in 1:J)", c("pslope[j]", "plasy[j]", "puasy[j]", "ped[j]", "passym[j]"))
  pof <- character(length(fix))
  pof[!isfix] <- sapply(1:length(pnl), function(i) if (pnl[i]) poj[!isfix][i] else po[!isfix][i])
  pof[isfix] <- po[isfix]
  ppc <- paste(pof[1], " <- ", pp[1], if (pnlr[1]) " + rnslope", ";", pof[2], " <- ", if (!isfix[2]) "1/(1 + exp(-1*(", pp[2], if (pnlr[2]) " + rnlasy", if (!isfix[2]) ")))", ";", pof[3], " <- ", if (!isfix[3]) "1/(1 + exp(-1*(", pp[3], if (pnlr[3]) " + rnuasy", if (!isfix[3]) ")))", ";", pof[4], " <- ", pp[4], if (pnlr[4]) " + rned", ";", pof[5], " <- ", pp[5], if (pnlr[5]) " + rnassym", ";")  
  
  modp <- paste(c("rnslope ~ normal(0, sigma_slope);", "rnlasy ~ normal(0, sigma_lasy);", "rnuasy ~ normal(0, sigma_uasy);", "rned ~ normal(0, sigma_ed);", "rnassym ~ normal(0, sigma_assym);")[pnlr], collapse=" ") 
  
  
  #### model statements
  mod <- paste(c(paste("slope ~", priors$b, ";"), 
                 paste("lasy ~", priors$c, ";"), 
                 paste("uasy ~", priors$d, ";"), 
                 paste("ed ~", priors$e, ";"), 
                 paste("assym ~", priors$f, ";"))[!isfix], collapse=" ") 
  
  mody <- paste("y ~ bernoulli(mu);")
  
  modr <- paste(c("rslope ~ normal(0, sigma_slope);", "rlasy ~ normal(0, sigma_lasy);", "ruasy ~ normal(0, sigma_uasy);", "red ~ normal(0, sigma_ed);", "rassym ~ normal(0, sigma_assym);")[pnlr], collapse=" ") 
  modrsig <- paste(c(paste("sigmasq_slope ~", priors$sb, ";"), 
                     paste("sigmasq_lasy ~", priors$sc, ";"), 
                     paste("sigmasq_uasy ~", priors$sd, ";"), 
                     paste("sigmasq_ed ~", priors$se, ";"), 
                     paste("sigmasq_assym ~",priors$sf, ";"))[pnlr], collapse=" ") 
  
  # new random effects for derived ed
  moded <- paste(c("rnslope ~ normal(0, sigma_slope);", 
                   "rnlasy ~ normal(0, sigma_lasy);", 
                   "rnuasy ~ normal(0, sigma_uasy);", 
                   "rned ~ normal(0, sigma_ed);", 
                   "rnassym ~ normal(0, sigma_assym);")[pnlr], collapse=" ") 
  
  stancode <- paste("data {", 
                    "int<lower=0> N; int y[N];",
                    "real<lower=0> x[N];", 
                    if (!is.null(idc)) dJ,
                    if (!is.null(random)) dK,
                    dpc,
                    "} ",
                    "parameters {",
                    parac,
                    if (!is.null(random)) rparac, 
                    if (!is.null(random)) rnparac,  
                    "} ",
                    "transformed parameters {",
                    "real mu[N];",                   
                    trans,                      
                    "} ",
                    "model {",
                    mod,
                    mody,                     
                    if (!is.null(random)) modr,
                    if (!is.null(random)) modrsig,    
                    if (!is.null(random)) modp,                    
                    "}",
                    "generated quantities {",
                    "real residuals[N];", 
                    "real log_lik[N];",
                    popc,
                    "for (i in 1:N){",
                    "residuals[i] <- y[i] - mu[i];",
                    "log_lik[i] <- bernoulli_log(y[i], mu[i]);",
                    "}",
                    ppc,                                     
                    "}",
                    sep="")
  
  eval(parse(text="assign('stancode', stancode, envir=.GlobalEnv)"))
  fit <- stan(model_code = 'stancode', data = stan_dat, ...)
  
  out <- list()
  out$call <- callDetail
  out$data <- stan_dat
  out$model <- stancode
  out$stan <- fit
  out$fixedpars <- trap
  out$pars <- c("pslope", "plasy", "puasy", "ped", "passym")
  if (!is.null(random)) out$random <- c("slope", "lasy", "uasy", "ed", "assym")[pnlr] 
  out$curves <- list(pars=pnl, J=J, names=curvenames)
  out$fixed <- fix
  out$fct <- fct
  class(out) <- "standrc"
  return(out)
}
