standrc_priors <- function(b="normal(pb, 100)",
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
                           pf=NULL){
  list(b=b, c=c, d=d, e=e, f=f, sy=sy, sb=sb, sc=sc, sd=sd, se=se, sf=sf, pb=pb, pc=pc, pd=pd, pe=pe, pf=pf)  
}