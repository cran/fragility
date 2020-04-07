pval.logOR <- function(f0 = 0, f1 = 0, e0, n0, e1, n1,
  alternative = "two.sided", OR.null = 1){
  e0.m <- e0 + f0
  ne0.m <- n0 - e0.m
  e1.m <- e1 + f1
  ne1.m <- n1 - e1.m
  if((abs(e0.m) < 1e-6 & abs(e1.m) < 1e-6) |
    (abs(ne0.m) < 1e-6 & abs(ne1.m) < 1e-6)){
    pval <- 1
  }else{
    if(abs(e0.m) < 1e-6 | abs(ne0.m) < 1e-6 |
      abs(e1.m) < 1e-6 | abs(ne1.m) < 1e-6){
      e0.m <- e0.m + 0.5
      ne0.m <- ne0.m + 0.5
      e1.m <- e1.m + 0.5
      ne1.m <- ne1.m + 0.5
    }
    logOR <- log(e1.m/ne1.m) - log(e0.m/ne0.m)
    se.logOR <- sqrt(1/e0.m + 1/ne0.m + 1/e1.m + 1/ne1.m)
    if(alternative == "two.sided"){
      pval <- 2*pnorm(-abs(logOR - log(OR.null))/se.logOR)
    }
    if(alternative == "one.sided"){
      pval <- pnorm(-abs(logOR - log(OR.null))/se.logOR)
    }
  }
  return(pval)
}