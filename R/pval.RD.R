pval.RD <- function(f0 = 0, f1 = 0, e0, n0, e1, n1,
  alternative = "two.sided", RD.null = 0){
  e0.m <- e0 + f0
  ne0.m <- n0 - e0.m
  e1.m <- e1 + f1
  ne1.m <- n1 - e1.m
  if((abs(e0.m) < 1e-6 & abs(e1.m) < 1e-6) |
    (abs(ne0.m) < 1e-6 & abs(ne1.m) < 1e-6)){
    pval <- 1
  }else{
    RD <- e1.m/(e1.m + ne1.m) - e0.m/(e0.m + ne0.m)
    se.RD <- sqrt(e0.m/(e0.m + ne0.m)^3*ne0.m +
      e1.m/(e1.m + ne1.m)^3*ne1.m)
    if(alternative == "two.sided"){
      pval <- 2*pnorm(-abs(RD - RD.null)/se.RD)
    }
    if(alternative == "one.sided"){
      pval <- pnorm(-abs(RD - RD.null)/se.RD)
    }
  }
  return(pval)
}