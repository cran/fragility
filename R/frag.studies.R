frag.studies <- function(e0, n0, e1, n1, data, methods,
  modify0 = "both", modify1 = "both", alpha = 0.05,
  alternative = "two.sided", OR = 1, RR = 1, RD = 0){
  if(!missing(data)){
    e0 <- eval(substitute(e0), data, parent.frame())
    n0 <- eval(substitute(n0), data, parent.frame())
    e1 <- eval(substitute(e1), data, parent.frame())
    n1 <- eval(substitute(n1), data, parent.frame())
  }
  if(any(e0 < 0) | any(n0 < 0) | any(e1 < 0) | any(n1 < 0)){
    stop("Event counts and sample sizes must be nonnegative.")
  }
  if(any(e0 > n0) | any(e1 > n1)) stop("Event counts must not be larger than sample sizes.")
  if(any(e0%%1 != 0) | any(n0%%1 != 0) | any(e1%%1 != 0) | any(n1%%1 != 0)){
    message("Some event counts and/or sample sizes are not integer(s); they are rounded.")
    e0 <- round(e0)
    n0 <- round(n0)
    e1 <- round(e1)
    n1 <- round(n1)
  }
  if(length(n0) != length(e0) | length(e1) != length(e0) | length(n1) != length(e0)){
    stop("e0, n0, e1, and n1 do not have the same length.")
  }
  if(length(e0) == 1){
    stop("The fragility of a single study can be assessed using frag.study().")
  }
  if(missing(methods)){
    methods <- c("Fisher", "chisq", "OR", "RR", "RD")
  }
  if(any(!is.element(methods, c("Fisher", "chisq", "OR", "RR", "RD")))){
    stop("methods must be \"Fisher\", \"chisq\", \"OR\", \"RR\", and/or \"RD\".")
  }
  if(length(modify0) != 1 | length(modify1) != 1){
    stop("Only one type of event status modification\n    can be specified for modify0/modify1.")
  }
  if(any(!is.element(c(modify0, modify1), c("increase", "decrease", "both", "none")))){
    stop("modify0/modify1 must be \"increase\", \"decrease\", \"both\", or \"none\".")
  }
  if(alpha <= 0 | alpha >= 1) stop("alpha must be between 0 and 1.")
  if(length(alternative) != 1){
    stop("Only one type of alternative hypothesis\n    can be specified for alternative.")
  }
  if(any(!is.element(alternative, c("two.sided", "one.sided")))){
    stop("alternative must be \"two.sided\" or \"one.sided\".")
  }
  if(OR < 0) stop("The null value of OR must be nonnegative.")
  if(RR < 0) stop("The null value of RR must be nonnegative.")
  if(RD < -1 | RD > 1) stop("The null value of RD must be between -1 and 1.")

  out <- list(methods = methods, alpha = alpha,
    alternative = alternative, null = c("OR" = OR, "RR" = RR, "RD" = RD),
    modify0 = modify0, modify1 = modify1)

  pval <- FI <- FQ <- matrix(NA, length(e0), length(methods))
  for(sid in 1:length(e0)){
    temp.out <- frag.study.iter(e0 = e0[sid], n0 = n0[sid],
      e1 = e1[sid], n1 = n1[sid], methods, modify0, modify1, alpha,
      alternative, OR, RR, RD, allcase = FALSE)
    pval[sid,] <- temp.out$pval
    FI[sid,] <- temp.out$FI
    FQ[sid,] <- temp.out$FQ
  }

  colnames(pval) <- colnames(FI) <- colnames(FQ) <- methods
  if(modify0 == "none" & modify1 == "none"){
    out <- c(out, list(pval = pval))
  }else{
    out <- c(out, list(pval = pval, FI = FI, FQ = FQ))
  }
  class(out) <- c("frag.multi", "frag.studies")
  return(out)
}