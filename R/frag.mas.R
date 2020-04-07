frag.mas <- function(e0, n0, e1, n1, ma.id, data, measure = "OR",
  alpha = 0.05, mod.dir = "both", OR = 1, RR = 1, RD = 0,
  method = "DL", test = "z", ...){
  if(!missing(data)){
    e0 <- eval(substitute(e0), data, parent.frame())
    n0 <- eval(substitute(n0), data, parent.frame())
    e1 <- eval(substitute(e1), data, parent.frame())
    n1 <- eval(substitute(n1), data, parent.frame())
    ma.id <- eval(substitute(ma.id), data, parent.frame())
  }
  if(any(e0 < 0) | any(n0 < 0) | any(e1 < 0) | any(n1 < 0)){
    stop("Event counts and sample sizes must be nonnegative.")
  }
  if(any(e0 > n0) | any(e1 > n1)){
    stop("Event counts must not be larger than sample sizes.")
  }
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
  if(any(!is.element(measure, c("OR", "RR", "RD")))){
    stop("measure must be \"OR\", \"RR\", or \"RD\".")
  }
  if(length(measure) != 1){
    stop("Only one measure can be specified.")
  }
  if(alpha <= 0 | alpha >= 1) stop("alpha must be between 0 and 1.")
  if(any(!is.element(mod.dir, c("both", "one", "left", "right")))){
    stop("mod.dir must be \"both\", \"one\", \"left\", or \"right\".")
  }
  if(length(mod.dir) != 1){
    stop("Only one choice can be specified for mod.dir.")
  }
  if(OR <= 0) stop("The null value of OR must be positive.")
  if(RR <= 0) stop("The null value of RR must be positive.")
  if(RD < -1 | RD > 1) stop("The null value of RD must be between -1 and 1.")

  if(measure == "OR") null.val <- log(OR)
  if(measure == "RR") null.val <- log(RR)
  if(measure == "RD") null.val <- RD

  out <- list(measure = measure, alpha = alpha,
    null = null.val, mod.dir = mod.dir)

  ma.ids <- unique(ma.id)
  n.ma <- length(ma.ids)
  ci.ori <- matrix(NA, n.ma, 2)
  colnames(ci.ori) <- c("LB", "UB")
  est.ori <- pval.ori <- FI <- FQ <- rep(NA, n.ma)
  for(i in 1:n.ma){
    temp.e0 <- e0[ma.id == ma.ids[i]]
    temp.n0 <- n0[ma.id == ma.ids[i]]
    temp.e1 <- e1[ma.id == ma.ids[i]]
    temp.n1 <- n1[ma.id == ma.ids[i]]
    temp.out <- frag.ma(e0 = temp.e0, n0 = temp.n0, e1 = temp.e1, n1 = temp.n1,
      measure = measure, alpha = alpha, mod.dir = mod.dir, OR = OR, RR = RR, RD = RD,
      method = method, test = test, ...)
    est.ori[i] <- temp.out$est.ori
    ci.ori[i,] <- temp.out$ci.ori
    pval.ori[i] <- temp.out$pval
    FI[i] <- temp.out$FI
    FQ[i] <- temp.out$FQ
  }

  out <- c(out, list(est.ori = est.ori, ci.ori = ci.ori, pval.ori = pval.ori, FI = FI, FQ = FQ))
  class(out) <- c("frag.mas", "frag.multi")
  return(out)
}