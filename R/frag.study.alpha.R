frag.study.alpha <- function(e0, n0, e1, n1, data, methods,
  modify0 = "both", modify1 = "both",
  alpha.from = 0.005, alpha.to = 0.05, alpha.breaks = 100,
  alternative = "two.sided", OR = 1, RR = 1, RD = 0){
  if(!missing(data)){
    e0 <- eval(substitute(e0), data, parent.frame())
    n0 <- eval(substitute(n0), data, parent.frame())
    e1 <- eval(substitute(e1), data, parent.frame())
    n1 <- eval(substitute(n1), data, parent.frame())
  }
  if(e0 < 0 | n0 < 0 | e1 < 0 | n1 < 0) stop("Event counts and sample sizes must be nonnegative.")
  if(e0 > n0 | e1 > n1) stop("Event counts must not be larger than sample sizes.")
  if(e0%%1 != 0 | n0%%1 != 0 | e1%%1 != 0 | n1%%1 != 0){
    message("Some event counts and/or sample sizes are not integer(s); they are rounded.")
    e0 <- round(e0)
    n0 <- round(n0)
    e1 <- round(e1)
    n1 <- round(n1)
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
  if(alpha.from >= alpha.to) stop("alpha.from must be less than alpha.to.")
  if(alpha.from <= 0 | alpha.to > 1) stop("alpha.from and alpha.to must be between 0 and 1.")
  if(alpha.breaks%%1 != 0){
    message("alpha.breaks is not an integer; it is rounded.")
    alpha.breaks <- round(alpha.breaks)
  }
  if(alpha.breaks <= 0){
    stop("alpha.breaks must be a positive integer.")
  }
  if(length(alternative) != 1){
    stop("Only one type of alternative hypothesis\n    can be specified for alternative.")
  }
  if(any(!is.element(alternative, c("two.sided", "one.sided")))){
    stop("alternative must be \"two.sided\" or \"one.sided\".")
  }
  if(OR < 0) stop("The null value of OR must be nonnegative.")
  if(RR < 0) stop("The null value of RR must be nonnegative.")
  if(RD < -1 | RD > 1) stop("The null value of RD must be between -1 and 1.")

  ne0 <- n0 - e0
  ne1 <- n1 - e1
  ori.data <- matrix(c(e0, ne0, e1, ne1), 2, 2, byrow = TRUE)
  rownames(ori.data) <- c("group 0", "group 1")
  colnames(ori.data) <- c("event", "no event")

  if(modify0 == "increase") f0.range <- c(0, ne0)
  if(modify0 == "decrease") f0.range <- c(-e0, 0)
  if(modify0 == "both") f0.range <- c(-e0, ne0)
  if(modify0 == "none") f0.range <- c(0, 0)
  if(modify1 == "increase") f1.range <- c(0, ne1)
  if(modify1 == "decrease") f1.range <- c(-e1, 0)
  if(modify1 == "both") f1.range <- c(-e1, ne1)
  if(modify1 == "none") f1.range <- c(0, 0)
  f0.mods <- f0.range[1]:f0.range[2]
  f1.mods <- f1.range[1]:f1.range[2]
  alphas <- seq(from = alpha.from, to = alpha.to,
    length.out = alpha.breaks)

  out <- list(data = ori.data, methods = methods, alphas = alphas,
    alternative = alternative, null = c("OR" = OR, "RR" = RR, "RD" = RD),
    modify0 = modify0, modify1 = modify1,
    f0.range = f0.range, f1.range = f1.range)

  temp.Fisher <- function(f0, f1){
    fisher.test(rbind(c(e0 + f0, n0 - e0 - f0),
      c(e1 + f1, n1 - e1 - f1)),
      alternative = "two.sided")$p.value
  }
  temp.chisq <- function(f0, f1){
    if((abs(e0 + f0) < 1e-6 & abs(e1 + f1) < 1e-6) |
      (abs(n0 - e0 - f0) < 1e-6 & abs(n1 - e1 - f1) < 1e-6)){
      temp.out <- 1
    }else{
      temp.out <- suppressWarnings(
        chisq.test(rbind(c(e0 + f0, n0 - e0 - f0),
        c(e1 + f1, n1 - e1 - f1)))$p.value)
    }
    return(temp.out)
  }
  temp.OR <- function(f0, f1){
    pval.logOR(f0 = f0, f1 = f1, e0 = e0, n0 = n0,
      e1 = e1, n1 = n1, alternative = alternative, OR.null = OR)
  }
  temp.RR <- function(f0, f1){
    pval.logRR(f0 = f0, f1 = f1, e0 = e0, n0 = n0,
      e1 = e1, n1 = n1, alternative = alternative, RR.null = RR)
  }
  temp.RD <- function(f0, f1){
    pval.RD(f0 = f0, f1 = f1, e0 = e0, n0 = n0,
      e1 = e1, n1 = n1, alternative = alternative, RD.null = RD)
  }

  if(modify0 == "none" & modify1 == "none"){
    pval <- NULL
    for(i in 1:length(methods)){
      temp.m <- methods[i]
      temp <- get(paste0("temp.", temp.m))
      pval.temp <- temp(0, 0)
      pval <- c(pval, pval.temp)
    }
    names(pval) <- methods
    out <- c(out, list(pval = pval))
    class(out) <- c("frag.alpha", "frag.study.alpha")
    return(out)
  }

  uniq.mods <- 1:(max(abs(f0.mods)) + max(abs(f1.mods)))

  pval <- FI <- FI.avg <- FQ <- FQ.avg <- NULL
  for(i in 1:length(methods)){
    temp.m <- methods[i]
    temp <- get(paste0("temp.", temp.m))
    temp.pval <- temp(0, 0)
    temp.signif <- alphas > temp.pval
    temp.FI <- temp.FQ <- temp.check <- rep(NA, alpha.breaks)
    for(f in uniq.mods){
      tt.mods <- rbind(cbind((-f):f, f - abs((-f):f)),
        cbind((-(f - 1)):(f - 1), -(f - abs((-(f - 1)):(f - 1)))))
      tt.mods <- tt.mods[is.element(tt.mods[,1], f0.mods) & is.element(tt.mods[,2], f1.mods),]
      if(is.vector(tt.mods)) tt.mods <- matrix(tt.mods, ncol = 2)
      for(j in 1:dim(tt.mods)[1]){
        tt.pval <- temp(f0 = tt.mods[j,1], f1 = tt.mods[j,2])
        if(any(temp.signif & tt.pval >= alphas & is.na(temp.check))){
          tt.idx.alpha <- which(temp.signif & tt.pval >= alphas & is.na(temp.check))
          temp.FI[tt.idx.alpha] <- f
          temp.FQ[tt.idx.alpha] <- f/(n0 + n1)
          temp.check[tt.idx.alpha] <- 1
        }
        if(any(!temp.signif & tt.pval < alphas & is.na(temp.check))){
          tt.idx.alpha <- which(!temp.signif & tt.pval < alphas & is.na(temp.check))
          temp.FI[tt.idx.alpha] <- f
          temp.FQ[tt.idx.alpha] <- f/(n0 + n1)
          temp.check[tt.idx.alpha] <- 1
        }
        if(all(!is.na(temp.check))) break
      }
      if(all(!is.na(temp.check))) break
    }
    if(any(is.na(temp.check))){
      temp.FI[is.na(temp.check)] <- temp.FQ[is.na(temp.check)] <- NA
    }
    pval <- c(pval, temp.pval)
    FI <- cbind(FI, temp.FI)
    FI.avg <- c(FI.avg, mean(temp.FI))
    FQ <- cbind(FQ, temp.FQ)
    FQ.avg <- c(FQ.avg, mean(temp.FQ))
  }

  FI <- cbind(alphas, FI)
  FQ <- cbind(alphas, FQ)
  names(pval) <- names(FI.avg) <- names(FQ.avg) <- methods
  colnames(FI) <- colnames(FQ) <- c("alpha", methods)
  out <- c(out, list(pval = pval, FI = FI, FI.avg = FI.avg,
    FQ = FQ, FQ.avg = FQ.avg))
  class(out) <- c("frag.alpha", "frag.study.alpha")
  return(out)
}