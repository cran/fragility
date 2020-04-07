frag.study.iter <- function(e0, n0, e1, n1, methods,
  modify0, modify1, alpha, alternative, OR, RR, RD, allcase){
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

  out <- list(data = ori.data, methods = methods, alpha = alpha,
    alternative = alternative, null = c("OR" = OR, "RR" = RR, "RD" = RD),
    modify0 = modify0, modify1 = modify1,
    f0.range = f0.range, f1.range = f1.range, allcase = allcase)

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

  if(all(f0.mods == 0) & all(f1.mods == 0)){
    pval <- NULL
    for(i in 1:length(methods)){
      temp.m <- methods[i]
      temp <- get(paste0("temp.", temp.m))
      pval.temp <- temp(0, 0)
      pval <- c(pval, pval.temp)
    }
    names(pval) <- methods
    out <- c(out, list(pval = pval))
    class(out) <- c("frag.study")
    return(out)
  }

  uniq.mods <- 1:(max(abs(f0.mods)) + max(abs(f1.mods)))

  pval <- FI <- FQ <- dir <- NULL
  mods <- rep(list(NULL), length(methods))
  for(i in 1:length(methods)){
    temp.m <- methods[i]
    temp <- get(paste0("temp.", temp.m))
    temp.pval <- temp(0, 0)
    temp.signif <- temp.pval < alpha
    signif.alter <- 0
    temp.mods <- NULL
    for(f in uniq.mods){
      tt.mods <- rbind(cbind((-f):f, f - abs((-f):f)),
        cbind((-(f - 1)):(f - 1), -(f - abs((-(f - 1)):(f - 1)))))
      tt.mods <- tt.mods[is.element(tt.mods[,1], f0.mods) & is.element(tt.mods[,2], f1.mods),]
      if(is.vector(tt.mods)) tt.mods <- matrix(tt.mods, ncol = 2)
      for(j in 1:dim(tt.mods)[1]){
        tt.pval <- temp(f0 = tt.mods[j,1], f1 = tt.mods[j,2])
        tt.signif <- tt.pval < alpha
        if(tt.signif != temp.signif){
          temp.FI <- f
          temp.FQ <- temp.FI/(n0 + n1)
          signif.alter <- 1
          temp.mods <- rbind(temp.mods, c(tt.mods[j,1], tt.mods[j,2]))
          if(!allcase) break
        }
      }
      if(signif.alter == 1) break
    }
    if(signif.alter == 0){
      temp.FI <- temp.FQ <- NA
      if(temp.signif){
        temp.dir <- "significance cannot be altered"
      }else{
        temp.dir <- "non-significance cannot be altered"
      }
    }else{
      if(temp.signif){
        temp.dir <- "significance altered to non-significance"
      }else{
        temp.dir <- "non-significance altered to significance"
      }
      colnames(temp.mods) <- c("group 0", "group 1")
    }
    pval <- c(pval, temp.pval)
    FI <- c(FI, temp.FI)
    FQ <- c(FQ, temp.FQ)
    dir <- c(dir, temp.dir)
    if(!is.na(temp.FI)){
      mods[[i]] <- temp.mods
    }else{
      mods[[i]] <- NA
    }
  }

  dir <- as.matrix(dir)
  names(pval) <- names(FI) <- names(FQ) <-
    rownames(dir) <- names(mods) <- methods
  out <- c(out, list(pval = pval, FI = FI, FQ = FQ,
    dir = dir, mods = mods))
  class(out) <- c("frag.study")
  return(out)
}