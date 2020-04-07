frag.study.all <- function(e0, n0, e1, n1, methods,
  modify0, modify1, alpha, alternative, OR, RR, RD){
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
  tot.mods <- outer(abs(f0.mods), abs(f1.mods), "+")

  out <- list(data = ori.data, methods = methods, alpha = alpha,
    alternative = alternative, null = c("OR" = OR, "RR" = RR, "RD" = RD),
    modify0 = modify0, modify1 = modify1,
    f0.range = f0.range, f1.range = f1.range, tot.mods = tot.mods)

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
      pval.temp <-  temp(0, 0)
      pval <- c(pval, pval.temp)
    }
    names(pval) <- methods
    out <- c(out, list(pval = pval))
    class(out) <- c("frag.study", "frag.study.all")
    return(out)
  }

  pval <- NULL
  FI <- FQ <- dir <- NULL
  FI0 <- FQ0 <- dir0 <- NULL
  FI1 <- FQ1 <- dir1 <- NULL
  pvals <- mods <- mods0 <- mods1 <- rep(list(NULL), length(methods))
  for(i in 1:length(methods)){
    temp.m <- methods[i]
    temp <- get(paste0("temp.", temp.m))
    temp <- Vectorize(temp)
    temp.pvals <- outer(f0.mods, f1.mods, temp)
    temp.pval <- temp.pvals[f0.mods == 0, f1.mods == 0]
    temp.FI <- temp.FQ <- temp.dir <- temp.mods <- NA
    temp.FI0 <- temp.FQ0 <- temp.dir0 <- temp.mods0 <- NA
    temp.FI1 <- temp.FQ1 <- temp.dir1 <- temp.mods1 <- NA
    if(temp.pval < alpha){
      if(any(temp.pvals >= alpha)){
        temp.FI <- min(tot.mods[temp.pvals >= alpha])
        temp.FQ <- temp.FI/(n0 + n1)
        temp.dir <- "significance altered to non-significance"
        temp.mods <- which(tot.mods == temp.FI &
           temp.pvals >= alpha, arr.ind = TRUE)
        temp.mods <- cbind(f0.mods[temp.mods[,1]],
          f1.mods[temp.mods[,2]])
        colnames(temp.mods) <- c("group 0", "group 1")
        if(modify0 != "none" & modify1 != "none"){
          temp.pvals0 <- temp.pvals[, f1.mods == 0]
          tot.mods0 <- tot.mods[, f1.mods == 0]
          temp.pvals1 <- temp.pvals[f0.mods == 0,]
          tot.mods1 <- tot.mods[f0.mods == 0,]
          if(any(temp.pvals0 >= alpha)){
            temp.FI0 <- min(tot.mods0[temp.pvals0 >= alpha])
            temp.FQ0 <- temp.FI0/(n0 + n1)
            temp.dir0 <- "significance altered to non-significance"
            temp.mods0 <- f0.mods[tot.mods0 == temp.FI0 &
              temp.pvals0 >= alpha]
          }else{
            temp.FI0 <- temp.FQ0 <- NA
            temp.dir0 <- "significance cannot be altered"
          }
          if(any(temp.pvals1 >= alpha)){
            temp.FI1 <- min(tot.mods1[temp.pvals1 >= alpha])
            temp.FQ1 <- temp.FI1/(n0 + n1)
            temp.dir1 <- "significance altered to non-significance"
            temp.mods1 <- f1.mods[tot.mods1 == temp.FI1 &
              temp.pvals1 >= alpha]
          }else{
            temp.FI1 <- temp.FQ1 <- NA
            temp.dir1 <- "significance cannot be altered"
          }
        }
      }else{
        temp.FI <- temp.FQ <- NA
        temp.dir <- "significance cannot be altered"
      }
    }else{
      if(any(temp.pvals < alpha)){
        temp.FI <- min(tot.mods[temp.pvals < alpha])
        temp.FQ <- temp.FI/(n0 + n1)
        temp.dir <- "non-significance altered to significance"
        temp.mods <- which(tot.mods == temp.FI &
           temp.pvals < alpha, arr.ind = TRUE)
        temp.mods <- cbind(f0.mods[temp.mods[,1]],
          f1.mods[temp.mods[,2]])
        colnames(temp.mods) <- c("group 0", "group 1")
        if(modify0 != "none" & modify1 != "none"){
          temp.pvals0 <- temp.pvals[, f1.mods == 0]
          tot.mods0 <- tot.mods[, f1.mods == 0]
          temp.pvals1 <- temp.pvals[f0.mods == 0,]
          tot.mods1 <- tot.mods[f0.mods == 0,]
          if(any(temp.pvals0 < alpha)){
            temp.FI0 <- min(tot.mods0[temp.pvals0 < alpha])
            temp.FQ0 <- temp.FI0/(n0 + n1)
            temp.dir0 <- "non-significance altered to significance"
            temp.mods0 <- f0.mods[tot.mods0 == temp.FI0 &
              temp.pvals0 < alpha]
          }else{
            temp.FI0 <- temp.FQ0 <- NA
            temp.dir0 <- "non-significance cannot be altered"
          }
          if(any(temp.pvals1 < alpha)){
            temp.FI1 <- min(tot.mods1[temp.pvals1 < alpha])
            temp.FQ1 <- temp.FI1/(n0 + n1)
            temp.dir1 <- "non-significance altered to significance"
            temp.mods1 <- f1.mods[tot.mods1 == temp.FI1 &
              temp.pvals1 < alpha]
          }else{
            temp.FI1 <- temp.FQ1 <- NA
            temp.dir1 <- "non-significance cannot be altered"
          }
        }
      }else{
        temp.FI <- temp.FQ <- NA
        temp.dir <- "non-significance cannot be altered"
      }
    }

    pval <- c(pval, temp.pval)
    FI <- c(FI, temp.FI)
    FQ <- c(FQ, temp.FQ)
    dir <- c(dir, temp.dir)
    mods[[i]] <- temp.mods
    pvals[[i]] <- temp.pvals
    if(modify0 != "none" & modify1 != "none"){
      FI0 <- c(FI0, temp.FI0)
      FQ0 <- c(FQ0, temp.FQ0)
      dir0 <- c(dir0, temp.dir0)
      mods0[[i]] <- temp.mods0
      FI1 <- c(FI1, temp.FI1)
      FQ1 <- c(FQ1, temp.FQ1)
      dir1 <- c(dir1, temp.dir1)
      mods1[[i]] <- temp.mods1
    }
  }

  dir <- as.matrix(dir)
  names(pval) <- names(FI) <- names(FQ) <-
    rownames(dir) <- names(mods) <- names(pvals) <- methods
  out <- c(out, list(pval = pval, FI = FI, FQ = FQ,
    dir = dir, mods = mods, pvals = pvals))
  if(modify0 != "none" & modify1 != "none"){
    dir0 <- as.matrix(dir0)
    names(FI0) <- names(FQ0) <-
      rownames(dir0) <- names(mods0) <- methods
    out <- c(out, list(FI0 = FI0, FQ0 = FQ0,
      dir0 = dir0, mods0 = mods0))
    dir1 <- as.matrix(dir1)
    names(FI1) <- names(FQ1) <-
      rownames(dir1) <- names(mods1) <- methods
    out <- c(out, list(FI1 = FI1, FQ1 = FQ1,
      dir1 = dir1, mods1 = mods1))
  }

  class(out) <- c("frag.study", "frag.study.all")
  return(out)
}