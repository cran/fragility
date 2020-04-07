plot.frag.alpha <- function(x, method, fragility = "FI",
  percentage = TRUE, xlab, ylab, xlim, ylim, cex.pts,
  col.line, col.pval, col.sig, lty.pval, lwd, lwd.pval,
  pch, pch.na, tid1, tid2, FQ.nma = FALSE, ...){
  if(!inherits(x, "frag.alpha")){
    stop("The input must be an object of class \"frag.alpha\".")
  }
  if(inherits(x, "frag.study.alpha")){
    if(missing(method)){
      method <- x$methods[1]
    }
    if(!is.element(method, x$methods)){
      stop("method must be used to produce the object of class \"frag.alpha\".")
    }
    if(length(method) != 1){
      stop("Only one method can be specified for plotting each time.")
    }
  }
  if(!is.element(fragility, c("FI", "FQ")) | length(fragility) != 1){
    stop("fragility must be either \"FI\" or \"FQ\".")
  }
  if(missing(col.sig)) col.sig <- c("forestgreen", "firebrick")
  if(length(col.sig) != 2){
    stop("col.sig must be a vector of two elements of colors.")
  }
  if(inherits(x, "frag.nma.alpha")){
    tids <- rownames(x$pval.ori)
    if(missing(tid1) & missing(tid2)){
      tid.f <- strsplit(x$tid.f, split = " vs. ")
      for(i in 1:length(tid.f)){
        if(tid.f[[i]][1] != tid.f[[i]][2]){
          tid1 <- tid.f[[i]][1]
          tid2 <- tid.f[[i]][2]
          break
        }
      }
    }
    if((missing(tid1) & !missing(tid2)) | (!missing(tid1) & missing(tid2))){
      stop("Both tid1 and tid2 should be specified.")
    }
    if(length(tid1) != 1 | length(tid2) != 1){
      stop("Both tid1 and tid2 should have a length of one.")
    }
    if(!is.element(tid1, tids) | !is.element(tid2, tids)){
      stop("Both tid1 and tid2 should treatment IDs.")
    }
    if(tid1 == tid2) stop("tid1 and tid2 are the same.")
  }

  alphas <- x$alphas
  if(inherits(x, "frag.study.alpha")){
    pval <- x$pval[method]
    if(fragility == "FI"){
      ypts <- x$FI[, method]
    }else{
      ypts <- x$FQ[, method]
    }
  }
  if(inherits(x, "frag.ma.alpha")){
    pval <- x$pval.ori
    if(fragility == "FI"){
      ypts <- x$FI
    }else{
      ypts <- x$FQ
    }
  }
  if(inherits(x, "frag.nma.alpha")){
    pval <- x$pval.ori[tid1, tid2]
    if(fragility == "FI"){
      ypts <- x$FI[[paste(tid1, "vs.", tid2)]]
    }else{
      if(!FQ.nma){
        ypts <- x$FQ[[paste(tid1, "vs.", tid2)]]
      }else{
        ypts <- x$FQ.nma[[paste(tid1, "vs.", tid2)]]
      }
    }
  }
  sig2nonsig <- pval < alphas
  if(all(is.na(ypts))) stop("All fragility measures are NA.")
  na.idx <- is.na(ypts)
  if(length(na.idx) > 0) ypts[na.idx] <- max(ypts[!na.idx])
  if(fragility == "FQ" & percentage){
    ypts <- ypts * 100
  }

  if(missing(xlab)) xlab <- "Significance level"
  if(missing(xlim)) xlim <- c(min(alphas), max(alphas))
  if(missing(ylab)){
    if(fragility == "FI") ylab <- "Fragility index"
    if(fragility == "FQ"){
      ylab <- "Fragility quotient"
      if(percentage) ylab <- paste0(ylab, ", %")
    }
  }
  if(missing(ylim)){
    if(min(ypts) < max(ypts)){
      ylim <- c(min(ypts), max(ypts))
    }else{
      ylim <- NULL
    }
  }
  if(missing(pch)) pch <- 16
  if(missing(pch.na)) pch.na <- 1
  if(missing(col.line)) col.line <- "gray50"
  if(missing(col.pval)) col.pval <- "gray50"
  if(missing(lwd)) lwd <- 1
  if(missing(lwd.pval)) lwd.pval <- 1
  if(missing(cex.pts)) cex.pts <- 0.5
  if(missing(lty.pval)) lty.pval <- 2

  pchs <- rep(pch, length(alphas))
  if(length(na.idx) > 0) pchs[na.idx] <- pch.na
  cols <- rep(NA, length(alphas))
  if(any(sig2nonsig)) cols[sig2nonsig] <- col.sig[2]
  if(any(!sig2nonsig)) cols[!sig2nonsig] <- col.sig[1]

  plot.default(mean(alphas), mean(ypts),
    xlab = xlab, xlim = xlim, ylab = ylab, ylim = ylim, type = "n", ...)
  lines(alphas, ypts, col = col.line, lwd = lwd)
  points(alphas, ypts, cex = cex.pts, pch = pchs, col = cols)
  if(pval >= min(alphas) & pval <= max(alphas)){
    abline(v = pval, lty = lty.pval, col = col.pval, lwd = lwd.pval)
  }
}