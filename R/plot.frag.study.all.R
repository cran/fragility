plot.frag.study.all <- function(x, method, modify0, modify1, trun,
  xlab, ylab, xlim, ylim, cex.pts, cex.legend.pval, cex.legend.title,
  col.ori, col.ori.hl, col.f.hl, col.sig, lty.ori, lwd.ori,
  pch, pch.ori, pch.ori.hl, pch.f, pch.f.hl, pch.trun,
  adjust.legend, adjust.seg, legend.pvals, ...){
  if(!inherits(x, "frag.study.all")){
    stop("The input must be an object of class \"frag.study.all\".")
  }
  if(missing(method)){
    method <- x$methods[1]
  }
  if(!is.element(method, x$methods)){
    stop("method must be used to produce the object of class \"frag.study.all\".")
  }
  if(length(method) != 1){
    stop("Only one method can be specified for plotting each time.")
  }
  if(missing(modify0)){
    if(any(x$f0.range != 0)) modify0 <- TRUE else modify0 <- FALSE
  }
  if(missing(modify1)){
    if(any(x$f1.range != 0)) modify1 <- TRUE else modify1 <- FALSE
  }

  pvals <- x$pvals[[method]]
  f0.range <- x$f0.range
  f1.range <- x$f1.range
  alpha <- x$alpha
  tot.mods <- x$tot.mods
  mods <- x$mods[[method]]

  if(modify0 + modify1 == 0){
    message("No event status modified; no plot is generated.")
    return(invisible())
  }

  if(missing(trun)) trun <- 10
  if(missing(cex.pts)) cex.pts <- 0.5
  if(missing(pch.ori)) pch.ori <- 15
  if(missing(pch.f)) pch.f <- 17
  if(missing(lty.ori)) lty.ori <- 2
  if(missing(col.ori)) col.ori <- "gray50"
  if(missing(lwd.ori)) lwd.ori <- 1

  if(modify0 + modify1 == 1){
    if(modify0){
      if(all(f0.range == 0)){
        stop("Group 0 is not modified; no plot is generated.")
        return(invisible())
      }
      pvals.pts <- pvals[,f1.range[1]:f1.range[2] == 0]
      xpts <- f0.range[1]:f0.range[2]
      xlimit <- f0.range
      if(missing(xlab)){
        xlab <- "Number of modified events in group 0"
      }
      if(x$modify1 == "none"){
        mods.single <- x$mods[[method]][1]
      }else{
        mods.single <- x$mods0[[method]]
      }
    }
    if(modify1){
      if(all(f1.range == 0)){
        stop("Group 1 is not modified; no plot is generated.")
      }
      pvals.pts <- pvals[f0.range[1]:f0.range[2] == 0,]
      xpts <- f1.range[1]:f1.range[2]
      xlimit <- f1.range
      if(missing(xlab)){
        xlab <- "Number of modified events in group 1"
      }
      if(x$modify0 == "none"){
        mods.single <- x$mods[[method]][2]
      }else{
        mods.single <- x$mods1[[method]]
      }
    }
    trun.pts.idx <- -log10(pvals.pts) > trun
    if(any(trun.pts.idx)){
       pvals.pts[trun.pts.idx] <- 10^(-trun)
    }
    if(missing(ylab)) ylab <- expression(-log[10](italic(P)))
    if(missing(xlim)) xlim <- xlimit
    if(missing(ylim)) ylim <- c(0, max(-log10(pvals.pts)))
    if(missing(col.sig)) col.sig <- c(adjustcolor("forestgreen", alpha.f = 0.1),
      adjustcolor("firebrick", alpha.f = 0.1))
    if(length(col.sig) != 2){
      stop("col.sig must be a vector of two elements of colors.")
    }
    plot.default(mean(xlim), -log10(alpha),
      xlab = xlab, xlim = xlim, ylab = ylab, ylim = ylim, type = "n", ...)
    rect(xleft = xlimit[1] - 100*sum(x$data), ybottom = -1,
      xright = xlimit[2] + 100*sum(x$data), ytop = -log10(alpha),
      col = col.sig[1], border = "white")
    if(any(pvals.pts < alpha)){
      rect(xleft = xlimit[1] - 100*sum(x$data), ybottom = -log10(alpha),
        xright = xlimit[2] + 100*sum(x$data), ytop = max(-log10(pvals.pts))*10 + 10,
        col = col.sig[2], border = "white")
    }
    abline(v = 0, lty = lty.ori, col = col.ori, lwd = lwd.ori)
    if(missing(pch)) pch <- 1
    if(missing(pch.trun)) pch.trun <- 3
    pchs <- rep(pch, length(xpts))
    pchs[xpts == 0] <- pch.ori
    pchs[which(xpts %in% mods.single)] <- pch.f
    if(any(trun.pts.idx)){
      pchs[trun.pts.idx] <- pch.trun
    }
    cexs <- rep(cex.pts, length(xpts))
    cexs[xpts == 0] <- 2 * cex.pts
    cexs[which(xpts %in% mods.single)] <- 2 * cex.pts
    points(xpts, -log10(pvals.pts), cex = cexs, pch = pchs)
  }

  if(modify0 + modify1 == 2){
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    if(all(f0.range == 0) | all(f1.range == 0) ){
      stop("At least one group is not modified; no plot is generated.")
    }
    if(missing(adjust.legend)) adjust.legend <- 1
    layout(matrix(c(1, 2), ncol = 2, byrow = TRUE),
      widths = c(10, adjust.legend))
    xpts <- f0.range[1]:f0.range[2]
    ypts <- f1.range[1]:f1.range[2]
    xlimit <- f0.range
    if(missing(xlab)) xlab <- "Number of modified events in group 0"
    if(missing(xlim)) xlim <- xlimit
    ylimit <- f1.range
    if(missing(ylab)) ylab <- "Number of modified events in group 1"
    if(missing(ylim)) ylim <- ylimit
    if(missing(col.sig)) col.sig <- c("forestgreen", "firebrick")
    if(length(col.sig) != 2){
      stop("col.sig must be a vector of two elements of colors.")
    }
    if(missing(pch)) pch <- 16
    if(missing(pch.ori.hl)) pch.ori.hl <- 0
    if(missing(pch.f.hl)) pch.f.hl <- 2
    if(missing(col.ori.hl)) col.ori.hl <- "black"
    if(missing(col.f.hl)) col.f.hl <- "black"
    plot.default(mean(xlimit), mean(ylimit),
      xlab = xlab, xlim = xlim, ylab = ylab, ylim = ylim, type = "n", ...)
    ori.idx <- tot.mods == 0
    FI.idx <- matrix(0, dim(pvals)[1], dim(pvals)[2])
    for(i in 1:dim(mods)[1]){
      FI.idx[xpts == mods[i,1], ypts == mods[i,2]] <- 1
    }
    signif <- c(pvals) < alpha
    neg.log10.pvals <- -log10(c(pvals))
    neg.log10.pvals.trun <- neg.log10.pvals
    neg.log10.pvals.trun[neg.log10.pvals.trun > trun] <- trun
    pvals.col.trans <- 0.1 +
      abs(neg.log10.pvals.trun + log10(alpha))/
      abs((c(signif)*trun + log10(alpha)))*(1 - 0.1)
    pvals.col <- ifelse(c(signif), col.sig[2], col.sig[1])
    adjustcolor.vec <- Vectorize(adjustcolor)
    pchs <- rep(pch, length(c(pvals)))
    cexs <- rep(cex.pts, length(c(pvals)))
    pchs[c(ori.idx) == 1] <- pch.ori
    pchs[c(FI.idx) == 1] <- pch.f
    cexs[c(ori.idx) == 1 | c(FI.idx) == 1] <- 2 * cex.pts
    points(x = rep(xpts, length(ypts)),
      y = rep(ypts, each = length(xpts)),
      col = adjustcolor.vec(pvals.col, alpha.f = pvals.col.trans),
      cex = cexs, pch = pchs)
    abline(h = 0, lty = lty.ori, col = col.ori, lwd = lwd.ori)
    abline(v = 0, lty = lty.ori, col = col.ori, lwd = lwd.ori)
    points(x = rep(xpts, length(ypts))[c(ori.idx) == 1],
      y = rep(ypts, each = length(xpts))[c(ori.idx) == 1],
      col = col.ori.hl, cex = 2 * cex.pts, pch = pch.ori.hl)
    points(x = rep(xpts, length(ypts))[c(FI.idx) == 1],
      y = rep(ypts, each = length(xpts))[c(FI.idx) == 1],
      col = col.f.hl, cex = 2 * cex.pts, pch = pch.f.hl)

    mar.temp <- op$mar
    mar.temp[2] <- 0
    par(mar = mar.temp)
    plot.default(1, type = "n", xlim = c(0, trun/10*adjust.legend),
      ylim = c(0, trun), axes = FALSE, frame.plot = FALSE,
      xlab = "", ylab = "")
    if(missing(adjust.seg)) adjust.seg <- 10
    if(adjust.seg %% 1 != 0) adjust.seg <- round(adjust.seg)
    if(missing(cex.legend.pval)) cex.legend.pval <- 0.6
    if(missing(legend.pvals)) legend.pvals <- NULL
    if(missing(cex.legend.title)) cex.legend.title <- 1
    seg2 <- round((trun + log10(alpha))/(-log10(alpha)/adjust.seg))
    ybs <- c(seq(0, -log10(alpha) * (1 - 1/adjust.seg),
      length.out = adjust.seg), seq(-log10(alpha),
      trun - (trun + log10(alpha))/seg2, length.out = seg2))
    yus <- c(seq(-log10(alpha)/adjust.seg, -log10(alpha),
      length.out = adjust.seg), seq(-log10(alpha) +
      (trun + log10(alpha))/seg2, trun, length.out = seg2))
    cols <- c(adjustcolor.vec(col.sig[1], alpha.f =
      seq(1, 0.1, length.out = adjust.seg)),
      adjustcolor.vec(col.sig[2], alpha.f =
      seq(0.1, 1, length.out = seg2)))
    rect(xleft = 0, ybottom = ybs, xright = trun/10*adjust.legend,
      ytop = yus, col = cols, border = NA)
    abline(h = -log10(alpha), col = "white")
    abline(h = trun, col = "white")
    mtext(text = expression(italic(P)~value), side = 3, line = 0,
      cex = cex.legend.title)
    mtext(text = format(alpha, scientific = FALSE),
      at = -log10(alpha), side = 4, line = 0,
      cex = cex.legend.pval)
    mtext(text = bquote(""<=10^-.(trun)), at = trun, side = 4, line = 0,
      cex = cex.legend.pval)
    if(!is.null(legend.pvals)){
      mtext(text = legend.pvals, at = -log10(legend.pvals),
        side = 4, line = 0, cex = cex.legend.pval)
      abline(h = -log10(legend.pvals), col = "white")
    }
  }
}