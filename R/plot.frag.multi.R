plot.frag.multi <- function(x, method, dir = "both", fragility = "FI",
  percentage = TRUE, max.f = NULL, bar, names.arg, space = 0,
  breaks, freq, reverse = FALSE, xlab, ylab, main = NULL,
  cex.marker, col.border, col.sig, trun.marker = TRUE, ...){
  if(!inherits(x, "frag.multi")){
    stop("The input must be an object of class \"frag.multi\".")
  }
  if(inherits(x, "frag.studies")){
    if(missing(method)){
      method <- x$methods[1]
    }
    if(!is.element(method, x$methods)){
      stop("method must be used to produce the object of class \"frag.multi\".")
    }
    if(length(method) != 1){
      stop("Only one method can be specified for plotting each time.")
    }
  }
  if(!is.element(dir, c("both", "sig2nonsig", "nonsig2sig")) |
    length(dir) != 1){
    stop("dir must be \"both\", \"sig2nonsig\", or \"nonsig2sig\".")
  }
  if(!is.element(fragility, c("FI", "FQ")) | length(fragility) != 1){
    stop("fragility must be either \"FI\" or \"FQ\".")
  }
  if(!is.null(max.f)){
    if(length(max.f) != 1) stop("max.f must be a single value.")
    if(fragility == "FI" & max.f%%1 != 0) max.f <- floor(max.f)
  }
  if(missing(bar)){
    if(fragility == "FI") bar <- TRUE else bar <- FALSE
  }
  if(bar & fragility == "FQ"){
    stop("Bar plot cannot be generated for fragility quotient.")
  }
  if(missing(breaks)) breaks <- "Sturges"
  if(missing(freq)) freq <- NULL
  if(missing(cex.marker)) cex.marker <- 1
  if(missing(col.sig)) col.sig <- c(adjustcolor("forestgreen", alpha.f = 0.1),
    adjustcolor("firebrick", alpha.f = 0.1))
  if(length(col.sig) != 2){
    stop("col.sig must be a vector of two elements of colors.")
  }

  if(inherits(x, "frag.studies")){
    if(fragility == "FI"){
      pts <- x$FI[, method]
    }else{
      pts <- x$FQ[, method]
      if(percentage) pts <- pts * 100
    }
    sig2nonsig <- x$pval[,method] < x$alpha
  }
  if(inherits(x, "frag.mas")){
    if(fragility == "FI"){
      pts <- x$FI
    }else{
      pts <- x$FQ
      if(percentage) pts <- pts * 100
    }
    sig2nonsig <- x$pval.ori < x$alpha
  }
  if(any(is.na(pts))){
    message("Note that fragility measures of NA are not plotted.")
    sig2nonsig <- sig2nonsig[!is.na(pts)]
    pts <- pts[!is.na(pts)]
  }
  if(length(pts) == 0) stop("All fragility measures are NA.")

  if(all(sig2nonsig)){
    if(dir == "both") dir <- "sig2nonsig"
    if(dir == "nonsig2sig"){
      stop("No study is modified from non-significance to significance.")
    }
  }
  if(all(!sig2nonsig)){
    if(dir == "both") dir <- "nonsig2sig"
    if(dir == "sig2nonsig"){
      stop("No study is modified from significance to non-significance.")
    }
  }
  if(dir == "sig2nonsig") pts <- pts[sig2nonsig]
  if(dir == "nonsig2sig") pts <- pts[!sig2nonsig]

  trun.idx <- trun.lab <- NULL
  if(!is.null(max.f)){
    if(all(pts >= max.f)){
      stop("max.f must be larger than some fragility index/quotient.")
    }
    if(all(pts <= max.f)){
      max.f <- max(pts)
    }
    if(any(pts > max.f)){
      trun.idx <- which(pts > max.f)
      trun.lab <- paste0(">", max.f)
    }
  }
  if(is.null(max.f)){
    max.f <- max(pts)
  }
  if(length(trun.idx) != 0){
    if(bar){
      pts.trun <- max.f + 1
      pts[trun.idx] <- pts.trun
    }
    if(!bar){
      temp.hist <- hist(pts[-trun.idx], breaks = breaks, plot = FALSE)
      pts.trun <- max.f + mean(diff(temp.hist$breaks))
      pts[trun.idx] <- pts.trun
    }
  }

  if(bar){
    xbar <- xbar.name <- min(pts):max.f
    if(length(trun.idx) != 0){
      xbar <- c(xbar, pts.trun)
      xbar.name <- c(xbar.name, trun.lab)
    }
    bar.pts <- matrix(NA, nrow = ifelse(dir == "both", 2, 1),
      ncol = length(xbar))
    for(bar.idx in 1:length(xbar)){
      if(dir == "both"){
        bar.pts[1, bar.idx] <- sum(pts[!sig2nonsig] == xbar[bar.idx])
        bar.pts[2, bar.idx] <- sum(pts[sig2nonsig] == xbar[bar.idx])
      }else{
        bar.pts[1, bar.idx] <- sum(pts == xbar[bar.idx])
      }
    }
    if(dir == "both"){
      col.bp <- col.sig
      if(reverse){
        bar.pts <- bar.pts[c(2,1),]
        col.bp <- col.bp[c(2,1)]
      }
    }else{
      col.bp <- ifelse(dir == "sig2nonsig", col.sig[2], col.sig[1])
    }
    if(missing(xlab)) xlab <- "Fragility index"
    if(missing(ylab)) ylab <- "Frequency"
    if(missing(names.arg)) names.arg <- xbar.name
    if(missing(col.border)) col.border <- par("fg")
    barplot(bar.pts, names.arg = names.arg, space = space, col = col.bp,
      xlab = xlab, ylab = ylab, border = col.border, main = main, ...)
  }

  if(!bar){
    if(missing(xlab)){
      if(fragility == "FI") xlab <- "Fragility index"
      if(fragility == "FQ"){
        xlab <- "Fragility quotient"
        if(percentage) xlab <- paste0(xlab, ", %")
      }
    }
    if(missing(col.border)) col.border <- NULL
    if(dir == "both"){
      col.hist <- col.sig
      if(reverse){
        col.hist <- col.hist[c(2,1)]
        sig2nonsig <- !sig2nonsig
      }
      hp <- hist(pts, breaks = breaks, freq = freq, col = col.hist[2],
        xlab = xlab, main = main, border = col.border, ...)
      if(is.null(freq)) freq <- hp$equidist
      pts.add <- pts[!sig2nonsig]
      if(freq){
        hist(pts.add, breaks = hp$breaks, freq = freq,
          col = "white", lty = "blank", add = TRUE)
        hist(pts.add, breaks = hp$breaks, freq = freq,
          col = col.hist[1], border = col.border, add = TRUE)
      }else{
        hp.add <- hist(pts.add, breaks = hp$breaks, plot = FALSE)
        for(k in 1:length(hp.add$counts)){
          if(hp.add$counts[k] > 0){
            rect(xleft = hp.add$breaks[k], ybottom = 0,
              xright = hp.add$breaks[k + 1],
              ytop = hp$density[k] * hp.add$counts[k]/hp$counts[k],
              col = "white", border = NA)
            rect(xleft = hp.add$breaks[k], ybottom = 0,
              xright = hp.add$breaks[k + 1],
              ytop = hp$density[k] * hp.add$counts[k]/hp$counts[k],
              col = col.hist[1], border = col.border)
          }
        }
      }
    }else{
      col.hist <- ifelse(dir == "sig2nonsig", col.sig[2], col.sig[1])
      hp <- hist(pts, breaks = breaks, freq = freq, col = col.hist,
        xlab = xlab, border = col.border, main = main, ...)
      if(is.null(freq)) freq <- hp$equidist
    }
    if(trun.marker & !is.null(trun.lab)){
      last <- length(hp$mids)
      text(x = hp$mids[last], y = ifelse(freq, hp$counts[last], hp$density[last]),
        label = trun.lab, pos = 3, cex = cex.marker)
    }
  }
}