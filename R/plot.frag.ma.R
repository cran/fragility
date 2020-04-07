plot.frag.ma <- function(x, xlab, ylab, xlim, ylim, ybreaks = NULL,
  study.marker = TRUE, cex.marker, offset.marker, col.line, lwd,
  legend, x.legend, y.legend, cex.legend, ...){
  if(!inherits(x, "frag.ma")){
    stop("The input must be an object of class \"frag.ma\".")
  }
  if(!is.null(ybreaks)){
    if(is.vector(ybreaks)){
      if(length(ybreaks) != 2) stop("The length of ybreaks should be 2.")
      if(ybreaks[1] >= ybreaks[2]) stop("ybreaks[1] should be smaller than ybreaks[2].")
    }else{
      stop("ybreaks should be a vector.")
    }
  }

  if(missing(col.line)) col.line <- c("blue", "red")
  if(length(col.line) != 2){
    stop("col.line must be a vector of two elements of colors.")
  }

  xpts <- 0:length(x$sid.iter)
  cts <- c(sum(x$data$e0), sum(x$data$e1), sum(x$data.mod$e0), sum(x$data.mod$e1))
  ypts <- matrix(NA, nrow = length(xpts), ncol = 2)
  g0.mod <- c(0, x$g0.mod.iter)
  g0.mod.cum <- cumsum(g0.mod)
  g1.mod <- c(0, x$g1.mod.iter)
  g1.mod.cum <- cumsum(g1.mod)
  ypts <- cbind(rep(sum(x$data$e0), length(xpts)) + g0.mod.cum,
    rep(sum(x$data$e1), length(xpts)) + g1.mod.cum)
  if(!is.null(ybreaks)){
    if(ybreaks[1] <= min(max(ypts[,1]), max(ypts[,2])) | ybreaks[2] >= max(min(ypts[,1]), min(ypts[,2]))){
      stop("ybreaks should be between the total event counts in groups 0 and 1.")
    }
    ypts.ori <- ypts
    if(all(ypts[,1] < ypts[,2])){
      ypts[,2] <- ypts[,2] - diff(ybreaks)
      g.breaks <- 1
    }else{
      ypts[,1] <- ypts[,1] - diff(ybreaks)
      g.breaks <- 0
    }
  }

  if(missing(xlab)) xlab <- "Iteration of event status modification"
  if(missing(xlim)) xlim <- c(min(xpts), max(xpts))
  if(missing(ylab)) ylab <- "Group-specific total event count"
  if(missing(ylim)) ylim <- c(min(ypts), max(ypts))

  if(missing(lwd)) lwd <- 1
  if(missing(cex.marker)) cex.marker <- 0.8
  if(missing(offset.marker)) offset.marker <- 0.2
  if(study.marker){
    sid.marker <- as.character(x$sid.iter)
    if(length(sid.marker) > 1){
      for(i in 2:length(sid.marker)){
        if(x$sid.iter[i] == x$sid.iter[i-1] &
          (all(x$g0.mod.iter[(i-1):i] == 0) | all(x$g1.mod.iter[(i-1):i] == 0))) sid.marker[i] <- "*"
      }
    }
  }

  if(is.null(ybreaks)){
    plot.default(mean(xpts), mean(ypts[,1]),
      xlab = xlab, xlim = xlim, ylab = ylab, ylim = ylim, type = "n", xaxt = "n", ...)
  }else{
    plot.default(mean(xpts), mean(ypts[,1]),
      xlab = xlab, xlim = xlim, ylab = ylab, ylim = ylim, type = "n", xaxt = "n", yaxt = "n", ...)
    yt <- axTicks(2)
    yt.lab <- yt
    yt.lab[yt.lab > ybreaks[1]] <- yt.lab[yt.lab > ybreaks[1]] + diff(ybreaks)
    axis(2, at = yt, labels = yt.lab)
    axis.break(axis = 2, breakpos = ybreaks[1])
  }
  xt <- axTicks(1)
  xt <- xt[xt%%1 == 0]
  axis(1, at = xt)
  for(i in 1:(length(xpts) - 1)){
    lines(c(xpts[i], xpts[i+1]), rep(ypts[i,1], 2), col = col.line[1], lwd = lwd)
    lines(c(xpts[i], xpts[i+1]), rep(ypts[i,2], 2), col = col.line[2], lwd = lwd)
    if(g0.mod[i+1] != 0){
      lines(rep(xpts[i+1], 2), c(ypts[i,1], ypts[i+1,1]), col = col.line[1], lwd = lwd)
      if(study.marker) text(x = xpts[i+1], y = ypts[i+1,1], labels = sid.marker[i],
                         pos = ifelse(g0.mod[i+1] > 0, 3, 1), cex = cex.marker, offset = offset.marker)
    }
    if(g1.mod[i+1] != 0){
      lines(rep(xpts[i+1], 2), c(ypts[i,2], ypts[i+1,2]), col = col.line[2], lwd = lwd)
      if(study.marker) text(x = xpts[i+1], y = ypts[i+1,2], labels = sid.marker[i],
                         pos = ifelse(g1.mod[i+1] > 0, 3, 1), cex = cex.marker, offset = offset.marker)
    }
  }
  if(missing(legend)) legend <- paste("Group", 0:1)
  if(missing(x.legend)) x.legend <- "right"
  if(missing(y.legend)) y.legend <- NULL
  if(missing(cex.legend)) cex.legend <- 1
  legend(x.legend, y.legend, legend = legend, col = col.line, lwd = lwd)
}