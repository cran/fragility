plot.frag.nma <- function(x, tid1, tid2, xlab, ylab, xlim, ylim, ybreaks = NULL,
  study.marker = TRUE, cex.marker, offset.marker, col.line, lwd,
  legend, x.legend, y.legend, cex.legend, ...){
  if(!inherits(x, "frag.nma")){
    stop("The input must be an object of class \"frag.nma\".")
  }
  tids <- rownames(x$FI)
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

  xpts <- 0:length(x$sid.iter[[paste(tid1, "vs.", tid2)]])
  data.ori <- x$data
  data.mod <- x$data.mod[[paste(tid1, "vs.", tid2)]]
  cts <- c(sum(data.ori$e[data.ori$tid == tid1]), sum(data.ori$e[data.ori$tid == tid2]),
    sum(data.mod$e[data.mod$tid == tid1]), sum(data.mod$e[data.mod$tid == tid2]))
  ypts <- matrix(NA, nrow = length(xpts), ncol = 2)
  tid1.mod <- c(0, x$tid1.mod.iter[[paste(tid1, "vs.", tid2)]])
  tid1.mod.cum <- cumsum(tid1.mod)
  tid2.mod <- c(0, x$tid2.mod.iter[[paste(tid1, "vs.", tid2)]])
  tid2.mod.cum <- cumsum(tid2.mod)
  ypts <- cbind(rep(sum(data.ori$e[data.ori$tid == tid1]), length(xpts)) + tid1.mod.cum,
    rep(sum(data.ori$e[data.ori$tid == tid2]), length(xpts)) + tid2.mod.cum)
  if(!is.null(ybreaks)){
    if(ybreaks[1] <= min(max(ypts[,1]), max(ypts[,2])) | ybreaks[2] >= max(min(ypts[,1]), min(ypts[,2]))){
      stop("ybreaks should be between the total event counts in the two groups.")
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
    sid.marker <- as.character(x$sid.iter[[paste(tid1, "vs.", tid2)]])
    if(length(sid.marker) > 1){
      for(i in 2:length(sid.marker)){
        if(x$sid.iter[[paste(tid1, "vs.", tid2)]][i] == x$sid.iter[[paste(tid1, "vs.", tid2)]][i-1] &
          (all(x$tid1.mod.iter[[paste(tid1, "vs.", tid2)]][(i-1):i] == 0) |
          all(x$tid2.mod.iter[[paste(tid1, "vs.", tid2)]][(i-1):i] == 0))) sid.marker[i] <- "*"
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
    if(tid1.mod[i+1] != 0){
      lines(rep(xpts[i+1], 2), c(ypts[i,1], ypts[i+1,1]), col = col.line[1], lwd = lwd)
      if(study.marker) text(x = xpts[i+1], y = ypts[i+1,1], labels = sid.marker[i],
                         pos = ifelse(tid1.mod[i+1] > 0, 3, 1), cex = cex.marker, offset = offset.marker)
    }
    if(tid2.mod[i+1] != 0){
      lines(rep(xpts[i+1], 2), c(ypts[i,2], ypts[i+1,2]), col = col.line[2], lwd = lwd)
      if(study.marker) text(x = xpts[i+1], y = ypts[i+1,2], labels = sid.marker[i],
                         pos = ifelse(tid2.mod[i+1] > 0, 3, 1), cex = cex.marker, offset = offset.marker)
    }
  }
  if(tid1 < tid2) od <- 1:2 else od <- 2:1
  if(missing(legend)) legend <- paste("Treatment", c(tid1, tid2)[od])
  if(missing(x.legend)) x.legend <- "right"
  if(missing(y.legend)) y.legend <- NULL
  if(missing(cex.legend)) cex.legend <- 1
  legend(x.legend, y.legend, legend = legend, col = col.line[od], lwd = lwd)
}