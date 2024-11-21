frag.nma.alpha <- function(sid, tid, e, n, data, measure = "OR", random = TRUE,
  alpha.from = 0.005, alpha.to = 0.05, alpha.breaks = 10,
  mod.dir = "both", tid1.f, tid2.f, OR = 1, RR = 1, RD = 0,
  incr, allincr, addincr, allstudies, ...){
  if(!missing(data)){
    sid <- eval(substitute(sid), data, parent.frame())
    tid <- eval(substitute(tid), data, parent.frame())
    e <- eval(substitute(e), data, parent.frame())
    n <- eval(substitute(n), data, parent.frame())
  }
  if(any(e < 0) | any(n < 0)){
    stop("Event counts and sample sizes must be nonnegative.")
  }
  if(any(e > n)){
    stop("Event counts must not be larger than sample sizes.")
  }
  if(any(e%%1 != 0) | any(n%%1 != 0)){
    message("Some event counts and/or sample sizes are not integer(s); they are rounded.")
    e <- round(e)
    n <- round(n)
  }
  if(length(n) != length(e) | length(sid) != length(e) | length(tid) != length(e)){
    stop("sid, tid, e, and n do not have the same length.")
  }
  if(length(unique(sid)) == 1){
    stop("A network meta-analysis should contain more than one study.")
  }
  if(any(!is.element(measure, c("OR", "RR", "RD")))){
    stop("measure must be \"OR\", \"RR\", or \"RD\".")
  }
  if(length(measure) != 1){
    stop("Only one measure can be specified.")
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
  if(any(!is.element(mod.dir, c("both", "one", "left", "right")))){
    stop("mod.dir must be \"both\", \"one\", \"left\", or \"right\".")
  }
  if(OR <= 0) stop("The null value of OR must be positive.")
  if(RR <= 0) stop("The null value of RR must be positive.")
  if(RD < -1 | RD > 1) stop("The null value of RD must be between -1 and 1.")

  if(missing(incr)) incr <- 0.5
  if(missing(allincr)) allincr <- FALSE
  if(missing(addincr)) addincr <- FALSE
  if(missing(allstudies)) allstudies <- FALSE
  model <- ifelse(random, "random", "fixed")

  ori.data <- data.frame(sid = sid, tid = tid, e = e, n = n)
  if(measure == "OR") null.val <- log(OR)
  if(measure == "RR") null.val <- log(RR)
  if(measure == "RD") null.val <- RD
  alphas <- seq(from = alpha.from, to = alpha.to,
    length.out = alpha.breaks)

  pmeta.ori <- meta::pairwise(treat = tid, event = e, n = n, studlab = sid, data = ori.data, sm = measure,
    incr = incr, allincr = allincr, addincr = addincr, allstudies = allstudies)
  rslt.ori <- netmeta(TE = pmeta.ori$TE, seTE = pmeta.ori$seTE, treat1 = pmeta.ori$treat1, treat2 = pmeta.ori$treat2, studlab = pmeta.ori$studlab,
    sm = measure, level.ma = 0.95, ...)
  est.ori <- rslt.ori[[paste0("TE.", model)]]
  se.ori <- rslt.ori[[paste0("seTE.", model)]]
  pval.ori <- rslt.ori[[paste0("pval.", model)]]

  if(is.vector(mod.dir)){
    if(any(!is.element(mod.dir, c("both", "one", "left", "right")))){
      stop("mod.dir should be \"both\", \"one\", \"left\", or \"right\".")
    }
    if(length(mod.dir) != 1){
      stop("mod.dir should be a character string or a matrix.")
    }
    if(is.element(mod.dir, c("both", "left", "right"))){
      mod.dir <- matrix(mod.dir, dim(est.ori)[1], dim(est.ori)[2])
    }else{
      mod.dir <- matrix("left", dim(est.ori)[1], dim(est.ori)[2])
      for(i in 1:dim(est.ori)[1]){
        for(j in 1:dim(est.ori)[2]){
          if(est.ori[i, j] >= null.val) mod.dir[i, j] <- "right"
        }
      }
    }
  }else{
    if(is.matrix(mod.dir)){
      if(all(dim(mod.dir) == dim(est.ori))) stop("the dimension of mod.dir is inconsistent with treatments.")
      for(i in dim(mod.dir)[1]){
        for(j in 1:dim(mod.dir)[2]){
          if(!is.element(mod.dir[i, j], c("both", "one", "left", "right"))){
            stop("each element of mod.dir should be \"both\", \"one\", \"left\", or \"right\".")
          }
        }
      }
    }else{
      stop("mod.dir should be a character string or a matrix.")
    }
  }
  if((missing(tid1.f) & !missing(tid2.f)) | (!missing(tid1.f) & missing(tid2.f))){
    stop("Both tid1.f and tid2.f should be specified.")
  }
  if(missing(tid1.f) & missing(tid2.f)){
    tid1.f <- rep(rslt.ori$trts, length(rslt.ori$trts))
    tid2.f <- rep(rslt.ori$trts, each = length(rslt.ori$trts))
  }else{
    if(any(!is.element(tid1.f, rslt.ori$trts)) | any(!is.element(tid2.f, rslt.ori$trts))){
      stop("tid1.f and tid2.f should be treatment IDs.")
    }
    if(length(tid1.f) != length(tid2.f)){
      stop("tid1.f and tid2.f do not have the same length.")
    }
  }
  rownames(mod.dir) <- colnames(mod.dir) <- rslt.ori$trts

  out <- list(data = ori.data, measure = measure, alphas = alphas,
    null = null.val, est.ori = est.ori, se.ori = se.ori,
    pval.ori = pval.ori, mod.dir = mod.dir, tid.f = paste(tid1.f, "vs.", tid2.f))

  empt.list <- rep(list(NULL), length(tid1.f))
  names(empt.list) <- paste(tid1.f, "vs.", tid2.f)
  FI <- FQ <- FQ.nma <- FI.avg <- FQ.avg <- FQ.nma.avg <- empt.list

  judge.mod <- function(dat, incr2, tid1, tid2){
    dat1 <- dat[dat$tid == tid1,]
    dat2 <- dat[dat$tid == tid2,]
    if(incr2){
      sid1 <- dat1$sid[dat1$e > 0]
      sid2 <- dat2$sid[dat2$e < dat2$n]
    }else{
      sid1 <- dat1$sid[dat1$e < dat1$n]
      sid2 <- dat2$sid[dat2$e > 0]
    }
    out <- list(sid1 = sid1, sid2 = sid2)
    return(out)
  }

  for(k in 1:length(tid1.f)){
    tid1 <- tid1.f[k]
    tid2 <- tid2.f[k]
    if(tid1 == tid2) next
    if(k > 1){
      if(is.element(paste(tid2, "vs.", tid1), paste(tid1.f[1:(k - 1)], "vs.", tid2.f[1:(k - 1)])) &
         is.element(mod.dir[tid1, tid2], c("both", "one"))){
        FI[[paste(tid1, "vs.", tid2)]] <- FI[[paste(tid2, "vs.", tid1)]]
        FQ[[paste(tid1, "vs.", tid2)]] <- FQ[[paste(tid2, "vs.", tid1)]]
        FQ.nma[[paste(tid1, "vs.", tid2)]] <- FQ.nma[[paste(tid2, "vs.", tid1)]]
        FI.avg[[paste(tid1, "vs.", tid2)]] <- FI.avg[[paste(tid2, "vs.", tid1)]]
        FQ.avg[[paste(tid1, "vs.", tid2)]] <- FQ.avg[[paste(tid2, "vs.", tid1)]]
        FQ.nma.avg[[paste(tid1, "vs.", tid2)]] <- FQ.nma.avg[[paste(tid2, "vs.", tid1)]]
        next
      }
    }

    FI[[paste(tid1, "vs.", tid2)]] <- FQ[[paste(tid1, "vs.", tid2)]] <-
      FQ.nma[[paste(tid1, "vs.", tid2)]] <- rep(NA, alpha.breaks)

    signif.alpha <- (alphas > pval.ori[tid1, tid2])

    if(any(signif.alpha)){
      alphas.signif <- alphas[signif.alpha]
      FI.signif <- FQ.signif <- FQ.nma.signif <- rep(NA, sum(signif.alpha))
      less <- (est.ori[tid1, tid2] < null.val)
      for(a in 1:length(alphas.signif)){
        sid.iter <- NULL
        data.temp <- ori.data
        signif <- TRUE
        judge.temp <- judge.mod(dat = data.temp, incr2 = ifelse(less, FALSE, TRUE), tid1, tid2)
        sids.temp <- c(judge.temp$sid1, judge.temp$sid2)
        grp.idx.temp <- NULL
        if(length(judge.temp$sid1) > 0) grp.idx.temp <- c(grp.idx.temp, rep(tid1, length(judge.temp$sid1)))
        if(length(judge.temp$sid2) > 0) grp.idx.temp <- c(grp.idx.temp, rep(tid2, length(judge.temp$sid2)))
        moremod <- (length(sids.temp) > 0)
        while(signif & moremod){
          sid.temp <- tid1.mod.temp <- tid2.mod.temp <-
             ci.lb.temp <- ci.ub.temp <- NULL
          data.tt <- NULL
          for(i in 1:length(sids.temp)){
            data.tt[[i]] <- data.temp
            if(grp.idx.temp[i] == tid1){
              data.tt[[i]]$e[data.tt[[i]]$sid == sids.temp[i] & data.tt[[i]]$tid == tid1] <-
                data.tt[[i]]$e[data.tt[[i]]$sid == sids.temp[i] & data.tt[[i]]$tid == tid1] + ifelse(less, 1, -1)
              tid1.mod.temp <- c(tid1.mod.temp, ifelse(less, 1, -1))
              tid2.mod.temp <- c(tid2.mod.temp, 0)
            }else{
              data.tt[[i]]$e[data.tt[[i]]$sid == sids.temp[i] & data.tt[[i]]$tid == tid2] <-
                data.tt[[i]]$e[data.tt[[i]]$sid == sids.temp[i] & data.tt[[i]]$tid == tid2] + ifelse(less, -1, 1)
              tid1.mod.temp <- c(tid1.mod.temp, 0)
              tid2.mod.temp <- c(tid2.mod.temp, ifelse(less, -1, 1))
            }
            pmeta.temp <- meta::pairwise(treat = tid, event = e, n = n, studlab = sid, data = data.tt[[i]], sm = measure,
              incr = incr, allincr = allincr, addincr = addincr, allstudies = allstudies)
            out.temp <- netmeta(TE = pmeta.temp$TE, seTE = pmeta.temp$seTE, treat1 = pmeta.temp$treat1, treat2 = pmeta.temp$treat2, studlab = pmeta.temp$studlab,
              sm = measure, level.ma = 1 - alphas.signif[a], ...)
            sid.temp <- c(sid.temp, sids.temp[i])
            ci.lb.temp <- c(ci.lb.temp, out.temp[[paste0("lower.", model)]][tid1, tid2])
            ci.ub.temp <- c(ci.ub.temp, out.temp[[paste0("upper.", model)]][tid1, tid2])
            if((less & (out.temp[[paste0("upper.", model)]][tid1, tid2] >= null.val)) |
              (!less & (out.temp[[paste0("lower.", model)]][tid1, tid2] <= null.val))){
              signif <- FALSE
              break
            }
          }
          if(less){
            idx.iter <- which(ci.ub.temp == max(ci.ub.temp))[1]
          }else{
            idx.iter <- which(ci.lb.temp == min(ci.lb.temp))[1]
          }
          sid.iter <- c(sid.iter, sid.temp[idx.iter])
          data.temp <- data.tt[[idx.iter]]
          if(signif){
            judge.temp <- judge.mod(dat = data.temp, incr2 = ifelse(less, FALSE, TRUE), tid1, tid2)
            sids.temp <- c(judge.temp$sid1, judge.temp$sid2)
            grp.idx.temp <- NULL
            if(length(judge.temp$sid1) > 0) grp.idx.temp <- c(grp.idx.temp, rep(tid1, length(judge.temp$sid1)))
            if(length(judge.temp$sid2) > 0) grp.idx.temp <- c(grp.idx.temp, rep(tid2, length(judge.temp$sid2)))
            moremod <- (length(sids.temp) > 0)
          }
        }
        if(!signif){
          FI.signif[a] <- length(sid.iter)
          FQ.signif[a] <- length(sid.iter)/sum(ori.data$n[ori.data$tid == tid1 | ori.data$tid == tid2])
          FQ.nma.signif[a] <- length(sid.iter)/sum(ori.data$n)
        }
      }
    }

    if(any(!signif.alpha)){
      alphas.nonsignif <- alphas[!signif.alpha]
      FI.nonsignif <- FQ.nonsignif <- FQ.nma.nonsignif <- rep(NA, sum(!signif.alpha))
      left <- right <- TRUE
      if(mod.dir[tid1, tid2] == "left") right <- FALSE
      if(mod.dir[tid1, tid2] == "right") left <- FALSE

      left.fcn <- function(max.iter = Inf, a, ...){
        FI.nonsignif.left <- FQ.nonsignif.left <- FQ.nma.nonsignif.left <- NA
        sid.iter.left <- NULL
        data.temp <- ori.data
        signif <- FALSE
        judge.temp <- judge.mod(dat = data.temp, incr2 = TRUE, tid1, tid2)
        sids.temp <- c(judge.temp$sid1, judge.temp$sid2)
        grp.idx.temp <- NULL
        if(length(judge.temp$sid1) > 0) grp.idx.temp <- c(grp.idx.temp, rep(tid1, length(judge.temp$sid1)))
        if(length(judge.temp$sid2) > 0) grp.idx.temp <- c(grp.idx.temp, rep(tid2, length(judge.temp$sid2)))
        moremod <- (length(sids.temp) > 0)
        while(!signif & moremod & length(sid.iter.left) <= max.iter){
          sid.temp <- tid1.mod.temp <- tid2.mod.temp <-
           ci.lb.temp <- ci.ub.temp <- NULL
          data.tt <- NULL
          for(i in 1:length(sids.temp)){
            data.tt[[i]] <- data.temp
            if(grp.idx.temp[i] == tid1){
              data.tt[[i]]$e[data.tt[[i]]$sid == sids.temp[i] & data.tt[[i]]$tid == tid1] <-
                data.tt[[i]]$e[data.tt[[i]]$sid == sids.temp[i] & data.tt[[i]]$tid == tid1] - 1
              tid1.mod.temp <- c(tid1.mod.temp, -1)
              tid2.mod.temp <- c(tid2.mod.temp, 0)
            }else{
              data.tt[[i]]$e[data.tt[[i]]$sid == sids.temp[i] & data.tt[[i]]$tid == tid2] <-
                data.tt[[i]]$e[data.tt[[i]]$sid == sids.temp[i] & data.tt[[i]]$tid == tid2] + 1
              tid1.mod.temp <- c(tid1.mod.temp, 0)
              tid2.mod.temp <- c(tid2.mod.temp, 1)
            }
            pmeta.temp <- meta::pairwise(treat = tid, event = e, n = n, studlab = sid, data = data.tt[[i]], sm = measure,
              incr = incr, allincr = allincr, addincr = addincr, allstudies = allstudies)
            out.temp <- netmeta(TE = pmeta.temp$TE, seTE = pmeta.temp$seTE, treat1 = pmeta.temp$treat1, treat2 = pmeta.temp$treat2, studlab = pmeta.temp$studlab,
              sm = measure, level.ma = 1 - alphas.nonsignif[a], ...)
            sid.temp <- c(sid.temp, sids.temp[i])
            ci.lb.temp <- c(ci.lb.temp, out.temp[[paste0("lower.", model)]][tid1, tid2])
            ci.ub.temp <- c(ci.ub.temp, out.temp[[paste0("upper.", model)]][tid1, tid2])
            if(out.temp[[paste0("upper.", model)]][tid1, tid2] < null.val){
              signif <- TRUE
              break
            }
          }
          idx.iter <- which(ci.ub.temp == min(ci.ub.temp))[1]
          sid.iter.left <- c(sid.iter.left, sid.temp[idx.iter])
          data.temp <- data.tt[[idx.iter]]
          if(!signif){
            judge.temp <- judge.mod(dat = data.temp, incr2 = TRUE, tid1, tid2)
            sids.temp <- c(judge.temp$sid1, judge.temp$sid2)
            grp.idx.temp <- NULL
            if(length(judge.temp$sid1) > 0) grp.idx.temp <- c(grp.idx.temp, rep(tid1, length(judge.temp$sid1)))
            if(length(judge.temp$sid2) > 0) grp.idx.temp <- c(grp.idx.temp, rep(tid2, length(judge.temp$sid2)))
            moremod <- (length(sids.temp) > 0)
          }
        }
        if(signif){
          FI.nonsignif.left <- length(sid.iter.left)
          FQ.nonsignif.left <- length(sid.iter.left)/sum(ori.data$n[ori.data$tid == tid1 | ori.data$tid == tid2])
          FQ.nma.nonsignif.left <- length(sid.iter.left)/sum(ori.data$n)
        }
        out.left <- list(FI.nonsignif.left = FI.nonsignif.left, FQ.nonsignif.left = FQ.nonsignif.left,
          FQ.nma.nonsignif.left = FQ.nma.nonsignif.left)
        return(out.left)
      }

      right.fcn <- function(max.iter = Inf, a, ...){
        FI.nonsignif.right <- FQ.nonsignif.right <- FQ.nma.nonsignif.right <- NA
        sid.iter.right <- NULL
        data.temp <- ori.data
        signif <- FALSE
        judge.temp <- judge.mod(dat = data.temp, incr2 = FALSE, tid1, tid2)
        sids.temp <- c(judge.temp$sid1, judge.temp$sid2)
        grp.idx.temp <- NULL
        if(length(judge.temp$sid1) > 0) grp.idx.temp <- c(grp.idx.temp, rep(tid1, length(judge.temp$sid1)))
        if(length(judge.temp$sid2) > 0) grp.idx.temp <- c(grp.idx.temp, rep(tid2, length(judge.temp$sid2)))
        moremod <- (length(sids.temp) > 0)
        while(!signif & moremod & length(sid.iter.right) <= max.iter){
          sid.temp <- tid1.mod.temp <- tid2.mod.temp <-
             ci.lb.temp <- ci.ub.temp <- NULL
          data.tt <- NULL
          for(i in 1:length(sids.temp)){
            data.tt[[i]] <- data.temp
            if(grp.idx.temp[i] == tid1){
            data.tt[[i]]$e[data.tt[[i]]$sid == sids.temp[i] & data.tt[[i]]$tid == tid1] <-
              data.tt[[i]]$e[data.tt[[i]]$sid == sids.temp[i] & data.tt[[i]]$tid == tid1] + 1
            tid1.mod.temp <- c(tid1.mod.temp, 1)
            tid2.mod.temp <- c(tid2.mod.temp, 0)
            }else{
              data.tt[[i]]$e[data.tt[[i]]$sid == sids.temp[i] & data.tt[[i]]$tid == tid2] <-
                data.tt[[i]]$e[data.tt[[i]]$sid == sids.temp[i] & data.tt[[i]]$tid == tid2] - 1
              tid1.mod.temp <- c(tid1.mod.temp, 0)
              tid2.mod.temp <- c(tid2.mod.temp, -1)
            }
            pmeta.temp <- meta::pairwise(treat = tid, event = e, n = n, studlab = sid, data = data.tt[[i]], sm = measure,
              incr = incr, allincr = allincr, addincr = addincr, allstudies = allstudies)
            out.temp <- netmeta(TE = pmeta.temp$TE, seTE = pmeta.temp$seTE, treat1 = pmeta.temp$treat1, treat2 = pmeta.temp$treat2, studlab = pmeta.temp$studlab,
              sm = measure, level.ma = 1 - alphas.nonsignif[a], ...)
            sid.temp <- c(sid.temp, sids.temp[i])
            ci.lb.temp <- c(ci.lb.temp, out.temp[[paste0("lower.", model)]][tid1, tid2])
            ci.ub.temp <- c(ci.ub.temp, out.temp[[paste0("upper.", model)]][tid1, tid2])
            if(out.temp[[paste0("lower.", model)]][tid1, tid2] > null.val){
              signif <- TRUE
              break
            }
          }
          idx.iter <- which(ci.lb.temp == max(ci.lb.temp))[1]
          sid.iter.right <- c(sid.iter.right, sid.temp[idx.iter])
          data.temp <- data.tt[[idx.iter]]
          if(!signif){
            judge.temp <- judge.mod(dat = data.temp, incr2 = FALSE, tid1, tid2)
            sids.temp <- c(judge.temp$sid1, judge.temp$sid2)
            grp.idx.temp <- NULL
            if(length(judge.temp$sid1) > 0) grp.idx.temp <- c(grp.idx.temp, rep(tid1, length(judge.temp$sid1)))
            if(length(judge.temp$sid2) > 0) grp.idx.temp <- c(grp.idx.temp, rep(tid2, length(judge.temp$sid2)))
            moremod <- (length(sids.temp) > 0)
          }
        }
        if(signif){
          FI.nonsignif.right <- length(sid.iter.right)
          FQ.nonsignif.right <- length(sid.iter.right)/sum(ori.data$n[ori.data$tid == tid1 | ori.data$tid == tid2])
          FQ.nma.nonsignif.right <- length(sid.iter.right)/sum(ori.data$n)
        }
        out.right <- list(FI.nonsignif.right = FI.nonsignif.right, FQ.nonsignif.right = FQ.nonsignif.right,
          FQ.nma.nonsignif.right = FQ.nma.nonsignif.right)
        return(out.right)
      }

      if(left & right){
        for(a in 1:length(alphas.nonsignif)){
          if(est.ori[tid1, tid2] >= null.val){
            out.right <- right.fcn(max.iter = Inf, a, ...)
            out.left <- left.fcn(max.iter = ifelse(!is.na(max(out.right$FI.nonsignif.right)), max(out.right$FI.nonsignif.right) + 1, Inf), a, ...)
          }else{
            out.left <- left.fcn(max.iter = Inf, a, ...)
            out.right <- right.fcn(max.iter = ifelse(!is.na(max(out.left$FI.nonsignif.left)), max(out.left$FI.nonsignif.left) + 1, Inf), a, ...)
          }
          if(!is.na(out.left$FI.nonsignif.left) & !is.na(out.right$FI.nonsignif.right)){
            lr <- (out.left$FI.nonsignif.left < out.right$FI.nonsignif.right)
          }
          if(is.na(out.left$FI.nonsignif.left)){
            lr <- FALSE
          }
          if(is.na(out.right$FI.nonsignif.right)){
            lr <- TRUE
          }
          FI.nonsignif[a] <- ifelse(lr, out.left$FI.nonsignif.left, out.right$FI.nonsignif.right)
          FQ.nonsignif[a] <- ifelse(lr, out.left$FQ.nonsignif.left, out.right$FQ.nonsignif.right)
          FQ.nma.nonsignif[a] <- ifelse(lr, out.left$FQ.nma.nonsignif.left, out.right$FQ.nma.nonsignif.right)
        }
      }else{
        for(a in 1:length(alphas.nonsignif)){
          if(left){
            out.left <- left.fcn(max.iter = Inf, a, ...)
            FI.nonsignif[a] <- out.left$FI.nonsignif.left
            FQ.nonsignif[a] <- out.left$FQ.nonsignif.left
            FQ.nma.nonsignif[a] <- out.left$FQ.nma.nonsignif.left
          }else{
            out.right <- right.fcn(max.iter = Inf, a, ...)
            FI.nonsignif[a] <- out.right$FI.nonsignif.right
            FQ.nonsignif[a] <- out.right$FQ.nonsignif.right
            FQ.nma.nonsignif[a] <- out.left$FQ.nma.nonsignif.right
          }
        }
      }
    }

    FI.temp <- FQ.temp <- FQ.nma.temp <- rep(NA, alpha.breaks)
    if(any(signif.alpha) & any(!signif.alpha)){
      FI.temp[signif.alpha] <- FI.signif
      FQ.temp[signif.alpha] <- FQ.signif
      FQ.nma.temp[signif.alpha] <- FQ.nma.signif
      FI.temp[!signif.alpha] <- FI.nonsignif
      FQ.temp[!signif.alpha] <- FQ.nonsignif
      FQ.nma.temp[!signif.alpha] <- FQ.nma.nonsignif
    }
    if(all(signif.alpha)){
      FI.temp <- FI.signif
      FQ.temp <- FQ.signif
      FQ.nma.temp <- FQ.nma.signif
    }
    if(all(!signif.alpha)){
      FI.temp <- FI.nonsignif
      FQ.temp <- FQ.nonsignif
      FQ.nma.temp <- FQ.nma.nonsignif
    }

    FI[[paste(tid1, "vs.", tid2)]] <- FI.temp
    FQ[[paste(tid1, "vs.", tid2)]] <- FQ.temp
    FQ.nma[[paste(tid1, "vs.", tid2)]] <- FQ.nma.temp
    FI.avg[[paste(tid1, "vs.", tid2)]] <- mean(FI.temp)
    FQ.avg[[paste(tid1, "vs.", tid2)]] <- mean(FQ.temp)
    FQ.nma.avg[[paste(tid1, "vs.", tid2)]] <- mean(FQ.nma.temp)
  }

  out <- c(out, list(FI = FI, FI.avg = FI.avg, FQ = FQ, FQ.avg = FQ.avg, FQ.nma = FQ.nma, FQ.nma.avg = FQ.nma.avg))
  class(out) <- c("frag.alpha", "frag.nma.alpha")
  return(out)
}