frag.nma <- function(sid, tid, e, n, data, measure = "OR", random = TRUE,
  alpha = 0.05, mod.dir = "both", tid1.f, tid2.f, OR = 1, RR = 1, RD = 0,
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
  if(alpha <= 0 | alpha >= 1) stop("alpha must be between 0 and 1.")
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

  pmeta.ori <- meta::pairwise(treat = tid, event = e, n = n, studlab = sid, data = ori.data, sm = measure,
    incr = incr, allincr = allincr, addincr = addincr, allstudies = allstudies)
  rslt.ori <- netmeta(TE = pmeta.ori$TE, seTE = pmeta.ori$seTE, treat1 = pmeta.ori$treat1, treat2 = pmeta.ori$treat2, studlab = pmeta.ori$studlab,
    sm = measure, level.ma = 1 - alpha, ...)
  est.ori <- rslt.ori[[paste0("TE.", model)]]
  ci.lb.ori <- rslt.ori[[paste0("lower.", model)]]
  ci.ub.ori <- rslt.ori[[paste0("upper.", model)]]
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

  out <- list(data = ori.data, measure = measure, alpha = alpha,
    null = null.val, est.ori = est.ori,
    ci.lb.ori = ci.lb.ori, ci.ub.ori = ci.ub.ori, pval.ori = pval.ori,
    mod.dir = mod.dir, tid.f = paste(tid1.f, "vs.", tid2.f))

  empt.list <- rep(list(NULL), length(tid1.f))
  names(empt.list) <- paste(tid1.f, "vs.", tid2.f)
  sid.iter <- tid1.mod.iter <- tid2.mod.iter <-
    est.iter <- ci.lb.iter <- ci.ub.iter <- data.mod <- empt.list
  na.mtx <- matrix(NA, dim(est.ori)[1], dim(est.ori)[2])
  rownames(na.mtx) <- colnames(na.mtx) <- rslt.ori$trts
  FI <- FQ <- FQ.nma <- dir <- na.mtx

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
        FI[tid1, tid2] <- FI[tid2, tid1]
        FQ[tid1, tid2] <- FQ[tid2, tid1]
        FQ.nma[tid1, tid2] <- FQ.nma[tid2, tid1]
        dir[tid1, tid2] <- dir[tid2, tid1]
        sid.iter[[paste(tid1, "vs.", tid2)]] <- sid.iter[[paste(tid2, "vs.", tid1)]]
        tid1.mod.iter[[paste(tid1, "vs.", tid2)]] <- tid2.mod.iter[[paste(tid2, "vs.", tid1)]]
        tid2.mod.iter[[paste(tid1, "vs.", tid2)]] <- tid1.mod.iter[[paste(tid2, "vs.", tid1)]]
        est.iter[[paste(tid1, "vs.", tid2)]] <- -est.iter[[paste(tid2, "vs.", tid1)]]
        ci.lb.iter[[paste(tid1, "vs.", tid2)]] <- -ci.ub.iter[[paste(tid2, "vs.", tid1)]]
        ci.ub.iter[[paste(tid1, "vs.", tid2)]] <- -ci.lb.iter[[paste(tid2, "vs.", tid1)]]
        data.mod[[paste(tid1, "vs.", tid2)]] <- data.mod[[paste(tid2, "vs.", tid1)]]
        next
      }
    }

    if(ci.ub.ori[tid1, tid2] < null.val | ci.lb.ori[tid1, tid2] > null.val){
      less <- (ci.ub.ori[tid1, tid2] < null.val)
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
           est.temp <- ci.lb.temp <- ci.ub.temp <- NULL
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
            sm = measure, level.ma = 1 - alpha, ...)
          sid.temp <- c(sid.temp, sids.temp[i])
          est.temp <- c(est.temp, out.temp[[paste0("TE.", model)]][tid1, tid2])
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
        sid.iter[[paste(tid1, "vs.", tid2)]] <- c(sid.iter[[paste(tid1, "vs.", tid2)]], sid.temp[idx.iter])
        tid1.mod.iter[[paste(tid1, "vs.", tid2)]] <- c(tid1.mod.iter[[paste(tid1, "vs.", tid2)]], tid1.mod.temp[idx.iter])
        tid2.mod.iter[[paste(tid1, "vs.", tid2)]] <- c(tid2.mod.iter[[paste(tid1, "vs.", tid2)]], tid2.mod.temp[idx.iter])
        est.iter[[paste(tid1, "vs.", tid2)]] <- c(est.iter[[paste(tid1, "vs.", tid2)]], est.temp[idx.iter])
        ci.lb.iter[[paste(tid1, "vs.", tid2)]] <- c(ci.lb.iter[[paste(tid1, "vs.", tid2)]], ci.lb.temp[idx.iter])
        ci.ub.iter[[paste(tid1, "vs.", tid2)]] <- c(ci.ub.iter[[paste(tid1, "vs.", tid2)]], ci.ub.temp[idx.iter])
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
      data.mod[[paste(tid1, "vs.", tid2)]] <- data.temp
      if(signif){
        FI[tid1, tid2] <- FQ[tid1, tid2] <- FQ.nma[tid1, tid2] <- NA
        dir[tid1, tid2] <- "significance cannot be altered"
      }else{
        FI[tid1, tid2] <- length(sid.iter[[paste(tid1, "vs.", tid2)]])
        FQ[tid1, tid2] <- FI[tid1, tid2]/sum(ori.data$n[ori.data$tid == tid1 | ori.data$tid == tid2])
        FQ.nma[tid1, tid2] <- FI[tid1, tid2]/sum(ori.data$n)
        dir[tid1, tid2] <- "significance altered to non-significance"
      }
    }

    if(ci.ub.ori[tid1, tid2] >= null.val & ci.lb.ori[tid1, tid2] <= null.val){
      left <- right <- TRUE
      if(mod.dir[tid1, tid2] == "left") right <- FALSE
      if(mod.dir[tid1, tid2] == "right") left <- FALSE

      left.fcn <- function(max.iter = Inf, ...){
        sid.iter.left <- tid1.mod.iter.left <- tid2.mod.iter.left <-
          est.iter.left <- ci.lb.iter.left <- ci.ub.iter.left <- NULL
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
            est.temp <- ci.lb.temp <- ci.ub.temp <- NULL
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
              sm = measure, level.ma = 1 - alpha, ...)
            sid.temp <- c(sid.temp, sids.temp[i])
            est.temp <- c(est.temp, out.temp[[paste0("TE.", model)]][tid1, tid2])
            ci.lb.temp <- c(ci.lb.temp, out.temp[[paste0("lower.", model)]][tid1, tid2])
            ci.ub.temp <- c(ci.ub.temp, out.temp[[paste0("upper.", model)]][tid1, tid2])
            if(out.temp[[paste0("upper.", model)]][tid1, tid2] < null.val){
              signif <- TRUE
              break
            }
          }
          idx.iter <- which(ci.ub.temp == min(ci.ub.temp))[1]
          sid.iter.left <- c(sid.iter.left, sid.temp[idx.iter])
          tid1.mod.iter.left <- c(tid1.mod.iter.left, tid1.mod.temp[idx.iter])
          tid2.mod.iter.left <- c(tid2.mod.iter.left, tid2.mod.temp[idx.iter])
          est.iter.left <- c(est.iter.left, est.temp[idx.iter])
          ci.lb.iter.left <- c(ci.lb.iter.left, ci.lb.temp[idx.iter])
          ci.ub.iter.left <- c(ci.ub.iter.left, ci.ub.temp[idx.iter])
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
        data.mod.left <- data.temp
        if(!signif){
          FI.left <- FQ.left <- FQ.nma.left <- NA
          dir.left <- "non-significance cannot be altered"
        }else{
          FI.left <- length(sid.iter.left)
          FQ.left <- FI.left/sum(ori.data$n[ori.data$tid == tid1 | ori.data$tid == tid2])
          FQ.nma.left <- FI.left/sum(ori.data$n)
          dir.left <- "non-significance altered to significance"
        }
        out.left <- list(FI.left = FI.left, FQ.left = FQ.left, FQ.nma.left = FQ.nma.left, dir.left = dir.left,
          sid.iter.left = sid.iter.left, tid1.mod.iter.left = tid1.mod.iter.left, tid2.mod.iter.left = tid2.mod.iter.left,
          est.iter.left = est.iter.left, ci.lb.iter.left = ci.lb.iter.left, ci.ub.iter.left = ci.ub.iter.left,
          data.mod.left = data.mod.left)
        return(out.left)
      }

      right.fcn <- function(max.iter = Inf, ...){
        sid.iter.right <- tid1.mod.iter.right <- tid2.mod.iter.right <-
          est.iter.right <- ci.lb.iter.right <- ci.ub.iter.right <- NULL
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
             est.temp <- ci.lb.temp <- ci.ub.temp <- NULL
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
              sm = measure, level.ma = 1 - alpha, ...)
            sid.temp <- c(sid.temp, sids.temp[i])
            est.temp <- c(est.temp, out.temp[[paste0("TE.", model)]][tid1, tid2])
            ci.lb.temp <- c(ci.lb.temp, out.temp[[paste0("lower.", model)]][tid1, tid2])
            ci.ub.temp <- c(ci.ub.temp, out.temp[[paste0("upper.", model)]][tid1, tid2])
            if(out.temp[[paste0("lower.", model)]][tid1, tid2] > null.val){
              signif <- TRUE
              break
            }
          }
          idx.iter <- which(ci.lb.temp == max(ci.lb.temp))[1]
          sid.iter.right <- c(sid.iter.right, sid.temp[idx.iter])
          tid1.mod.iter.right <- c(tid1.mod.iter.right, tid1.mod.temp[idx.iter])
          tid2.mod.iter.right <- c(tid2.mod.iter.right, tid2.mod.temp[idx.iter])
          est.iter.right <- c(est.iter.right, est.temp[idx.iter])
          ci.lb.iter.right <- c(ci.lb.iter.right, ci.lb.temp[idx.iter])
          ci.ub.iter.right <- c(ci.ub.iter.right, ci.ub.temp[idx.iter])
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
        data.mod.right <- data.temp
        if(!signif){
          FI.right <- FQ.right <- FQ.nma.right <- NA
          dir.right <- "non-significance cannot be altered"
        }else{
          FI.right <- length(sid.iter.right)
          FQ.right <- FI.right/sum(ori.data$n[ori.data$tid == tid1 | ori.data$tid == tid2])
          FQ.nma.right <- FI.right/sum(ori.data$n)
          dir.right <- "non-significance altered to significance"
        }  
        out.right <- list(FI.right = FI.right, FQ.right = FQ.right, FQ.nma.right = FQ.nma.right, dir.right = dir.right,
          sid.iter.right = sid.iter.right, tid1.mod.iter.right = tid1.mod.iter.right, tid2.mod.iter.right = tid2.mod.iter.right,
          est.iter.right = est.iter.right, ci.lb.iter.right = ci.lb.iter.right, ci.ub.iter.right = ci.ub.iter.right,
          data.mod.right = data.mod.right)
        return(out.right)
      }

      if(left & right){
        if(est.ori[tid1, tid2] >= null.val){
          out.right <- right.fcn(max.iter = Inf, ...)
          out.left <- left.fcn(max.iter = ifelse(!is.na(out.right$FI.right), out.right$FI.right + 1, Inf), ...)
        }else{
          out.left <- left.fcn(max.iter = Inf, ...)
          out.right <- right.fcn(max.iter = ifelse(!is.na(out.left$FI.left), out.left$FI.left + 1, Inf), ...)
        }
        if(!is.na(out.left$FI.left) & !is.na(out.right$FI.right)){
          lr <- (out.left$FI.left < out.right$FI.right)
        }
        if(is.na(out.left$FI.left)){
          lr <- FALSE
        }
        if(is.na(out.right$FI.right)){
          lr <- TRUE
        }
        if(lr){
          FI[tid1, tid2] <- out.left$FI.left
          FQ[tid1, tid2] <- out.left$FQ.left
          FQ.nma[tid1, tid2] <- out.left$FQ.nma.left
          sid.iter[[paste(tid1, "vs.", tid2)]] <- out.left$sid.iter.left
          tid1.mod.iter[[paste(tid1, "vs.", tid2)]] <- out.left$tid1.mod.iter.left
          tid2.mod.iter[[paste(tid1, "vs.", tid2)]] <- out.left$tid2.mod.iter.left
          est.iter[[paste(tid1, "vs.", tid2)]] <- out.left$est.iter.left
          ci.lb.iter[[paste(tid1, "vs.", tid2)]] <- out.left$ci.lb.iter.left
          ci.ub.iter[[paste(tid1, "vs.", tid2)]] <- out.left$ci.ub.iter.left
          data.mod[[paste(tid1, "vs.", tid2)]] <- out.left$data.mod.left
        }else{
          FI[tid1, tid2] <- out.right$FI.right
          FQ[tid1, tid2] <- out.right$FQ.right
          FQ.nma[tid1, tid2] <- out.right$FQ.nma.right
          sid.iter[[paste(tid1, "vs.", tid2)]] <- out.right$sid.iter.right
          tid1.mod.iter[[paste(tid1, "vs.", tid2)]] <- out.right$tid1.mod.iter.right
          tid2.mod.iter[[paste(tid1, "vs.", tid2)]] <- out.right$tid2.mod.iter.right
          est.iter[[paste(tid1, "vs.", tid2)]] <- out.right$est.iter.right
          ci.lb.iter[[paste(tid1, "vs.", tid2)]] <- out.right$ci.lb.iter.right
          ci.ub.iter[[paste(tid1, "vs.", tid2)]] <- out.right$ci.ub.iter.right
          data.mod[[paste(tid1, "vs.", tid2)]] <- out.right$data.mod.right
        }
        if(!is.na(FI[tid1, tid2])){
          dir[tid1, tid2] <- "non-significance altered to significance"
        }else{
          dir[tid1, tid2] <- "non-significance cannot be altered"
        }
      }else{
        if(left){
          out.left <- left.fcn(max.iter = Inf, ...)
          FI[tid1, tid2] <- out.left$FI.left
          FQ[tid1, tid2] <- out.left$FQ.left
          FQ.nma[tid1, tid2] <- out.left$FQ.nma.left
          dir[tid1, tid2] <- out.left$dir.left
          sid.iter[[paste(tid1, "vs.", tid2)]] <- out.left$sid.iter.left
          tid1.mod.iter[[paste(tid1, "vs.", tid2)]] <- out.left$tid1.mod.iter.left
          tid2.mod.iter[[paste(tid1, "vs.", tid2)]] <- out.left$tid2.mod.iter.left
          est.iter[[paste(tid1, "vs.", tid2)]] <- out.left$est.iter.left
          ci.lb.iter[[paste(tid1, "vs.", tid2)]] <- out.left$ci.lb.iter.left
          ci.ub.iter[[paste(tid1, "vs.", tid2)]] <- out.left$ci.ub.iter.left
          data.mod[[paste(tid1, "vs.", tid2)]] <- out.left$data.mod.left
        }else{
          out.right <- right.fcn(max.iter = Inf, ...)
          FI[tid1, tid2] <- out.right$FI.right
          FQ[tid1, tid2] <- out.right$FQ.right
          FQ.nma[tid1, tid2] <- out.right$FQ.nma.right
          dir[tid1, tid2] <- out.right$dir.right
          sid.iter[[paste(tid1, "vs.", tid2)]] <- out.right$sid.iter.right
          tid1.mod.iter[[paste(tid1, "vs.", tid2)]] <- out.right$tid1.mod.iter.right
          tid2.mod.iter[[paste(tid1, "vs.", tid2)]] <- out.right$tid2.mod.iter.right
          est.iter[[paste(tid1, "vs.", tid2)]] <- out.right$est.iter.right
          ci.lb.iter[[paste(tid1, "vs.", tid2)]] <- out.right$ci.lb.iter.right
          ci.ub.iter[[paste(tid1, "vs.", tid2)]] <- out.right$ci.ub.iter.right
          data.mod[[paste(tid1, "vs.", tid2)]] <- out.right$data.mod.right
        }
      }
    }
  }

  out <- c(out, list(FI = FI, FQ = FQ, FQ.nma = FQ.nma, dir = dir,
    sid.iter = sid.iter, tid1.mod.iter = tid1.mod.iter, tid2.mod.iter = tid2.mod.iter,
    est.iter = est.iter, ci.lb.iter = ci.lb.iter, ci.ub.iter = ci.ub.iter, data.mod = data.mod))
  class(out) <- c("frag.nma")
  return(out)
}