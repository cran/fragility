frag.ma <- function(e0, n0, e1, n1, data, measure = "OR",
  alpha = 0.05, mod.dir = "both", OR = 1, RR = 1, RD = 0,
  method = "DL", test = "z", ...){
  if(!missing(data)){
    e0 <- eval(substitute(e0), data, parent.frame())
    n0 <- eval(substitute(n0), data, parent.frame())
    e1 <- eval(substitute(e1), data, parent.frame())
    n1 <- eval(substitute(n1), data, parent.frame())
  }
  if(any(e0 < 0) | any(n0 < 0) | any(e1 < 0) | any(n1 < 0)){
    stop("Event counts and sample sizes must be nonnegative.")
  }
  if(any(e0 > n0) | any(e1 > n1)){
    stop("Event counts must not be larger than sample sizes.")
  }
  if(any(e0%%1 != 0) | any(n0%%1 != 0) | any(e1%%1 != 0) | any(n1%%1 != 0)){
    message("Some event counts and/or sample sizes are not integer(s); they are rounded.")
    e0 <- round(e0)
    n0 <- round(n0)
    e1 <- round(e1)
    n1 <- round(n1)
  }
  if(length(n0) != length(e0) | length(e1) != length(e0) | length(n1) != length(e0)){
    stop("e0, n0, e1, and n1 do not have the same length.")
  }
  if(length(e0) == 1){
    stop("The fragility of a single study can be assessed using frag.study().")
  }
  if(any(!is.element(measure, c("OR", "RR", "RD")))){
    stop("measure must be \"OR\", \"RR\", or \"RD\".")
  }
  if(length(measure) != 1){
    stop("Only one measure can be specified.")
  }
  if(alpha <= 0 | alpha >= 1) stop("alpha must be between 0 and 1.")
  if(any(!is.element(mod.dir, c("both", "one", "left", "right")))){
    stop("mod.dir must be \"both\", \"one\", \"left\", or \"right\".")
  }
  if(length(mod.dir) != 1){
    stop("Only one choice can be specified for mod.dir.")
  }
  if(OR <= 0) stop("The null value of OR must be positive.")
  if(RR <= 0) stop("The null value of RR must be positive.")
  if(RD < -1 | RD > 1) stop("The null value of RD must be between -1 and 1.")

  ne0 <- n0 - e0
  ne1 <- n1 - e1
  ori.data <- data.frame(e0 = e0, ne0 = ne0, n0 = n0,
    e1 = e1, ne1 = ne1, n1 = n1)
  if(measure == "OR") null.val <- log(OR)
  if(measure == "RR") null.val <- log(RR)
  if(measure == "RD") null.val <- RD

  rslt.ori <- rma.uni(ai = e1, bi = ne1, ci = e0, di = ne0,
    measure = measure, level = (1 - alpha) * 100, method = method, test = test, ...)
  est.ori <- as.numeric(rslt.ori$beta)
  ci.lb.ori <- rslt.ori$ci.lb
  ci.ub.ori <- rslt.ori$ci.ub
  pval.ori <- rslt.ori$pval
  if(mod.dir == "one"){
    if(est.ori >= null.val) mod.dir <- "right" else mod.dir <- "left"
  }
  out <- list(data = ori.data, measure = measure, alpha = alpha,
    null = null.val, est.ori = est.ori,
    ci.ori = c("LB" = ci.lb.ori, "UB" = ci.ub.ori), pval.ori = pval.ori, mod.dir = mod.dir)

  judge.mod <- function(dat, incr0){
    if(incr0){
      sid0 <- which(dat$e0 < dat$n0)
      sid1 <- which(dat$e1 > 0)
    }else{
      sid0 <- which(dat$e0 > 0)
      sid1 <- which(dat$e1 < dat$n1)
    }
    out <- list(sid0 = sid0, sid1 = sid1)
    return(out)
  }

  if(ci.ub.ori < null.val | ci.lb.ori > null.val){
    sid.iter <- g0.mod.iter <- g1.mod.iter <-
      est.iter <- ci.lb.iter <- ci.ub.iter <- NULL
    less <- (ci.ub.ori < null.val)
    data.temp <- ori.data
    signif <- TRUE
    judge.temp <- judge.mod(dat = data.temp, incr0 = ifelse(less, FALSE, TRUE))
    sids.temp <- c(judge.temp$sid0, judge.temp$sid1)
    grp.idx.temp <- NULL
    if(length(judge.temp$sid0) > 0) grp.idx.temp <- c(grp.idx.temp, rep(0, length(judge.temp$sid0)))
    if(length(judge.temp$sid1) > 0) grp.idx.temp <- c(grp.idx.temp, rep(1, length(judge.temp$sid1)))
    moremod <- (length(sids.temp) > 0)
    while(signif & moremod){
      sid.temp <- g0.mod.temp <- g1.mod.temp <-
         est.temp <- ci.lb.temp <- ci.ub.temp <- NULL
      for(i in 1:length(sids.temp)){
        mod.temp <- rep(0, dim(data.temp)[1])
        if(grp.idx.temp[i] == 0){
          mod.temp[sids.temp[i]] <- ifelse(less, -1, 1)
          out.temp <- rma.uni(ai = data.temp$e1, bi = data.temp$ne1,
            ci = data.temp$e0 + mod.temp, di = data.temp$ne0 - mod.temp,
            measure = measure, level = (1 - alpha) * 100,
            method = method, test = test, ...)
          g0.mod.temp <- c(g0.mod.temp, ifelse(less, -1, 1))
          g1.mod.temp <- c(g1.mod.temp, 0)
        }else{
          mod.temp[sids.temp[i]] <- ifelse(less, 1, -1)
          out.temp <- rma.uni(ai = data.temp$e1 + mod.temp, bi = data.temp$ne1 - mod.temp,
            ci = data.temp$e0, di = data.temp$ne0,
            measure = measure, level = (1 - alpha) * 100,
            method = method, test = test, ...)
          g0.mod.temp <- c(g0.mod.temp, 0)
          g1.mod.temp <- c(g1.mod.temp, ifelse(less, 1, -1))
        }
        sid.temp <- c(sid.temp, sids.temp[i])
        est.temp <- c(est.temp, as.numeric(out.temp$beta))
        ci.lb.temp <- c(ci.lb.temp, out.temp$ci.lb)
        ci.ub.temp <- c(ci.ub.temp, out.temp$ci.ub)
        if((less & (out.temp$ci.ub >= null.val)) | (!less & (out.temp$ci.lb <= null.val))){
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
      g0.mod.iter <- c(g0.mod.iter, g0.mod.temp[idx.iter])
      g1.mod.iter <- c(g1.mod.iter, g1.mod.temp[idx.iter])
      est.iter <- c(est.iter, est.temp[idx.iter])
      ci.lb.iter <- c(ci.lb.iter, ci.lb.temp[idx.iter])
      ci.ub.iter <- c(ci.ub.iter, ci.ub.temp[idx.iter])
      data.temp$e0[sid.temp[idx.iter]] <- data.temp$e0[sid.temp[idx.iter]] + g0.mod.temp[idx.iter]
      data.temp$ne0[sid.temp[idx.iter]] <- data.temp$ne0[sid.temp[idx.iter]] - g0.mod.temp[idx.iter]
      data.temp$e1[sid.temp[idx.iter]] <- data.temp$e1[sid.temp[idx.iter]] + g1.mod.temp[idx.iter]
      data.temp$ne1[sid.temp[idx.iter]] <- data.temp$ne1[sid.temp[idx.iter]] - g1.mod.temp[idx.iter]
      if(signif){
        judge.temp <- judge.mod(dat = data.temp, incr0 = ifelse(less, FALSE, TRUE))
        sids.temp <- c(judge.temp$sid0, judge.temp$sid1)
        grp.idx.temp <- NULL
        if(length(judge.temp$sid0) > 0) grp.idx.temp <- c(grp.idx.temp, rep(0, length(judge.temp$sid0)))
        if(length(judge.temp$sid1) > 0) grp.idx.temp <- c(grp.idx.temp, rep(1, length(judge.temp$sid1)))
        moremod <- (length(sids.temp) > 0)
      }
    }
    data.mod <- data.temp
    if(signif){
      FI <- FQ <- NA
      dir <- "significance cannot be altered"
    }else{
      FI <- length(sid.iter)
      FQ <- FI/sum(ori.data$n0 + ori.data$n1)
      dir <- "significance altered to non-significance"
    }
  }

  if(ci.ub.ori >= null.val & ci.lb.ori <= null.val){
    left <- right <- TRUE
    if(mod.dir == "left") right <- FALSE
    if(mod.dir == "right") left <- FALSE

    left.fcn <- function(max.iter = Inf, ...){
      sid.iter.left <- g0.mod.iter.left <- g1.mod.iter.left <-
        est.iter.left <- ci.lb.iter.left <- ci.ub.iter.left <- NULL
      data.temp <- ori.data
      signif <- FALSE
      judge.temp <- judge.mod(dat = data.temp, incr0 = TRUE)
      sids.temp <- c(judge.temp$sid0, judge.temp$sid1)
      grp.idx.temp <- NULL
      if(length(judge.temp$sid0) > 0) grp.idx.temp <- c(grp.idx.temp, rep(0, length(judge.temp$sid0)))
      if(length(judge.temp$sid1) > 0) grp.idx.temp <- c(grp.idx.temp, rep(1, length(judge.temp$sid1)))
      moremod <- (length(sids.temp) > 0)
      while(!signif & moremod & length(sid.iter.left) <= max.iter){
        sid.temp <- g0.mod.temp <- g1.mod.temp <-
           est.temp <- ci.lb.temp <- ci.ub.temp <- NULL
        for(i in 1:length(sids.temp)){
          mod.temp <- rep(0, dim(data.temp)[1])
          if(grp.idx.temp[i] == 0){
            mod.temp[sids.temp[i]] <- 1
            out.temp <- rma.uni(ai = data.temp$e1, bi = data.temp$ne1,
              ci = data.temp$e0 + mod.temp, di = data.temp$ne0 - mod.temp,
              measure = measure, level = (1 - alpha) * 100,
              method = method, test = test, ...)
            g0.mod.temp <- c(g0.mod.temp, 1)
            g1.mod.temp <- c(g1.mod.temp, 0)
          }else{
            mod.temp[sids.temp[i]] <- -1
            out.temp <- rma.uni(ai = data.temp$e1 + mod.temp, bi = data.temp$ne1 - mod.temp,
              ci = data.temp$e0, di = data.temp$ne0,
              measure = measure, level = (1 - alpha) * 100,
              method = method, test = test, ...)
            g0.mod.temp <- c(g0.mod.temp, 0)
            g1.mod.temp <- c(g1.mod.temp, -1)
          }
          sid.temp <- c(sid.temp, sids.temp[i])
          est.temp <- c(est.temp, as.numeric(out.temp$beta))
          ci.lb.temp <- c(ci.lb.temp, out.temp$ci.lb)
          ci.ub.temp <- c(ci.ub.temp, out.temp$ci.ub)
          if(out.temp$ci.ub < null.val){
            signif <- TRUE
            break
          }
        }
        idx.iter <- which(ci.ub.temp == min(ci.ub.temp))[1]
        sid.iter.left <- c(sid.iter.left, sid.temp[idx.iter])
        g0.mod.iter.left <- c(g0.mod.iter.left, g0.mod.temp[idx.iter])
        g1.mod.iter.left <- c(g1.mod.iter.left, g1.mod.temp[idx.iter])
        est.iter.left <- c(est.iter.left, est.temp[idx.iter])
        ci.lb.iter.left <- c(ci.lb.iter.left, ci.lb.temp[idx.iter])
        ci.ub.iter.left <- c(ci.ub.iter.left, ci.ub.temp[idx.iter])
        data.temp$e0[sid.temp[idx.iter]] <- data.temp$e0[sid.temp[idx.iter]] + g0.mod.temp[idx.iter]
        data.temp$ne0[sid.temp[idx.iter]] <- data.temp$ne0[sid.temp[idx.iter]] - g0.mod.temp[idx.iter]
        data.temp$e1[sid.temp[idx.iter]] <- data.temp$e1[sid.temp[idx.iter]] + g1.mod.temp[idx.iter]
        data.temp$ne1[sid.temp[idx.iter]] <- data.temp$ne1[sid.temp[idx.iter]] - g1.mod.temp[idx.iter]
        if(!signif){
          judge.temp <- judge.mod(dat = data.temp, incr0 = TRUE)
          sids.temp <- c(judge.temp$sid0, judge.temp$sid1)
          grp.idx.temp <- NULL
          if(length(judge.temp$sid0) > 0) grp.idx.temp <- c(grp.idx.temp, rep(0, length(judge.temp$sid0)))
          if(length(judge.temp$sid1) > 0) grp.idx.temp <- c(grp.idx.temp, rep(1, length(judge.temp$sid1)))
          moremod <- (length(sids.temp) > 0)
        }
      }
      data.mod.left <- data.temp
      if(!signif){
        FI.left <- FQ.left <- NA
        dir.left <- "non-significance cannot be altered"
      }else{
        FI.left <- length(sid.iter.left)
        FQ.left <- FI.left/sum(ori.data$n0 + ori.data$n1)
        dir.left <- "non-significance altered to significance"
      }
      out.left <- list(FI.left = FI.left, FQ.left = FQ.left, dir.left = dir.left,
        sid.iter.left = sid.iter.left, g0.mod.iter.left = g0.mod.iter.left, g1.mod.iter.left = g1.mod.iter.left,
        est.iter.left = est.iter.left, ci.lb.iter.left = ci.lb.iter.left, ci.ub.iter.left = ci.ub.iter.left,
        data.mod.left = data.mod.left)
      return(out.left)
    }

    right.fcn <- function(max.iter = Inf, ...){
      sid.iter.right <- g0.mod.iter.right <- g1.mod.iter.right <-
        est.iter.right <- ci.lb.iter.right <- ci.ub.iter.right <- NULL
      data.temp <- ori.data
      signif <- FALSE
      judge.temp <- judge.mod(dat = data.temp, incr0 = FALSE)
      sids.temp <- c(judge.temp$sid0, judge.temp$sid1)
      grp.idx.temp <- NULL
      if(length(judge.temp$sid0) > 0) grp.idx.temp <- c(grp.idx.temp, rep(0, length(judge.temp$sid0)))
      if(length(judge.temp$sid1) > 0) grp.idx.temp <- c(grp.idx.temp, rep(1, length(judge.temp$sid1)))
      moremod <- (length(sids.temp) > 0)
      while(!signif & moremod & length(sid.iter.right) <= max.iter){
        sid.temp <- g0.mod.temp <- g1.mod.temp <-
           est.temp <- ci.lb.temp <- ci.ub.temp <- NULL
        for(i in 1:length(sids.temp)){
          mod.temp <- rep(0, dim(data.temp)[1])
          if(grp.idx.temp[i] == 0){
            mod.temp[sids.temp[i]] <- -1
            out.temp <- rma.uni(ai = data.temp$e1, bi = data.temp$ne1,
              ci = data.temp$e0 + mod.temp, di = data.temp$ne0 - mod.temp,
              measure = measure, level = (1 - alpha) * 100,
              method = method, test = test, ...)
            g0.mod.temp <- c(g0.mod.temp, -1)
            g1.mod.temp <- c(g1.mod.temp, 0)
          }else{
            mod.temp[sids.temp[i]] <- 1
            out.temp <- rma.uni(ai = data.temp$e1 + mod.temp, bi = data.temp$ne1 - mod.temp,
              ci = data.temp$e0, di = data.temp$ne0,
              measure = measure, level = (1 - alpha) * 100,
              method = method, test = test, ...)
            g0.mod.temp <- c(g0.mod.temp, 0)
            g1.mod.temp <- c(g1.mod.temp, 1)
          }
          sid.temp <- c(sid.temp, sids.temp[i])
          est.temp <- c(est.temp, as.numeric(out.temp$beta))
          ci.lb.temp <- c(ci.lb.temp, out.temp$ci.lb)
          ci.ub.temp <- c(ci.ub.temp, out.temp$ci.ub)
          if(out.temp$ci.lb > null.val){
            signif <- TRUE
            break
          }
        }
        idx.iter <- which(ci.lb.temp == max(ci.lb.temp))[1]
        sid.iter.right <- c(sid.iter.right, sid.temp[idx.iter])
        g0.mod.iter.right <- c(g0.mod.iter.right, g0.mod.temp[idx.iter])
        g1.mod.iter.right <- c(g1.mod.iter.right, g1.mod.temp[idx.iter])
        est.iter.right <- c(est.iter.right, est.temp[idx.iter])
        ci.lb.iter.right <- c(ci.lb.iter.right, ci.lb.temp[idx.iter])
        ci.ub.iter.right <- c(ci.ub.iter.right, ci.ub.temp[idx.iter])
        data.temp$e0[sid.temp[idx.iter]] <- data.temp$e0[sid.temp[idx.iter]] + g0.mod.temp[idx.iter]
        data.temp$ne0[sid.temp[idx.iter]] <- data.temp$ne0[sid.temp[idx.iter]] - g0.mod.temp[idx.iter]
        data.temp$e1[sid.temp[idx.iter]] <- data.temp$e1[sid.temp[idx.iter]] + g1.mod.temp[idx.iter]
        data.temp$ne1[sid.temp[idx.iter]] <- data.temp$ne1[sid.temp[idx.iter]] - g1.mod.temp[idx.iter]
        if(!signif){
          judge.temp <- judge.mod(dat = data.temp, incr0 = FALSE)
          sids.temp <- c(judge.temp$sid0, judge.temp$sid1)
          grp.idx.temp <- NULL
          if(length(judge.temp$sid0) > 0) grp.idx.temp <- c(grp.idx.temp, rep(0, length(judge.temp$sid0)))
          if(length(judge.temp$sid1) > 0) grp.idx.temp <- c(grp.idx.temp, rep(1, length(judge.temp$sid1)))
          moremod <- (length(sids.temp) > 0)
        }
      }
      data.mod.right <- data.temp
      if(!signif){
        FI.right <- FQ.right <- NA
        dir.right <- "non-significance cannot be altered"
      }else{
        FI.right <- length(sid.iter.right)
        FQ.right <- FI.right/sum(ori.data$n0 + ori.data$n1)
        dir.right <- "non-significance altered to significance"
      }
      out.right <- list(FI.right = FI.right, FQ.right = FQ.right, dir.right = dir.right,
        sid.iter.right = sid.iter.right, g0.mod.iter.right = g0.mod.iter.right, g1.mod.iter.right = g1.mod.iter.right,
        est.iter.right = est.iter.right, ci.lb.iter.right = ci.lb.iter.right, ci.ub.iter.right = ci.ub.iter.right,
        data.mod.right = data.mod.right)
      return(out.right)
    }

    if(left & right){
      if(est.ori >= null.val){
        out.right <- right.fcn(max.iter = Inf, ...)
        out.left <- left.fcn(max.iter = ifelse(!is.na(out.right$FI.right), out.right$FI.right + 1, Inf), ...)
      }else{
        out.left <- left.fcn(max.iter = Inf, ...)
        out.right <- right.fcn(max.iter = ifelse(!is.na(out.left$FI.left), out.left$FI.left + 1, Inf) , ...)
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
        FI <- out.left$FI.left
        FQ <- out.left$FQ.left
        sid.iter <- out.left$sid.iter.left
        g0.mod.iter <- out.left$g0.mod.iter.left
        g1.mod.iter <- out.left$g1.mod.iter.left
        est.iter <- out.left$est.iter.left
        ci.lb.iter <- out.left$ci.lb.iter.left
        ci.ub.iter <- out.left$ci.ub.iter.left
        data.mod <- out.left$data.mod.left
      }else{
        FI <- out.right$FI.right
        FQ <- out.right$FQ.right
        sid.iter <- out.right$sid.iter.right
        g0.mod.iter <- out.right$g0.mod.iter.right
        g1.mod.iter <- out.right$g1.mod.iter.right
        est.iter <- out.right$est.iter.right
        ci.lb.iter <- out.right$ci.lb.iter.right
        ci.ub.iter <- out.right$ci.ub.iter.right
        data.mod <- out.right$data.mod.right
      }
      if(!is.na(FI)){
        dir <- "non-significance altered to significance"
      }else{
        dir <- "non-significance cannot be altered"
      }
    }else{
      if(left){
        out.left <- left.fcn(max.iter = Inf, ...)
        FI <- out.left$FI.left
        FQ <- out.left$FQ.left
        dir <- out.left$dir.left
        sid.iter <- out.left$sid.iter.left
        g0.mod.iter <- out.left$g0.mod.iter.left
        g1.mod.iter <- out.left$g1.mod.iter.left
        est.iter <- out.left$est.iter.left
        ci.lb.iter <- out.left$ci.lb.iter.left
        ci.ub.iter <- out.left$ci.ub.iter.left
        data.mod <- out.left$data.mod.left
      }else{
        out.right <- right.fcn(max.iter = Inf, ...)
        FI <- out.right$FI.right
        FQ <- out.right$FQ.right
        dir <- out.right$dir.right
        sid.iter <- out.right$sid.iter.right
        g0.mod.iter <- out.right$g0.mod.iter.right
        g1.mod.iter <- out.right$g1.mod.iter.right
        est.iter <- out.right$est.iter.right
        ci.lb.iter <- out.right$ci.lb.iter.right
        ci.ub.iter <- out.right$ci.ub.iter.right
        data.mod <- out.right$data.mod.right
      }
    }
  }

  ci.iter <- cbind(ci.lb.iter, ci.ub.iter)
  colnames(ci.iter) <- c("LB", "UB")
  out <- c(out, list(FI = FI, FQ = FQ, dir = dir,
    sid.iter = sid.iter, g0.mod.iter = g0.mod.iter, g1.mod.iter = g1.mod.iter,
    est.iter = est.iter, ci.iter = ci.iter, data.mod = data.mod))
  class(out) <- c("frag.ma")
  return(out)
}