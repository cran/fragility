frag.ma.alpha <- function(e0, n0, e1, n1, data, measure = "OR",
  alpha.from = 0.005, alpha.to = 0.05, alpha.breaks = 100,
  mod.dir = "both", OR = 1, RR = 1, RD = 0,
  method = "DL", test = "z", drop00 = FALSE, ...){
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
  alphas <- seq(from = alpha.from, to = alpha.to,
    length.out = alpha.breaks)

  dz.ori <- FALSE
  if(drop00 & all((e1 == 0 & e0 == 0) | (ne1 == 0 & ne0 == 0))){
    dz.ori <- TRUE
  }
  if(dz.ori){
    rslt.ori <- data.frame(beta = 0, se = Inf, pval = 1, test = test)
  }else{
    rslt.ori <- rma.uni(ai = e1, bi = ne1, ci = e0, di = ne0,
      measure = measure, level = 95,
      method = method, test = test, drop00 = drop00, ...)
  }
  est.ori <- as.numeric(rslt.ori$beta)
  se.ori <- rslt.ori$se
  pval.ori <- rslt.ori$pval
  if(mod.dir == "one"){
    if(est.ori >= null.val) mod.dir <- "right" else mod.dir <- "left"
  }
  out <- list(data = ori.data, measure = measure, alphas = alphas,
    null = null.val, est.ori = est.ori, se.ori = se.ori, test = rslt.ori$test,
    pval.ori = pval.ori, mod.dir = mod.dir)

  FI <- FQ <- rep(NA, alpha.breaks)
  signif.alpha <- (alphas > pval.ori)

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

  if(any(signif.alpha)){
    alphas.signif <- alphas[signif.alpha]
    FI.signif <- FQ.signif <- rep(NA, sum(signif.alpha))
    less <- (est.ori < null.val)
    for(a in 1:length(alphas.signif)){
      sid.iter <- NULL
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
           ci.lb.temp <- ci.ub.temp <- NULL
        for(i in 1:length(sids.temp)){
          mod.temp <- rep(0, dim(data.temp)[1])
          if(grp.idx.temp[i] == 0){
            mod.temp[sids.temp[i]] <- ifelse(less, -1, 1)
            ai.temp <- data.temp$e1
            bi.temp <- data.temp$ne1
            ci.temp <- data.temp$e0 + mod.temp
            di.temp <- data.temp$ne0 - mod.temp
            dz.temp <- FALSE
            if(drop00 & all((ai.temp == 0 & ci.temp == 0) | (bi.temp == 0 & di.temp == 0))){
              dz.temp <- TRUE
            }
            if(dz.temp){
              out.temp <- data.frame(beta = 0, ci.lb = -Inf, ci.ub = Inf)
            }else{
              out.temp <- rma.uni(ai = ai.temp, bi = bi.temp, ci = ci.temp, di = di.temp,
                measure = measure, level = (1 - alphas.signif[a]) * 100,
                method = method, test = test, drop00 = drop00, ...)
            }
            g0.mod.temp <- c(g0.mod.temp, ifelse(less, -1, 1))
            g1.mod.temp <- c(g1.mod.temp, 0)
          }else{
            mod.temp[sids.temp[i]] <- ifelse(less, 1, -1)
            ai.temp <- data.temp$e1 + mod.temp
            bi.temp <- data.temp$ne1 - mod.temp
            ci.temp <- data.temp$e0
            di.temp <- data.temp$ne0
            dz.temp <- FALSE
            if(drop00 & all((ai.temp == 0 & ci.temp == 0) | (bi.temp == 0 & di.temp == 0))){
              dz.temp <- TRUE
            }
            if(dz.temp){
              out.temp <- data.frame(beta = 0, ci.lb = -Inf, ci.ub = Inf)
            }else{
              out.temp <- rma.uni(ai = ai.temp, bi = bi.temp, ci = ci.temp, di = di.temp,
                measure = measure, level = (1 - alphas.signif[a]) * 100,
                method = method, test = test, drop00 = drop00, ...)
            }
            g0.mod.temp <- c(g0.mod.temp, 0)
            g1.mod.temp <- c(g1.mod.temp, ifelse(less, 1, -1))
          }
          sid.temp <- c(sid.temp, sids.temp[i])
          ci.lb.temp <- c(ci.lb.temp, out.temp$ci.lb)
          ci.ub.temp <- c(ci.ub.temp, out.temp$ci.ub)
          if((less & (out.temp$ci.ub >= null.val)) |
            (!less & (out.temp$ci.lb <= null.val))){
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
      if(!signif){
        FI.signif[a] <- length(sid.iter)
        FQ.signif[a] <- length(sid.iter)/sum(ori.data$n0 + ori.data$n1)
      }
    }
  }

  if(any(!signif.alpha)){
    alphas.nonsignif <- alphas[!signif.alpha]
    FI.nonsignif <- FQ.nonsignif <- rep(NA, sum(!signif.alpha))
    left <- right <- TRUE
    if(mod.dir == "left") right <- FALSE
    if(mod.dir == "right") left <- FALSE

    left.fcn <- function(max.iter = Inf, a, ...){
      FI.nonsignif.left <- FQ.nonsignif.left <- NA
      sid.iter.left <- NULL
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
         ci.lb.temp <- ci.ub.temp <- NULL
        for(i in 1:length(sids.temp)){
          mod.temp <- rep(0, dim(data.temp)[1])
          if(grp.idx.temp[i] == 0){
            mod.temp[sids.temp[i]] <- 1
            ai.temp <- data.temp$e1
            bi.temp <- data.temp$ne1
            ci.temp <- data.temp$e0 + mod.temp
            di.temp <- data.temp$ne0 - mod.temp
            dz.temp <- FALSE
            if(drop00 & all((ai.temp == 0 & ci.temp == 0) | (bi.temp == 0 & di.temp == 0))){
              dz.temp <- TRUE
            }
            if(dz.temp){
              out.temp <- data.frame(beta = 0, ci.lb = -Inf, ci.ub = Inf)
            }else{
              out.temp <- rma.uni(ai = ai.temp, bi = bi.temp, ci = ci.temp, di = di.temp,
                measure = measure, level = (1 - alphas.nonsignif[a]) * 100,
                method = method, test = test, drop00 = drop00, ...)
            }
            g0.mod.temp <- c(g0.mod.temp, 1)
            g1.mod.temp <- c(g1.mod.temp, 0)
          }else{
            mod.temp[sids.temp[i]] <- -1
            ai.temp <- data.temp$e1 + mod.temp
            bi.temp <- data.temp$ne1 - mod.temp
            ci.temp <- data.temp$e0
            di.temp <- data.temp$ne0
            dz.temp <- FALSE
            if(drop00 & all((ai.temp == 0 & ci.temp == 0) | (bi.temp == 0 & di.temp == 0))){
              dz.temp <- TRUE
            }
            if(dz.temp){
              out.temp <- data.frame(beta = 0, ci.lb = -Inf, ci.ub = Inf)
            }else{
              out.temp <- rma.uni(ai = ai.temp, bi = bi.temp, ci = ci.temp, di = di.temp,
                measure = measure, level = (1 - alphas.nonsignif[a]) * 100,
                method = method, test = test, drop00 = drop00, ...)
            }
            g0.mod.temp <- c(g0.mod.temp, 0)
            g1.mod.temp <- c(g1.mod.temp, -1)
          }
          sid.temp <- c(sid.temp, sids.temp[i])
          ci.lb.temp <- c(ci.lb.temp, out.temp$ci.lb)
          ci.ub.temp <- c(ci.ub.temp, out.temp$ci.ub)
          if(out.temp$ci.ub < null.val){
            signif <- TRUE
            break
          }
        }
        idx.iter <- which(ci.ub.temp == min(ci.ub.temp))[1]
        sid.iter.left <- c(sid.iter.left, sid.temp[idx.iter])
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
      if(signif){
        FI.nonsignif.left <- length(sid.iter.left)
        FQ.nonsignif.left <- length(sid.iter.left)/sum(ori.data$n0 + ori.data$n1)
      }
      out.left <- list(FI.nonsignif.left = FI.nonsignif.left, FQ.nonsignif.left = FQ.nonsignif.left)
      return(out.left)
    }

    right.fcn <- function(max.iter = Inf, a, ...){
      FI.nonsignif.right <- FQ.nonsignif.right <- NA
      sid.iter.right <- NULL
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
           ci.lb.temp <- ci.ub.temp <- NULL
        for(i in 1:length(sids.temp)){
          mod.temp <- rep(0, dim(data.temp)[1])
          if(grp.idx.temp[i] == 0){
            mod.temp[sids.temp[i]] <- -1
            ai.temp <- data.temp$e1
            bi.temp <- data.temp$ne1
            ci.temp <- data.temp$e0 + mod.temp
            di.temp <- data.temp$ne0 - mod.temp
            dz.temp <- FALSE
            if(drop00 & all((ai.temp == 0 & ci.temp == 0) | (bi.temp == 0 & di.temp == 0))){
              dz.temp <- TRUE
            }
            if(dz.temp){
              out.temp <- data.frame(beta = 0, ci.lb = -Inf, ci.ub = Inf)
            }else{
              out.temp <- rma.uni(ai = ai.temp, bi = bi.temp, ci = ci.temp, di = di.temp,
                measure = measure, level = (1 - alphas.nonsignif[a]) * 100,
                method = method, test = test, drop00 = drop00, ...)
            }
            g0.mod.temp <- c(g0.mod.temp, -1)
            g1.mod.temp <- c(g1.mod.temp, 0)
          }else{
            mod.temp[sids.temp[i]] <- 1
            ai.temp <- data.temp$e1 + mod.temp
            bi.temp <- data.temp$ne1 - mod.temp
            ci.temp <- data.temp$e0
            di.temp <- data.temp$ne0
            dz.temp <- FALSE
            if(drop00 & all((ai.temp == 0 & ci.temp == 0) | (bi.temp == 0 & di.temp == 0))){
              dz.temp <- TRUE
            }
            if(dz.temp){
              out.temp <- data.frame(beta = 0, ci.lb = -Inf, ci.ub = Inf)
            }else{
              out.temp <- rma.uni(ai = ai.temp, bi = bi.temp, ci = ci.temp, di = di.temp,
                measure = measure, level = (1 - alphas.nonsignif[a]) * 100,
                method = method, test = test, drop00 = drop00, ...)
            }
            g0.mod.temp <- c(g0.mod.temp, 0)
            g1.mod.temp <- c(g1.mod.temp, 1)
          }
          sid.temp <- c(sid.temp, sids.temp[i])
          ci.lb.temp <- c(ci.lb.temp, out.temp$ci.lb)
          ci.ub.temp <- c(ci.ub.temp, out.temp$ci.ub)
          if(out.temp$ci.lb > null.val){
            signif <- TRUE
            break
          }
        }
        idx.iter <- which(ci.lb.temp == max(ci.lb.temp))[1]
        sid.iter.right <- c(sid.iter.right, sid.temp[idx.iter])
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
      if(signif){
        FI.nonsignif.right <- length(sid.iter.right)
        FQ.nonsignif.right <- length(sid.iter.right)/sum(ori.data$n0 + ori.data$n1)
      }
      out.right <- list(FI.nonsignif.right = FI.nonsignif.right, FQ.nonsignif.right = FQ.nonsignif.right)
      return(out.right)
    }

    if(left & right){
      for(a in 1:length(alphas.nonsignif)){
        if(est.ori >= null.val){
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
      }
    }else{
      for(a in 1:length(alphas.nonsignif)){
        if(left){
          out.left <- left.fcn(max.iter = Inf, a, ...)
          FI.nonsignif[a] <- out.left$FI.nonsignif.left
          FQ.nonsignif[a] <- out.left$FQ.nonsignif.left
        }else{
          out.right <- right.fcn(max.iter = Inf, a, ...)
          FI.nonsignif[a] <- out.right$FI.nonsignif.right
          FQ.nonsignif[a] <- out.right$FQ.nonsignif.right
        }
      }
    }
  }

  if(any(signif.alpha) & any(!signif.alpha)){
    FI[signif.alpha] <- FI.signif
    FQ[signif.alpha] <- FQ.signif
    FI[!signif.alpha] <- FI.nonsignif
    FQ[!signif.alpha] <- FQ.nonsignif
  }
  if(all(signif.alpha)){
    FI <- FI.signif
    FQ <- FQ.signif
  }
  if(all(!signif.alpha)){
    FI <- FI.nonsignif
    FQ <- FQ.nonsignif
  }
  FI.avg <- mean(FI)
  FQ.avg <- mean(FQ)
  out <- c(out, list(FI = FI, FI.avg = FI.avg, FQ = FQ, FQ.avg = FQ.avg))
  class(out) <- c("frag.alpha", "frag.ma.alpha")
  return(out)
}