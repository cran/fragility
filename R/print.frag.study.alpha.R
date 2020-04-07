print.frag.study.alpha <- function(x, ...){
  if(!inherits(x, "frag.study.alpha")){
    stop("The input must be an object of \"frag.study.alpha\".")
  }

  cat("___________________________________\n")
  cat(paste0("Original data:\n"))
  print(x$data)

  if(all(x$f0.range == 0) & all(x$f1.range == 0)){
    cat("No event modification occurs in both groups 0 and 1\n")
  }
  if(any(x$f0.range != 0) & all(x$f1.range == 0)){
    cat("Range of event modification in group 0:\n")
    plu1 <- ifelse(abs(x$f0.range[1]) == 1, "", "s")
    plu2 <- ifelse(abs(x$f0.range[2]) == 1, "", "s")
    if(all(x$f0.range != 0)){
      cat(paste0("  up to ", abs(x$f0.range[1]), " event", plu1, " modified to be non-event", plu1, ";\n",
        "  up to ", abs(x$f0.range[2]), " non-event", plu2, " modified to be event", plu2, "\n"))
    }
    if(x$f0.range[1] == 0){
      cat(paste0("  up to ", abs(x$f0.range[2]), " non-event", plu2, " modified to be event", plu2, "\n"))
    }
    if(x$f0.range[2] == 0){
      cat(paste0("  up to ", abs(x$f0.range[1]), " event", plu1, " modified to be non-event", plu1, "\n"))
    }
    cat("No event modification occurs in group 1\n")
  }
  if(all(x$f0.range == 0) & any(x$f1.range != 0)){
    cat("No event modification occurs in group 0\n")
    cat("Range of event modification in group 1:\n")
    plu1 <- ifelse(abs(x$f1.range[1]) == 1, "", "s")
    plu2 <- ifelse(abs(x$f1.range[2]) == 1, "", "s")
    if(all(x$f1.range != 0)){
      cat(paste0("  up to ", abs(x$f1.range[1]), " event", plu1, " modified to be non-event", plu1, ";\n",
        "  up to ", abs(x$f1.range[2]), " non-event", plu2, " modified to be event", plu2, "\n"))
    }
    if(x$f1.range[1] == 0){
      cat(paste0("  up to ", abs(x$f1.range[2]), " non-event", plu2, " modified to be event", plu2, "\n"))
    }
    if(x$f1.range[2] == 0){
      cat(paste0("  up to ", abs(x$f1.range[1]), " event", plu1, " modified to be non-event", plu1, "\n"))
    }
  }
  if(any(x$f0.range != 0) & any(x$f1.range != 0)){
    cat("Range of event modification in group 0:\n")
    plu1 <- ifelse(abs(x$f0.range[1]) == 1, "", "s")
    plu2 <- ifelse(abs(x$f0.range[2]) == 1, "", "s")
    if(all(x$f0.range != 0)){
      cat(paste0("  up to ", abs(x$f0.range[1]), " event", plu1, " modified to be non-event", plu1, ";\n",
        "  up to ", abs(x$f0.range[2]), " non-event", plu2, " modified to be event", plu2, "\n"))
    }
    if(x$f0.range[1] == 0){
      cat(paste0("  up to ", abs(x$f0.range[2]), " non-event", plu2, " modified to be event", plu2, "\n"))
    }
    if(x$f0.range[2] == 0){
      cat(paste0("  up to ", abs(x$f0.range[1]), " event", plu1, " modified to be non-event", plu1, "\n"))
    }
    cat("Range of event modification in group 1:\n")
    plu1 <- ifelse(abs(x$f1.range[1]) == 1, "", "s")
    plu2 <- ifelse(abs(x$f1.range[2]) == 1, "", "s")
    if(all(x$f1.range != 0)){
      cat(paste0("  up to ", abs(x$f1.range[1]), " event", plu1, " modified to be non-event", plu1, ";\n",
        "  up to ", abs(x$f1.range[2]), " non-event", plu2, " modified to be event", plu2, "\n"))
    }
    if(x$f1.range[1] == 0){
      cat(paste0("  up to ", abs(x$f1.range[2]), " non-event", plu2, " modified to be event", plu2, "\n"))
    }
    if(x$f1.range[2] == 0){
      cat(paste0("  up to ", abs(x$f1.range[1]), " event", plu1, " modified to be non-event", plu1, "\n"))
    }
  }

  methods.name <- function(m){
    if(m == "Fisher") return("Fisher's exact test")
    if(m == "chisq") return("chi-squared test")
    if(m == "OR") return("odds ratio")
    if(m == "RR") return("relative risk")
    if(m == "RD") return("risk difference")
  }
  methods.name <- Vectorize(methods.name)
  cat("___________________________________\n")
  cat(paste0("Significance level varies from ", min(x$alphas),
    " to ", max(x$alphas), "\n"))
  if(any(is.element(c("OR", "RR", "RD"), x$methods))){
    null.val <- NULL
    if(is.element("OR", x$methods)){
      null.val <- c(null.val, paste0("OR = ", x$null["OR"]))
    }
    if(is.element("RR", x$methods)){
      null.val <- c(null.val, paste0("RR = ", x$null["RR"]))
    }
    if(is.element("RD", x$methods)){
      null.val <- c(null.val, paste0("RD = ", x$null["RD"]))
    }
    null.val <- paste(null.val, collapse = ", ")
    cat(paste0("Null hypothesis: ", null.val, "\n"))
  }

  if(x$alternative == "one.sided"){
    alt <- NULL
    for(i in 1:length(x$methods)){
      if(is.element(x$methods[i], c("Fisher", "chisq"))){
        alt <- c(alt, " (two-sided) ")
      }
      if(is.element(x$methods[i], c("OR", "RR", "RD"))){
        alt <- c(alt, " (one-sided) ")
      }
    }
  }else{
    alt <- rep(" ", length(x$methods))
  }
  cat(paste0("p-value", ifelse(x$alternative == "two.sided",
    " (two-sided)", ""), ":\n",
    paste(paste0("  ", format(round(x$pval, 3), nsmall = 3),
    alt, "based on ", methods.name(x$methods)), collapse = "\n"), "\n"))

  if(any(x$f0.range != 0) | any(x$f1.range != 0)){
    cat("___________________________________\n")
    cat("Fragility index (FI) and fragility quotient (FQ):\n")
    for(i in 1:length(x$methods)){
      temp.m <- x$methods[i]
      cat(paste0("Based on ", methods.name(temp.m), ",\n"))
      if(!is.na(x$FI.avg[temp.m])){
        cat(paste0("  Average FI = ", format(round(x$FI.avg[temp.m], 2), nsmall = 2),
          " (min = ", min(x$FI[,temp.m]),
          ", max = ", max(x$FI[,temp.m]), ");\n",
          "  Average FQ = ",
          format(round(100*x$FQ.avg[temp.m], 1), nsmall = 1), "%",
          " (min = ", format(round(100*min(x$FQ[,temp.m]), 1), nsmall = 1),
          "%, max = ", format(round(100*max(x$FQ[,temp.m]), 1), nsmall = 1),
          "%)\n"))

      }else{
        cat(paste0("Average FI and FQ = NA\n"))
      }
    }
  }
}