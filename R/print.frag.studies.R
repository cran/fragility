print.frag.studies <- function(x, ...){
  if(!inherits(x, "frag.studies")){
    stop("The input must be an object of \"frag.studies\".")
  }

  nstudy <- dim(x$pval)[1]
  cat("The input dataset contains", nstudy, "studies\n")

  methods.name <- function(m){
    if(m == "Fisher") return("Fisher's exact test")
    if(m == "chisq") return("chi-squared test")
    if(m == "OR") return("odds ratio")
    if(m == "RR") return("relative risk")
    if(m == "RD") return("risk difference")
  }
  methods.name <- Vectorize(methods.name)
  cat(paste0("Significance level = ", x$alpha, "\n"))
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
    " (two-sided)", ""), " is based on:\n",
    paste(paste0("  ", methods.name(x$methods)), alt, collapse = "\n")))

  if(x$modify0 != "none" | x$modify1 != "none"){
    cat("\n\nFragility index (FI) and fragility quotient (FQ):\n")
    for(i in 1:length(x$methods)){
      temp.m <- x$methods[i]
      study.sig <- sum(x$pval[,temp.m] < x$alpha)
      plu.sig <- ifelse(study.sig <= 1, "study", "studies")
      plus.sig <- ifelse(study.sig <= 1, "s", "")
      study.nonsig <- sum(x$pval[,temp.m] >= x$alpha)
      plu.nonsig <- ifelse(study.nonsig <= 1, "study", "studies")
      plus.nonsig <- ifelse(study.nonsig <= 1, "s", "")
      cat(paste0("Based on ", methods.name(temp.m), ",\n"))

      cat(paste0("  ", study.sig, " ", plu.sig, " yield", plus.sig,
        " significance"))
      if(study.sig == 0){
        cat(";\n")
      }
      if(study.sig < 5 & study.sig > 0){
        cat(" with\n")
        FI.sig <- x$FI[x$pval[,temp.m] < x$alpha, temp.m]
        FQ.sig <- x$FQ[x$pval[,temp.m] < x$alpha, temp.m]
        FQ.sig <- paste0(format(round(100*FQ.sig, 1), nsmall = 1), "%")
        FQ.sig <- gsub(" ", "", FQ.sig)
        FQ.sig[FQ.sig == "NA%"] <- NA
        cat(paste0("    FI = ", paste(FI.sig, collapse = ", "), " and\n"))
        cat(paste0("    FQ = ", paste(FQ.sig, collapse = ", "), ";\n"))
      }
      if(study.sig >= 5){
        FI.sig <- x$FI[x$pval[,temp.m] < x$alpha, temp.m]
        FQ.sig <- x$FQ[x$pval[,temp.m] < x$alpha, temp.m]
        FI.sig.smry <- quantile(FI.sig, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, type = 3)
        FQ.sig.smry <- quantile(FQ.sig, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, type = 3)
        FI.sig.smry <- format(round(FI.sig.smry), nsmall = 0)
        FI.sig.smry <- gsub(" ", "", FI.sig.smry)
        FQ.sig.smry <- paste0(format(round(100*FQ.sig.smry, 1), nsmall = 1), "%")
        FQ.sig.smry <- gsub(" ", "", FQ.sig.smry)
        cat(paste0(" with\n    median FI = ",
          FI.sig.smry[3], ", range ", FI.sig.smry[1], "-", FI.sig.smry[5],
          ", IQR ", FI.sig.smry[2], "-", FI.sig.smry[4], " and\n"))
        cat(paste0("    median FQ = ",
          FQ.sig.smry[3], ", range ", FQ.sig.smry[1], "-", FQ.sig.smry[5], 
          ", IQR ", FQ.sig.smry[2], "-", FQ.sig.smry[4]))
        if(all(!is.na(FI.sig))){
          cat(";\n")
        }else{
          cat(paste0("\n    while ", sum(is.na(FI.sig)),
            ifelse(sum(is.na(FI.sig)) <= 1, " study has ", " studies have "), "FI = FQ = NA;\n"))
        }
      }

      cat(paste0("  ", study.nonsig, " ", plu.nonsig, " yield", plus.nonsig,
        " non-significance"))
      if(study.nonsig == 0){
        cat("\n")
      }
      if(study.nonsig < 5 & study.nonsig > 0){
        cat(" with\n")
        FI.nonsig <- x$FI[x$pval[,temp.m] >= x$alpha, temp.m]
        FQ.nonsig <- x$FQ[x$pval[,temp.m] >= x$alpha, temp.m]
        FQ.nonsig <- paste0(format(round(100*FQ.nonsig, 1), nsmall = 1), "%")
        FQ.nonsig <- gsub(" ", "", FQ.nonsig)
        FQ.nonsig[FQ.nonsig == "NA%"] <- NA
        cat(paste0("    FI = ", paste(FI.nonsig, collapse = ", "), " and\n"))
        cat(paste0("    FQ = ", paste(FQ.nonsig, collapse = ", "), ";\n"))
      }
      if(study.nonsig >= 5){
        FI.nonsig <- x$FI[x$pval[,temp.m] >= x$alpha, temp.m]
        FQ.nonsig <- x$FQ[x$pval[,temp.m] >= x$alpha, temp.m]
        FI.nonsig.smry <- quantile(FI.nonsig, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, type = 3)
        FQ.nonsig.smry <- quantile(FQ.nonsig, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, type = 3)
        FI.nonsig.smry <- format(round(FI.nonsig.smry), nsmall = 0)
        FI.nonsig.smry <- gsub(" ", "", FI.nonsig.smry)
        FQ.nonsig.smry <- paste0(format(round(100*FQ.nonsig.smry, 1), nsmall = 1), "%")
        FQ.nonsig.smry <- gsub(" ", "", FQ.nonsig.smry)
        cat(paste0(" with\n    median FI = ",
          FI.nonsig.smry[3], ", range ", FI.nonsig.smry[1], "-", FI.nonsig.smry[5],
          ", IQR ", FI.nonsig.smry[2], "-", FI.nonsig.smry[4], " and\n"))
        cat(paste0("    median FQ = ",
          FQ.nonsig.smry[3], ", range ", FQ.nonsig.smry[1], "-", FQ.nonsig.smry[5], 
          ", IQR ", FQ.nonsig.smry[2], "-", FQ.nonsig.smry[4]))
        if(all(!is.na(FI.nonsig))){
          cat(";\n")
        }else{
          cat(paste0("\n    while ", sum(is.na(FI.nonsig)),
            ifelse(sum(is.na(FI.nonsig)) <= 1, " study has ", " studies have "), "FI = FQ = NA;\n"))
        }
      }

      cat("  overall, among all studies,\n")
      if(nstudy < 5){
        FI.overall <- x$FI[, temp.m]
        FQ.overall <- x$FQ[, temp.m]
        FQ.overall <- paste0(format(round(100*FQ.overall, 1), nsmall = 1), "%")
        FQ.overall <- gsub(" ", "", FQ.overall)
        FQ.overall[FQ.overall == "NA%"] <- NA
        cat(paste0("    FI = ", paste(FI.overall, collapse = ", "), " and\n"))
        cat(paste0("    FQ = ", paste(FQ.overall, collapse = ", "), "\n"))
      }
      if(nstudy >= 5){
        FI.overall <- x$FI[, temp.m]
        FQ.overall <- x$FQ[, temp.m]
        FI.overall.smry <- quantile(FI.overall, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, type = 3)
        FQ.overall.smry <- quantile(FQ.overall, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE, type = 3)
        FI.overall.smry <- format(round(FI.overall.smry), nsmall = 0)
        FI.overall.smry <- gsub(" ", "", FI.overall.smry)
        FQ.overall.smry <- paste0(format(round(100*FQ.overall.smry, 1), nsmall = 1), "%")
        FQ.overall.smry <- gsub(" ", "", FQ.overall.smry)
        cat(paste0("    median FI = ",
          FI.overall.smry[3], ", range ", FI.overall.smry[1], "-", FI.overall.smry[5],
          ", IQR ", FI.overall.smry[2], "-", FI.overall.smry[4], " and\n"))
        cat(paste0("    median FQ = ",
          FQ.overall.smry[3], ", range ", FQ.overall.smry[1], "-", FQ.overall.smry[5], 
          ", IQR ", FQ.overall.smry[2], "-", FQ.overall.smry[4]))
        if(all(!is.na(FI.overall))){
          cat("\n")
        }else{
          cat(paste0("\n    while ", sum(is.na(FI.overall)),
            ifelse(sum(is.na(FI.overall)) <= 1, " study has ", " studies have "), "FI = FQ = NA\n"))
        }
      }
    }
  }else{
    cat("No event modification occurs in both groups 0 and 1\n")
  }
}