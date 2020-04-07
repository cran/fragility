print.frag.mas <- function(x, ...){
  if(!inherits(x, "frag.mas")){
    stop("The input must be an object of \"frag.mas\".")
  }

  n.ma <- length(x$FI)
  cat("The input dataset contains", n.ma, "meta-analyses\n")

  cat(paste0("Significance level = ", x$alpha, "\n"))
  cat(paste0("The effect size is ", x$measure, ifelse(is.element(x$measure, c("OR", "RR")), " (on a logarithmic scale)", ""), "\n"))
  cat(paste0("The null value of is ", x$null, "\n"))

  cat("\nFragility index (FI) and fragility quotient (FQ):\n")

  ma.sig <- sum(x$pval.ori < x$alpha)
  plu.sig <- ifelse(ma.sig <= 1, "meta-analysis", "meta-analyses")
  plus.sig <- ifelse(ma.sig <= 1, "s", "")
  ma.nonsig <- sum(x$pval.ori >= x$alpha)
  plu.nonsig <- ifelse(ma.nonsig <= 1, "meta-analysis", "meta-analyses")
  plus.nonsig <- ifelse(ma.nonsig <= 1, "s", "")

  cat(paste0("  ", ma.sig, " ", plu.sig, " yield", plus.sig, " significance"))
  if(ma.sig == 0){
    cat(";\n")
  }
  if(ma.sig < 5 & ma.sig > 0){
    cat(" with\n")
    FI.sig <- x$FI[x$pval.ori < x$alpha]
    FQ.sig <- x$FQ[x$pval.ori < x$alpha]
    FQ.sig <- paste0(format(round(100*FQ.sig, 1), nsmall = 1), "%")
    FQ.sig <- gsub(" ", "", FQ.sig)
    FQ.sig[FQ.sig == "NA%"] <- NA
    cat(paste0("    FI = ", paste(FI.sig, collapse = ", "), " and\n"))
    cat(paste0("    FQ = ", paste(FQ.sig, collapse = ", "), ";\n"))
  }
  if(ma.sig >= 5){
    FI.sig <- x$FI[x$pval.ori < x$alpha]
    FQ.sig <- x$FQ[x$pval.ori < x$alpha]
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
        ifelse(sum(is.na(FI.sig)) <= 1, " meta-analysis has ", " meta-analyses have "), "FI = FQ = NA;\n"))
    }
  }

  cat(paste0("  ", ma.nonsig, " ", plu.nonsig, " yield", plus.nonsig, " non-significance"))
  if(ma.nonsig == 0){
    cat("\n")
  }
  if(ma.nonsig < 5 & ma.nonsig > 0){
    cat(" with\n")
    FI.nonsig <- x$FI[x$pval.ori >= x$alpha]
    FQ.nonsig <- x$FQ[x$pval.ori >= x$alpha]
    FQ.nonsig <- paste0(format(round(100*FQ.nonsig, 1), nsmall = 1), "%")
    FQ.nonsig <- gsub(" ", "", FQ.nonsig)
    FQ.nonsig[FQ.nonsig == "NA%"] <- NA
    cat(paste0("    FI = ", paste(FI.nonsig, collapse = ", "), " and\n"))
    cat(paste0("    FQ = ", paste(FQ.nonsig, collapse = ", "), ";\n"))
  }
  if(ma.nonsig >= 5){
    FI.nonsig <- x$FI[x$pval.ori >= x$alpha]
    FQ.nonsig <- x$FQ[x$pval.ori >= x$alpha]
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
        ifelse(sum(is.na(FI.nonsig)) <= 1, " meta-analysis has ", " meta-analyses have "), "FI = FQ = NA;\n"))
    }
  }

  cat("  overall, among all meta-analyses,\n")
  if(n.ma < 5){
    FI.overall <- x$FI
    FQ.overall <- x$FQ
    FQ.overall <- paste0(format(round(100*FQ.overall, 1), nsmall = 1), "%")
    FQ.overall <- gsub(" ", "", FQ.overall)
    FQ.overall[FQ.overall == "NA%"] <- NA
    cat(paste0("    FI = ", paste(FI.overall, collapse = ", "), " and\n"))
    cat(paste0("    FQ = ", paste(FQ.overall, collapse = ", "), "\n"))
  }
  if(n.ma >= 5){
    FI.overall <- x$FI
    FQ.overall <- x$FQ
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
        ifelse(sum(is.na(FI.overall)) <= 1, " meta-analysis has ", " meta-analyses have "), "FI = FQ = NA\n"))
    }
  }
}