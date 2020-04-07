print.frag.ma <- function(x, ...){
  if(!inherits(x, "frag.ma")){
    stop("The input must be an object of \"frag.ma\".")
  }

  cat(paste0("Original meta-analysis contains\n"))
  cat(paste0("  ", dim(x$data)[1], " studies;\n"))
  cat(paste0("  ", format(round(sum(x$data$e0)), scientific = FALSE, big.mark = ","), " total events and ",
    format(round(sum(x$data$n0)), scientific = FALSE, big.mark = ","), " total sample sizes in group 0;\n"))
  cat(paste0("  ", format(round(sum(x$data$e1)), scientific = FALSE, big.mark = ","), " total events and ",
    format(round(sum(x$data$n1)), scientific = FALSE, big.mark = ","), " total sample sizes in group 1\n"))
  cat(paste0("Significance level = ", x$alpha, "\n"))
  cat(paste0("The effect size is ", x$measure, ifelse(is.element(x$measure, c("OR", "RR")), " (on a logarithmic scale)", ""), "\n"))
  cat(paste0("The null value of is ", x$null, "\n"))
  cat(paste0("The estimated overall effect size is\n"))
  cat(paste0("  ", format(round(x$est.ori, 3), nsmall = 3), " with CI (",
    format(round(x$ci.ori[1], 3), nsmall = 3), ", ", format(round(x$ci.ori[2], 3), nsmall = 3), ") and p-value ",
    format(round(x$pval.ori, 3), nsmall = 3), "\n"))
  if(!is.na(x$FI)){
    cat(paste0("Fragility index (FI) = ", x$FI, " and fragility quotient (FQ) = ", format(round(100*x$FQ, 1), nsmall = 1), "%\n"))
    cat(paste0("  for ", x$dir, "\n"))
  }else{
    cat(paste0("FI = FQ = NA, i.e.,\n  ", x$dir, "\n"))
  }
}