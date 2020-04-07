print.frag.ma.alpha <- function(x, ...){
  if(!inherits(x, "frag.ma.alpha")){
    stop("The input must be an object of \"frag.ma.alpha\".")
  }

  cat(paste0("Original meta-analysis contains\n"))
  cat(paste0("  ", dim(x$data)[1], " studies;\n"))
  cat(paste0("  ", format(round(sum(x$data$e0)), scientific = FALSE, big.mark = ","), " total events and ",
    format(round(sum(x$data$n0)), scientific = FALSE, big.mark = ","), " total sample sizes in group 0;\n"))
  cat(paste0("  ", format(round(sum(x$data$e1)), scientific = FALSE, big.mark = ","), " total events and ",
    format(round(sum(x$data$n1)), scientific = FALSE, big.mark = ","), " total sample sizes in group 1\n"))
  cat(paste0("Significance level varies from ", min(x$alphas), " to ", max(x$alphas), "\n"))
  cat(paste0("The effect size is ", x$measure, ifelse(is.element(x$measure, c("OR", "RR")), " (on a logarithmic scale)", ""), "\n"))
  cat(paste0("The null value of is ", x$null, "\n"))
  cat(paste0("The estimated overall effect size is\n"))
  cat(paste0("  ", format(round(x$est.ori, 3), nsmall = 3), " with p-value ",
    format(round(x$pval.ori, 3), nsmall = 3), "\n"))
  if(!is.na(x$FI.avg)){
    cat(paste0("Average fragility index (FI) = ", format(round(x$FI.avg, 2), nsmall = 2),
      " (min = ", min(x$FI), ", max = ", max(x$FI), ")\n",
      "Average fragility quotient (FQ) = ", format(round(100*x$FQ.avg, 1), nsmall = 1), "%",
      " (min = ", format(round(100*min(x$FQ), 1), nsmall = 1), "%, max = ", format(round(100*max(x$FQ), 1), nsmall = 1), "%)\n"))
  }else{
    cat(paste0("Average fragility index (FI) and fragility quotient (FQ) = NA\n"))
  }
}