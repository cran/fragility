print.frag.nma <- function(x, ...){
  if(!inherits(x, "frag.nma")){
    stop("The input must be an object of \"frag.nma\".")
  }

  cat(paste0("Original network meta-analysis (NMA) contains\n"))
  cat(paste0("  ", length(unique(x$data$sid)), " studies and ", length(unique(x$data$tid)), " treatments\n"))
  cat(paste0("Significance level = ", x$alpha, "\n"))
  cat(paste0("The effect size is ", x$measure, ifelse(is.element(x$measure, c("OR", "RR")), " (on a logarithmic scale)", ""), "\n"))
  cat(paste0("The null value of is ", x$null, "\n"))
  cat("Fragility index (FI):\n")
  print(x$FI)
  cat("Fragility quotient (FQ), based on the associated comparison:\n")
  print(x$FQ)
  cat("Fragility quotient (FQ), based on the total sample size in the NMA:\n")
  print(x$FQ.nma)
  cat("See the manual for details to retrieve more information.\n")
}