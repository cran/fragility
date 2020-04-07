print.frag.nma.alpha <- function(x, ...){
  if(!inherits(x, "frag.nma.alpha")){
    stop("The input must be an object of \"frag.nma.alpha\".")
  }

  cat(paste0("Original network meta-analysis (NMA) contains\n"))
  cat(paste0("  ", length(unique(x$data$sid)), " studies and ", length(unique(x$data$tid)), " treatments\n"))
  cat(paste0("Significance level varies from ", min(x$alphas), " to ", max(x$alphas), "\n"))
  cat(paste0("The effect size is ", x$measure, ifelse(is.element(x$measure, c("OR", "RR")), " (on a logarithmic scale)", ""), "\n"))
  cat(paste0("The null value of is ", x$null, "\n"))
  cat("Average fragility index (FI):\n")
  for(i in 1:length(x$FI.avg)){
    temp.comp <- names(x$FI.avg)[i]
    temp.trts <- unlist(strsplit(temp.comp, split = " vs. "))
    if(temp.trts[1] != temp.trts[2]){
      cat(paste0("  ", x$FI.avg[[i]], " for comparison ", temp.comp, "\n"))
    }
  }
  cat("Average fragility quotient (FQ), based on the associated comparison:\n")
  for(i in 1:length(x$FQ.avg)){
    temp.comp <- names(x$FQ.avg)[i]
    temp.trts <- unlist(strsplit(temp.comp, split = " vs. "))
    if(temp.trts[1] != temp.trts[2]){
      cat(paste0("  ", format(round(x$FQ.avg[[i]], 10), nsmall = 10), " for comparison ", temp.comp, "\n"))
    }
  }
  cat("Average fragility quotient (FQ), based on the total sample size in the NMA:\n")
  for(i in 1:length(x$FQ.nma.avg)){
    temp.comp <- names(x$FQ.nma.avg)[i]
    temp.trts <- unlist(strsplit(temp.comp, split = " vs. "))
    if(temp.trts[1] != temp.trts[2]){
      cat(paste0("  ", format(round(x$FQ.nma.avg[[i]], 10), nsmall = 10), " for comparison ", temp.comp, "\n"))
    }
  }
}