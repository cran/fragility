\name{dat.copd}
\alias{dat.copd}
\docType{data}
\title{
Network Meta-Analysis of Chronic Obstructive Pulmonary Disease
}
\description{
An illustrative example of network meta-analysis with a binary outcome.
}
\usage{data("dat.copd")}
\format{
  A data frame with 3 studies and 4 treatments, with the following 4 variables.
  \describe{
    \item{\code{sid}}{a numeric vector indicating study IDs.}
    \item{\code{tid}}{a numeric vector indicating treatment IDs.}
    \item{\code{e}}{a numeric vector indicating event counts.}
    \item{\code{n}}{a numeric vector indicating sample sizes.}
  }
}
\details{
The dataset was reported in Table 1 of Woods et al.\sspace(2010), and was also used as an illustrative example in the \CRANpkg{netmeta}. It investigates the effects of four treatments on chronic obstructive pulmonary disease (COPD), including 1) placebo; 2) fluticasone; 3) salmeterol; and 4) salmeterol fluticasone combination (SFC). Of note, Woods et al.\sspace(2010) originally used this dataset to illustrate the usage of network meta-analysis of count and hazard ratio statistics; this package only focuses on the event counts to show the relevant functions' usage for binary outcomes.
}
\source{
Woods BS, Hawkins N, Scott DA (2010). "Network meta-analysis on the log-hazard scale, combining count and hazard ratio statistics accounting for multi-arm trials: a tutorial." \emph{BMC Medical Research Methodology}, \bold{10}, 54. <\doi{10.1186/1471-2288-10-54}>
}
\keyword{datasets}
