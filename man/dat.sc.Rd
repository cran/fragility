\name{dat.sc}
\alias{dat.sc}
\docType{data}
\title{
Network Meta-Analysis of Smoking Cessation
}
\description{
An illustrative example of network meta-analysis with a binary outcome.
}
\usage{data("dat.sc")}
\format{
  A data frame containing 24 studies and 4 treatments, with the following 4 variables.
  \describe{
    \item{\code{sid}}{a numeric vector indicating study IDs.}
    \item{\code{tid}}{a numeric vector indicating treatment IDs.}
    \item{\code{e}}{a numeric vector indicating event counts.}
    \item{\code{n}}{a numeric vector indicating sample sizes.}
  }
}
\details{
The dataset was originally reported by Hasselblad (1998) and reanalyzed by Lu and Ades (2006) using formal network meta-analysis methods. It investigates the effects of four treatments on smoking cessation, including 1) no contact; 2) self-help; 3) individual counseling; and 4) group counseling. The outcome is binary, which is smoking cessation status (yes or no) of a participant after treatment.
}
\source{
Hasselblad V (1998). "Meta-analysis of multitreatment studies." \emph{Medical Decision Making}, \bold{18}(1), 37--43. <\doi{10.1177/0272989X9801800110}>

Lu G, Ades AE (2006). "Assessing evidence inconsistency in mixed treatment comparisons." \emph{Journal of the American Statistical Association}, \bold{101}(474), 447--59. <\doi{10.1198/016214505000001302}>
}
\keyword{datasets}
