\name{dat.ad}
\alias{dat.ad}
\docType{data}
\title{
Randomized Controlled Trials of Antidepressant Drugs
}
\description{
An illustrative example of clinical studies with binary outcomes.
}
\usage{data("dat.ad")}
\format{
  A data frame with 347 randomized controlled trials of antidepressant drugs, which contains the following 4 variables.
  \describe{
    \item{\code{e0}}{a numeric vector indicating event counts in treatment group 0.}
    \item{\code{n0}}{a numeric vector indicating sample sizes in treatment group 0.}
    \item{\code{e1}}{a numeric vector indicating event counts in treatment group 1.}
    \item{\code{n1}}{a numeric vector indicating sample sizes in treatment group 1.}
  }
}
\details{
The binary outcome is acceptability (dropout due to any cause). The 347 trials are selected from Cipriani et al.\sspace(2018) after excluding trials with missing data or with more than two treatment arms.
}
\source{
Cipriani A, Furukawa TA, Salanti G, Chaimani A, Atkinson LZ, Ogawa Y, Leucht S, Ruhe HG, Turner EH, Higgins JPT, Egger M, Takeshima N, Hayasaka Y, Imai H, Shinohara K, Tajika A, Ioannidis JPA, Geddes JR (2018). "Comparative efficacy and acceptability of 21 antidepressant drugs for the acute treatment of adults with major depressive disorder: a systematic review and network meta-analysis." \emph{Lancet}, \bold{391}(10128), 1357--66. <\doi{10.1016/S0140-6736(17)32802-7}>
}
\keyword{datasets}
