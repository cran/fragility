\name{dat.ns}
\alias{dat.ns}
\docType{data}
\title{
Meta-analyses on Nutrition Support
}
\description{
An illustrative example of meta-analyses with binary outcomes.
}
\usage{data("dat.ns")}
\format{
  A data frame with 564 meta-analyses on nutrition support, which contains the following 5 variables.
  \describe{
    \item{\code{ma.id}}{a numeric vector indicating meta-analysis indexes.}
    \item{\code{e0}}{a numeric vector indicating event counts in treatment group 0.}
    \item{\code{n0}}{a numeric vector indicating sample sizes in treatment group 0.}
    \item{\code{e1}}{a numeric vector indicating event counts in treatment group 1.}
    \item{\code{n1}}{a numeric vector indicating sample sizes in treatment group 1.}
  }
}
\details{
The meta-analyses are extracted from Feinberg et al.\sspace(2017) after excluding meta-analyses with non-binary outcomes. The selected 564 meta-analyses focus on different binary outcomes.
}
\source{
Feinberg J, Nielsen EE, Korang SK, Halberg Engell K, Nielsen MS, Zhang K, Didriksen M, Lund L, Lindahl N, Hallum S, Liang N, Xiong W, Yang X, Brunsgaard P, Garioud A, Safi S, Lindschou J, Kondrup J, Gluud C, Jakobsen JC (2017). "Nutrition support in hospitalised adults at nutritional risk." \emph{Cochrane Database of Systematic Reviews}, \bold{5}, Art.{\sspace}No.: CD011598. <\doi{10.1002/14651858.CD011598.pub2}>
}
\keyword{datasets}
