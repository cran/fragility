\name{print.frag.study}
\alias{print.frag.study}
\title{
Print Method for "frag.study" Objects
}
\description{
Prints information about the fragility of an individual study with a binary outcome.
}
\usage{
\method{print}{frag.study}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{"frag.study"} produced by the function \code{\link{frag.study}}.}
  \item{\dots}{other arguments.}
}
\details{
  The output consists of three parts: 1) the information about the original data, where the 2\eqn{\times}2 table and the ranges of event status modifications in both groups are presented; 2) the information about the significance, including the pre-specified significance level, null value(s) (if the odds ratio, relative risk, and/or risk difference are used), and the p-value(s) with the associated method and the alternative hypothesis; and 3) the information about the fragility, including the fragility index and fragility quotient based on each method considered, the direction of significance change, and the corresponding minimal event status modification(s) for altering significance or non-significance. If both groups have no event status modifications, the third part will not be printed because the fragility needs to be assessed based on certain event status modifications.
}
\value{
None.
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{frag.study}}
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%}
\keyword{print}
