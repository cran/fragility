\name{print.frag.study.alpha}
\alias{print.frag.study.alpha}
\title{
Print Method for "frag.study.alpha" Objects
}
\description{
Prints information about the fragility of an individual study with a binary outcome at different statistical significance levels.
}
\usage{
\method{print}{frag.study.alpha}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{"frag.study.alpha"} produced by the function \code{\link{frag.study.alpha}}.}
  \item{\dots}{other arguments.}
}
\details{
  The output consists of three parts: 1) the information about the original data, where the 2\eqn{\times}2 table and the ranges of event status modifications in both groups are presented; 2) the information about the significance, including the pre-specified significance level, null value(s) (if the odds ratio, relative risk, and/or risk difference are used), and the p-value(s) with the associated method and the alternative hypothesis; and 3) the information about the fragility, including the averages, minimums, and maximums of fragility index and fragility quotient based on each method considered across different statistical significance levels. If both groups have no event status modifications, the third part will not be printed because the fragility needs to be assessed based on certain event status modifications.
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
\code{\link{frag.study.alpha}}
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%}
\keyword{print}
