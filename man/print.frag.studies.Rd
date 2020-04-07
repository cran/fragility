\name{print.frag.studies}
\alias{print.frag.studies}
\title{
Print Method for "frag.studies" Objects
}
\description{
Prints information about the fragility of multiple individual studies with binary outcomes.
}
\usage{
\method{print}{frag.studies}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{"frag.studies"} produced by the function \code{\link{frag.studies}}.}
  \item{\dots}{other arguments.}
}
\details{
  The output consists of two parts: 1) the information about the number of studies, the pre-specified significance level, null value(s) (if the odds ratio, relative risk, and/or risk difference are used), the method(s) used to produce p-value(s), and the alternative hypothesis; and 2) the information about the fragility of the multiple studies, including the median, range (from minimum to maximum), and interquartile range (IQR, from first quartile to third quartile) of fragility index and those of fragility quotient based on each method considered. The second part about the fragility presents the information grouped by the direction of significance change (significance altered to non-significance and non-significance altered to significance) as well as the overall (ungrouped) information. If both groups have no event status modifications, the second part will not be printed, because the fragility needs to be assessed based on certain event status modifications.
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
\code{\link{frag.studies}}
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%}
\keyword{print}
