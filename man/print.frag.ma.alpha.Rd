\name{print.frag.ma.alpha}
\alias{print.frag.ma.alpha}
\title{
Print Method for "frag.ma.alpha" Objects
}
\description{
Prints information about the fragility of a meta-analysis with a binary outcome at different statistical significance levels.
}
\usage{
\method{print}{frag.ma.alpha}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{"frag.ma.alpha"} produced by the function \code{\link{frag.ma.alpha}}.}
  \item{\dots}{other arguments.}
}
\details{
  The output includes the characteristics of the meta-analysis (the number of studies, events, and samples), the estimated overall effect size, the information about significance (the range of considered significance levels, the null value, and the p-value), and the averages, minimums, and maximums of fragility index and fragility quotient across different statistical significance levels.
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
\code{\link{frag.ma.alpha}}
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%}
\keyword{print}