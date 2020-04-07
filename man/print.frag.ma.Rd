\name{print.frag.ma}
\alias{print.frag.ma}
\title{
Print Method for "frag.ma" Objects
}
\description{
Prints information about the fragility of a meta-analysis with a binary outcome.
}
\usage{
\method{print}{frag.ma}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{"frag.ma"} produced by the function \code{\link{frag.ma}}.}
  \item{\dots}{other arguments.}
}
\details{
  The output includes the characteristics of the meta-analysis (the number of studies, events, and samples), the estimated overall effect size, the information about significance (the pre-specified significance level, the null value, the confidence interval, and the p-value), and the fragility index and quotient.
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
\code{\link{frag.ma}}
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%}
\keyword{print}