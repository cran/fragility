\name{print.frag.mas}
\alias{print.frag.mas}
\title{
Print Method for "frag.mas" Objects
}
\description{
Prints information about the fragility of multiple meta-analyses with binary outcomes.
}
\usage{
\method{print}{frag.mas}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{"frag.mas"} produced by the function \code{\link{frag.mas}}.}
  \item{\dots}{other arguments.}
}
\details{
  The output includes the number of meta-analyses, the pre-specified significance level, the null value of effect size, and the descriptive information about the fragility of the multiple meta-analyses, such as the median, range (from minimum to maximum), and interquartile range (IQR, from the first quartile to the third quartile) of fragility index and those of fragility quotient.
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
\code{\link{frag.mas}}
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%}
\keyword{print}