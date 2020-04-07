\name{print.frag.nma.alpha}
\alias{print.frag.nma.alpha}
\title{
Print Method for "frag.nma.alpha" Objects
}
\description{
Prints information about the fragility of treatment comparison(s) in a network meta-analysis with a binary outcome at different statistical significance levels.
}
\usage{
\method{print}{frag.nma.alpha}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{"frag.nma.alpha"} produced by the function \code{\link{frag.nma.alpha}}.}
  \item{\dots}{other arguments.}
}
\details{
  The output includes the basic characteristics of the network meta-analysis (the numbers of studies and treatments), the effect size, and the brief information about the averages of fragility indexes and fragility quotients across different statistical significance levels. More information can be retrieved from the list object \code{x}. For example, the fragility indexes at different significance levels can be retrieved from \code{x$FI}.
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
\code{\link{frag.nma.alpha}}
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%}
\keyword{print}