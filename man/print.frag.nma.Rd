\name{print.frag.nma}
\alias{print.frag.nma}
\title{
Print Method for "frag.nma" Objects
}
\description{
Prints information about the fragility of a network meta-analysis with a binary outcome.
}
\usage{
\method{print}{frag.nma}(x, ...)
}
\arguments{
  \item{x}{an object of class \code{"frag.nma"} produced by the function \code{\link{frag.nma}}.}
  \item{\dots}{other arguments.}
}
\details{
  Due to space limit, the output only includes the basic characteristics of the network meta-analysis (the numbers of studies and treatments), the effect size, and the brief information about fragility index and quotient. More information can be retrieved from the list object \code{x}. For example, the estimated effect sizes can be retrieved from \code{x$est.ori}, and the directions of significance change corresponding to the fragility index/quotient can be retrieved from \code{x$dir}.
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
\code{\link{frag.nma}}
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%}
\keyword{print}