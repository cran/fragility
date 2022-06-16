\name{plot.frag.alpha}
\alias{plot.frag.alpha}
\title{
Plot Method for "frag.alpha" Objects
}
\description{
Visualizes the fragility at different statistical significance levels.
}
\usage{
\method{plot}{frag.alpha}(x, method, fragility = "FI", percentage = TRUE,
     xlab, ylab, xlim, ylim, cex.pts, col.line, col.pval, col.sig,
     lty.pval, lwd, lwd.pval, pch, pch.na, tid1, tid2, FQ.nma = FALSE, ...)
}
\arguments{
  \item{x}{an object of class \code{"frag.alpha"} produced by the function \code{\link{frag.study.alpha}}, \code{\link{frag.ma.alpha}}, or \code{\link{frag.nma.alpha}}.}
  \item{method}{a character string indicating the method used to calculate the p-value when \code{x} is of class \code{"frag.study.alpha"} produced by \code{\link{frag.study.alpha}}. It must be an element of \code{x$methods}; the default is the first method in \code{x$methods}.}
  \item{fragility}{a character string of either \code{"FI"} (fragility index, the default) or \code{"FQ"} (fragility quotient), indicating which fragility measure will be presented in the plot.}
  \item{percentage}{a logical value indicating whether the fragility quotient is presented in percentage (\code{TRUE}, the default) or not (\code{FALSE}). It is only used when the argument \code{fragility} = \code{"FQ"}.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for the y axis.}
  \item{xlim}{the x limits \code{c(x1, x2)} of the plot.}
  \item{ylim}{the y limits \code{c(y1, y2)} of the plot.}
  \item{cex.pts}{the size of points in the plot (the default is 0.5).}
  \item{col.line}{the color of the line segments that join the points (the default is \code{"gray50"}).}
  \item{col.pval}{the color of the vertical line depicting the p-value of the original data (the default is \code{"gray50"}). It is only used when the p-value is between the lower and upper limits of the considered significance levels.}
  \item{col.sig}{a vector of the two colors for two directions of significance change. The default includes \code{"forestgreen"} (for non-significance altered to significance) and \code{"firebrick"} (for significance altered to non-significance).}
  \item{lty.pval}{the type of the vertical line depicting the p-value of the original data. The default is 2 (dashed). It is only used when the p-value is between the lower and upper limits of the considered significance levels.}
  \item{lwd}{the width of the line segments that join the points (the default is 1).}
  \item{lwd.pval}{the width of the vertical line depicting the p-value of the original data. The default is 1. It is only used when the p-value is between the lower and upper limits of the considered significance levels.}
  \item{pch}{the symbol of points in the plot. The default is 16 (filled circle).}
  \item{pch.na}{the symbol of points depicting fragility index/quotient = \code{NA} (if any, in cases that significance or non-significance cannot be altered). The default is 1 (unfilled circle).}
  \item{tid1}{the first treatment of the comparison of interest; if not specified, the first comparison in \code{x$tid.f} will be used. It is only used when \code{x} is of class \code{"frag.nma.alpha"} for a network meta-analysis.}
  \item{tid2}{the second treatment of the comparison of interest. Only one comparison can be specified by \code{tid1} and \code{tid2} at one time for plotting. It is only used when \code{x} is of class \code{"frag.nma.alpha"}.}
  \item{FQ.nma}{a logical value indicating whether the fragility quotient in a network meta-analysis is the total sample size associated with the specific treatment comparison (\code{FALSE}, the default) or the total sample size in the whole network meta-analysis across all treatment groups (\code{TRUE}). It is only used when \code{x} is of class \code{"frag.nma.alpha"}.}
  \item{\dots}{other arguments that can be passed to \code{\link[graphics]{plot.default}}.}
}
\details{
A scatter plot is generated, where the x axis presents the different significance levels, and the y axis presents the corresponding fragility indexes or fragility quotients.
}
\value{
None.
}
\references{
Lin L (2021). "Factors that impact fragility index and their visualizations." \emph{Journal of Evaluation in Clinical Practice}, \bold{27}(2), 356--64. <\doi{10.1111/jep.13428}>

Lin L, Chu H (2022). "Assessing and visualizing fragility of clinical results with binary outcomes in R using the fragility package." \emph{PLOS ONE}, \bold{17}(6), e0268754. <\doi{10.1371/journal.pone.0268754}>
}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{frag.study.alpha}}, \code{\link{frag.ma.alpha}}, and \code{\link{frag.nma.alpha}} for assessing fragility of an individual study, a meta-analysis, and a network meta-analysis at multiple significance levels, respectively.
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%}
\keyword{plot}
