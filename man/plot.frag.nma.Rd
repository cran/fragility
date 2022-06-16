\name{plot.frag.nma}
\alias{plot.frag.nma}
\title{
Plot Method for "frag.nma" Objects
}
\description{
Visualizes the iterative process of event status modifications for deriving the fragility index/quotient of a network meta-analysis with a binary outcome.
}
\usage{
\method{plot}{frag.nma}(x, tid1, tid2, xlab, ylab, xlim, ylim, ybreaks = NULL,
     study.marker = TRUE, cex.marker, offset.marker, col.line, lwd,
     legend, x.legend, y.legend, cex.legend, ...)
}
\arguments{
  \item{x}{an object of class \code{"frag.nma"} produced by the function \code{\link{frag.nma}}.}
  \item{tid1}{the first treatment of the comparison of interest; if not specified, the first comparison in \code{x$tid.f} will be used.}
  \item{tid2}{the second treatment of the comparison of interest. Only one comparison can be specified by \code{tid1} and \code{tid2} at one time for plotting.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for the y axis.}
  \item{xlim}{the x limits \code{c(x1, x2)} of the plot.}
  \item{ylim}{the y limits \code{c(y1, y2)} of the plot.}
  \item{ybreaks}{a vector of the two numeric values for breaking the y axis or \code{NULL} (the default, not breaking the y axis). The two values should be between the total event counts in the two treatment groups specified by \code{tid1} and \code{tid2}, and the first value should be less than the second value. It is useful when the total event counts in the two groups differ a lot, leading to a large distance between the two lines representing the event status modifications in the two groups.}
  \item{study.marker}{a logical value indicating whether study labels involved in event status modifications are presented in the plot (the default is \code{TRUE}). When using the default, an asterisk represents that the study with an event status modification remains the same as in the previous iteration.}
  \item{cex.marker}{the size of study labels (when \code{study.marker} = \code{TRUE}) with the default 0.8.}
  \item{offset.marker}{a numeric value controlling the distance between study labels from the lines (the default is 0.2). It is passed to the function \code{\link[graphics]{text}}.}
  \item{col.line}{a vector of the two colors for two lines representing the event status modifications in the two treatment groups. The default includes \code{"blue"} (for \code{tid1}) and \code{"red"} (for \code{tid2}).}
  \item{lwd}{the width of the lines (the default is 1).}
  \item{legend}{the legend of plot. The default gives the two groups.}
  \item{x.legend}{the x co-ordinate or a keyword, such as \code{"right"} (the default), to be used to position the legend. It is passed to \code{\link[graphics]{legend}}.}
  \item{y.legend}{the y co-ordinate to be used to position the legend (the default is \code{NULL}).}
  \item{cex.legend}{the size of legend text (the default is 1).}
  \item{\dots}{other arguments that can be passed to \code{\link[graphics]{plot.default}}.}
}
\details{
A plot showing the iterative process of event status modifications is generated for a specific treatment comparison. The x axis presents the iterations, and the y axis gives the group-specific total event counts.
}
\value{
None.
}
\references{
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
\code{\link{frag.nma}}
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%}
\keyword{plot}