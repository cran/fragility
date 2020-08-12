\name{plot.frag.study.all}
\alias{plot.frag.study.all}
\title{
Plot Method for "frag.study.all" Objects
}
\description{
Visualizes the fragility of an individual study with a binary outcome.
}
\usage{
\method{plot}{frag.study.all}(x, method, modify0, modify1, trun,
     xlab, ylab, xlim, ylim, cex.pts, cex.legend.pval, cex.legend.title,
     col.ori, col.ori.hl, col.f.hl, col.sig, lty.ori, lwd.ori,
     pch, pch.ori, pch.ori.hl, pch.f, pch.f.hl, pch.trun,
     adjust.legend, adjust.seg, legend.pvals, \dots)
}
\arguments{
  \item{x}{an object of class \code{"frag.study.all"}. It is produced by the function \code{\link{frag.study}} with its argument \code{all} = \code{TRUE}.}
  \item{method}{a character string (\code{"Fisher"}, \code{"chisq"}, \code{"OR"}, \code{"RR"}, or \code{"RD"}) indicating the method used to calculate the p-value. It must be an element of \code{x$methods}; the default is the first method in \code{x$methods}.}
  \item{modify0}{a logical value indicating whether event status is modified in group 0 for plotting. The default is \code{TRUE} if the range of event status modifications stored in the object \code{x} (\code{x$f0.range}) is not 0; otherwise, the default is \code{FALSE}.}
  \item{modify1}{a logical value indicating whether event status is modified in group 1 for plotting. The default is \code{TRUE} if the range of event status modifications stored in the object \code{x} (\code{x$f1.range}) is not 0; otherwise, the default is \code{FALSE}.}
  \item{trun}{a positive numeric value indicating truncation of p-value (on a base-10 logarithmic scale); p-values smaller than this threshold (10^\code{-trun}) are truncated. The default is 10.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for the y axis.}
  \item{xlim}{the x limits \code{c(x1, x2)} of the plot.}
  \item{ylim}{the y limits \code{c(y1, y2)} of the plot.}
  \item{cex.pts}{the size of points in the plot (the default is 0.5).}
  \item{cex.legend.pval}{the text size of p-values in the legend (the default is 0.6). It is only used when both arguments \code{modify0} and \code{modify1} are \code{TRUE}.}
  \item{cex.legend.title}{the size of the legend title (the default is 1). It is only used when both arguments \code{modify0} and \code{modify1} are \code{TRUE}.}
  \item{col.ori}{the color of the line(s) depicting the original data (without event status modifications). The default is \code{"gray50"}.}
  \item{col.ori.hl}{the color of the point depicting the original data (without event status modifications). The default is \code{"black"}. It is only used when both arguments \code{modify0} and \code{modify1} are \code{TRUE}.}
  \item{col.f.hl}{the color of the point(s) for highlighting the minimal event status modifications for altering statistical significance or non-significance. The default is \code{"black"}. It is only used when both arguments \code{modify0} and \code{modify1} are \code{TRUE}.}
  \item{col.sig}{a vector of two colors for non-significant and significant results, accordingly. The default includes \code{"forestgreen"} (for non-significant results) and \code{"firebrick"} (for significant results). When both arguments \code{modify0} and \code{modify1} are \code{TRUE}, they serve as colors of points representing p-values based on different event status modifications, with opacity automatically adjusted for depicting the p-value magnitudes (based on the function \code{\link[grDevices]{adjustcolor}} with different opacity levels for its argument \code{alpha.f}). When only one of the arguments \code{modify0} and \code{modify1} is \code{TRUE}, they serve as the background colors for distinguishing non-significance and significance, with the default opacity \code{alpha.f = 0.1} in \code{\link[grDevices]{adjustcolor}}.}
  \item{lty.ori}{the type of the line(s) depicting the original data (without event status modifications). The default is 2 (dashed).}
  \item{lwd.ori}{the width of the line(s) depicting the original data (without event status modifications). The default is 1.}
  \item{pch}{the symbol of the points in the plot. The default is 16 (filled circle) when both arguments \code{modify0} and \code{modify1} are \code{TRUE}, and is 1 (circle) when only one of the arguments \code{modify0} and \code{modify1} is \code{TRUE}.}
  \item{pch.ori}{the symbol of the point depicting the original data (without event status modifications). The default is 15 (filled square).}
  \item{pch.ori.hl}{the symbol of the point for highlighting the original data (without event status modifications). The default is 0 (square). It is only used when both arguments \code{modify0} and \code{modify1} are \code{TRUE}.}
  \item{pch.f}{the symbol of the point depicting the original data (without event status modifications). The default is 15 (filled square).}
  \item{pch.f.hl}{the symbol of the point(s) for highlighting the minimal event status modifications for altering statistical significance or non-significance. The default is 2 (triangle point up). It is only used when both arguments \code{modify0} and \code{modify1} are \code{TRUE}.}
  \item{pch.trun}{the symbol of the point(s) depicting truncated p-values. The default is 3 (plus). It is only used when only one of the arguments \code{modify0} and \code{modify1} is \code{TRUE}.}
  \item{adjust.legend}{a positive numeric value for adjusting the width of the legend of p-values. The default is 1.}
  \item{adjust.seg}{a positive integer for adjusting the number of segments in the legend of p-values. The default is 10.}
  \item{legend.pvals}{a numeric value or a vector of numeric values that give additional p-values (e.g., 0.005 and 0.1) to be shown in the legend. The default is \code{NULL}.}
  \item{\dots}{other arguments that can be passed to \code{\link[graphics]{plot.default}}.}
}
\details{
When both arguments \code{modify0} and \code{modify1} are \code{TRUE}, the generated plot presents p-values (with different colors representing their magnitudes) based on all possible event status modifications. The modifications in group 0 are presented on the x axis, and those in group 1 are presented on the y axis. When only one of the arguments \code{modify0} and \code{modify1} is \code{TRUE}, a scatter plot is generated, which presents p-values (on a base-10 logarithmic scale) on the y axis against event status modifications in group 0 (if \code{modify0} = \code{TRUE}) or group 1 (if \code{modify1} = \code{TRUE}) on the x axis.
}
\value{
None.
}
\references{
Lin L. "Factors that impact fragility index and their visualizations." \emph{Journal of Evaluation in Clinical Practice}, in press. <\doi{10.1111/jep.13428}>
}
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
\keyword{plot}
