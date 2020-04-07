\name{plot.frag.multi}
\alias{plot.frag.multi}
\title{
Plot Method for "frag.multi" Objects
}
\description{
Visualizes the fragility of multiple studies or multiple meta-analyses with binary outcomes.
}
\usage{
\method{plot}{frag.multi}(x, method, dir = "both", fragility = "FI", percentage = TRUE,
     max.f = NULL, bar, names.arg, space = 0, breaks, freq, reverse = FALSE,
     xlab, ylab, main = NULL, cex.marker, col.border, col.sig,
     trun.marker = TRUE, ...)
}
\arguments{
  \item{x}{an object of class \code{"frag.multi"} produced by the function \code{\link{frag.studies}} or \code{\link{frag.mas}}.}
  \item{method}{a character string indicating the method used to calculate the p-value when \code{x} is of class \code{"frag.studies"} produced by \code{\link{frag.studies}}. It must be an element of \code{x$methods}; the default is the first method in \code{x$methods}.}
  \item{dir}{a character string indicating the type of fragility measures to be considered in the plot with a certain direction of the significance change. It can be one of \code{"sig2nonsig"} (significance altered to non-significance), \code{"nonsig2sig"} (non-significance altered to significance), and \code{"both"} (both directions, the default).}
  \item{fragility}{a character string of either \code{"FI"} (fragility index, the default) or \code{"FQ"} (fragility quotient), indicating which fragility measure will be presented in the plot.}
  \item{percentage}{a logical value indicating whether the fragility quotient is presented in percentage (\code{TRUE}, the default) or not (\code{FALSE}). It is only used when the argument \code{fragility} = \code{"FQ"}.}
  \item{max.f}{a positive numeric value indicating the maximum fragility index/quotient to be presented in the plot. Fragility indexes/quotients larger than \code{max.f} (if any) will be truncated. The default is \code{NULL} (i.e., no truncation).}
  \item{bar}{a logical value indicating whether a bar plot (\code{TRUE}) or a histogram (\code{FALSE}) will be generated. The bar plot is only available for fragility index (the argument \code{fragility} = \code{"FI"}); the default is \code{bar} = \code{TRUE} in this case. For fragility quotient (\code{fragility} = \code{"FQ"}), \code{bar} is always set to \code{FALSE}.}
  \item{names.arg}{a vector of names to be plotted below each bar. It is only used in the bar plot (\code{bar} = \code{TRUE}) for fragility index (\code{fragility} = \code{"FI"}). The default includes values increasing from the minimum fragility index to the maximum fragility index (or \code{max.f}) by 1.}
  \item{space}{the amount of space (as a fraction of the average bar width) left before each bar. The default is 0. It is only used in the bar plot (\code{bar} = \code{TRUE}) for fragility index (\code{fragility} = \code{"FI"}).}
  \item{breaks}{the breaks on the x axis in the histogram. It is only used when \code{bar} = \code{FALSE}. It is passed to the function \code{\link[graphics]{hist}} with the same default value.}
  \item{freq}{a logical value indicating whether the y axis in the histogram presents frequencies (\code{TRUE}) or densities (\code{FALSE}). It is only used when \code{bar} = \code{FALSE}. It is passed to the function \code{\link[graphics]{hist}} with the same default value (\code{TRUE} if and only if breaks are equidistant; otherwise, \code{FALSE}).}
  \item{reverse}{a logical value indicating whether the groups of two directions of significance change. It is only used when \code{dir = "both"}. The default is \code{FALSE}, where the part corresponding to non-significance altered to significance is stacked at the bottom and that corresponding to significance altered to non-significance is stacked at the top in the bar plot or histogram.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for the y axis.}
  \item{main}{a main title for the plot.}
  \item{cex.marker}{the size of the text about the truncated fragility measures when the argument \code{trun.marker} = \code{TRUE} and \code{bar} = \code{FALSE}. The default is 1.}
  \item{col.border}{the color to be used for the border of the bars in the bar plot or the histogram.}
  \item{col.sig}{a vector of two colors for the two directions of significance change. The default includes \code{"forestgreen"} (for non-significance altered to significance) and \code{"firebrick"} (for significance altered to non-significance), with the opacity setting to \code{alpha.f} = 0.1 in \code{\link[grDevices]{adjustcolor}}.}
  \item{trun.marker}{a logical value indicating whether a text is placed at place of truncated fragility measures (the default is \code{TRUE}) in the histogram (\code{bar} = \code{FALSE}). The text gives information about the truncation.}
  \item{\dots}{other arguments that can be passed to \code{\link[graphics]{barplot}} (when \code{bar} = \code{TRUE}) or \code{\link[graphics]{hist}} (when \code{bar} = \code{FALSE}).}
}
\details{
A bar plot or a histogram is generated. In the bar plot, the x axis presents the values of fragility indexes among the multiple studies or multiple meta-analyses, and the y axis presents the corresponding frequencies (counts). In the histogram, the x axis presents the intervals of fragility indexes or fragility quotients, and the y axis presents the corresponding frequencies (counts) or densities.
}
\value{
None.
}
\references{
Lin L. "Factors that impact fragility index and their visualizations." Unpublished manuscript.
}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{frag.studies}} and \code{\link{frag.mas}} for assessing fragility of multiple studies and multiple meta-analyses, respectively.
}
%%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%}
\keyword{plot}
