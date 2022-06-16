\name{frag.studies}
\alias{frag.studies}
\title{
Assessing Fragility of Multiple Individual Studies With Binary Outcomes
}
\description{
Produces fragility indexes or fragility quotients for altering statistical significance of multiple individual studies with binary outcomes.
}
\usage{
frag.studies(e0, n0, e1, n1, data, methods,
             modify0 = "both", modify1 = "both", alpha = 0.05,
             alternative = "two.sided", OR = 1, RR = 1, RD = 0)
}
\arguments{
  \item{e0}{a numeric vector or the corresponding column name in the argument \code{data}, indicating event counts in multiple studies' group 0.}
  \item{n0}{a numeric vector or the corresponding column name in the argument \code{data}, indicating sample sizes in multiple studies' group 0.}
  \item{e1}{a numeric vector or the corresponding column name in the argument \code{data}, indicating event counts in multiple studies' group 1.}
  \item{n1}{a numeric vector or the corresponding column name in the argument \code{data}, indicating sample sizes in multiple studies' group 1.}
  \item{data}{an optional data frame containing the dataset of multiple individual studies with binary outcomes. If \code{data} is specified, the previous arguments, \code{e0}, \code{n0}, \code{e1}, and \code{n1}, should be specified as their corresponding column names in \code{data}.}
  \item{methods}{a character string or a vector of character strings indicating the method(s) used to calculate the p-value(s) of the association between treatment and outcome. It must consist of \code{"Fisher"} (Fisher's exact test), \code{"chisq"} (chi-squared test), \code{"OR"} (odds ratio), \code{"RR"} (relative risk), and/or \code{"RD"} (risk difference). The default is the vector consisting of all aforementioned five methods.}
  \item{modify0}{a character string indicating the event status modifications in group 0. It should be one of \code{"increase"} (increasing event numbers), \code{"decrease"} (decreasing event numbers), \code{"both"} (the default, modifying event numbers in both directions), and \code{"none"} (no modification).}
  \item{modify1}{a character string indicating the event status modifications in group 1. It has the same usage with \code{modify0}, with the default being \code{"both"}.}
  \item{alpha}{a numeric value between 0 and 1, indicating the statistical significance level (the default is 0.05).}
  \item{alternative}{a character string of either \code{"two.sided"} (the default) or \code{"one.sided"}, indicating the alternative hypothesis. It is only used for \code{"OR"}, \code{"RR"}, and \code{"RD"} in the argument \code{methods}.}
  \item{OR}{a numeric positive value indicating the value of odds ratio under the null hypothesis (the default is 1). It is used only if the argument \code{methods} includes \code{"OR"}.}
  \item{RR}{a numeric positive value indicating the value of relative risk under the null hypothesis (the default is 1). It is used only if the argument \code{methods} includes \code{"RR"}.}
  \item{RD}{a numeric value between \eqn{-}1 indicating the value of risk difference under the null hypothesis (the default is 0). It is used only if the argument \code{methods} includes \code{"RD"}.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  An object of classes \code{"frag.multi"} and \code{"frag.studies"}. The object is a list containing the following components:
  \item{methods}{methods used to calculate the p-value(s) of the association between treatment and outcome.}
  \item{alpha}{pre-specified statistical significance level.}
  \item{alternative}{alternative hypothesis.}
  \item{null}{values of odds ratio, relative risk, and risk difference under the null hypothesis; they are only used if \code{"OR"}, \code{"RR"}, and \code{"RD"} are included in the argument \code{methods}, respectively.}
  \item{modify0}{type of event status modifications in group 0.}
  \item{modify1}{type of event status modifications in group 1.}
  \item{pval}{p-value(s) produced by the method(s) included in the argument \code{methods} for each study.}
  \item{FI}{fragility indexes based on the method(s) among all studies.}
  \item{FQ}{fragility quotients based on the method(s) among all studies.}
  Of note, if both arguments \code{modify0} and \code{modify0} are \code{"none"} (i.e., no modification), the items \code{FI} and \code{FQ} in the above list are not produced, because fragility cannot be assessed without event status modifications.
}
\references{
Lin L, Chu H (2022). "Assessing and visualizing fragility of clinical results with binary outcomes in R using the fragility package." \emph{PLOS ONE}, \bold{17}(6), e0268754. <\doi{10.1371/journal.pone.0268754}>

Walsh M, Srinathan SK, McAuley DF, Mrkobrada M, Levine O, Ribic C, Molnar AO, Dattani ND, Burke A, Guyatt G, Thabane L, Walter SD, Pogue J, Devereaux PJ (2014). "The statistical significance of randomized controlled trial results is frequently fragile: a case for a Fragility Index." \emph{Journal of Clinical Epidemiology}, \bold{67}(6), 622--8. <\doi{10.1016/j.jclinepi.2013.10.019}>
}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{frag.study}} for assessing fragility of an individual study.
}
\examples{
## Load datasets of trials on antidepressant drugs
data(dat.ad)

## Assess fragility of first 50 trials on antidepressant drugs
out1 <- frag.studies(e0, n0, e1, n1, data = dat.ad[1:50,], methods = "OR")
out1

## Visualize the results using the bar plot of fragility indexes
plot(out1)

## Additional arguments that can be passed to barplot()
plot(out1, cex.name = 0.7, mgp = c(2, 0.5, 0),
  main = "Bar plot of fragility index")

## Reverse the bars representing two directions of significance change
plot(out1, reverse = TRUE)

## Add space between bars, suppress borders, and change colors of bars
plot(out1, space = 0.5, col.border = NA, col.sig = c("blue", "red"),
  cex.name = 0.7, mgp = c(2, 0.5, 0))
## More options, such as legends, can be added to the plot
legend("toprigh", fill = c("blue", "red"), border = NA,
  legend = c("Non-significance altered to significance",
  "Significance altered to non-significance"))

## Visualize the results using the histogram of fragility indexes
plot(out1, bar = FALSE)

## Present density, instead of frequencies (counts), in the histogram
plot(out1, bar = FALSE, freq = FALSE)

## Visualize the results using the histogram of fragility quotients
plot(out1, fragility = "FQ")

\donttest{
## Assess fragility of the complete datasets of 347 trials;
##  it takes longer time
out2 <- frag.studies(e0, n0, e1, n1, data = dat.ad)
out2

## Generate the bar plot for each method
plot(out2, method = "Fisher", cex.name = 0.7, mgp = c(2, 0.5, 0),
  ylim = c(0, 60), main = "Fisher's exact Test")
plot(out2, method = "chisq", cex.name = 0.7, mgp = c(2, 0.5, 0),
  ylim = c(0, 60), main = "Chi-squared test")
plot(out2, method = "OR", cex.name = 0.7, mgp = c(2, 0.5, 0),
  ylim = c(0, 60), main = "Odds ratio")
plot(out2, method = "RR", cex.name = 0.7, mgp = c(2, 0.5, 0),
  ylim = c(0, 60), main = "Relative risk")
plot(out2, method = "RD", cex.name = 0.7, mgp = c(2, 0.5, 0),
  ylim = c(0, 60), main = "Risk difference")

## Restrict the bar plot to studies with
##  significance altered to non-significance
plot(out2, dir = "sig2nonsig", method = "Fisher", cex.name = 0.7,
  mgp = c(2, 0.5, 0), main = "Fisher's exact Test")

## Restrict the bar plot to studies with
##  non-significance altered to significance
plot(out2, dir = "nonsig2sig", method = "Fisher", cex.name = 0.7,
  mgp = c(2, 0.5, 0), main = "Fisher's exact Test")

## Truncate the fragility index at 16
plot(out2, method = "Fisher", max.f = 16, cex.name = 0.7,
  mgp = c(2, 0.5, 0))

## A marker for the truncation will be presented in the histogram
plot(out2, method = "Fisher", bar = FALSE, max.f = 16,
  mgp = c(2, 0.5, 0))

## Visualize the fragility quotients of all trials
plot(out2, method = "Fisher", fragility = "FQ", max.f = 16,
  mgp = c(2, 0.5, 0))

## Use more breaks in the histogram
plot(out2, method = "Fisher", fragility = "FQ", max.f = 16,
  breaks = 20, mgp = c(2, 0.5, 0))
}
}
\keyword{fragility}
