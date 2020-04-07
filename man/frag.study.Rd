\name{frag.study}
\alias{frag.study}
\title{
Assessing Fragility of an Individual Study With a Binary Outcome
}
\description{
Produces results, such as fragility index, fragility quotient, and event status modifications for altering statistical significance, of an individual study with a binary outcome.
}
\usage{
frag.study(e0, n0, e1, n1, data, all = FALSE, methods,
           modify0 = "both", modify1 = "both", alpha = 0.05,
           alternative = "two.sided", OR = 1, RR = 1, RD = 0,
           allcase = TRUE)
}
\arguments{
  \item{e0}{a numeric value or the corresponding column name in the argument \code{data}, indicating event count in group 0.}
  \item{n0}{a numeric value or the corresponding column name in the argument \code{data}, indicating sample size in group 0.}
  \item{e1}{a numeric value or the corresponding column name in the argument \code{data}, indicating event count in group 1.}
  \item{n1}{a numeric value or the corresponding column name in the argument \code{data}, indicating sample size in group 1.}
  \item{data}{an optional data frame containing the dataset of the individual study with a binary outcome. If \code{data} is specified, the previous arguments, \code{e0}, \code{n0}, \code{e1}, and \code{n1}, should be specified as their corresponding column names in \code{data}.}
  \item{all}{a logical value indicating whether all eligible event status modifications will be considered (the default is \code{FALSE}).}
  \item{methods}{a character string or a vector of character strings indicating the method(s) used to calculate the p-value(s) of the association between treatment and outcome. It must consist of \code{"Fisher"} (Fisher's exact test), \code{"chisq"} (chi-squared test), \code{"OR"} (odds ratio), \code{"RR"} (relative risk), and/or \code{"RD"} (risk difference). The default is the vector consisting of all aforementioned five methods.}
  \item{modify0}{a character string indicating the event status modifications in group 0. It should be one of \code{"increase"} (increasing event numbers), \code{"decrease"} (decreasing event numbers), \code{"both"} (the default, modifying event numbers in both directions), and \code{"none"} (no modification).}
  \item{modify1}{a character string indicating the event status modifications in group 1. It has the same usage with \code{modify0}, with the default being \code{"both"}.}
  \item{alpha}{a numeric value between 0 and 1, indicating the statistical significance level (the default is 0.05).}
  \item{alternative}{a character string of either \code{"two.sided"} (the default) or \code{"one.sided"}, indicating the alternative hypothesis. It is only used for \code{"OR"}, \code{"RR"}, and \code{"RD"} in the argument \code{methods}.}
  \item{OR}{a numeric positive value indicating the value of odds ratio under the null hypothesis (the default is 1). It is used only if the argument \code{methods} includes \code{"OR"}.}
  \item{RR}{a numeric positive value indicating the value of relative risk under the null hypothesis (the default is 1). It is used only if the argument \code{methods} includes \code{"RR"}.}
  \item{RD}{a numeric value between \eqn{-}1 indicating the value of risk difference under the null hypothesis (the default is 0). It is used only if the argument \code{methods} includes \code{"RD"}.}
  \item{allcase}{a logical value indicating whether obtaining all cases of minimal event status modifications for altering statistical significance or non-significance (the default is \code{TRUE}).}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  An object of class \code{"frag.study"} (if the argument \code{all} = \code{FALSE}) or of both classes \code{"frag.study"} and \code{"frag.study.all"} (\code{all} = \code{TRUE}). The object is a list containing the following components:
  \item{data}{original data in the form of 2\eqn{\times}2 table.}
  \item{methods}{methods used to calculate the p-value(s) of the association between treatment and outcome.}
  \item{alpha}{pre-specified statistical significance level.}
  \item{alternative}{alternative hypothesis.}
  \item{null}{values of odds ratio, relative risk, and risk difference under the null hypothesis; they are only used if \code{"OR"}, \code{"RR"}, and \code{"RD"} are included in the argument \code{methods}, respectively.}
  \item{modify0}{type of event status modifications in group 0.}
  \item{modify1}{type of event status modifications in group 1.}
  \item{f0.range}{range of event status modifications in group 0.}
  \item{f1.range}{range of event status modifications in group 1.}
  \item{allcase}{a logical value indicating whether all cases of minimal event status modifications for altering statistical significance are obtained; it is only included in the object of class \code{"frag.study"}.}
  \item{tot.mods}{a matrix containing all possible event status modifications with rows representing modifications in group 0 and columns representing those in group 1; it is determined by \code{f0.range} and \code{f1.range}, and is only included in the object of class \code{"frag.study.all"}.}
  \item{pval}{p-value(s) produced by the method(s) included in the argument \code{methods}.}
  \item{FI}{fragility index(es) based on the method(s).}
  \item{FQ}{fragility quotient(s) based on the method(s).}
  \item{dir}{direction of the change of statistical significance.}
  \item{mods}{minimal event status modification(s) that alter statistical significance or non-significance.}
  \item{pvals}{p-values based on all possible event status modifications (as in \code{tot.mods}); it is only included in the object of class \code{"frag.study.all"}.}
  \item{FI0}{fragility index(es) based on the method(s) when event status modifications are restricted to group 0; it is only included in the object of class \code{"frag.study.all"} if both arguments \code{modify0} and \code{modify1} are not \code{"none"}.}
  \item{FQ0}{fragility quotient(s) based on the method(s) when event status modifications are restricted to group 0; it is only included in the object of class \code{"frag.study.all"} if both arguments \code{modify0} and \code{modify1} are not \code{"none"}.}
  \item{dir0}{direction of the change of statistical significance when event status modifications are restricted to group 0; it is only included in the object of class \code{"frag.study.all"} if both arguments \code{modify0} and \code{modify1} are not \code{"none"}.}
  \item{mods0}{minimal event status modification(s) that alter statistical significance or non-significance when event status modifications are restricted to group 0; it is only included in the object of class \code{"frag.study.all"} if both arguments \code{modify0} and \code{modify1} are not \code{"none"}.}
  \item{FI1}{fragility index(es) based on the method(s) when event status modifications are restricted to group 1; it is only included in the object of class \code{"frag.study.all"} if both arguments \code{modify0} and \code{modify1} are not \code{"none"}.}
  \item{FQ1}{fragility quotient(s) based on the method(s) when event status modifications are restricted to group 1; it is only included in the object of class \code{"frag.study.all"} if both arguments \code{modify0} and \code{modify1} are not \code{"none"}.}
  \item{dir1}{direction of the change of statistical significance when event status modifications are restricted to group 1; it is only included in the object of class \code{"frag.study.all"} if both arguments \code{modify0} and \code{modify1} are not \code{"none"}.}
  \item{mods1}{minimal event status modification(s) that alter statistical significance or non-significance when event status modifications are restricted to group 1; it is only included in the object of class \code{"frag.study.all"} if both arguments \code{modify0} and \code{modify1} are not \code{"none"}.}
  Of note, if both arguments \code{modify0} and \code{modify0} are \code{"none"} (i.e., no modification), the items following \code{pval} in the above list are not produced, because fragility cannot be assessed without event status modifications.
}
\references{
Walsh M, Srinathan SK, McAuley DF, Mrkobrada M, Levine O, Ribic C, Molnar AO, Dattani ND, Burke A, Guyatt G, Thabane L, Walter SD, Pogue J, Devereaux PJ (2014). "The statistical significance of randomized controlled trial results is frequently fragile: a case for a Fragility Index." \emph{Journal of Clinical Epidemiology}, \bold{67}(6), 622--8. <\doi{10.1016/j.jclinepi.2013.10.019}>
}
%%\author{
%%  ~~who you are~~
%%}
\note{
  If the significance or non-significance cannot be altered based on the given range of event status modifications, the fragility index are fragility quotient are defined to be \code{NA} (not available).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{frag.study.alpha}} for assessing fragility of an individual study at multiple significance levels; \code{\link{frag.studies}} for assessing fragility of multiple studies; \code{\link{frag.ma}} for assessing fragility of a meta-analysis; and \code{\link{frag.nma}} for assessing fragility of a network meta-analysis.
}
\examples{
## Load datasets of trials on antidepressant drugs
data(dat.ad)

## Assess fragility of trial 13
out1 <- frag.study(e0 = e0, n0 = n0, e1 = e1, n1 = n1, data = dat.ad[13,])
out1

## Consider all possible event status modifications
##  and focus on odds ratio
out2 <- frag.study(e0 = e0, n0 = n0, e1 = e1, n1 = n1, data = dat.ad[13,],
  all = TRUE, methods = "OR")
out2

## Visualize the fragility
\donttest{
plot(out2)
}

## Restrict the modifications to treatment group 0
plot(out2, modify1 = FALSE)

## Change the truncation threshold of p-value to 10^-30
plot(out2, modify1 = FALSE, trun = 30)

## Restrict the modifications to treatment group 1
plot(out2, modify0 = FALSE)

## Assess fragility of trial 13
##  by only increasing event counts in treatment group 0
out3 <- frag.study(e0 = e0, n0 = n0, e1 = e1, n1 = n1, data = dat.ad[13,],
  all = TRUE, modify0 = "increase", modify1 = "none")
out3

## The default plot is for the first method (i.e., Fisher's exact test)
##  if the argument method is not specified
plot(out3)

## Generate plots based on other methods
plot(out3, method = "chisq")
plot(out3, method = "OR")
plot(out3, method = "RR")
plot(out3, method = "RD")

## Assess fragility of trial 13
##  by only modifying event status in treatment group 1
out4 <- frag.study(e0 = e0, n0 = n0, e1 = e1, n1 = n1, data = dat.ad[13,],
  all = TRUE, modify0 = "none", modify1 = "both")
out4
plot(out4)

## Change background colors in the plot
plot(out4, col.sig = c(adjustcolor("blue", alpha.f = 0.1),
  adjustcolor("green", alpha.f = 0.1)))
## Adjust the plot by specifying additional arguments of plot.default()
plot(out4, mgp = c(2, 1, 0), main = "Trial 13", font.main = 1)

## An example of non-significance altered to significance
out5 <- frag.study(e0 = e0, n0 = n0, e1 = e1, n1 = n1, data = dat.ad[21,],
  all = TRUE)
out5
plot(out5, method = "Fisher", legend.pvals = 0.001,
  main = paste("Fragility index =", out5$FI["Fisher"], "in trial 21"))

## Change the statistical significance level to 0.005
out6 <- frag.study(e0 = e0, n0 = n0, e1 = e1, n1 = n1, data = dat.ad[21,],
  all = TRUE, alpha = 0.005)
out6
plot(out6)

## Change the null value of odds ratio to 2
out7 <- frag.study(e0 = e0, n0 = n0, e1 = e1, n1 = n1, data = dat.ad[21,],
  all = TRUE, methods = c("OR"), alpha = 0.05, OR = 2)
out7
plot(out7)
}
\keyword{fragility}
