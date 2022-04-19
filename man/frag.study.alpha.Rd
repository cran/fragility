\name{frag.study.alpha}
\alias{frag.study.alpha}
\title{
Assessing Fragility of an Individual Study at Different Significance Levels
}
\description{
Produces fragility index or fragility quotient for altering statistical significance of an individual study with a binary outcome at different significance levels.
}
\usage{
frag.study.alpha(e0, n0, e1, n1, data, methods,
                 modify0 = "both", modify1 = "both",
                 alpha.from = 0.005, alpha.to = 0.05, alpha.breaks = 100,
                 alternative = "two.sided", OR = 1, RR = 1, RD = 0)
}
\arguments{
  \item{e0}{a numeric value or the corresponding column name in the argument \code{data}, indicating event count in group 0.}
  \item{n0}{a numeric value or the corresponding column name in the argument \code{data}, indicating sample size in group 0.}
  \item{e1}{a numeric value or the corresponding column name in the argument \code{data}, indicating event count in group 1.}
  \item{n1}{a numeric value or the corresponding column name in the argument \code{data}, indicating sample size in group 1.}
  \item{data}{an optional data frame containing the dataset of the individual study with a binary outcome. If \code{data} is specified, the previous arguments, \code{e0}, \code{n0}, \code{e1}, and \code{n1}, should be specified as their corresponding column names in \code{data}.}
  \item{methods}{a character string or a vector of character strings indicating the method(s) used to calculate the p-value(s) of the association between treatment and outcome. It must consist of \code{"Fisher"} (Fisher's exact test), \code{"chisq"} (chi-squared test), \code{"OR"} (odds ratio), \code{"RR"} (relative risk), and/or \code{"RD"} (risk difference). The default is the vector consisting of all aforementioned five methods.}
  \item{modify0}{a character string indicating the event status modifications in group 0. It should be one of \code{"increase"} (increasing event numbers), \code{"decrease"} (decreasing event numbers), \code{"both"} (the default, modifying event numbers in both directions), and \code{"none"} (no modification).}
  \item{modify1}{a character string indicating the event status modifications in group 1. It has the same usage with \code{modify0}, with the default being \code{"both"}.}
  \item{alpha.from}{a numeric value between 0 and 1, indicating the smallest value of the statistical significance levels to be considered (the default is 0.005). It should be smaller than the argument \code{alpha.to}.}
  \item{alpha.to}{a numeric value between 0 and 1, indicating the largest value of the statistical significance levels to be considered (the default is 0.05). It should be larger than the argument \code{alpha.from}.}
  \item{alpha.breaks}{a positive integer indicating the number of statistical significance levels to be considered (the default is 100). The candidate significance levels are thus equally-spaced between \code{alpha.from} and \code{alpha.to}.}
  \item{alternative}{a character string of either \code{"two.sided"} (the default) or \code{"one.sided"}, indicating the alternative hypothesis. It is only used for \code{"OR"}, \code{"RR"}, and \code{"RD"} in the argument \code{methods}.}
  \item{OR}{a numeric positive value indicating the value of odds ratio under the null hypothesis (the default is 1). It is used only if the argument \code{methods} includes \code{"OR"}.}
  \item{RR}{a numeric positive value indicating the value of relative risk under the null hypothesis (the default is 1). It is used only if the argument \code{methods} includes \code{"RR"}.}
  \item{RD}{a numeric value between \eqn{-}1 indicating the value of risk difference under the null hypothesis (the default is 0). It is used only if the argument \code{methods} includes \code{"RD"}.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  An object of classes \code{"frag.alpha"} and \code{"frag.study.alpha"}. The object is a list containing the following components:
  \item{data}{original data in the form of 2\eqn{\times}2 table.}
  \item{methods}{methods used to calculate the p-value(s) of the association between treatment and outcome.}
  \item{alphas}{different statistical significance levels.}
  \item{alternative}{alternative hypothesis.}
  \item{null}{values of odds ratio, relative risk, and risk difference under the null hypothesis; they are only used if \code{"OR"}, \code{"RR"}, and \code{"RD"} are included in the argument \code{methods}, respectively.}
  \item{modify0}{type of event status modifications in group 0.}
  \item{modify1}{type of event status modifications in group 1.}
  \item{f0.range}{range of event status modifications in group 0.}
  \item{f1.range}{range of event status modifications in group 0.}
  \item{pval}{p-value(s) produced by the method(s) included in the argument \code{methods}.}
  \item{FI}{fragility indexes based on the method(s) at each statistical significance level in \code{alphas}.}
  \item{FI.avg}{average fragility index based on each method.}
  \item{FQ}{fragility quotients based on the method(s) at each statistical significance level in \code{alphas}.}
  \item{FQ.avg}{average fragility quotient based on each method.}
  Of note, if both arguments \code{modify0} and \code{modify0} are \code{"none"} (i.e., no modification), the items following \code{pval} in the above list are not produced, because fragility cannot be assessed without event status modifications.
}
\references{
Benjamin DJ, Berger JO, Johannesson M, Nosek BA, Wagenmakers EJ, Berk R, Bollen KA, Brembs B, Brown L, Camerer C, Cesarini D, Chambers CD, Clyde M, Cook TD, De Boeck P, Dienes Z, Dreber A, Easwaran K, Efferson C, Fehr E, Fidler F, Field AP, Forster M, George EI, Gonzalez R, Goodman S, Green E, Green DP, Greenwald AG, Hadfield JD, Hedges LV, Held L, Ho TH, Hoijtink H, Hruschka DJ, Imai K, Imbens G, Ioannidis JPA, Jeon M, Jones JH, Kirchler M, Laibson D, List J, Little R, Lupia A, Machery E, Maxwell SE, McCarthy M, Moore DA, Morgan SL, Munafo M, Nakagawa S, Nyhan B, Parker TH, Pericchi L, Perugini M, Rouder J, Rousseau J, Savalei V, Schonbrodt FD, Sellke T, Sinclair B, Tingley D, Van Zandt T, Vazire S, Watts DJ, Winship C, Wolpert RL, Xie Y, Young C, Zinman J, Johnson VE (2018). "Redefine statistical significance." \emph{Nature Human Behaviour}, \bold{2}, 6--10. <\doi{10.1038/s41562-017-0189-z}>

Ioannidis JPA (2018). "The proposal to lower P value thresholds to .005." \emph{JAMA}, \bold{319}(14), 1429--30. <\doi{10.1001/jama.2018.1536}>
}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{frag.study}} for assessing fragility of an individual study at a specific significance level.
}
\examples{
## Load datasets of trials on antidepressant drugs
data(dat.ad)

## Assess fragility of trial 13 at multiple significance levels
out1 <- frag.study.alpha(e0, n0, e1, n1, data = dat.ad[13,])
out1

## Generate plots to show the change of fragility measure,
##  where the default plot is for the first method
##  (i.e., Fisher's exact test) if the argument method is not specified
plot(out1)

## The plot for odds ratio
plot(out1, method = "OR")

## Change the choices of significance levels
out2 <- frag.study.alpha(e0, n0, e1, n1, data = dat.ad[13,],
  methods = "OR", alpha.from = 0.001, alpha.to = 0.1, alpha.breaks = 500)
out2
plot(out2)

## Generate the plot with the x-axis on the log scale
plot(out2, log = "x")

## Generate the plot for fragility quotient (in percentage)
plot(out2, fragility = "FQ", log = "x")

## Generate the plot for fragility quotient (in decimal form)
plot(out2, fragility = "FQ", percentage = FALSE, log = "x")

## Change the size of points
plot(out2, cex.pts = 1, log = "x")

## Change the colors
plot(out2, col.line = "blue", col.pval = "orange",
  col.sig = c("gray", "black"), log = "x")

## Change the line type
plot(out2, lty.pval = 3, log = "x")

## Change the line width
plot(out2, lwd = 2, lwd.pval = 2, log = "x")

## Change the point symbol
plot(out2, pch = 1, log = "x")
}
\keyword{fragility}
