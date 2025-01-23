\name{frag.ma}
\alias{frag.ma}
\title{
Assessing the Fragility of a Meta-Analysis With a Binary Outcome
}
\description{
Produces results, such as fragility index, fragility quotient, and event status modifications, for altering the statistical significance of a meta-analysis with a binary outcome.
}
\usage{
frag.ma(e0, n0, e1, n1, data, measure = "OR", alpha = 0.05,
        mod.dir = "both", OR = 1, RR = 1, RD = 0,
        method = "DL", test = "z", drop00 = FALSE, ...)
}
\arguments{
  \item{e0}{a numeric vector or the corresponding column name in the argument \code{data}, indicating event counts in group 0 among the studies in the meta-analysis.}
  \item{n0}{a numeric vector or the corresponding column name in the argument \code{data}, indicating sample sizes in group 0 among the studies in the meta-analysis.}
  \item{e1}{a numeric vector or the corresponding column name in the argument \code{data}, indicating event counts in group 1 among the studies in the meta-analysis.}
  \item{n1}{a numeric vector or the corresponding column name in the argument \code{data}, indicating sample sizes in group 1 among the studies in the meta-analysis.}
  \item{data}{an optional data frame containing the dataset of the collected studies in the meta-analysis with a binary outcome. If \code{data} is specified, the previous arguments, \code{e0}, \code{n0}, \code{e1}, and \code{n1}, should be specified as their corresponding column names in \code{data}.}
  \item{measure}{a character string indicating the measure of treatment effect (i.e., effect size) for the binary outcome. It should be one of \code{"OR"} (odds ratio, the default), \code{"RR"} (relative risk), and \code{"RD"} (risk difference).}
  \item{alpha}{a numeric value between 0 and 1, indicating the statistical significance level (the default is 0.05). It determines the confidence level \eqn{(1-\alpha)\times100\%} of the confidence interval when deriving the fragility index/quotient.}
  \item{mod.dir}{a character string indicating the direction of the confidence interval change due to event status modifications when the original confidence interval covers the null value (i.e., non-significance altered to significance). It is not used when significance is altered to non-significance. It should be one of \code{"left"} (moving to the left side of the null value), \code{"right"} (moving to the right side of the null value), \code{"one"} (based on the direction of the point estimate of the original overall effect size), and \code{"both"} (both directions, the default).}
  \item{OR}{a numeric positive value indicating the value of odds ratio under the null hypothesis (the default is 1). It is used only if the argument \code{measure} is \code{"OR"}.}
  \item{RR}{a numeric positive value indicating the value of relative risk under the null hypothesis (the default is 1). It is used only if the argument \code{measure} is \code{"RR"}.}
  \item{RD}{a numeric value between \eqn{-}1 indicating the value of risk difference under the null hypothesis (the default is 0). It is used only if the argument \code{measure} is \code{"RD"}.}
  \item{method}{a character string specifying the method for performing meta-analysis. It is passed to \code{\link[metafor]{rma.uni}}. Options include the fixed-effect model by specifying \code{"FE"} and the random-effects model by specifying \code{"DL"} (DerSimonian--Laird estimator), \code{"HE"} (Hedges estimator), \code{"HS"} (Hunter--Schmidt estimator), \code{"SJ"} (Sidik--Jonkman estimator), \code{"ML"} (maximum likelihood estimator), \code{"REML"} (restricted maximum likelihood estimator), \code{"EB"} (empirical Bayes estimator), \code{"PM"} (Paule--Mandel estimator), or \code{"GENQ"} (generalized \eqn{W}-statistic estimator) for different methods to estimate the heterogeneity variance. See details in the manual of the package \CRANpkg{metafor}.}
  \item{test}{a character string specifying how confidence intervals are derived. It is passed to \code{\link[metafor]{rma.uni}}. Options include \code{"z"}, \code{"t"}, \code{"knha"} (Hartung--Knapp--Sidik--Jonkman method), and \code{"adhoc"}; the first indicates Wald-type confidence intervals based on the standard normal distribution (the default), while the latter three yield confidence intervals based on the t-distribution. See details in the manual of the package \CRANpkg{metafor}.}
  \item{drop00}{a logical value indicating whether studies with no events (or only events) in both groups should be dropped when calculating the observed effect sizes or outcomes. It is passed to \code{\link[metafor]{rma.uni}}.}
  \item{\dots}{other arguments that can be passed to \code{\link[metafor]{rma.uni}}. See details in the manual of the package \CRANpkg{metafor}.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  An object of class \code{"frag.ma"}. The object is a list containing the following components:
  \item{data}{original data in the form of a data frame with six columns: event counts, non-event counts, and sample sizes in group 0 and those in group 1 among all studies in the meta-analysis.}
  \item{measure}{measure of treatment effect (i.e., effect size).}
  \item{alpha}{pre-specified statistical significance level.}
  \item{null}{value of odds ratio, relative risk, or risk difference (specified by \code{measure}) under the null hypothesis. Note that odds ratio and relative risk are presented on a natural logarithmic scale in all output values.}
  \item{est.ori}{the point estimate of the overall effect size based on the original meta-analysis.}
  \item{ci.ori}{the confidence interval of the overall effect size based on the original meta-analysis.}
  \item{pval.ori}{the p-value of the overall effect size based on the original meta-analysis.}
  \item{mod.dir}{the direction of the confidence interval change due to event status modifications when the original confidence interval covers the null value (i.e., non-significance altered to significance).}
  \item{FI}{fragility index.}
  \item{FQ}{fragility quotient.}
  \item{dir}{direction of the change of statistical significance.}
  \item{sid.iter}{study index(es) with event status modifications during the iterative algorithm.}
  \item{g0.mod.iter}{event status modification(s) in group 0 during the iterative algorithm.}
  \item{g1.mod.iter}{event status modification(s) in group 1 during the iterative algorithm.}
  \item{est.iter}{point estimate(s) of the overall effect size during the iterative algorithm.}
  \item{ci.iter}{confidence interval(s) of the overall effect size during the iterative algorithm.}
  \item{data.mod}{data with event status modification(s) that alter the significance or non-significance.}
}
\references{
Atal I, Porcher R, Boutron I, Ravaud P (2019). "The statistical significance of meta-analyses is frequently fragile: definition of a fragility index for meta-analyses." \emph{Journal of Clinical Epidemiology}, \bold{111}, 32--40. <\doi{10.1016/j.jclinepi.2019.03.012}>

Lin L, Chu H (2022). "Assessing and visualizing fragility of clinical results with binary outcomes in R using the fragility package." \emph{PLOS ONE}, \bold{17}(6), e0268754. <\doi{10.1371/journal.pone.0268754}>

Lin L, Xing A, Chu H, Murad MH, Xu C, Baer BR, Wells MT, Sanchez-Ramos L (2023). "Assessing the robustness of results from clinical trials and meta-analyses with the fragility index." \emph{American Journal of Obstetrics and Gynecology}, \bold{228}(3), 276--282. <\doi{10.1016/j.ajog.2022.08.053}>

Wang Z, Xing X, Mun EY, Wu C, Lin L (2025). "The role of double-zero-event studies in evidence synthesis: Evaluating robustness using the fragility index." \emph{Journal of Evaluation in Clinical Practice}, \bold{31}(1), e14301. <\doi{10.1111/jep.14301}>
}
%%\author{
%%  ~~who you are~~
%%}
\note{
  If the significance or non-significance cannot be altered based on the given range of event status modifications, the fragility index and fragility quotient are defined to be \code{NA} (not available).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{frag.ma.alpha}} for assessing fragility of a meta-analysis at multiple significance levels; \code{\link{frag.mas}} for assessing fragility of multiple meta-analyses; \code{\link{frag.study}} for assessing fragility of an individual study; and \code{\link{frag.nma}} for assessing fragility of a network meta-analysis.
}
\examples{
## Load datasets of meta-analyses on nutrition support
data(dat.ns)

## Assess fragility of the 11th meta-analysis
out1 <- frag.ma(e0, n0, e1, n1, data = dat.ns[dat.ns$ma.id == 11,])
out1

## Visualize the process of event status modifications
##  for altering the non-significance to significance
plot(out1)

## Assess fragility based on relative risk
out2 <- frag.ma(e0, n0, e1, n1, data = dat.ns[dat.ns$ma.id == 11,],
  measure = "RR")
out2

## Assess fragility based on risk difference
out3 <- frag.ma(e0, n0, e1, n1, data = dat.ns[dat.ns$ma.id == 11,],
  measure = "RD")
out3

## Use the REML method to estimate the heterogeneity variance
out4 <- frag.ma(e0, n0, e1, n1, data = dat.ns[dat.ns$ma.id == 11,],
  method = "REML")
out4

## Use the Hartung-Knapp-Sidik-Jonkman to derive
##  the confidence interval in meta-analysis
out5 <- frag.ma(e0, n0, e1, n1, data = dat.ns[dat.ns$ma.id == 11,],
  method = "REML", test = "knha")
out5

## Change the significance level to 0.005
out6 <- frag.ma(e0, n0, e1, n1, data = dat.ns[dat.ns$ma.id == 11,],
  alpha = 0.005)
out6
plot(out6)

## Reduce the distance between the two lines by breaking the y-axis
plot(out6, ybreaks = c(5, 9))

## Suppress the markers of study IDs around the lines
plot(out6, ybreaks = c(5, 9), study.marker = FALSE)

## Change the size of markers of study IDs
plot(out6, ybreaks = c(5, 9), cex.marker = 0.5)

## Change the distance between the markers of study IDs and the lines
plot(out6, ybreaks = c(5, 9), offset.marker = 0)

## Change the colors and width of the lines
plot(out6, ybreaks = c(5, 9), col.line = c("gray", "orange"), lwd = 2)

## Change the location of the legend
plot(out6, ybreaks = c(5, 9), x.legend = "topleft")

\donttest{
## Assess fragility of the first meta-analysis;
##  it takes longer time
out7 <- frag.ma(e0, n0, e1, n1, data = dat.ns[dat.ns$ma.id == 1,])
out7
plot(out7)
plot(out7, ybreaks = c(840, 880))

## Achieve significance by moving the CI based on the direction of the 
##  original overall effect size estimate
out8 <- frag.ma(e0, n0, e1, n1, data = dat.ns[dat.ns$ma.id == 1,],
  mod.dir = "one")
out8

## Achieve significance by moving the CI to the left side
out9 <- frag.ma(e0, n0, e1, n1, data = dat.ns[dat.ns$ma.id == 1,],
  mod.dir = "left")
out9

## Achieve significance by moving the CI to the right side;
##  it takes longer time than that moving the CI to the left side
##  as the original overall effect size estimate was less than the null
out10 <- frag.ma(e0, n0, e1, n1, data = dat.ns[dat.ns$ma.id == 1,],
  mod.dir = "right")
out10
dev.new(width = 12, height = 6)
plot(out10, cex.marker = 0.5)
}
}
\keyword{fragility}
