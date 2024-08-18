\name{frag.mas}
\alias{frag.mas}
\title{
Assessing Fragility of Multiple Meta-Analyses With Binary Outcomes
}
\description{
Produces fragility indexes or fragility quotients for altering statistical significance of multiple meta-analyses with binary outcomes.
}
\usage{
frag.mas(e0, n0, e1, n1, ma.id, data, measure = "OR", alpha = 0.05,
         mod.dir = "both", OR = 1, RR = 1, RD = 0,
         method = "DL", test = "z", drop00 = FALSE, ...)
}
\arguments{
  \item{e0}{a numeric vector or the corresponding column name in the argument \code{data}, indicating event counts in group 0 among the studies in all meta-analyses.}
  \item{n0}{a numeric vector or the corresponding column name in the argument \code{data}, indicating sample sizes in group 0 among the studies in all meta-analyses.}
  \item{e1}{a numeric vector or the corresponding column name in the argument \code{data}, indicating event counts in group 1 among the studies in all meta-analyses.}
  \item{n1}{a numeric vector or the corresponding column name in the argument \code{data}, indicating sample sizes in group 1 among the studies in all meta-analyses.}
  \item{ma.id}{a numeric vector or the corresponding column name in the argument \code{data}, indicating the index of each meta-analysis.}
  \item{data}{an optional data frame containing the dataset of the multiple meta-analyses with binary outcomes. If \code{data} is specified, the previous arguments, \code{e0}, \code{n0}, \code{e1}, \code{n1}, and \code{ma.id}, should be specified as their corresponding column names in \code{data}.}
  \item{measure}{a character string indicating the measure of treatment effect (i.e., effect size) for the binary outcome. It should be one of \code{"OR"} (odds ratio, the default), \code{"RR"} (relative risk), and \code{"RD"} (risk difference).}
  \item{alpha}{a numeric value between 0 and 1, indicating the statistical significance level (the default is 0.05). It determines the confidence level \eqn{(1-\alpha)\times100\%} of the confidence interval when deriving the fragility index/quotient.}
  \item{mod.dir}{a character string indicating the direction of the confidence interval change due to event status modifications when the original confidence interval covers the null value (i.e., non-significance altered to significance). It is not used when significance is altered to non-significance. It should be one of \code{"left"} (moving to the left side of the null value), \code{"right"} (moving to the right side of the null value), \code{"one"} (based on the direction of the point estimate of the original overall effect size), and \code{"both"} (both directions, the default).}
  \item{OR}{a numeric positive value indicating the value of odds ratio under the null hypothesis (the default is 1). It is used only if the argument \code{measure} is \code{"OR"}.}
  \item{RR}{a numeric positive value indicating the value of relative risk under the null hypothesis (the default is 1). It is used only if the argument \code{measure} is \code{"RR"}.}
  \item{RD}{a numeric value between \eqn{-}1 indicating the value of risk difference under the null hypothesis (the default is 0). It is used only if the argument \code{measure} is \code{"RD"}.}
  \item{method}{a character string specifying the method for performing meta-analysis. It is the same with the argument in \code{\link{frag.ma}}, and is passed to \code{\link[metafor]{rma.uni}}. See details in the manual of the package \CRANpkg{metafor}.}
  \item{test}{a character string specifying how confidence intervals are derived. It is the same with the argument in \code{\link{frag.ma}}, and is passed to \code{\link[metafor]{rma.uni}}. See details in the manual of the package \CRANpkg{metafor}.}
  \item{drop00}{a logical value indicating whether studies with no events (or only events) in both groups should be dropped when calculating the observed effect sizes or outcomes. It is passed to \code{\link[metafor]{rma.uni}}.}
  \item{\dots}{other arguments that can be passed to \code{\link[metafor]{rma.uni}}. See details in the manual of the package \CRANpkg{metafor}.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  An object of classes \code{"frag.mas"} and \code{"frag.multi"}. The object is a list containing the following components:
  \item{measure}{measure of treatment effect (i.e., effect size).}
  \item{alpha}{pre-specified statistical significance level.}
  \item{null}{value of odds ratio, relative risk, or risk difference (specified by \code{measure}) under the null hypothesis. Note that odds ratio and relative risk are presented on a natural logarithmic scale in all output values.}
  \item{mod.dir}{the direction of the confidence interval change due to event status modifications when the original confidence interval covers the null value (i.e., non-significance altered to significance).}
  \item{est.ori}{the point estimates of the overall effect sizes based on the original datasets of the multiple meta-analyses.}
  \item{ci.ori}{the confidence intervals of the overall effect sizes based on the original datasets of the multiple meta-analyses.}
  \item{pval.ori}{the p-values of the overall effect sizes based on the original datasets of the multiple meta-analyses.}
  \item{FI}{fragility indexes of the multiple meta-analyses.}
  \item{FQ}{fragility quotients of the multiple meta-analyses.}
}
\references{
Atal I, Porcher R, Boutron I, Ravaud P (2019). "The statistical significance of meta-analyses is frequently fragile: definition of a fragility index for meta-analyses." \emph{Journal of Clinical Epidemiology}, \bold{111}, 32--40. <\doi{10.1016/j.jclinepi.2019.03.012}>

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
\code{\link{frag.ma}} for assessing fragility of a meta-analysis; and \code{\link{frag.studies}} for assessing fragility of multiple individual studies.
}
\examples{
## Load datasets of meta-analyses on nutrition support
data(dat.ns)

## Assess fragility of the 11th, 12th, and 13th meta-analyses
out1 <- frag.mas(e0, n0, e1, n1, ma.id,
  data = dat.ns[dat.ns$ma.id >= 11 & dat.ns$ma.id <= 13,], mod.dir = "one")
out1

## Generate bar plot or histogram to show fragility measures' distribution;
##  see more options of usage in examples of frag.studies()
plot(out1)

\donttest{
## Assess fragility of all meta-analyses
##  in the Cochrane systematic review on nutrition support;
##  it takes longer time
out2 <- frag.mas(e0, n0, e1, n1, ma.id, data = dat.ns)
out2
plot(out2)
plot(out2, max.f = 40, cex.name = 0.6, mgp = c(2, 0.5, 0))
plot(out2, fragility = "FQ", breaks = 20)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{fragility}
