\name{frag.ma.alpha}
\alias{frag.ma.alpha}
\title{
Assessing Fragility of a Meta-Analysis at Different Significance Levels
}
\description{
Produces fragility index or fragility quotient for altering statistical significance of a meta-analysis with a binary outcome at different significance levels.
}
\usage{
frag.ma.alpha(e0, n0, e1, n1, data, measure = "OR",
              alpha.from = 0.005, alpha.to = 0.05, alpha.breaks = 100,
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
  \item{alpha.from}{a numeric value between 0 and 1, indicating the smallest value of the statistical significance levels to be considered (the default is 0.005). It should be smaller than the argument \code{alpha.to}.}
  \item{alpha.to}{a numeric value between 0 and 1, indicating the largest value of the statistical significance levels to be considered (the default is 0.05). It should be larger than the argument \code{alpha.from}.}
  \item{alpha.breaks}{a positive integer indicating the number of statistical significance levels to be considered (the default is 100). The candidate significance levels are thus equally-spaced between \code{alpha.from} and \code{alpha.to}.}
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
  An object of classes \code{"frag.alpha"} and \code{"frag.ma.alpha"}. The object is a list containing the following components:
  \item{data}{original data in the form of a data frame with six columns: event counts, non-event counts, and sample sizes in group 0 and those in group 1 among all studies in the meta-analysis.}
  \item{measure}{measure of treatment effect (i.e., effect size).}
  \item{alphas}{different statistical significance levels.}
  \item{null}{value of odds ratio, relative risk, or risk difference (specified by \code{measure}) under the null hypothesis. Note that odds ratio and relative risk are presented on a natural logarithmic scale in all output values.}
  \item{est.ori}{the point estimate of the overall effect size based on the original meta-analysis.}
  \item{se.ori}{the standard error of the overall effect size based on the original meta-analysis.}
  \item{test}{the method used to derive confidence intervals by the function \code{\link[metafor]{rma.uni}}. It is one of \code{"z"}, \code{"t"}, \code{"knha"} (Hartung--Knapp--Sidik--Jonkman method), and \code{"adhoc"}; the first indicates Wald-type confidence intervals based on the standard normal distribution, while the latter three yield confidence intervals based on the t-distribution.}
  \item{pval.ori}{the p-value of the overall effect size based on the original meta-analysis.}
  \item{mod.dir}{the direction of the confidence interval change due to event status modifications when the original confidence interval covers the null value (i.e., non-significance altered to significance).}
  \item{FI}{fragility indexes at different statistical significance levels in \code{alphas}.}
  \item{FI.avg}{average fragility index.}
  \item{FQ}{fragility quotients at different statistical significance levels in \code{alphas}.}
  \item{FQ.avg}{average fragility quotient.}
}
\references{
Benjamin DJ, Berger JO, Johannesson M, Nosek BA, Wagenmakers EJ, Berk R, Bollen KA, Brembs B, Brown L, Camerer C, Cesarini D, Chambers CD, Clyde M, Cook TD, De Boeck P, Dienes Z, Dreber A, Easwaran K, Efferson C, Fehr E, Fidler F, Field AP, Forster M, George EI, Gonzalez R, Goodman S, Green E, Green DP, Greenwald AG, Hadfield JD, Hedges LV, Held L, Ho TH, Hoijtink H, Hruschka DJ, Imai K, Imbens G, Ioannidis JPA, Jeon M, Jones JH, Kirchler M, Laibson D, List J, Little R, Lupia A, Machery E, Maxwell SE, McCarthy M, Moore DA, Morgan SL, Munafo M, Nakagawa S, Nyhan B, Parker TH, Pericchi L, Perugini M, Rouder J, Rousseau J, Savalei V, Schonbrodt FD, Sellke T, Sinclair B, Tingley D, Van Zandt T, Vazire S, Watts DJ, Winship C, Wolpert RL, Xie Y, Young C, Zinman J, Johnson VE (2018). "Redefine statistical significance." \emph{Nature Human Behaviour}, \bold{2}, 6--10. <\doi{10.1038/s41562-017-0189-z}>

Ioannidis JPA (2018). "The proposal to lower P value thresholds to .005." \emph{JAMA}, \bold{319}(14), 1429--30. <\doi{10.1001/jama.2018.1536}>

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
\code{\link{frag.ma}} for assessing fragility of a meta-analysis at a specific significance level; and \code{\link{frag.study.alpha}} for assessing fragility of an individual study at multiple significance levels
}
\examples{
## Load datasets of meta-analyses on nutrition support
data(dat.ns)

## Assess fragility of the 11th meta-analysis
##  at multiple significance levels
out1 <- frag.ma.alpha(e0, n0, e1, n1, data = dat.ns[dat.ns$ma.id == 11,],
  alpha.from = 0.01, alpha.to = 0.05, alpha.breaks = 5)
out1

## Generate plot to show fragility measures against significance levels;
##  see more options of usage in examples of frag.study.alpha()
plot(out1)
plot(out1, fragility = "FQ")

\donttest{
## Assess fragility of the 6th meta-analysis
##  at multiple significance levels;
##  it takes longer time
out2 <- frag.ma.alpha(e0, n0, e1, n1, data = dat.ns[dat.ns$ma.id == 6,])
out2
plot(out2)
}
}
\keyword{fragility}