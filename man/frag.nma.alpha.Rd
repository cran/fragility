\name{frag.nma.alpha}
\alias{frag.nma.alpha}
\title{
Assessing Fragility of a Network Meta-Analysis at Different Significance Levels
}
\description{
Produces fragility index or fragility quotient for altering statistical significance of treatment comparison(s) in a network meta-analysis with a binary outcome at different significance levels.
}
\usage{
frag.nma.alpha(sid, tid, e, n, data, measure = "OR", random = TRUE,
               alpha.from = 0.005, alpha.to = 0.05, alpha.breaks = 10,
               mod.dir = "both", tid1.f, tid2.f, OR = 1, RR = 1, RD = 0,
               incr, allincr, addincr, allstudies, ...)
}
\arguments{
  \item{sid}{a numeric vector or the corresponding column name in the argument \code{data}, indicating the study IDs in the network meta-analysis.}
  \item{tid}{a numeric vector or the corresponding column name in the argument \code{data}, indicating the treatment IDs in the network meta-analysis.}
  \item{e}{a numeric vector or the corresponding column name in the argument \code{data}, indicating the event count in each study's each treatment group in the network meta-analysis.}
  \item{n}{a numeric vector or the corresponding column name in the argument \code{data}, indicating the sample size in each study's each treatment group in the network meta-analysis.}
  \item{data}{an optional data frame containing the dataset of the collected studies on multiple treatment comparisons in the network meta-analysis with a binary outcome. If \code{data} is specified, the previous arguments, \code{sid}, \code{tid}, \code{e}, and \code{n}, should be specified as their corresponding column names in \code{data}.}
  \item{measure}{a character string indicating the measure of treatment effect (i.e., effect size) for the binary outcome. It should be one of \code{"OR"} (odds ratio, the default), \code{"RR"} (relative risk), and \code{"RD"} (risk difference).}
  \item{random}{a logical value indicating whether the network meta-analysis is performed under the fixed-effects setting (\code{FALSE}) or the random-effects setting (\code{TRUE}, the default).}
  \item{alpha.from}{a numeric value between 0 and 1, indicating the smallest value of the statistical significance levels to be considered (the default is 0.005). It should be smaller than the argument \code{alpha.to}.}
  \item{alpha.to}{a numeric value between 0 and 1, indicating the largest value of the statistical significance levels to be considered (the default is 0.05). It should be larger than the argument \code{alpha.from}.}
  \item{alpha.breaks}{a positive integer indicating the number of statistical significance levels to be considered (the default is 100). The candidate significance levels are thus equally-spaced between \code{alpha.from} and \code{alpha.to}.}
  \item{mod.dir}{a character string or a square matrix of character strings indicating the direction of the confidence interval change for each treatment comparison due to event status modifications when the original confidence interval covers the null value (i.e., non-significance altered to significance). It is not used when significance is altered to non-significance. It should be one of \code{"left"} (moving to the left side of the null value), \code{"right"} (moving to the right side of the null value), \code{"one"} (based on the direction of the point estimate of the original overall effect size), and \code{"both"} (both directions, the default). When it is specified as a matrix, the dimension should be the number of treatments, and the element (h, k) specifies the direction for treatment h compared with treatment k.}
  \item{tid1.f}{the first treatment(s) of comparison(s) of interest; if not specified, the fragility of all comparisons will be assessed.}
  \item{tid2.f}{the second treatment(s) of comparison(s) of interest.}
  \item{OR}{a numeric positive value indicating the value of odds ratio under the null hypothesis (the default is 1). It is used only if the argument \code{measure} is \code{"OR"}.}
  \item{RR}{a numeric positive value indicating the value of relative risk under the null hypothesis (the default is 1). It is used only if the argument \code{measure} is \code{"RR"}.}
  \item{RD}{a numeric value between \eqn{-}1 indicating the value of risk difference under the null hypothesis (the default is 0). It is used only if the argument \code{measure} is \code{"RD"}.}
  \item{incr}{A numerical value which is added to each cell frequency for studies with a zero cell count. It is passed to the function \code{\link[netmeta]{pairwise}} in the package \CRANpkg{netmeta}.}
  \item{allincr}{A logical value indicating whether \code{incr} is added to each cell frequency of all studies if at least one study has a zero cell count. If \code{FALSE} (the default), \code{incr} is added only to each cell frequency of studies with a zero cell count. It is passed to the function \code{\link[netmeta]{pairwise}} in the package \CRANpkg{netmeta}.}
  \item{addincr}{A logical value indicating whether \code{incr} is added to each cell frequency of all studies irrespective of zero cell counts (the default is \code{FALSE}). It is passed to the function \code{\link[netmeta]{pairwise}} in the package \CRANpkg{netmeta}.}
  \item{allstudies}{A logical value indicating whether studies with zero or all events in two treatment groups are to be included in the network meta-analysis (the default is \code{FALSE}). It is used only for \code{measure} = \code{"OR"} and \code{"RR"}. It is passed to the function \code{\link[netmeta]{pairwise}} in the package \CRANpkg{netmeta}.}
  \item{\dots}{other arguments that can be passed to \code{\link[netmeta]{netmeta}}. See details in the manual of the package \CRANpkg{netmeta}.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  An object of classes \code{"frag.alpha"} and \code{"frag.nma.alpha"}. The object is a list containing the following components:
  \item{data}{original data of the network meta-analysis.}
  \item{measure}{measure of treatment effect (i.e., effect size).}
  \item{alphas}{different statistical significance levels.}
  \item{null}{value of odds ratio, relative risk, or risk difference (specified by \code{measure}) under the null hypothesis. Note that odds ratio and relative risk are presented on a natural logarithmic scale in all output values.}
  \item{est.ori}{a matrix containing the point estimate of the overall effect size for each pair of treatment comparison based on the original network meta-analysis. Each element presents the treatment corresponding to the row compared with the treatment corresponding to the column; the similar interpretation applies to other outputs in the form of matrix.}
  \item{se.ori}{a matrix containing the standard error of the overall effect size for each pair of treatment comparison based on the original network meta-analysis.}
  \item{pval.ori}{the p-value of the overall effect size for each pair of treatment comparison based on the original network meta-analysis.}
  \item{mod.dir}{the direction of the confidence interval change due to event status modifications when the original confidence interval covers the null value (i.e., non-significance altered to significance).}
  \item{tid.f}{treatment comparisons whose fragility is assessed.}
  \item{FI}{fragility indexes at different statistical significance levels in \code{alphas}.}
  \item{FI.avg}{average fragility indexes.}
  \item{FQ}{fragility quotients at different statistical significance levels in \code{alphas}, calculated as fragility indexes divided by the total sample sizes associated with the specific treatment comparisons.}
  \item{FQ.avg}{average fragility quotients based on the total sample sizes associated with the specific treatment comparisons.}
  \item{FQ.nma}{fragility quotients at different statistical significance levels in \code{alphas}, calculated as fragility indexes divided by the total sample size in the whole network meta-analysis (across all treatment groups).}
  \item{FQ.nma.avg}{average fragility quotients based on the total sample size in the whole network meta-analysis (across all treatment groups).}
}
\references{
Benjamin DJ, Berger JO, Johannesson M, Nosek BA, Wagenmakers EJ, Berk R, Bollen KA, Brembs B, Brown L, Camerer C, Cesarini D, Chambers CD, Clyde M, Cook TD, De Boeck P, Dienes Z, Dreber A, Easwaran K, Efferson C, Fehr E, Fidler F, Field AP, Forster M, George EI, Gonzalez R, Goodman S, Green E, Green DP, Greenwald AG, Hadfield JD, Hedges LV, Held L, Ho TH, Hoijtink H, Hruschka DJ, Imai K, Imbens G, Ioannidis JPA, Jeon M, Jones JH, Kirchler M, Laibson D, List J, Little R, Lupia A, Machery E, Maxwell SE, McCarthy M, Moore DA, Morgan SL, Munafo M, Nakagawa S, Nyhan B, Parker TH, Pericchi L, Perugini M, Rouder J, Rousseau J, Savalei V, Schonbrodt FD, Sellke T, Sinclair B, Tingley D, Van Zandt T, Vazire S, Watts DJ, Winship C, Wolpert RL, Xie Y, Young C, Zinman J, Johnson VE (2018). "Redefine statistical significance." \emph{Nature Human Behaviour}, \bold{2}, 6--10. <\doi{10.1038/s41562-017-0189-z}>

Ioannidis JPA (2018). "The proposal to lower P value thresholds to .005." \emph{JAMA}, \bold{319}(14), 1429--30. <\doi{10.1001/jama.2018.1536}>
}
%%\author{
%%  ~~who you are~~
%%}
\note{
When the network meta-analysis contains many studies and many treatments, the assess of fragility may be very computationally demanding. In such cases, one may focus on certain treatment comparisons that are of primary interest by specifying the arguments \code{tid1.f} and \code{tid2.f}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{frag.nma}}
}
\examples{
\donttest{
## Load datasets of network meta-analysis of
##  chronic obstructive pulmonary disease (COPD)
data(dat.copd)

## Assess fragility of network  meta-analysis of COPD
##  at multiple significance levels
out1 <- frag.nma.alpha(sid, tid, e, n, data = dat.copd)
out1

## Generate plot to show fragility measures against significance levels;
##  see more options of usage in examples of frag.study.alpha()
plot(out1, tid1 = 2, tid2 = 1)
plot(out1, tid1 = 3, tid2 = 1)
plot(out1, tid1 = 4, tid2 = 1)
plot(out1, tid1 = 3, tid2 = 2)
plot(out1, tid1 = 4, tid2 = 2)
plot(out1, tid1 = 4, tid2 = 3)
}
}
\keyword{fragility}