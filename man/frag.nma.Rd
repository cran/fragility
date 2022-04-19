\name{frag.nma}
\alias{frag.nma}
\title{
Assessing Fragility of a Network Meta-Analysis With a Binary Outcome
}
\description{
Produces results, such as fragility index, fragility quotient, and event status modifications for altering statistical significance, of treatment comparison(s) in a network meta-analysis with a binary outcome.
}
\usage{
frag.nma(sid, tid, e, n, data, measure = "OR", random = TRUE, alpha = 0.05,
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
  \item{alpha}{a numeric value between 0 and 1, indicating the statistical significance level (the default is 0.05). It determines the confidence level \eqn{(1-\alpha)\times100\%} of the confidence interval when deriving the fragility index/quotient.}
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
  An object of class \code{"frag.nma"}. The object is a list containing the following components:
  \item{data}{original data of the network meta-analysis.}
  \item{measure}{measure of treatment effect (i.e., effect size).}
  \item{alpha}{pre-specified statistical significance level.}
  \item{null}{value of odds ratio, relative risk, or risk difference (specified by \code{measure}) under the null hypothesis. Note that odds ratio and relative risk are presented on a natural logarithmic scale in all output values.}
  \item{est.ori}{a matrix containing the point estimate of the overall effect size for each pair of treatment comparisons based on the original network meta-analysis. Each element presents the treatment corresponding to the row compared with the treatment corresponding to the column; a similar interpretation applies to other outputs in the form of a matrix.}
  \item{ci.lb.ori}{the lower bound of the confidence interval of the overall effect size for each pair of treatment comparisons based on the original network meta-analysis.}
  \item{ci.lb.ori}{the upper bound of the confidence interval of the overall effect size for each pair of treatment comparisons based on the original network meta-analysis.}
  \item{pval.ori}{the p-value of the overall effect size for each pair of treatment comparisons based on the original network meta-analysis.}
  \item{mod.dir}{the direction of the confidence interval change due to event status modifications when the original confidence interval covers the null value (i.e., non-significance altered to significance).}
  \item{tid.f}{treatment comparisons whose fragility is assessed.}
  \item{FI}{fragility index.}
  \item{FQ}{fragility quotient, calculated as fragility index divided by the total sample size associated with the specific treatment comparison.}
  \item{FQ.nma}{fragility quotient, calculated as fragility index divided by the total sample size in the whole network meta-analysis (across all treatment groups).}
  \item{dir}{direction of the change of statistical significance.}
  \item{sid.iter}{a list of study index(es) with event status modifications during the iterative algorithm for specific pairs of treatment comparisons.}
  \item{tid1.mod.iter}{a list of event status modification(s) in the first treatment group during the iterative algorithm for specific pairs of treatment comparisons.}
  \item{tid2.mod.iter}{a list of event status modification(s) in the second treatment group during the iterative algorithm for specific pairs of treatment comparisons.}
  \item{est.iter}{a list of point estimate(s) of the overall effect size during the iterative algorithm for specific pairs of treatment comparisons.}
  \item{ci.lb.iter}{a list of lower bound(s) of confidence interval(s) of the overall effect size during the iterative algorithm for specific pairs of treatment comparisons.}
  \item{ci.ub.iter}{a list of upper bound(s) of confidence interval(s) of the overall effect size during the iterative algorithm for specific pairs of treatment comparisons.}
  \item{data.mod}{a list of data with event status modification(s) that alter the significance or non-significance for specific pairs of treatment comparisons.}
}
\references{
Schwarzer G, Carpenter JR, Rucker G (2015). \emph{Meta-Analysis with R}. Cham, Switzerland: Springer. <\doi{10.1007/978-3-319-21416-0}>

Xing A, Chu H, Lin L (2020). "Fragility index of network meta-analysis with application to smoking cessation data." \emph{Journal of Clinical Epidemiology}, \bold{127}, 29--39. <\doi{10.1016/j.jclinepi.2020.07.003}>
}
%%\author{
%%  ~~who you are~~
%%}
\note{
When the network meta-analysis contains many studies and many treatments, the assessment of fragility may be very computationally demanding. In such cases, one may focus on certain treatment comparisons that are of primary interest by specifying the arguments \code{tid1.f} and \code{tid2.f}. Also, if the significance or non-significance cannot be altered based on the given range of event status modifications, the fragility index and fragility quotient are defined to be \code{NA} (not available).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{frag.nma.alpha}} for assessing fragility of a network meta-analysis at multiple significance levels; \code{\link{frag.study}} for assessing fragility of an individual study; and \code{\link{frag.ma}} for assessing fragility of a meta-analysis.
}
\examples{
\donttest{
## Load datasets of network meta-analysis of
##  chronic obstructive pulmonary disease (COPD)
data(dat.copd)

## Assess fragility of network meta-analysis of COPD
out1 <- frag.nma(sid, tid, e, n, data = dat.copd)
out1

## Retrieve the estimated log odds ratio in the original dataset
out1$est.ori
## Retrieve the fragility index of each treatment comparison
out1$FI
## Retrieve the direction of significance change for comparison 2 vs. 3
out1$dir["2", "3"]
## Retrieve the study IDs with iterative event status modifications
##  for comparison 2 vs. 3
out1$sid.iter[["2 vs. 3"]]
## Retrieve the event status modifications
##  in the first group (i.e., treatment 2)
out1$tid1.mod.iter[["2 vs. 3"]]
## Retrieve the event status modifications
##  in the second group (i.e., treatment 3)
out1$tid2.mod.iter[["2 vs. 3"]]

## Visualize the process of event status modifications
##  for altering the non-significance to significance for each comparison;
##  see more options of usage in examples of frag.ma()
plot(out1, tid1 = 2, tid2 = 1)
plot(out1, tid1 = 3, tid2 = 1)
plot(out1, tid1 = 4, tid2 = 1)
plot(out1, tid1 = 3, tid2 = 2)
plot(out1, tid1 = 4, tid2 = 2)
plot(out1, tid1 = 4, tid2 = 3)

## Assess fragility of treatment comparisons 2 vs. 1 and 3 vs. 1 only
out2 <- frag.nma(sid, tid, e, n, data = dat.copd,
  tid1.f = c(2, 3), tid2.f = c(1, 1))
out2

## The default plot presents the process of event status modification
##  for the first comparison, i.e., 2 vs. 1
plot(out2)

## Load datasets of network meta-analysis of smoking cessation
data(dat.sc)

## Assess fragility of network meta-analysis of smoking cessation
out3 <- frag.nma(sid, tid, e, n, data = dat.sc)
out3

## Visualize the process of event status modifications
plot(out3, tid1 = 2, tid2 = 1, ybreaks = c(170, 595),
  x.legend = "topright")
plot(out3, tid1 = 3, tid2 = 1, ybreaks = c(635, 1200),
  x.legend = "bottomright")
plot(out3, tid1 = 4, tid2 = 1, ybreaks = c(105, 600))
plot(out3, tid1 = 3, tid2 = 2, ybreaks = c(160, 1205))
plot(out3, tid1 = 4, tid2 = 2, ybreaks = c(110, 140),
  x.legend = "topright")
plot(out3, tid1 = 4, tid2 = 3, ybreaks = c(130, 1205),
  x.legend = "bottomright")
}
}
\keyword{fragility}