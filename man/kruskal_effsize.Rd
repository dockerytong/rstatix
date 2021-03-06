% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kruskal_effesize.R
\name{kruskal_effsize}
\alias{kruskal_effsize}
\title{Kruskal-Wallis Effect Size}
\usage{
kruskal_effsize(
  data,
  formula,
  ci = FALSE,
  conf.level = 0.95,
  ci.type = "perc",
  nboot = 1000
)
}
\arguments{
\item{data}{a data.frame containing the variables in the formula.}

\item{formula}{a formula of the form \code{x ~ group} where \code{x} is a
numeric variable giving the data values and \code{group} is a factor with
one or multiple levels giving the corresponding groups. For example,
\code{formula = TP53 ~ cancer_group}.}

\item{ci}{If TRUE, returns confidence intervals by bootstrap. May be slow.}

\item{conf.level}{The level for the confidence interval.}

\item{ci.type}{The type of confidence interval to use. Can be any of "norm",
"basic", "perc", or "bca". Passed to \code{boot::boot.ci}.}

\item{nboot}{The number of replications to use for bootstrap.}
}
\value{
return a data frame with some of the following columns: \itemize{
 \item \code{.y.}: the y variable used in the test. \item \code{n}: Sample
 counts. \item \code{effsize}: estimate of the effect size. \item
 \code{magnitude}: magnitude of effect size. \item \code{conf.low,conf.high}:
 lower and upper bound of the effect size confidence interval.}
}
\description{
Compute the effect size for Kruskal-Wallis test as the eta
 squared based on the H-statistic: \code{eta2[H] = (H - k + 1)/(n - k)};
 where \code{H} is the value obtained in the Kruskal-Wallis test; \code{k} is
 the number of groups; \code{n} is the total number of observations.


 The eta-squared estimate assumes values from 0 to 1 and multiplied by 100%
 indicates the percentage of variance in the dependent variable explained by
 the independent variable. The interpretation values commonly in published
 litterature are: \code{0.01- < 0.06} (small effect), \code{0.06 - < 0.14}
 (moderate effect) and \code{>= 0.14} (large effect).

Confidence intervals are calculated by bootstap.
}
\examples{
# Load data
#:::::::::::::::::::::::::::::::::::::::
data("ToothGrowth")
df <- ToothGrowth

# Kruskal-wallis rank sum test
#:::::::::::::::::::::::::::::::::::::::::
df \%>\% kruskal_effsize(len ~ dose)

# Grouped data
df \%>\%
  group_by(supp) \%>\%
  kruskal_effsize(len ~ dose)
}
\references{
Maciej Tomczak and Ewa Tomczak. The need to report effect size
 estimates revisited. An overview of some recommended measures of effect
 size. Trends in Sport Sciences. 2014; 1(21):19-25.

 http://imaging.mrc-cbu.cam.ac.uk/statswiki/FAQ/effectSize

 http://www.psy.gla.ac.uk/~steve/best/effect.html
}
