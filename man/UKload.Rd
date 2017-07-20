\name{UKload}
\alias{UKload}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{UK electricity load data}
\description{
Dataset on UK electricity demand, taken from the national grid (\url{http://www2.nationalgrid.com/}).
} %- end description

\usage{
data(UKload)
}
%- maybe also `usage' for other objects documented here.
\format{ \code{UKload} contains the following variables:
\describe{
\item{NetDemand}{net electricity demand between 11:30am and 12am.}
\item{wM}{instantaneous temperature, averaged over several English cities.}
\item{wM_s95}{exponential smooth of \code{wM}, that is \code{wM_s95[i] = a*wM[i] + (1-a)*wM_s95[i]} with \code{a=0.95}}.
\item{Posan}{periodic index in \code{[0, 1]} indicating the position along the year.}
\item{Dow}{factor variable indicating the day of the week.}
\item{Trend}{progressive counter, useful for defining the long term trend.}
\item{NetDemand.48}{lagged version of \code{NetDemand}, that is \code{NetDemand.48[i] = NetDemand[i-2]}.}
\item{Holy}{binary variable indicating holidays.}
\item{Year}{should be obvious.}
\item{Date}{should be obvious.}
}

}

\details{ See Fasiolo et al. (2016) for details.}

\value{matrix of replicate data series
}
\references{
Fasiolo, M., Goude, Y., Nedellec, R. and Wood, S. N. (2016). Fast calibrated additive quantile regression. Available at \url{https://github.com/mfasiolo/qgam/blob/master/draft_qgam.pdf}.
}
\examples{
  library(qgam)
  data(UKload)
  plot(UKload$NetDemand, type = 'l')
}
