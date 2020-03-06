\name{AUDem}
\alias{AUDem}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Australian electricity demand data}
\description{
Data set on electricity demand from Sidney, Australia. The data has been downloaded from \url{https://www.ausgrid.com.au}, and it originally contained electricity demand from 300 customers, at 30min resolution. We discarded 53
customers because their demand was too irregular, and we integrated the demand data with temperature data from the
National Climatic Data Center, covering the same period. 
} %- end description

\usage{
data(AUDem)
}
%- maybe also `usage' for other objects documented here.
\format{\code{AUDem} is a list, where \code{AUDem$meanDem} is a \code{data.frame} containing the following variables:
\describe{
\item{doy}{the day of the year, from 1 to 365;}
\item{tod}{the time of day, ranging from 18 to 22, where 18 indicates the period from 17:00 to
17:30, 18.5 the period from 17:30 to 18:00 and so on;}
\item{dem}{the demand (in KW) during a 30min period, averaged over the 247 households;}
\item{dow}{factor variable indicating the day of the week;}
\item{temp}{the external temperature at Sidney airport, in degrees Celsius;}
\item{date}{local date and time;}
\item{dem48}{the lagged mean demand, that is the average demand (dem) during the same
30min period of the previous day;}
}
The second element is \code{AUDem$qDem48} which is a matrix with as many rows as \code{AUDem$meanDem}. Each rows contains 20 equally spaced empirical quantiles of the lagged individual electricity demand of the 247 customers.
}

\value{A list where \code{AUDem$meanDem} is a data.frame and \code{AUDem$qDem48} a matrix.
}
\examples{
library(qgam)
data(AUDem)
  
# Mean demand over the period
plot(AUDem$meanDem$dem, type = 'l')
  
# 20 quantiles of individual demand over 5 days
matplot(seq(0.01, 0.99, length.out = 20), 
        t(AUDem$qDem48[c(1, 50, 75, 100, 250), ]), 
        type = 'l', 
        ylab = "Electricity demand (KW)",
        xlab = expression("Probability level " * "(p)"), 
        lty = 1)
}
