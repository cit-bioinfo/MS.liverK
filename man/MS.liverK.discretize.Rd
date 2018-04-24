\name{MS.liverK.discretize}
\alias{MS.liverK.discretize}
\title{
discretizes the signatures output from \code{MS.liverK.subtypes}
}
\description{
Given the results obtained with \code{liverCancerSbtypes()}, this function discretizes all continuous outputs, on a scale from 1 to 5. The lowest value for each signature is set to 1 and the hisghest to 5.
}
\usage{
MS.liverK.discretize(subtypes)
}

\arguments{
  \item{subtypes}{
  unsup element from \code{MS.liverK.subtypes()}
}
}

\value{
a new subtypes list with the same organization, where the value from signatures has been discretized with 5 levels.
}
\author{
F Petitprez
}
\note{
This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact CITR@ligue-cancer.net
}
\keyword{ MS.liverK.discretize }
