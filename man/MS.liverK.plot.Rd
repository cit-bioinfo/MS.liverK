\name{MS.liverK.plot}
\alias{MS.liverK.plot}
\title{
plots the results of \code{MS.liverK.subtypes}
}
\description{
Given the results obtained with \code{MS.liverK.subtypes()}, this function plots a segment plot and heatmaps representing the repartition and correspondance of classes.
}
\usage{
MS.liverK.plot(subtypes,col=NULL)
}

\arguments{
  \item{subtypes}{
  unsup element from \code{MS.liverK.subtypes()}
}
  \item{col}{
  named vector containing the colors to use for each class. Names must contain each of \code{"A","B","subgroup A","subgroup B","S1","S2","S3","G1","G2","G3","G4","G5","G6","CTNNB1","Inflammation","Polysomy chr7","Proliferation","Unannotated","+","-"}, for Lee's, Roessler's, Hoshida's, Chiang's, Boyault's and EPCAM.AFP classifications.
}
}

\value{
a plot of all results from subtypes. Returns \code{NULL}.
}
\author{
F Petitprez
}
\note{
This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact CITR@ligue-cancer.net
}
\keyword{ MS.liverK.plot }
