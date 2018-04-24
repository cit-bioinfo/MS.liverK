\name{MS.liverK.convert}
\alias{MS.liverK.convert}
\title{
converts a dataset by mapping it on a platform used to classify samples
}
\description{
This function maps each gene of the used-provided dataset on the scale used by GSE20238 by applying a linear regression on deciles for each gene.
}
\usage{
MS.liverK.convert(userData)
}

\arguments{
  \item{userData}{
  A gene expression profile matrix, with genes in line and samples in colmumns, aggregated by genes. The genes, specified as HUGO gene symbols, must be put as row names.
}
}

\value{
a data frame containing the mapped gene expressions, that can be used as input for \code{MS.liverK.subtypes}.
}
\author{
F Petitprez
}
\note{
This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact CITR@ligue-cancer.net
}
\keyword{ MS.liverK.convert }
