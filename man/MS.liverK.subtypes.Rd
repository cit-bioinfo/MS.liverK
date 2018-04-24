\name{MS.liverK.subtypes}
\alias{MS.liverK.subtypes}
\title{
What are the subtypes of the provided cancer samples ?
}
\description{
Given a set of probes measures for several samples and the corresponding genes symbols, returns the molecular subtype of the samples according to various supervised or unsupervised classifications.
}
\usage{
MS.liverK.subtypes(probesData,probesSymbols = NULL, multipleProbeSep = " /// ",
    signatureChoice = c("Down","-Down","Up","Up-Down")[4], allHCC = TRUE,
    PrognosticSignatures = FALSE)
}

\arguments{
  \item{probesData}{
  data frame containing the expression data. Probes must be in lines and samples in columns. Probes names must be specified as rownames.
}
  \item{probesSymbols}{
  data frame containing the correspondance between probes (first column) and gene symbols (second column). They must be provided as characters, and not as factors. If expression data is already aggregated by gene and that gene symbols are set as row names, keep NULL (default).
}
  \item{multipleProbeSep}{
  The separator used for probes linked to several gene symbols.
}
  \item{signatureChoice}{
  One of "Down", "-Down", "Up", "Up-Down" (default), stipulating which score function is used to compute the signatures score : either the mean transcription of down-regulated genes, its opposite, the mean transcription of up-regulated genes or the difference between the two latter.
}
  \item{allHCC}{
  boolean. Are all samples in the cohort HCC ? Default to TRUE. If FALSE, a new category is added to the Boyault classification, corresponding to non tumoral or HCA samples.
}
  \item{PrognosticSignatures}{
  boolean. Return additional prognostic signatures scores ? Default to FALSE. If TRUE, additionalprognostic signatures are evaluated on the data, and are returned under prediction$signatures.
}
}

\details{
MS.liverK.subtypes function runs 6 unsupervised classifiers, designed and published by Lee et al., Hoshida et al., Chiang et al., Roessler et al. and Boyault et al. It also runs a series of supervised signatures on the data, and returns a continuous and a discontinuous score for them.

liverCancerSubtypes.testDataFormat checks whether the parameters are correctly formated and stops the execution if a problem is found. It is automatically run when calling liverCancerSubtypes.
}
\value{
List of 3:
 \item{molecularSubtypes}{List of 2: \code{prediction} contains the predicted subtypes for the 6 unsupervised classifications. \code{signatures} contains the continuous score for a series of gene signatures based on these subtypes.}
 \item{prognosis}{List of 2: \code{prediction} contains the two-class prognosis classification by Nault et Zucman-Rossi. \code{signatures} contains the continuous output for a series of prognosis genes signatures.}
 \item{biologicalPathwaysSignatures}{continuous output for biological pathways genes signatures}
}
\references{
Boyault, S., Rickman, D.S., de Reynies, A., Balabaud, C., Rebouissou, S., Jeannot, E., Herault, A., Saric, J., Belghiti, J., Franco, D., et al. (2007). Transcriptome classification of HCC is related to gene alterations and to new therapeutic targets. Hepatol. Baltim. Md 45, 42-52.

Chiang, D.Y., Villanueva, A., Hoshida, Y., Peix, J., Newell, P., Minguez, B., LeBlanc, A.C., Donovan, D.J., Thung, S.N., Sole, M., et al. (2008). Focal gains of VEGFA and molecular classification of hepatocellular carcinoma. Cancer Res. 68, 6779-6788.

Hoshida, Y., Nijman, S.M.B., Kobayashi, M., Chan, J.A., Brunet, J.-P., Chiang, D.Y., Villanueva, A., Newell, P., Ikeda, K., Hashimoto, M., et al. (2009). Integrative transcriptome analysis reveals common molecular subclasses of human hepatocellular carcinoma. Cancer Res. 69, 7385-7392.

Lee, J.-S., Chu, I.-S., Heo, J., Calvisi, D.F., Sun, Z., Roskams, T., Durnez, A., Demetris, A.J., and Thorgeirsson, S.S. (2004). Classification and prediction of survival in hepatocellular carcinoma by gene expression profiling. Hepatol. Baltim. Md 40, 667-676.

Nault, J.C., and Zucman Rossi, J. (2013). Molecular classification of hepatocellular adenomas. Int. J. Hepatol. 2013, 315947.

Roessler, S., Jia, H.-L., Budhu, A., Forgues, M., Ye, Q.-H., Lee, J.-S., Thorgeirsson, S.S., Sun, Z., Tang, Z.-Y., Qin, L.-X., et al. (2010). A unique metastasis gene signature enables prediction of tumor relapse in early-stage hepatocellular carcinoma patients. Cancer Res. 70, 10202-10212.

Subramanian, A., Tamayo, P., Mootha, V.K., Mukherjee, S., Ebert, B.L., Gillette, M.A., Paulovich, A., Pomeroy, S.L., Golub, T.R., Lander, E.S., et al. (2005). Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. Proc. Natl. Acad. Sci. U. S. A. 102, 15545-15550.
}
\author{
A de Reynies, F Petitprez
}
\note{
This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact CITR@ligue-cancer.net
}
%%\examples{
%%}
\keyword{ MS.liverK.subtypes }
