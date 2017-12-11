\name{plot.BinQmet}
\alias{plot.BinQmet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Graphical representation of the results of the statistical analysis of binary Q-methodology data}
\description{Draws the graphs provided by the statistical analysis of binary Q-methodology data. Tunning parameters such as the colors representing the 'representative of the concept' and 'not representative of the concept' areas or the axis of MFA to be plotted can be specified. The function also allows to plot interactive outputs realized with plotly.}
\usage{
plot.BinQmet(res, choice = "all", interact = FALSE, col.pos.ratings = NULL, col.neg.ratings = NULL, col.clust.part = NULL, axis = c(1, 2), ext.dev.Rstudio = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
An object of class BinQmet returned by the function BinQmet.}
  \item{choice}{
A character element specifying the graphs to be plotted ('stim' for the concept representation mapping, 'part.seg' for the representation of the segmentation process of the participants (dendrogram + partitioning if a consolidation process has been implemented), 'part.mul' for the multidimensional representation of the structure of disagreement (PCA) among the panel of participants, 'all' for all representations. By default, all representations are plotted.)}
  \item{interact}{
A boolean specifying if the graphical outputs should be interactive (with plotly) or not.}
  \item{col.pos.ratings}{
A color element that can be specified by a character such as "blue" or by a hexadecimal code such as "#0D10D5". All hexadecimal codes can be found at \url{http://htmlcolorcodes.com}. By default, no color is specified and default color is used.}
  \item{col.neg.ratings}{
A color element that can be specified by a character such as "blue" or by a hexadecimal code such as "#0D10D5". All hexadecimal codes can be found at \url{http://htmlcolorcodes.com}. By default, no color is specified and default color is used.}
  \item{col.clust.part}{
A vector with as many color as clusters of participants. Colors can be specified by a character such as "blue" or by a hexadecimal code such as "#0D10D5". All hexadecimal codes can be found at \url{http://htmlcolorcodes.com}. By default, no colors are specified and default colors are used.}
  \item{axis}{
A length 2 numeric vector specifying the MFA components to be plotted. By default, the first 2 components are plotted.}
\item{ext.dev.Rstudio}{
A boolean specifying if the graphical outputs should be plotted in the Rstudio plot pane or not.}
}
\value{
Returns the graphical representations.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Margot Brard \email{margot.brard@agrocampus-ouest.fr}}
\seealso{
\code{\link{analyse.BinQmet}}
}
\examples{
# Small pedagogic data set
data(pedagdata)
res.pedag <- analyse.BinQmet(pedagdata, id.info.stim = 21 : ncol(pedagdata), type.info.stim = c(rep("cont", 4), "cat"), id.info.part = 9 : nrow(pedagdata), type.info.part = c(rep("cat", 2), "cont"))
plot.BinQmet(res.pedag, choice = "stim", col.pos.ratings = "orange")
plot.BinQmet(res.pedag, choice = "stim", interact = TRUE, col.pos.ratings = "orange")
plot.BinQmet(res.pedag, choice = "part.mul", interact = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
