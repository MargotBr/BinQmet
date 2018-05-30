\name{BinQmet-package}
\alias{BinQmet-package}
\alias{BinQmet}
\docType{package}
\title{
\packageTitle{BinQmet}
}
\description{
\packageDescription{BinQmet}
}

\author{
Margot Brard, Sebastien Le, David Causeur

Maintainer: Margot Brard \email{margot.brard@agrocampus-ouest.fr}
}

\keyword{ package }
\seealso{
\code{\link{AnalyseBinQmet}}
\code{\link{plot.AnalyseBinQmet}}
}
\examples{
# Small pedagogic data set
data(pedagdata)
res.pedag <- AnalyseBinQmet(pedagdata, id.info.stim = 21 : ncol(pedagdata), type.info.stim = c(rep("cont", 4), "cat"), id.info.part = 9 : nrow(pedagdata), type.info.part = c(rep("cat", 2), "cont"))
res.pedag
plot.AnalyseBinQmet(res.pedag, choice = "stim", interact = TRUE, col.pos.ratings = "orange")
plot.AnalyseBinQmet(res.pedag, choice = "part.mul", interact = TRUE)
}
