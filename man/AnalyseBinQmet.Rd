\name{AnalyseBinQmet}
\alias{AnalyseBinQmet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Statistical analysis of binary Q-method data}
\description{Performs the statistical analysis of binary Q-method data through a stimulus-oriented approach and a participant-oriented approach. The method is able to manage covariates providing informations about the participants and the stimuli to facilitate the interpretation of the results.}
\usage{
AnalyseBinQmet(dta, id.info.stim = NULL, type.info.stim = NULL, id.info.part = NULL, type.info.part = NULL, graph = TRUE, ext.dev.Rstudio = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dta}{
A numeric matrix or data frame of dimensions S x P (S=number of stimuli, R=number of participants).}
 \item{id.info.stim}{
A vector of integer elements composed of the identification of the columns containing the supplementary information (i.e. covariates) about the stimuli. This argument is optional and, by default, it is fixed to NULL, meaning that dta does not contain supplementary information about the stimuli.}
  \item{type.info.stim}{
A vector of character elements composed of the type of the covariates about the stimuli. This vector must be of the same length that id.info.stim. A continuous covariate is associated to 'cont' and a categorical covariate is associated to 'cat'. This argument is optional and, by default, it is fixed to NULL, meaning that dta does not contain supplementary information about the stimuli.}
  \item{id.info.part}{
A vector of integer elements composed of the identification of the lines containing the supplementary information (i.e. covariates) about the participants. This argument is optional and, by default, it is fixed to NULL, meaning that dta does not contain supplementary information about the participants.}
  \item{type.info.part}{
A vector of character elements composed of the type of the covariates about the participants. This vector must be of the same length that id.info.part. A continuous covariate is associated to 'cont' and a categorical covariate is associated to 'cat'. This argument is optional and, by default, it is fixed to NULL, meaning that dta does not contain supplementary information about the participants.}
  \item{graph}{
A boolean specifying if the graphical outputs should be plotted or not. By default, they are plotted.}
  \item{ext.dev.Rstudio}{
A boolean specifying if the graphical outputs should be plotted in the Rstudio plot pane or not.}
\item{...}{
Supplementary arguments of the AgreeClustBin function of the AgreeClust package.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{BinQmet.data}{A matrix of dimensions S x (R + 2) (S=number of stimuli, R=number of participants) corresponding to the multiple table necessary to analyse binary Q-method data. The first table is of dimensions S x R and contains information about the structure into groups of stimuli, and the second table is of dimensions S x 2 and contains information about the degree of association with the concept.}
\item{res.mfa}{All the results of the MFA realized during the stimulus-oriented statistical analysis.}
\item{concept.surface}{The response surface computed during the stimulus-oriented statistical analysis. The response (predicted \% of participants assessing each point as representative of the concept) is expressed according to the coordinates of the points provided by MFA.}
\item{res.AgreeClust}{All the results of the agreement-based clustering of the participants realized during the participant-oriented statistical analysis.}
}

\references{
%% ~put references to the literature/web site here ~
}
\author{Margot Brard \email{margot.brard@agrocampus-ouest.fr}}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{plot.AnalyseBinQmet}}
}
\examples{
# Small pedagogic data set
data(pedagdata)
res.pedag <- AnalyseBinQmet(pedagdata, id.info.stim = 21 : ncol(pedagdata), type.info.stim = c(rep("cont", 4), "cat"), id.info.part = 9 : nrow(pedagdata), type.info.part = c(rep("cat", 2), "cont"))
res.pedag
plot.AnalyseBinQmet(res.pedag, interact = TRUE)

# 'Perfumes' data set
data(perfumes)
res.perfumes <- AnalyseBinQmet(perfumes, id.info.stim = 121 : ncol(perfumes), type.info.stim = rep("cont", 27), id.info.part = 40 : nrow(perfumes), type.info.part = rep("cat", 2))
res.perfumes
plot.AnalyseBinQmet(res.perfumes, interact = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
