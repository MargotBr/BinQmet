The `BinQmet` package is dedicated to the statistical analysis of binary Q-method data. It first analyses data trough a stimulus-oriented approach. In this context, it considers a Multiple Factor Analysis that provides a multidimensional representation of the stimuli. On this representation, two stimuli are close (resp. distant) if they have been perceived as similarly representative (resp. differently representative) of the concept by the participants. This representation is supplemented by coloured areas representing the degree of representativeness of the concept on the factorial plane. The BinQmethod package then analyses data trough a participant-oriented approach. In this context, the structure of disagreement among the panel of participants is captured through the profiles of residuals of a no-latent class regression model adjusted on the entire set of binary ratings, and can be visualized by using exploratory data analysis tools. The disagreement between two participants is then quantify in a concise way through the Euclidean distance between their respective profiles of residuals, this disagreement index being used as a basis to construct a dendrogram representing the structure of disagreement among the panel. The proper number of disagreed clusters among the panel of participants is then chosen by implementing a sequential strategy to test the significance of each K-clusters structure of disagreement. Finally, this participant-oriented approach provides a segmentation of the participants, each cluster of participants being assimilated to a pattern of perception of the concept through the set of stimuli.

# <span style="color: #EA485C">Installing BinQmet</span>

To get the current development version from GitHub:

  ```{r eval=FALSE}
install.packages("devtools")
devtools::install_github("MargotBr/BinQmet", build_vignettes = TRUE)
library(BinQmet)
```

# <span style="color: #EA485C">Using BinQmet</span>

To get an overview of the functionalities of the package, read the corresponding vignette:

  ```{r eval=FALSE}
vignette("BinQmet")
```
# BinQmet
