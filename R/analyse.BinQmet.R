analyse.BinQmet <- function(dta, id.info.stim = NULL, type.info.stim = NULL, id.info.part = NULL, type.info.part = NULL, graph = TRUE, ext.dev.Rstudio = FALSE, ...) {

  options(warn = -1)

  # load packages
  suppressPackageStartupMessages(require(FactoMineR, quietly = TRUE))
  suppressPackageStartupMessages(require(reshape2, quietly = TRUE))
  suppressPackageStartupMessages(require(ggplot2, quietly = TRUE))
  suppressPackageStartupMessages(require(ggrepel, quietly = TRUE))
  suppressPackageStartupMessages(require(AgreeClust, quietly = TRUE))
  if(!suppressWarnings(suppressPackageStartupMessages(require("AgreeClust", quietly = TRUE)))){
    devtools::install_github("MargotBr/AgreeClust", build_vignettes = TRUE)
    suppressPackageStartupMessages(require("AgreeClust", character.only=TRUE))
  }

  # save the data set
  dta.sauv <- dta

  # remove external information about raters and stimuli
  if (!is.null(id.info.part)) {
    dta <- dta[-id.info.part,]
    dta <- droplevels(dta)
  }
  if (!is.null(id.info.stim)) {
    dta <- dta[, -id.info.stim]
    dta <- droplevels(dta)
  }

  # calculate the numbers of raters and stimuli
  nbrater <- ncol(dta)
  nbstim <- nrow(dta)

  # create a res object to save the results
  res <- list()

  # return the important arguments
  res[[1]] <- list(dta.sauv, id.info.part, type.info.part, id.info.stim, type.info.stim)
  names(res[[1]]) <- c("dta", "id.info.part", "type.info.part", "id.info.stim", "type.info.stim")

  # create the factorial map of the stimuli
  axes <- 1 : 2
  dta.quanti <- apply(apply(dta, 2, as.character), 2, as.numeric)
  sum.pos.ratings <- apply(dta.quanti, 1, sum)
  sum.neg.ratings <- ncol(dta) - sum.pos.ratings
  dta.binQ <- cbind.data.frame(dta, sum.pos.ratings, sum.neg.ratings)
  res[[2]] <- dta.binQ
  colnames(dta.binQ)[(ncol(dta) + 1) : ncol(dta.binQ)] <- c("1", "0")
  res.mfa <- MFA(dta.binQ, group = c(ncol(dta), 2), type <- c("n", "f"), name.group = c("groups", "association"), graph = FALSE)
  coord.stim <- as.data.frame(res.mfa$ind$coord[, axes])
  colnames(coord.stim) <- c("axeA", "axeB")
  res[[3]] <- res.mfa

  # project the positive ratings area
  col.pos.ratings <- "#EA485C"
  col.neg.ratings <- "#A9A9A9"
  resolution <- 200
  coord.and.categories <- cbind.data.frame(rownames(dta), res.mfa$ind$coord[rownames(dta), axes], dta)
  x1 <- coord.and.categories[, 2]
  x2 <- coord.and.categories[, 3] # coordonnee 2
  x12 <- scale(x1, center = TRUE, scale = FALSE)[, ] * scale(x2,center = TRUE, scale = FALSE)[, ] 	# interaction between 1st & 2nd coordinate
  XX <- cbind.data.frame(x1, x2, x12) # expliocative variables
  size.x1 <- diff(range(x1))
  size.x2 <- diff(range(x2))
  by <- max(size.x1, size.x2) / resolution
  f1 <- seq((min(x1) - size.x1 * 0.05), (max(x1) + size.x1 * 0.05), by) # coordinates fictive points
  f2 <- seq((min(x2) - size.x2 * 0.05), (max(x2) + size.x2 * 0.05), by) # coordinates fictive points
  grille <- expand.grid(f1 = f1, f2 = f2) # gcoordinates fictive points
  grille.x1 <- grille[, 1]
  grille.x2 <- grille[, 2]
  grille.x12 <- scale(grille.x1, center = TRUE, scale = FALSE)[, ] * scale(grille.x2, center = TRUE, scale = FALSE)[, ]
  grille.XX <- cbind.data.frame(grille.x1, grille.x2, grille.x12)
  colnames(grille.XX) <- colnames(XX)
  compute.nb.pos.by.rater <- function(i) {
    ratings <- coord.and.categories[, i + 3] # ratings of the rater
    dta.mod <- cbind.data.frame(XX, ratings)
    mod <- glm(ratings ~ x1 + x2 + x12, data = dta.mod, family = binomial) # logistic regression
    predict.prob.ratings.fictive <- predict.glm(object = mod, newdata = grille.XX, type = "response")
    predict.ratings.fictive <- as.numeric(predict.prob.ratings.fictive >= 0.5)
    nb.pos.by.rater <- matrix(predict.ratings.fictive, nrow = length(f1), ncol = length(f2))
    return(nb.pos.by.rater)
  }
  list.nb.pos.by.rater <- lapply(1 : nbrater, compute.nb.pos.by.rater)
  nb.pos.ratings <- Reduce("+", list.nb.pos.by.rater)
  dimnames(nb.pos.ratings) <- list(as.character(f1), as.character(f2))
  nb.pos.ratings <- nb.pos.ratings / nbrater * 100
  mat.surface <- melt(nb.pos.ratings)
  colnames(mat.surface) <- c("f1","f2","z")
  palette.col.neg <- colorRampPalette(c(col.neg.ratings, "white"))
  col.neg <- palette.col.neg(49)
  palette.col.pos <- colorRampPalette(c("white", col.pos.ratings))
  col.pos  <- palette.col.pos(50)
  col <- c(col.neg, col.pos)
  plot.stim <- ggplot(NULL) +
    labs(x = paste("Dim ", 1," - ", round(res.mfa$eig[axes[1], 2], 2), " %", sep = ""), y = paste("Dim ", 2, " - ", round(res.mfa$eig[axes[2], 2], 2), " %", sep = "")) +
    coord_fixed()+
    geom_raster(data = mat.surface, aes(f1, f2, fill = z)) +
    scale_fill_gradientn(colours = c(col.neg, "#FFFFFFFF", col.pos), name="% of participants who assessed \n this stimulus as representative \n of the concept", guide = "colorbar", limits = c(0, 100)) +
    geom_contour(data = mat.surface, aes(f1, f2, z = z), colour = "black") +
    geom_point(data = coord.stim, aes(axeA, axeB), size = 1) +
    geom_text_repel(data = coord.stim, aes(axeA, axeB, label = rownames(coord.stim)), size = 3) +
    ggtitle("Concept representation mapping") +
    theme(
      legend.key = element_rect(colour = "white", fill = "white"),
      panel.background = element_rect(fill = 'white', colour = "#444444"),
      panel.grid.major = element_line(colour = "white"),
      panel.grid.minor = element_line(colour = "white"),
      axis.text = element_text(colour = "#444444"),
      axis.ticks = element_line(colour = "#444444"),
      axis.title = element_text(colour = "#444444"),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12, color = "#444444")
    )
  if (graph == TRUE) {
    if ((Sys.getenv("RSTUDIO") == "1") == FALSE | ext.dev.Rstudio == TRUE) {
      dev.new(noRStudioGD = TRUE)
    }
    print(plot.stim)
  }

  # segment the panel of participants
  res.AgreeClust <- AgreeClustBin(dta = dta.sauv, id.info.rater = id.info.part, type.info.rater = type.info.part, id.info.stim = id.info.stim, type.info.stim = type.info.stim, graph = FALSE, ...)
  res[[4]] <- res.AgreeClust
  if (graph == TRUE) {
    plot.AgreeClust(res.AgreeClust, name.rater = "participant", ext.dev.Rstudio = ext.dev.Rstudio)
  }

  # return the results
  names(res) <- c("call", "BinQmet.data", "res.mfa", "res.AgreeClust")
  message("Analysis performed")
  options(warn = 0)
  class(res) <- c("BinQmet", "list ")
  return(res)

}
