plot.AnalyseBinQmet <- function(res, choice = "all", interact = FALSE, col.pos.ratings = NULL, col.neg.ratings = NULL, col.clust.part = NULL, axis = c(1, 2), ext.dev.Rstudio = FALSE, vignette = FALSE) {

  options(warn = -1)

  # load packages
  suppressPackageStartupMessages(require(grid, quietly = TRUE))
  suppressPackageStartupMessages(require(gridExtra, quietly = TRUE))
  suppressPackageStartupMessages(require(ggplot2, quietly = TRUE))
  suppressPackageStartupMessages(require(ggrepel, quietly = TRUE))
  suppressPackageStartupMessages(require(plotly, quietly = TRUE))
  suppressPackageStartupMessages(require(reshape2, quietly = TRUE))
  suppressPackageStartupMessages(require(doBy, quietly = TRUE))
  if(!suppressWarnings(suppressPackageStartupMessages(require("AgreeClust", quietly = TRUE)))){
    devtools::install_github("MargotBr/AgreeClust", build_vignettes = TRUE)
    suppressPackageStartupMessages(require("AgreeClust", character.only=TRUE))
  }

  # check the format of the arguments
  if (!inherits(res, "AnalyseBinQmet")) {
    stop("Non convenient data - res should be an AnalyseBinQmet object")
  }
  choice <- match.arg(choice, c("all", "stim", "part.seg", "part.mul"), several.ok = TRUE)
  mat.partition <- cbind.data.frame(res$res.AgreeClust$partition, names(res$res.AgreeClust$partition))
  colnames(mat.partition) <- c("Cluster", "Rater")
  mat.partition$Cluster <- as.factor(mat.partition$Cluster)
  nb.clust <- nlevels(mat.partition$Cluster)
  if (!is.null(col.clust.part)) {
    if (length(col.clust.part) < nb.clust) {
      stop("Non convenient specification of colors - col.clust should contain at least as many elements as clusters")
    }
  }
  if (length(axis) != 2 | length(unique(axis)) != 2 | class(axis) != "numeric") {
    stop("Non convenient specification of axis - axis should be a numeric vector of 2 different elements")
  }

  # calculate the numbers of raters and stimuli
  if (!is.null(res$call$id.info.part)) {
    dta <- res$call$dta[-res$call$id.info.part,]
    dta <- droplevels(dta)
  }
  if (!is.null(res$call$id.info.stim)) {
    dta <- dta[, -res$call$id.info.stim]
    dta <- droplevels(dta)
  }
  if(is.null(res$call$id.info.part) & is.null(res$call$id.info.stim)) {
    dta <- res$call$dta
    dta <- droplevels(dta)
  }
  nbrater <- ncol(dta)
  nbstim <- nrow(dta)

  # stimulus-oriented analysis
  if (!is.na(match("all", choice)) | !is.na(match("stim", choice))) {
    if (is.null(col.pos.ratings)) {
      col.pos.ratings <- "#EA485C"
    }
    if (is.null(col.neg.ratings)) {
      col.neg.ratings <- "#A9A9A9"
    }
    res.mfa <- res$res.mfa
    coord.stim <- as.data.frame(res.mfa$ind$coord[, axis])
    colnames(coord.stim) <- c("axeA", "axeB")
    resolution <- 200
    coord.and.categories <- cbind.data.frame(rownames(res$BinQmet.data), res.mfa$ind$coord[rownames(res$BinQmet.data), axis], res$BinQmet.data)
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
    if (interact == FALSE) {
      plot.stim <- ggplot(NULL) +
        labs(x = paste("Dim ", axis[1]," - ", round(res.mfa$eig[axis[1], 2], 2), " %", sep = ""), y = paste("Dim ", axis[2], " - ", round(res.mfa$eig[axis[2], 2], 2), " %", sep = "")) +
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
      if ((Sys.getenv("RSTUDIO") == "1") == FALSE | ext.dev.Rstudio == TRUE) {
        dev.new(noRStudioGD = TRUE)
      }
      print(plot.stim)
    } else if (interact == TRUE) {
      coord.and.categories <- coord.and.categories[1 : 3]
      colnames(coord.and.categories) <- c("Stimulus", "f1", "f2")
      col <- cbind.data.frame(seq(0, 1, by = 0.01), c(col.neg, "#FFFFFF", "#FFFFFF", col.pos))
      text.tooltip <- paste("Stimulus:", coord.and.categories$Stimulus)
      if (!is.null(res$call$id.info.part)) {
        dta <- res$call$dta[-res$call$id.info.part,]
        dta <- droplevels(dta)
      }
      if (!is.null(res$call$id.info.stim)) {
        dta <- dta[, -res$call$id.info.stim]
        dta <- droplevels(dta)
      }
      melted.data <- melt(as.matrix(dta))
      colnames(melted.data) <- c("Stimulus", "Participant", "Rating")
      mat.partition <- cbind.data.frame(names(res$res.AgreeClust$partition), res$res.AgreeClust$partition)
      colnames(mat.partition) <- c("Participant", "Cluster")
      melted.cluster <- merge(melted.data, mat.partition, by = "Participant")
      melted.cluster$Rating <- as.numeric(as.character(melted.cluster$Rating))
      sum.pos.ratings <- summaryBy(Rating ~ Stimulus : Cluster, data = melted.cluster, FUN = "sum")
      colnames(sum.pos.ratings) <- c("Stimulus", "Cluster", "NbPosRating")
      size.clust <- summary(as.factor(mat.partition$Cluster))
      for (i in 1 : nlevels(as.factor(mat.partition$Cluster))) {
        tab.clust <- sum.pos.ratings[which(sum.pos.ratings$Cluster == levels(as.factor(mat.partition$Cluster))[i]), ]
        tab.clust[, "NbPosRating"] <- tab.clust[, "NbPosRating"] / size.clust[which(names(size.clust) == levels(as.factor(mat.partition$Cluster))[i])] * 100
        text.tooltip <- paste(text.tooltip,
                              paste(paste0("<br>Positive ratings in Cluster ", levels(as.factor(mat.partition$Cluster))[i], ":"),
                                    paste0(round(tab.clust[order(tab.clust$Stimulus, coord.and.categories$Stimulus), "NbPosRating"], 1), "%")))
      }
      if (!is.null(res$call$id.info.stim)) {
        info.stim <- res$call$dta[, res$call$id.info.stim]
        for (i in 1 : length(res$call$id.info.stim)) {
          info <- info.stim[1 : length(coord.and.categories$Stimulus), i][order(info.stim[1 : length(coord.and.categories$Stimulus), i], coord.and.categories$Stimulus)]
          text.tooltip <- paste(text.tooltip, paste0("<br>", colnames(info.stim)[i], ": ", as.character(unlist(info))))
        }
      }
      plot.stim.interact <- plot_ly(mat.surface) %>%
        add_trace(mat.surface,
                  x = mat.surface$f1 ,
                  y = mat.surface$f2,
                  z = mat.surface$z,
                  type = "contour",
                  contours = list(showlabels = TRUE, start = 0, end = 100, coloring = "heatmap"),
                  colorscale = col,
                  hoverinfo = "none",
                  hoverlabel = list(bgcolor = "black", font = list(color = "white")),
                  colorbar = list(len = 1,
                                  lenmode = "fraction",
                                  title = "% of participants who assessed \n this stimulus as representative \n of the concept")) %>%
        add_trace(coord.and.categories,
                  x = coord.and.categories$f1 ,
                  y = coord.and.categories$f2,
                  name = "stimulus",
                  hoverinfo = 'text',
                  text = text.tooltip,
                  type = "scatter", mode = "markers",
                  marker = list(color = "black", size = 4, line = list(color = 'black', width = 2)),
                  showlegend = FALSE) %>%
        layout(
          titlefont = list(size = 14, color = "#444444"),
          title = "<b>Concept representation mapping<b>",
          xaxis = list(zerolinecolor = "white", scaleanchor = "y", showgrid = FALSE, title = paste("Dim ", axis[1], " - ", round(res.mfa$eig[axis[1],2],2), "%", sep=""), titlefont = list(color = "#444444"), tickfont = list(size = 8, color = "#444444"), showline = TRUE, mirror = "ticks", linecolor = "#444444", linewidth = 1),
          yaxis = list(zerolinecolor = "white", scaleanchor = "x", showgrid = FALSE, title = paste("Dim ", axis[2], " - ", round(res.mfa$eig[axis[2],2],2), "%", sep=""), titlefont = list(color = "#444444"), tickfont = list(size = 8, color = "#444444"), showline = TRUE, mirror = "ticks", linecolor = "#444444", linewidth = 1),
          legend = list(yanchor = "top")
        )
      if (vignette == FALSE) {
        if (ext.dev.Rstudio == TRUE) {
          old.viewer <- options()$viewer
          options(viewer = NULL)
          print(plot.stim.interact)
          options(viewer = old.viewer)
        } else {
          print(plot.stim.interact)
        }
      }
    }
  }

  # participant-oriented analysis
  if (!is.na(match("all", choice)) | !is.na(match("part.seg", choice))) {
    res.plotAgreeClust <- plot.AgreeClust(res$res.AgreeClust, choice = "seg", interact = interact, col.clust = col.clust.part, axis = axis, name.rater = "participant", ext.dev.Rstudio = ext.dev.Rstudio, vignette = vignette)
  }
  if (!is.na(match("all", choice)) | !is.na(match("part.mul", choice))) {
    res.plotAgreeClust <- plot.AgreeClust(res$res.AgreeClust, choice = "mul", interact = interact, col.clust = col.clust.part, axis = axis, name.rater = "participant", ext.dev.Rstudio = ext.dev.Rstudio, vignette = vignette)
  }

  # end the function
  options(warn = 0)
  if (vignette == TRUE & interact == TRUE) {
    message("Representations plotted")
    if (!is.na(match("stim", choice))) {
      return(plot.stim.interact)
    } else if (!is.na(match("part.seg", choice)) | !is.na(match("part.mul", choice))) {
      return(res.plotAgreeClust)
    } else if (!is.na(match("all", choice))) {
      return(list(plot.stim.interact, res.plotAgreeClust))
    }
  } else {
    return(message("Representations plotted"))
  }

}
