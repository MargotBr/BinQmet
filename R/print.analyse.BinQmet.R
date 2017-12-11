print.analyse.BinQmet <- function (res){

  if (!inherits(res, "BinQmet")){
    stop("Non convenient data - res should be a BinQmet object")
  }

  cat("** Results for the analysis of the binary Q-method (BinQmet) **\n")
  cat("\n")
  cat("The analysis was performed on", (nrow(res$call$dta) - length(res$call$id.info.part)),
      "stimuli assessed by", (ncol(res$call$dta) - length(res$call$id.info.stim)), "participants\n")
  cat("The results are available in the following objects:\n\n")
  res.desc <- array("", c(5, 2), list(1 : 5, c("name", "description")))
  res.desc[1, ] <- c("$call", "arguments used in the BinQmet function")
  res.desc[2, ] <- c("$BinQmet.data", "multiple table corresponding to a binary Q-method dataset")
  res.desc[3, ] <- c("$res.mfa", "MFA results")
  res.desc[4, ] <- c("$concept.surface", "response surface (i.e. strength of the representativeness of the concept)")
  res.desc[5, ] <- c("$res.AgreeClust", "results of the segmentation of the panel of participants (AgreeClust package)")
  print(res.desc)

}
