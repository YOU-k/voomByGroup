#sce_tmp as reference
sim_exp_lib <- function(sce=sce_tmp,nCells=Ncells){
  lib.sizes <- colSums(counts(sce))
  if (length(lib.sizes) > 5000) {
    message("NOTE: More than 5000 cells provided. ",
            "5000 sampled library sizes will be used to test normality.")
    lib.sizes.sampled <- sample(lib.sizes, 5000, replace = FALSE)
  } else {
    lib.sizes.sampled <- lib.sizes
  }
  
  norm.test <- shapiro.test(lib.sizes.sampled)
  lib.norm <- norm.test$p.value > 0.2
  
  if (lib.norm) {
    fit <- fitdistrplus::fitdist(lib.sizes, "norm")
    lib.loc <- unname(fit$estimate["mean"])
    lib.scale <- unname(fit$estimate["sd"])
    message("NOTE: Library sizes have been found to be normally ",
            "distributed instead of log-normal. You may want to check ",
            "this is correct.")
  } else {
    fit <- fitdistrplus::fitdist(lib.sizes, "lnorm")
    lib.loc <- unname(fit$estimate["meanlog"])
    lib.scale <- unname(fit$estimate["sdlog"])
  }
  
    if (lib.norm) {
      exp.lib.sizes <- rnorm(nCells, lib.loc, lib.scale)
      min.lib <- min(exp.lib.sizes[exp.lib.sizes > 0])
      exp.lib.sizes[exp.lib.sizes < 0] <- min.lib / 2
    } else {
      exp.lib.sizes <- rlnorm(nCells, lib.loc, lib.scale)
    }
  return(exp.lib.sizes)
}
