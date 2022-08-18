voomByGroup <- function (counts, group = NULL, design = NULL, lib.size = NULL, dynamic = NULL, normalize.method = "none", 
                         span = 0.5, save.plot = FALSE, print = TRUE, plot = c("none", "all", "separate", "combine"), 
                         col.lines = NULL, pos.legend = c("inside", "outside", "none"), 
                         fix.y.axis = FALSE, ...) 
  # 14 June 2017 (Last updated 6 May 2022)
  # Charity Law, Xueyi Dong and Yue You
{
  # Counts
  out <- list()
  if (is(counts, "DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    if(is.null(group))
      group <- counts$samples$group
    # if (is.null(design) && diff(range(as.numeric(counts$sample$group))) > 0) 
    #   design <- model.matrix(~group, data = counts$samples)
    if (is.null(lib.size)) 
      lib.size <- with(counts$samples, lib.size * norm.factors)
    counts <- counts$counts
  }
  else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts, "ExpressionSet"))
    if (isExpressionSet) {
      if (length(Biobase::fData(counts))) 
        out$genes <- Biobase::fData(counts)
      if (length(Biobase::pData(counts))) 
        out$targets <- Biobase::pData(counts)
      counts <- Biobase::exprs(counts)
    }
    else {
      counts <- as.matrix(counts)
    }
  }
  if (nrow(counts) < 2L) 
    stop("Need at least two genes to fit a mean-variance trend")
  # Library size
  if(is.null(lib.size))
    lib.size <- colSums(counts)
  # Group  
  if(is.null(group))
    group <- rep("Group1", ncol(counts))
  group <- as.factor(group)
  intgroup <- as.integer(group)
  levgroup <- levels(group)
  ngroups <- length(levgroup)
  # Design matrix  
  if (is.null(design)) {
    design <- matrix(1L, ncol(counts), 1)
    rownames(design) <- colnames(counts)
    colnames(design) <- "GrandMean"
  }
  # Dynamic
  if (is.null(dynamic)) {
    dynamic <- rep(FALSE, ngroups)
  }
  # voom by group  
  if(print)
    cat("Group:\n")
  E <- w <- counts
  xy <- line <- as.list(rep(NA, ngroups))
  names(xy) <- names(line) <- levgroup
  for (lev in 1L:ngroups) {
    if(print)
      cat(lev, levgroup[lev], "\n")
    i <- intgroup == lev
    countsi <- counts[, i]
    libsizei <- lib.size[i]
    designi <- design[i, , drop = FALSE]
    QR <- qr(designi)
    if(QR$rank<ncol(designi))
      designi <- designi[,QR$pivot[1L:QR$rank], drop = FALSE]
    if(ncol(designi)==ncol(countsi))
      designi <- matrix(1L, ncol(countsi), 1)
    voomi <- voom(counts = countsi, design = designi, lib.size = libsizei, normalize.method = normalize.method, 
                  span = span, plot = FALSE, save.plot = TRUE, ...)
    E[, i] <- voomi$E
    w[, i] <- voomi$weights
    xy[[lev]] <- voomi$voom.xy
    line[[lev]] <- voomi$voom.line
  }
  #voom overall
  if (TRUE %in% dynamic){
    voom_all <- voom(counts = counts, design = design, lib.size = lib.size, normalize.method = normalize.method, 
                     span = span, plot = FALSE, save.plot = TRUE, ...)
    E_all <- voom_all$E
    w_all <- voom_all$weights
    xy_all <- voom_all$voom.xy
    line_all <- voom_all$voom.line
    
    dge <- DGEList(counts)
    disp <- estimateCommonDisp(dge)
    disp_all <- disp$common
  }
  # Plot, can be "both", "none", "separate", or "combine"
  plot <- plot[1]
  if(plot!="none"){
    disp.group <- c()
    for (lev in levgroup) {
      dge.sub <- DGEList(counts[,group == lev])
      disp <- estimateCommonDisp(dge.sub)
      disp.group[lev] <- disp$common
    }
    if(plot %in% c("all", "separate")){
      if (fix.y.axis == TRUE) {
        yrange <- sapply(levgroup, function(lev){
          c(min(xy[[lev]]$y), max(xy[[lev]]$y))
        }, simplify = TRUE)
        yrange <- c(min(yrange[1,]) - 0.1, max(yrange[2,]) + 0.1)
      }
      for (lev in 1L:ngroups) {
        if (fix.y.axis == TRUE){
          plot(xy[[lev]], xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", pch = 16, cex = 0.25, ylim = yrange)
        } else {
          plot(xy[[lev]], xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", pch = 16, cex = 0.25)
        }
        title(paste("voom: Mean-variance trend,", levgroup[lev]))
        lines(line[[lev]], col = "red")
        legend("topleft", bty="n", paste("BCV:", round(sqrt(disp.group[lev]), 3)), text.col="red")
      }
    }
    
    
    if(plot %in% c("all", "combine")){
      if(is.null(col.lines))
        col.lines <- 1L:ngroups
      if(length(col.lines)<ngroups)
        col.lines <- rep(col.lines, ngroups)
      xrange <- unlist(lapply(line, `[[`, "x"))
      xrange <- c(min(xrange)-0.3, max(xrange)+0.3)
      yrange <- unlist(lapply(line, `[[`, "y"))
      yrange <- c(min(yrange)-0.1, max(yrange)+0.3)
      plot(1L,1L, type="n", ylim=yrange, xlim=xrange, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )")
      title("voom: Mean-variance trend")
      if (TRUE %in% dynamic){
        for (dy in which(dynamic)){
          line[[dy]] <- line_all 
          disp.group[dy] <- disp_all
          levgroup[dy] <- paste0(levgroup[dy]," (all)")
        }
        
      }
      for (lev in 1L:ngroups) 
        lines(line[[lev]], col=col.lines[lev], lwd=2)
      pos.legend <- pos.legend[1]
      disp.order <- order(disp.group, decreasing = TRUE)
      text.legend <- paste(levgroup, ", BCV: ", round(sqrt(disp.group), 3), sep="")
      if(pos.legend %in% c("inside", "outside")){
        if(pos.legend=="outside"){
          plot(1,1, type="n", yaxt="n", xaxt="n", ylab="", xlab="", frame.plot=FALSE)
          legend("topleft", text.col=col.lines[disp.order], text.legend[disp.order], bty="n")
        } else {
          legend("topright", text.col=col.lines[disp.order], text.legend[disp.order], bty="n")
        }
      }
    }
  }
  # Output  
  if (TRUE %in% dynamic){
    E[,intgroup %in% which(dynamic)] <- E_all[,intgroup %in% which(dynamic)]
    w[,intgroup %in% which(dynamic)] <- w_all[,intgroup %in% which(dynamic)]
  }
  out$E <- E
  out$weights <- w
  out$design <- design
  if(save.plot){
    out$voom.line <- line
    out$voom.xy <- xy
  }
  new("EList", out)
}
