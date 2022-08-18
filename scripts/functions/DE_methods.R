apply_voombygroup <- function(y,design,cond,contr.matrix,dynamic){
  source("/stornext/HPCScratch/home/you.y/simulation/functions/voomByGroup.R")
  voomByGroup(y, design = design,group=cond,plot="all",dynamic = dynamic) -> y
  fit <- lmFit(y, design)
  vfit <- contrasts.fit(fit, contrasts=contr.matrix)
  tfit <- treat(vfit, lfc=0)
  return(tfit)
}
  
apply_voomwithqualityweight <- function(y,design,contr.matrix){
  voomWithQualityWeights(y, design = design, plot=TRUE) -> y
  fit <- lmFit(y, design)
  vfit <- contrasts.fit(fit, contrasts=contr.matrix)
  tfit <- treat(vfit, lfc=0)
  return(tfit)
}

apply_voomwithqualityweight_block <- function(y,design,contr.matrix,cond){
  voomWithQualityWeights(y, design = design, plot=TRUE,var.group = cond) -> y
  fit <- lmFit(y, design)
  vfit <- contrasts.fit(fit, contrasts=contr.matrix)
  tfit <- treat(vfit, lfc=0)
  return(tfit)
}


apply_voom <- function(y, design,contr.matrix){
  voom(y,design = design,plot = TRUE) -> y
  fit <- lmFit(y, design)
  vfit <- contrasts.fit(fit, contrasts=contr.matrix)
  tfit <- treat(vfit, lfc=0)
  return(tfit)
}
  
apply_glmQLFit <- function(y,design,contr.matrix){
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  res <- list()
  for (i in 1:ncol(contr.matrix)){
    res[[i]]<- glmQLFTest(fit, contrast = contr.matrix[,i])
  }
  return(res)
}
  
apply_glmLRT <- function(y,design,contr.matrix){
  y <- estimateDisp(y, design)
  fit <- glmFit(y, design)
  res <- list()
  for (i in 1:ncol(contr.matrix)){
    res[[i]]<- glmLRT(fit, contrast = contr.matrix[,i])
  }
  return(res)
}



