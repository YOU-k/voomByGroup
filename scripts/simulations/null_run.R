library(limma)
library(edgeR)
library(SingleCellExperiment)

x <- read.table("../../data/itreg_pb.txt")
x <- as.matrix(x)

#get baseline prop
prop <- goodTuringProportions(x)
props1 <- list(c(1,1,1),
               c(1,1,1),
               c(1,1,1),
               c(1,1,1)) 

props <- props1

#get ncells from loop
nx=250
n1=round(sum(props[[1]]*nx))
n2=round(sum(props[[2]]*nx))
n3=round(sum(props[[3]]*nx))
n4=round(sum(props[[4]]*nx))

NCells=n1+n2+n3+n4

#read in reference data
readRDS("../../data/sce_tmp.rds") -> sce_tmp

#sim exp lib size.
set.seed(666)
source("../functions/sim_exp_lib.R")
sim_exp_lib(sce = sce_tmp,nCells = NCells) -> exp.lib.sizes
round(exp.lib.sizes) -> libsize

library(tidyverse)
df_all <- tibble()
for (r in 1:50) {
  nx=250
  n1=round(sum(props[[1]]*nx))
  n2=round(sum(props[[2]]*nx))
  n3=round(sum(props[[3]]*nx))
  n4=round(sum(props[[4]]*nx))
  
  
  NCells=n1+n2+n3+n4
  subject <- c(rep(1:3,props[[1]]*nx),rep(4:6,props[[2]]*nx),
               rep(7:9,props[[3]]*nx),rep(10:12,props[[4]]*nx))
  #define groups
  group <- factor(c(rep(1,n1),rep(2,n2),rep(3,n3),rep(4,n4)))
  table(subject,group)
  
  NGenes=10000
  
  is.expr <- rowSums(prop>1e-6)>=6
  prop <- prop[is.expr,]
  i <- seq(from=0,to=1,length=length(prop))
  f <- approxfun(i,sort(prop),rule=2)
  baselineprop <- f( (1:NGenes)/(NGenes+1) )
  baselineprop <- baselineprop/sum(baselineprop)
  
  i <- sample(1:NGenes,NCells)
  i1 <- i[1:50] #seprate two groups for simulate groups factor in further analysis #
  i2 <- i[51:100]
  i3 <- i[101:150]
  i4 <- i[151:200]
  fc <- 2
  
  ide_tb <- data.frame(ide=i,de=c(rep(c("de1","de2","de3","de4"),each=50),rep("no",(NCells-200))))
  baselineprop1 <- baselineprop2 <-  baselineprop3 <- baselineprop4 <- baselineprop #calculate baseline propration for each gene
  baselineprop1[i1] <- baselineprop1[i1]*fc #add fc as fold changes# 
  baselineprop2[i2] <- baselineprop2[i2]*fc #same as above#
  baselineprop3[i3] <- baselineprop3[i3]*fc
  baselineprop4[i4] <- baselineprop4[i4]*fc
  
  mu0.1 <- matrix(baselineprop1,NGenes,1) %*% matrix(libsize[1:n1],1,n1)
  mu0.2 <- matrix(baselineprop2,NGenes,1) %*% matrix(libsize[(n1+1):(n1+n2)],1,n2)
  mu0.3 <- matrix(baselineprop3,NGenes,1) %*% matrix(libsize[(n1+n2+1):(n1+n2+n3)],1,n3)
  mu0.4 <- matrix(baselineprop4,NGenes,1) %*% matrix(libsize[(n1+n2+n3+1):(n1+n2+n3+n4)],1,n4)
  
  cbind(mu0.1,mu0.2,mu0.3,mu0.4) -> mu
  
  
  source("../functions/simscRNAseq_group4.R")
  y <- simscRNAseq_group4(mu,phi=c(0.4,0.4,0.4,0.4),subject,correlation=0.1,group=group)
  
  counts2 <- y

  
  #create single cell object
  SingleCellExperiment(list(counts=counts2)) -> sce
  sce$cluster_id <- "iTreg"
  sce$sample_id <- subject
  sce$group_id <- group
  
  
  #aggregate to pseudobulk
  library(muscat)
  #sce <- sce[,colSums(counts(sce)>0) >50]
  pb <- aggregateData(sce)
  formula <- ~ 0+ group_id
  cd <- as.data.frame(colData(pb))
  design <- model.matrix(formula, cd)
  
  design
  
  contr.matrix <- makeContrasts(
    group1vsgroup2= group_id1- group_id2,
    group1vsgroup3= group_id1- group_id3,
    group1vsgroup4= group_id1- group_id4,
    group2vsgroup3= group_id2- group_id3,
    group2vsgroup4= group_id2- group_id4,
    group3vsgroup4= group_id3- group_id4,
    
    
    levels=colnames(design)
  )
  contr.matrix
  
  DGEList(assay(pb,"iTreg"), remove.zeros = FALSE) -> y
  cond <- pb$group_id[!colSums(y$counts)==0]
  y <- y[rowSums(y$counts)>30,]
  #keep.exprs <- filterByExpr(y, group=pb$group_id)
  #y <- y[keep.exprs,, keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  
  plotMDS(y,col=c(rep("black",3), rep("red",3),rep("green",3),rep("blue",3)), 
          labels= c(rep("Grp1",3), rep("Grp2",3),rep("Grp3",3),rep("Grp4",3)))
  
  
  
  #estimateDisp(y)
  #design <- model.matrix(formula, cd1)
  source("../functions/DE_methods.R")
  pdf("vbg_bcv.pdf",width = 4,height = 4)
  apply_voombygroup(y,design = design,cond = cond,contr.matrix = contr.matrix,dynamic = NULL) -> fit_vbg
  dev.off()
  
  apply_voomwithqualityweight(y,design = design,contr.matrix = contr.matrix) -> fit_vwq
  apply_voomwithqualityweight_block(y,design = design,contr.matrix = contr.matrix,cond=cond) -> fit_vwq_block
  apply_voom(y,design = design,contr.matrix = contr.matrix) -> fit_voom
  apply_glmLRT(y,design = design,contr.matrix = contr.matrix) -> res_lrt
  apply_glmQLFit(y,design = design,contr.matrix = contr.matrix) -> res_qlf
  
  
  for (d in 1:ncol(contr.matrix)){
    lapply(list(fit_vbg,fit_vwq,fit_vwq_block,fit_voom), 
           function(x){topTreat(x, coef=d, n=Inf)}) -> treat_tables
    lapply(treat_tables, 
           function(x){as.data.frame(x) %>% mutate(gene=rownames(x)) %>%
               mutate(de=ide_tb$de[match(gene,ide_tb$ide)])}) -> tmp_t
    bind_rows(tmp_t)-> voom_tmp
    voom_tmp$method <- rep(c("vbg","vwq","vwq_block","voom"),each=nrow(y))
    #edger methods
    res_lrt[[d]] -> res_lrt_d
    res_qlf[[d]] -> res_qlf_d
    lapply(list(res_lrt_d,res_qlf_d), 
           function(x){topTags(x,n=Inf)}) -> tag_tables 
    lapply(tag_tables, 
           function(x){as.data.frame(x) %>% mutate(gene=rownames(x)) %>%
               mutate(de=ide_tb$de[match(gene,ide_tb$ide)])}) -> tmp_e
    bind_rows(tmp_e)-> edger_tmp
    edger_tmp$method <- rep(c("lrt","qlf"),each=nrow(y))
    
    bind_rows(voom_tmp %>% mutate(Pvalue=P.Value,FDR=adj.P.Val),
              edger_tmp)[,c("gene","de","logFC","Pvalue","FDR","method")] -> tmp_all
    tmp_all$rep = r
    tmp_all$contr <- colnames(contr.matrix)[d]
    rbind(df_all,tmp_all) -> df_all
  }
  
}
#saveRDS(df_all,"../../data/df_all.rds")
