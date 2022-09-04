pbmc_pb <- readRDS("../../data/pbmc_pb.rds")
k="CD56dim CD16+ NK"
pbmc_pb -> pb_tmp

#remove sample (A_2, J_18 for cell number below 85.)
pb_tmp <- pb_tmp[,!pb_tmp$condition=="Severe"]


!colSums(assay(pb_tmp,k))==0 -> keep
assay(pb_tmp,k)[,keep] -> counts_k
counts_k[rowSums(counts_k)>50,]-> counts_k

colSums(counts_k)[order(colSums(counts_k))]
DGEList(counts_k) -> y
y <- calcNormFactors(y)
y$samples$group <- pb_tmp$group_id[keep]

cd <- y$samples$group
model.matrix(~0+cd) -> design

contr.matrix <- makeContrasts(
  AvsH = cdAsymptomatic - cdHC,
  MvsH = cdModerate - cdHC,
  MvsA = cdModerate- cdAsymptomatic,
  levels = colnames(design))

contr.matrix

cols=y$samples$group
table(cols)
cols[cols=="Asymptomatic"]<- "#1B9E77"
cols[cols=="HC"]<- "#7570B3"
cols[cols=="Moderate"]<- "#E7298A"
cols[cols=="Severe"]<- "#E7298A"

par(mfrow=c(1,1))
limma::plotMDS(y,col=cols)
source("../functions/DE_methods.R")
#voom
apply_voom(y,design,contr.matrix) -> voom_result
decideTests(voom_result) -> voom_de
topTreat(voom_result,coef = 1,n=Inf) -> t_voom

#voomwithqualityweight
apply_voomwithqualityweight(y,design,contr.matrix) -> vwq_result
decideTests(vwq_result) -> vwq_de
topTreat(vwq_result,coef = 1,n=Inf) -> t_vwq

#voomwithqualityweight_block
apply_voomwithqualityweight_block(y,design,contr.matrix,cond=cd) -> vwq_b_result
decideTests(vwq_b_result) -> vwq_b_de
topTreat(vwq_b_result,coef = 1,n=Inf) -> t_vwq_b


#voombygroup
apply_voombygroup(y,design=design,cd,contr.matrix,dynamic=NULL) -> vbg_result
decideTests(vbg_result) -> vbg_de
topTreat(vbg_result,coef = 1,n=Inf) -> t_vbg


#lrt
apply_glmLRT(y,design,contr.matrix) -> lrt_result
lapply(lrt_result, 
       function(x){decideTestsDGE(x)}) -> de_lrt 
topTags(lrt_result[[1]],n=Inf) -> t_lrt

#qlfit
apply_glmQLFit(y,design, contr.matrix) -> qlfit_result
lapply(qlfit_result, 
       function(x){decideTestsDGE(x)}) -> de_ql 
topTags(qlfit_result[[1]],n=Inf) -> t_ql

# combine results for comparisons
de_all <- data.frame(
  AvsH=c(voom_de[,1],vwq_de[,1],vwq_b_de[,1],vbg_de[,1],de_lrt[[1]],de_ql[[1]]),
  MvsH=c(voom_de[,2],vwq_de[,2],vwq_b_de[,2],vbg_de[,2],de_lrt[[2]],de_ql[[2]]),
  MvsA=c(voom_de[,3],vwq_de[,3],vwq_b_de[,3],vbg_de[,3],de_lrt[[3]],de_ql[[3]]),
  method=rep(c("voom","vwq","vwq_block","vbg","lrt","qlfit"),each=nrow(voom_result))
)
bind_rows(de_all %>% dplyr::group_by(method,AvsH) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(design="AvsH",sig=AvsH),
          de_all %>% dplyr::group_by(method,MvsH) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(design="MvsH",sig=MvsH),
          de_all %>% dplyr::group_by(method,MvsA) %>% dplyr::summarise(n=n()) %>% dplyr::mutate(design="MvsA",sig=MvsA)) %>% 
  dplyr::select(method,design,n,sig)-> de_summary
as.character(de_summary$sig) -> de_summary$sig

de_summary[!de_summary$sig=="0",] -> de_tmp
de_tmp[de_tmp$design=="AvsH",] -> de_tmp

de_tmp$method[de_tmp$method=="vbg"] <- "voomByGroup"
de_tmp$method[de_tmp$method=="vwq"] <- "VQW"
de_tmp$method[de_tmp$method=="vwq_block"] <- "VQWB"
de_tmp$method[de_tmp$method=="lrt"] <- "edgeR LRT"
de_tmp$method[de_tmp$method=="qlfit"] <- "edgeR QL"

factor(de_tmp$method,levels = c("edgeR LRT","edgeR QL","voom","VQW","VQWB","voomByGroup")) -> de_tmp$method
de_tmp$sig[de_tmp$sig=="1"] <- "Up"
de_tmp$sig[de_tmp$sig=="-1"] <- "Down"


#bar plot of DE genes
ggplot(de_tmp[!de_tmp$method=="VQW",])+geom_col(aes(x=method,y=n,fill=sig),position = "dodge") +
  facet_grid(.~design) +theme_bw() + 
  scale_fill_brewer(palette = "Paired") + labs(x="",y="Number of DE genes",fill="") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#venn plot of DE genes
vennDiagram(cbind(de_ql[[1]],voom_de[,1],vwq_b_de[,1],vbg_de[,1]),
            circle.col=c("turquoise", "salmon","orange","green"))

# GO analysis 
library(org.Hs.eg.db)
entrez.ids <- mapIds(org.Hs.eg.db, keys=rownames(vbg_de), 
                     column="ENTREZID", keytype="SYMBOL")

get_go <- function(genes) {
  is.de <- rownames(vbg_de) %in% genes
  go.out <- goana(unique(entrez.ids[is.de]), species="Hs", 
                  universe=unique(entrez.ids))
  # Only keeping biological process terms that are not overly general.
  go.out <- go.out[order(go.out$P.DE),]
  go.useful <- go.out[go.out$Ont=="BP" & go.out$N <= 100& go.out$P.DE<0.01,]
  return(go.useful)
}

rownames(vbg_de)[vbg_de[,1]==1]-> vbg_genes
get_go(vbg_genes) -> vbg_go
vbg_go[grepl("interfer",vbg_go$Term),]
head(vbg_go)


#find genes in GO
library(org.Hs.eg.db)
library(GO.db)
gamma_go_id="GO:0060333"


vbg_gamma <- c()
for (id in gamma_go_id) {
  allegs = get(id, org.Hs.egGO2ALLEGS)
  genes = unlist(mget(allegs,org.Hs.egSYMBOL))
  vbg_genes[vbg_genes%in% genes] -> vbg_tmp
  c(vbg_tmp,vbg_gamma) -> vbg_gamma
}
vbg_gamma<- unique(vbg_gamma)
vbg_gamma



y_tmp <- y[,y$samples$group%in% c("HC","Asymptomatic")]
lcpm <- cpm(y_tmp,log = TRUE)
library(pheatmap)
pheatmap(lcpm[c(vbg_gamma),order(y_tmp$samples$group)],
         annotation_col = data.frame(group=y_tmp$samples$group,row.names = rownames(y_tmp$samples)),
         cluster_cols = FALSE,scale = "row",
         cluster_rows = FALSE)



vbg_go[grepl("interfer",vbg_go$Term),] ->GOtable

GOtable$Term = factor(GOtable$Term,levels = rev(GOtable$Term))
ggplot(data=GOtable,aes(x=Term,y=-log10(P.DE),fill=-log10(P.DE)))+
  geom_bar(stat="identity",show.legend=F)+
  scale_fill_continuous(type = "viridis")+
  labs(y="-log10(P value)",x="GO Biological Processes")+
  theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  coord_flip()



#inferon gamma genes pval compare-wise comparisons
data.frame(vbg=t_vbg$adj.P.Val[match(rownames(t_vbg),rownames(t_vbg))],
           voom=t_voom$adj.P.Val[match(rownames(t_vbg),rownames(t_voom))],
           vwq_b=t_vwq_b$adj.P.Val[match(rownames(t_vbg),rownames(t_vwq_b))],
           lrt=t_lrt$table$FDR[match(rownames(t_vbg),rownames(t_lrt))],
           ql=t_ql$table$FDR[match(rownames(t_vbg),rownames(t_ql))]) -> df

method <-c ("vbg","voom","voomQWB","edgeR LRT","edgeR QL")
par(mfrow=c(2,3))
for (i in 2:5){
  plot(x=df[,1],y=df[,i],xlab="Adjusted P values (voomByGroup)", ylab="Adjusted P values",main=method[i],col="grey")
  points(df_gamma[,1], df_gamma[,2], col='blue',cex=2)
  abline(0,1,lwd=2,col="red")
}

