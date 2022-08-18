setwd("/stornext/HPCScratch/home/you.y/simulation/ngroup_sim/bcv_change/")
#0.05
library(RColorBrewer) 
brewer.pal(6,"Paired") -> col
names(col) <- c("LL","Low vs Low","LH","Low vs High","HH","High vs High")


df_summary_all <- tibble()
de_list <- list(c("de1","de2"),
                c("de1","de3"),
                c("de1","de4"),
                c("de2","de3"),
                c("de2","de4"),
                c("de3","de4"))
unique(df_all$contr) -> uct
for (ct in 1:length(uct)){
  df_all %>% dplyr::filter(FDR<0.05) %>% dplyr::filter(contr==uct[ct]) -> df_sub
  df_sub$true_positive_n <- 0
  df_sub$true_positive_n[df_sub$de %in% de_list[[ct]]] <- 1
  
  df_sub$false_positive_n <- 1
  df_sub$false_positive_n[df_sub$de %in% de_list[[ct]]] <- 0
  df_sub %>% group_by(method,rep) %>% 
    summarise(tpn=sum(true_positive_n),n=n(),
              fpn=sum(false_positive_n)) -> df_summary
  df_summary$TPR <- df_summary$tpn/100
  
  df_summary$p <- df_summary$fpn+df_summary$tpn
  df_summary$FDR <- df_summary$fpn/df_summary$p
  df_summary$design <- uct[ct]
  bind_rows(df_summary,df_summary_all) -> df_summary_all
}
table(df_summary_all$design)
df_summary_all %>% group_by(method,design) %>% summarise(n=sum(p),FDR=mean(FDR)) -> df_tmp

df_tmp$design[df_tmp$design=="group1vsgroup2"] <- "Low vs Low"
df_tmp$design[df_tmp$design=="group1vsgroup4"] <- "Low vs High"
df_tmp$design[df_tmp$design=="group1vsgroup3"] <- "Low vs High"
df_tmp$design[df_tmp$design=="group2vsgroup3"] <- "Low vs High"
df_tmp$design[df_tmp$design=="group2vsgroup4"] <- "Low vs High"
df_tmp$design[df_tmp$design=="group3vsgroup4"] <- "High vs High"


library(ggsci)
ggplot(df_tmp,aes(x=method,y=FDR,col=design)) + 
  geom_point(size=3.2,shape=1,stroke=1) +
  geom_hline(yintercept = 0.05,linetype="dashed")+
  theme_bw()+
  scale_color_manual(values=col)+
  labs(y="FDR (alpha=0.05)",x="",col="Comparison")+theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


