library(DHARMa)
library(igraph)
library(dplyr)
library(here)
library(ggplot2)
library(ggpubr)
library(glmmTMB)
library(sjPlot)

source(here("functions", "katz_function.R")) #Loading function to calculate Katz centrality

data<-read.table(here("data", "Cosmo_et_al_fulldataset.txt"), header=TRUE)
net<-read.table(here("data", "Cosmo_et_al_network.txt"), header=TRUE) #Importing data
net_names<-net[,1] #Subsetting names
net_data<-as.matrix(net[, 2:19]) #Subsetting only interaction matrix
rownames(net_data)<-net_names

net_weighted<-graph_from_incidence_matrix(net_data, weighted=TRUE)
net_binary<-graph_from_incidence_matrix(net_data, weighted=NULL)

adj.weight<-get.adjacency(net_weighted, attr="weight", sparse=FALSE)
adj.binary<-get.adjacency(net_binary, sparse=FALSE)

#Binary and weighted degree centrality

deg_bin<-rowSums(adj.binary)[1:123]
deg_w<-rowSums(adj.weight)[1:123]

#Closeness centrality

#Shortest paths

spaths_bin<-1/distances(net_binary)
diag(spaths_bin)<-0
spaths_w<-1/distances(net_weighted, weights=(1/E(net_weighted)$weight))
diag(spaths_w)<-0

#Closeness centrality measures

closeness_bin<-rowSums(spaths_bin)[1:123]
closeness_w<-rowSums(spaths_w)[1:123]

#Katz centrality

katz_bin<-katz(adj.binary)[1:123]
katz_w<-katz(adj.weight)[1:123]

#Statistical models

net_df<-data.frame(sample=rownames(net_data), deg_bin, deg_w, closeness_bin, closeness_w, katz_bin, katz_w)
net_df<-left_join(net_df, data)

#net_df$deg_bin<-scale(net_df$deg_bin)
#net_df$deg_w<-scale(net_df$deg_w)
#net_df$closeness_bin<-scale(net_df$closeness_bin)
#net_df$closeness_w<-scale(net_df$closeness_w)
#net_df$katz_bin<-scale(net_df$katz_bin)
#net_df$katz_w<-scale(net_df$katz_w)

#Degree centrality

model_deg<-glmmTMB(deg_bin~scale(among_structural_pd)+scale(within_structural_pd)+scale(within_comp_PC1)+
                         scale(among_comp_PC1)+scale(within_comp_PC2)+scale(among_comp_PC2)+
                         scale(plant_height_cm)+scale(nleaves)+stratum+season+(1|plant), data=net_df)

model_degw<-glmmTMB(deg_w~scale(among_structural_pd)+scale(within_structural_pd)+scale(within_comp_PC1)+
                     scale(among_comp_PC1)+scale(within_comp_PC2)+scale(among_comp_PC2)+
                     scale(plant_height_cm)+scale(nleaves)+stratum+season+(1|plant), data=net_df)
#Closeness centrality

model_close<-glmmTMB(closeness_bin~scale(among_structural_pd)+scale(within_structural_pd)+scale(within_comp_PC1)+
                     scale(among_comp_PC1)+scale(within_comp_PC2)+scale(among_comp_PC2)+
                     scale(plant_height_cm)+scale(nleaves)+stratum+season+(1|plant), data=net_df)

model_closew<-glmmTMB(closeness_w~scale(among_structural_pd)+scale(within_structural_pd)+scale(within_comp_PC1)+
                     scale(among_comp_PC1)+scale(within_comp_PC2)+scale(among_comp_PC2)+
                     scale(plant_height_cm)+scale(nleaves)+stratum+season+(1|plant), data=net_df)

#Katz centrality

model_katz<-glmmTMB(katz_bin~scale(among_structural_pd)+scale(within_structural_pd)+scale(within_comp_PC1)+
                     scale(among_comp_PC1)+scale(within_comp_PC2)+scale(among_comp_PC2)+
                     scale(plant_height_cm)+scale(nleaves)+stratum+season+(1|plant), data=net_df)

model_katzw<-glmmTMB(katz_w~scale(among_structural_pd)+scale(within_structural_pd)+scale(within_comp_PC1)+
                     scale(among_comp_PC1)+scale(within_comp_PC2)+scale(among_comp_PC2)+
                     scale(plant_height_cm)+scale(nleaves)+stratum+season+(1|plant), data=net_df)

#Getting model confidence intervals and coefficient estimates

ci_list<-list(
ci_deg<-confint(model_deg),
ci_degw<-confint(model_degw),
ci_close<-confint(model_close),
ci_closew<-confint(model_closew),
ci_katz<-confint(model_katz),
ci_katzw<-confint(model_katzw)
)

ci_df<-lapply(ci_list, FUN=function(x){
  x<-as.data.frame(x[2:7, ])
  df<-data.frame(variable=rownames(x), lowerci=x$`2.5 %`, upperci=x$`97.5 %`, coef=x$Estimate)
  
})

ci_df<-do.call(rbind, ci_df)

variable<-rep(c("Among-plant structural PD", "Within-plant structural PD", "Within-plant compositional PD (PC1)", 
"Among-plant compositional PD (PC1)", "Within-plant compositional PD (PC2)", 
"Among-plant compositional PD (PC2)"), 6)

models<-rep(c("Degree", 
              "Closeness centrality", "Katz centrality"),
              each=12)

sit<-rep(rep(c("Binary network", 
               "Weighted network"), each=6), 3)

ci_df$variable<-variable
ci_df$models<-models
ci_df$sit<-sit

ci_df$sit<-as.factor(ci_df$sit)
ci_df$models<-as.factor(ci_df$models)

ci_df_r <- ci_df %>%
  mutate(model_r = factor(models, levels=c("Degree", "Closeness centrality", "Katz centrality")))

p_ci<-ggplot(ci_df_r, aes(reorder(variable, -coef), coef)) +
  geom_linerange(
    aes(ymin = lowerci, ymax = upperci, color=sit),
    position = position_dodge(0.5), size = 0.6)+
  geom_point(aes(color=sit, shape=model_r), position = position_dodge(0.5), size=3.0)+
  scale_color_manual(values=c("tan3", "royalblue3"))+
  facet_grid(sit~model_r, scales="free")

p_ci<-ggpar(p_ci, legend="none", font.tickslab=18, font.x=20, font.legend=14, xlab="", ylab="Coefficient", orientation=c("horizontal"), ggtheme=theme_bw())+
  geom_hline(yintercept=0, linetype="dashed", color = "black")+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), strip.text = element_text(size = 16))

p_ci

summary(model_degw)

#save as pdf, 7x16in

# Model tables

tab_model(model_deg, model_close, model_katz)
tab_model(model_degw, model_closew, model_katzw)


