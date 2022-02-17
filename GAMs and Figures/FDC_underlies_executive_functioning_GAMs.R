#This script is investigating whether white matter microstructure (FDC) underlies improvement in executive functioning
#using Generelized Additive Models (GAMs). We used 2 factors for to describe executive functioning performance: Executive 
#Efficiency and Executive and Complex Reasoning Accuracy. All GAMs include sex, mean DWI framewise displacement and number of 
#DWI bad slices (index of scan quality) as covariates. Included at the end of the script are sensitivity analyses in which we 
#controlled for Total Brain Volume. 

######################
#### READ IN DATA ####
######################
source("/cbica/projects/pnc_fixel_cs/GAMs/scripts/FDC_development_GAMs.R")

########################
#### VISUALIZE DATA ####  ##Can skip this section##
########################
smooth_plot_EE <- function(yvar){
  require(ggplot2)
  ggplot(df_fdc, aes_(x=~F3_Executive_Efficiency,y=as.name(yvar))) + geom_point() + geom_smooth(method = "gam", formula = y ~s(x))
}
plot_EE_list<-lapply(names(df_fdc[c(2:15)]),smooth_plot_EE)
plot_EE<-plot_grid(plotlist=plot_EE_list)
plot_EE

#######################################################
#### BUILD GAM MODELS - EXECUTIVE FUNCTION EFFECTS ####
#######################################################
#1.Executive Efficiency 
gamModels_EE <- lapply(Components, function(x) {
  gam(substitute(i ~ s(Age) + F3_Executive_Efficiency + oSex + raw_num_bad_slices + mean_fd, list(i = as.name(x))), method="REML", data = df_fdc)
})

#2.Executive & Complex Reasoning Accuracy
gamModels_ECRA <- lapply(Components, function(x) {
  gam(substitute(i ~ s(Age) + F1_Exec_Comp_Res_Accuracy + oSex + raw_num_bad_slices + mean_fd, list(i = as.name(x))), method="REML", data = df_fdc)
})

#Look at model summaries
models_EE <- lapply(gamModels_EE, summary)
models_ECRA<-lapply(gamModels_ECRA, summary)
models_EE[[1]]

#Pull t-values and p-values
cog_tvalue_EE <- sapply(gamModels_EE, function(v) summary(v)$p.table[2,3])
cog_tvalue_EE <- as.data.frame(cog_tvalue_EE)
cog_pvalue_EE <- sapply(gamModels_EE, function(v) summary(v)$p.table[2,4])
cog_pvalue_EE <- as.data.frame(cog_pvalue_EE)
cog_pvalue_EE <- round(cog_pvalue_EE,3)
cog_pvalue_EE_fdr <- as.data.frame(p.adjust(cog_pvalue_EE[,1], method="fdr"))
cog_pvalue_EE_fdr <- round(cog_pvalue_EE_fdr,3)

cog_tvalue_ECRA <- sapply(gamModels_ECRA, function(v) summary(v)$p.table[2,3])
cog_tvalue_ECRA <- as.data.frame(cog_tvalue_ECRA)
cog_pvalue_ECRA <- sapply(gamModels_ECRA, function(v) summary(v)$p.table[2,4])
cog_pvalue_ECRA <- as.data.frame(cog_pvalue_ECRA)
cog_pvalue_ECRA <- round(cog_pvalue_ECRA,3)
cog_pvalue_ECRA_fdr <- as.data.frame(p.adjust(cog_pvalue_ECRA[,1], method="fdr"))
cog_pvalue_ECRA_fdr <- round(cog_pvalue_ECRA_fdr,3)

cog_stats<-cbind(bundles,cog_tvalue_EE,cog_pvalue_EE_fdr,cog_tvalue_ECRA,cog_pvalue_ECRA_fdr)
cog_stats<-cog_stats%>%
  dplyr::rename(t_EE=cog_tvalue_EE)%>%
  dplyr::rename(p_fdr_EE=`p.adjust(cog_pvalue_EE[, 1], method = "fdr")`) %>%
  dplyr::rename(t_ECRA=cog_tvalue_ECRA)%>%
  dplyr::rename(p_fdr_ECRA=`p.adjust(cog_pvalue_ECRA[, 1], method = "fdr")`)

#Pull Partial R2
df_fdc_cog_min<-df_fdc %>% #no missing data here on any variable, thus can use df_fdc_min data for reduced model testing
  dplyr::select(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,Age,oSex,mean_fd,raw_num_bad_slices,TBV,F1_Exec_Comp_Res_Accuracy,
                F3_Executive_Efficiency)
count(df_fdc_cog_min[rowSums(is.na(df_fdc_cog_min))==0,])

redmodel_cog <- lapply(Components, function(x) {
  gam(substitute(i ~ s(Age) + oSex + raw_num_bad_slices + mean_fd, list(i = as.name(x))), method="REML", data = df_fdc_cog_min)
})

partialR2_EE_V1 <- partialRsq(gamModels_EE[[1]],redmodel_cog[[1]])
partialR2_EE_V2 <- partialRsq(gamModels_EE[[2]],redmodel_cog[[2]])
partialR2_EE_V3 <- partialRsq(gamModels_EE[[3]],redmodel_cog[[3]])
partialR2_EE_V4 <- partialRsq(gamModels_EE[[4]],redmodel_cog[[4]])
partialR2_EE_V5 <- partialRsq(gamModels_EE[[5]],redmodel_cog[[5]])
partialR2_EE_V6 <- partialRsq(gamModels_EE[[6]],redmodel_cog[[6]])
partialR2_EE_V7 <- partialRsq(gamModels_EE[[7]],redmodel_cog[[7]])
partialR2_EE_V8 <- partialRsq(gamModels_EE[[8]],redmodel_cog[[8]])
partialR2_EE_V9 <- partialRsq(gamModels_EE[[9]],redmodel_cog[[9]])
partialR2_EE_V10 <- partialRsq(gamModels_EE[[10]],redmodel_cog[[10]])
partialR2_EE_V11 <- partialRsq(gamModels_EE[[11]],redmodel_cog[[11]])
partialR2_EE_V12 <- partialRsq(gamModels_EE[[12]],redmodel_cog[[12]])
partialR2_EE_V13 <- partialRsq(gamModels_EE[[13]],redmodel_cog[[13]])
partialR2_EE_V14 <- partialRsq(gamModels_EE[[14]],redmodel_cog[[14]])

partialR2_EE<-as.data.frame(cbind(partialR2_EE_V1[[1]],partialR2_EE_V2[[1]],partialR2_EE_V3[[1]],partialR2_EE_V4[[1]],partialR2_EE_V5[[1]],
                                   partialR2_EE_V6[[1]],partialR2_EE_V7[[1]],partialR2_EE_V8[[1]],partialR2_EE_V9[[1]],partialR2_EE_V10[[1]],
                                   partialR2_EE_V11[[1]],partialR2_EE_V12[[1]],partialR2_EE_V13[[1]],partialR2_EE_V14[[1]]))
partialR2_EE<-as.data.frame(t(partialR2_EE))
partialR2_EE<-partialR2_EE %>%
  dplyr::rename(partialR2_EE=V1)

partialR2_ECRA_V1 <- partialRsq(gamModels_ECRA[[1]],redmodel_cog[[1]])
partialR2_ECRA_V2 <- partialRsq(gamModels_ECRA[[2]],redmodel_cog[[2]])
partialR2_ECRA_V3 <- partialRsq(gamModels_ECRA[[3]],redmodel_cog[[3]])
partialR2_ECRA_V4 <- partialRsq(gamModels_ECRA[[4]],redmodel_cog[[4]])
partialR2_ECRA_V5 <- partialRsq(gamModels_ECRA[[5]],redmodel_cog[[5]])
partialR2_ECRA_V6 <- partialRsq(gamModels_ECRA[[6]],redmodel_cog[[6]])
partialR2_ECRA_V7 <- partialRsq(gamModels_ECRA[[7]],redmodel_cog[[7]])
partialR2_ECRA_V8 <- partialRsq(gamModels_ECRA[[8]],redmodel_cog[[8]])
partialR2_ECRA_V9 <- partialRsq(gamModels_ECRA[[9]],redmodel_cog[[9]])
partialR2_ECRA_V10 <- partialRsq(gamModels_ECRA[[10]],redmodel_cog[[10]])
partialR2_ECRA_V11 <- partialRsq(gamModels_ECRA[[11]],redmodel_cog[[11]])
partialR2_ECRA_V12 <- partialRsq(gamModels_ECRA[[12]],redmodel_cog[[12]])
partialR2_ECRA_V13 <- partialRsq(gamModels_ECRA[[13]],redmodel_cog[[13]])
partialR2_ECRA_V14 <- partialRsq(gamModels_ECRA[[14]],redmodel_cog[[14]])

partialR2_ECRA<-as.data.frame(cbind(partialR2_ECRA_V1[[1]],partialR2_ECRA_V2[[1]],partialR2_ECRA_V3[[1]],partialR2_ECRA_V4[[1]],partialR2_ECRA_V5[[1]],
                                  partialR2_ECRA_V6[[1]],partialR2_ECRA_V7[[1]],partialR2_ECRA_V8[[1]],partialR2_ECRA_V9[[1]],partialR2_ECRA_V10[[1]],
                                  partialR2_ECRA_V11[[1]],partialR2_ECRA_V12[[1]],partialR2_ECRA_V13[[1]],partialR2_ECRA_V14[[1]]))
partialR2_ECRA<-as.data.frame(t(partialR2_ECRA))
partialR2_ECRA<-partialR2_ECRA %>%
  dplyr::rename(partialR2_ECRA=V1)

#Join Partial R2 with t-stats and fdr p-values for both EE and ECRA
cog_stats<-cbind(cog_stats,partialR2_EE,partialR2_ECRA)
cog_stats$partialR2_EE<-as.numeric(as.character(cog_stats$partialR2_EE))
cog_stats<-cog_stats[order(-cog_stats$partialR2_EE),]

##############################################################################################
#### FIGURE 5B - BAR PLOT OF EXECUTIVE FUNCTIONING PARTIAL R2 FOR EACH COVARIANCE NETWORK ####
##############################################################################################
#Considering very similar results between EE and ECRA, in the paper we only presents results pertaining to the EE analysis (in Figure 5)
library(ggplot2)
cog_stats$bundles <- factor(cog_stats$bundles, levels = cog_stats$bundles)
figure5b<-ggplot(cog_stats, aes(x=bundles, y=partialR2_EE, fill=bundles)) + geom_bar(stat="identity") + scale_fill_manual("Processing Method", values = c("CC (middle)" = "#F77F85", "SLF" = "#008B45", "Splenium" = "#EE3B3B", "Fornix" = "#4169E1", "Sup. CST" = "#9CB9F5", "Inf. CST" = "#6E91EB", "U-fibers" = "#9CCBB3", "Rostrum" = "#FFB6C1", "SLF (parietal)" = "#D0E0D8", "Uncinate" = "#4EAB7C", "Sup. Cerebellum" = "#F7EAAA", "Vermix" = "#F3D370", "Int. capsule" = "#CAE1FF", "Middle CP" = "#EEB422")) + labs(x ="Covariance Networks", y="Partial R2 for Executive Efficiency") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(face="bold",size=18, angle = 50, hjust=1), axis.text.y = element_text(face="bold",size=20), axis.title = element_text(size=28), axis.line= element_line(colour = 'black', size = 1.5), legend.position = 'none', legend.text = element_text(size=20), legend.title = element_text(size = 24), plot.title = element_text(face="bold",size = 20)) + ggtitle("Association of executive efficiency with white matter microstructure covariance networks") +scale_y_continuous(breaks=c(0,0.02,0.04,0.06)) +
  theme(plot.title = element_text(hjust = 0.5))
figure5b
ggsave(plot = figure5b,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure5B_bargraph_partialR2_EE.png",device = "png",width = 200,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


########################################################
#### COVARIANCE NETWORKS PREDICT EXECUTIVE FUNCTION ####
########################################################
#Full model: EE or ECRA ~ Age + sex + covariance networks + father education + mother education + ethnicity
#Null model: EE or ECRA ~ Age + sex + father education + mother education + ethnicity
#Ftest to test the significant contribution of covariance networks

#Add demographic covariates
colnames(demo)
demo<-demo%>%
  dplyr::select(bblid,race,race2, ethnicity,fedu1,medu1)
df_fdc<-left_join(df_fdc,demo,by="bblid")

#Check and remove missings
sum(is.na(df_fdc$medu1))
df_fdc_demo<-df_fdc %>%
  filter(!is.na(medu1.y)) %>%
  filter(!is.na(fedu1.y)) %>%
  filter(!is.na(ethnicity.x))

#1.Full and Null Models for Executive Efficiency
FullModel_pred_EE<-lm(F3_Executive_Efficiency ~ Age + oSex + V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14 + ethnicity.x +fedu1.y+medu1.y,
                   data=df_fdc_demo)

NullModel_pred_EE<-lm(F3_Executive_Efficiency ~  Age + oSex + ethnicity.x + fedu1.y+medu1.y,data=df_fdc_demo)
ols_vif_tol(NullModel_pred_EE)

NMF_only_predict_EE<-lm(F3_Executive_Efficiency ~  V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14,data=df_fdc)

Full_pred_EE_summary<-summary(FullModel_pred_EE)
Null_pred_EE_summary<-summary(NullModel_pred_EE)
NMF_only_predict_EE_summary<-summary(NMF_only_predict_EE)
Full_pred_EE_summary
Null_pred_EE_summary

#1a.Compare Fullmodel and Nullmodel with an F test to test significant contribution of networks in predicting EE
Ftest_EE<-anova(FullModel_pred_EE,NullModel_pred_EE)
Ftest_EE # significant contribution of the covariance networks in predicting EE

#1b.Correlation between actual EE and predicted EE by NMF_only_predict_EE model
NMF_EE_Cor <- cor.test(predict(NMF_only_predict_EE), df_fdc$F3_Executive_Efficiency)
NMF_EE_Cor #the r value will be added in Figure 5C

#2.Full and Null Models for Executive & Complex Reasoning Accuracy
FullModel_pred_ECRA<-lm(F1_Exec_Comp_Res_Accuracy ~ Age + oSex + V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14 + ethnicity.x +fedu1.y+medu1.y,
                      data=df_fdc_demo)

NullModel_pred_ECRA<-lm(F1_Exec_Comp_Res_Accuracy ~  Age + oSex + ethnicity.x + fedu1.y+medu1.y,data=df_fdc_demo)

NMF_only_predict_ECRA<-lm(F1_Exec_Comp_Res_Accuracy ~  V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14,data=df_fdc)

Full_pred_ECRA_summary<-summary(FullModel_pred_ECRA)
Null_pred_ECRA_summary<-summary(NullModel_pred_ECRA)
NMF_only_predict_ECRA_summary<-summary(NMF_only_predict_ECRA)
Full_pred_ECRA_summary
Null_pred_ECRA_summary

#2a.Compare Fullmodel and Nullmodel with an F test to test significant contribution of networks in predicting ECRA
Ftest_ECRA<-anova(FullModel_pred_ECRA,NullModel_pred_ECRA)
Ftest_ECRA # significant contribution of the covariance networks in predicting EE

#2b.Correlation between actual EE and predicted EE by NMF_only_predict_EE model
NMF_ECRA_Cor <- cor.test(predict(NMF_only_predict_ECRA), df_fdc$F1_Exec_Comp_Res_Accuracy)
NMF_ECRA_Cor 


#########################################################################
#### FIGURE 5C - CORRELATION PLOT BETWEEN PREDICTED EE AND ACTUAL EE ####
#########################################################################
png(file='/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/Figures/Figure5C_NMF_predicted_vs_actual_EE.png',width=2000,height=2000,res=300)
par(mar=c(5,5,2,2), cex.axis=2, bty="l")
plot(predict(NMF_only_predict_EE), df_fdc$F3_Executive_Efficiency, ylab="Executive Efficiency (z)", xlab="Predicted Executive Efficiency", cex.lab=2, pch=19, col="darkorange1")
abline(a=0,b=1, col = "darkorange3", lwd = 4)
dev.off()


#######################################################################
#### FIGURE 5D - ASSOCIATION OF EE WITH FDC OF COVARIANCE NETWORKS ####
#######################################################################
#Create geom_smooth plot to show the above association for the three (3) most important covariance networks (based on partial R2)
#This part of the code was also taken from Adam Pines' 
#(https://github.com/PennBBL/multishell_diffusion/blob/master/PostProc/multishell_analyses.Rmd)

#Association between EE and Fornix (i.e.,V2)
EE_V2<-ggplot(df_fdc,aes(x=F3_Executive_Efficiency,y=V2)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'lm', formula = y~x, colour=('#4169E1'), fill = "#4169E1", alpha = .8) + xlab("Executive Efficiency (z)") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
EE_V2<- EE_V2+theme(text = element_text(size=this_font_size),
                    axis.text = element_text(size = this_font_size),
                    axis.title.y = element_text(size = this_font_size),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    legend.text = element_text(size = this_font_size),
                    legend.title = element_text(size = this_font_size),
                    axis.title = element_text(size = this_font_size),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "transparent",colour = NA),
                    plot.background = element_rect(fill = "transparent",colour = NA),
                    plot.margin = unit(c(0.2,0.2,0,0.2), "cm")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r = 39, b = 0, l = 0),size = 15))  #Top, left,Bottom, right
EE_V2
# save
ggsave(plot = EE_V2,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure5D_EE_V2.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)

#Association between EE and SLF (parietal) (i.e.,V10)
EE_V10<-ggplot(df_fdc,aes(x=F3_Executive_Efficiency,y=V10)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'lm', formula = y~x, colour=('#D0E0D8'), fill = "#D0E0D8", alpha = .8) + xlab("Executive Efficiency (z)") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
EE_V10<- EE_V10+theme(text = element_text(size=this_font_size),
                    axis.text = element_text(size = this_font_size),
                    axis.title.y = element_text(size = this_font_size),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    legend.text = element_text(size = this_font_size),
                    legend.title = element_text(size = this_font_size),
                    axis.title = element_text(size = this_font_size),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "transparent",colour = NA),
                    plot.background = element_rect(fill = "transparent",colour = NA),
                    plot.margin = unit(c(0.2,0.2,0,0.2), "cm")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r = 39, b = 0, l = 0),size = 15))  #Top, left,Bottom, right
EE_V10
# save
ggsave(plot = EE_V10,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure5D_EE_V10.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)

#Association between EE and SLF (i.e.,V5)
EE_V5<-ggplot(df_fdc,aes(x=F3_Executive_Efficiency,y=V5)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'lm', formula = y~x, colour=('#008B45'), fill = "#008B45", alpha = .8) + xlab("Executive Efficiency (z)") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
EE_V5<- EE_V5+theme(text = element_text(size=this_font_size),
                      axis.text = element_text(size = this_font_size),
                      axis.title.y = element_text(size = this_font_size),
                      axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      legend.text = element_text(size = this_font_size),
                      legend.title = element_text(size = this_font_size),
                      axis.title = element_text(size = this_font_size),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_rect(fill = "transparent",colour = NA),
                      plot.background = element_rect(fill = "transparent",colour = NA),
                      plot.margin = unit(c(0.2,0.2,0,0.2), "cm")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r = 39, b = 0, l = 0),size = 15))  #Top, left,Bottom, right
EE_V5
# save
ggsave(plot = EE_V5,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure5D_EE_V5.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


################################
#### SUPPLEMENTARY ANALYSES ####
################################

##############################################################
#### Controlling for Total Brain Volume (TBV) in the GAMs ####
##############################################################
#1.Executive Efficiency
gamModels_EE_TBV <- lapply(Components, function(x) {
  gam(substitute(i ~ s(Age) + F3_Executive_Efficiency + oSex + TBV + raw_num_bad_slices + mean_fd, list(i = as.name(x))), method="REML", data = df_fdc)
})

#2.Executive & Complex Reasoning Accuracy
gamModels_ECRA_TBV <- lapply(Components, function(x) {
  gam(substitute(i ~ s(Age) + F1_Exec_Comp_Res_Accuracy + oSex + TBV + raw_num_bad_slices + mean_fd, list(i = as.name(x))), method="REML", data = df_fdc)
})

#Look at model summaries
models_EE_TBV <- lapply(gamModels_EE_TBV, summary)
models_ECRA_TBV<-lapply(gamModels_ECRA_TBV, summary)
models_EE_TBV[[2]]

cog_tvalue_EE_TBV <- sapply(gamModels_EE_TBV, function(v) summary(v)$p.table[2,3])
cog_tvalue_EE_TBV <- as.data.frame(cog_tvalue_EE_TBV)
cog_pvalue_EE_TBV <- sapply(gamModels_EE_TBV, function(v) summary(v)$p.table[2,4])
cog_pvalue_EE_TBV <- as.data.frame(cog_pvalue_EE_TBV)
cog_pvalue_EE_TBV <- round(cog_pvalue_EE_TBV,3)
cog_pvalue_EE_TBV_fdr <- as.data.frame(p.adjust(cog_pvalue_EE_TBV[,1], method="fdr"))
cog_pvalue_EE_TBV_fdr <- round(cog_pvalue_EE_TBV_fdr,3)

cog_tvalue_ECRA_TBV <- sapply(gamModels_ECRA_TBV, function(v) summary(v)$p.table[2,3])
cog_tvalue_ECRA_TBV <- as.data.frame(cog_tvalue_ECRA_TBV)
cog_pvalue_ECRA_TBV <- sapply(gamModels_ECRA_TBV, function(v) summary(v)$p.table[2,4])
cog_pvalue_ECRA_TBV <- as.data.frame(cog_pvalue_ECRA_TBV)
cog_pvalue_ECRA_TBV <- round(cog_pvalue_ECRA_TBV,3)
cog_pvalue_ECRA_TBV_fdr <- as.data.frame(p.adjust(cog_pvalue_ECRA_TBV[,1], method="fdr"))
cog_pvalue_ECRA_TBV_fdr <- round(cog_pvalue_ECRA_TBV_fdr,3)

cog_stats_TBV<-cbind(bundles,cog_tvalue_EE_TBV,cog_pvalue_EE_TBV_fdr,cog_tvalue_ECRA_TBV,cog_pvalue_ECRA_TBV_fdr)
cog_stats_TBV<-cog_stats_TBV[order(-cog_tvalue_EE_TBV),]
cog_stats_TBV<-cog_stats_TBV%>%
  dplyr::rename(t_EE_TBV=cog_tvalue_EE_TBV)%>%
  dplyr::rename(p_fdr_EE_TBV=`p.adjust(cog_pvalue_EE_TBV[, 1], method = "fdr")`) %>%
  dplyr::rename(t_ECRA_TBV=cog_tvalue_ECRA_TBV)%>%
  dplyr::rename(p_fdr_ECRA_TBV=`p.adjust(cog_pvalue_ECRA_TBV[, 1], method = "fdr")`)

#Calculate Partial R2
redmodel_EE_TBV <- lapply(Components, function(x) {
  gam(substitute(i ~ s(Age) + oSex + raw_num_bad_slices + mean_fd + TBV , list(i = as.name(x))), method="REML", data = df_fdc)
})

partialR2_EE_TBV_V1 <- partialRsq(gamModels_EE_TBV[[1]],redmodel_EE_TBV[[1]])
partialR2_EE_TBV_V2 <- partialRsq(gamModels_EE_TBV[[2]],redmodel_EE_TBV[[2]])
partialR2_EE_TBV_V3 <- partialRsq(gamModels_EE_TBV[[3]],redmodel_EE_TBV[[3]])
partialR2_EE_TBV_V4 <- partialRsq(gamModels_EE_TBV[[4]],redmodel_EE_TBV[[4]])
partialR2_EE_TBV_V5 <- partialRsq(gamModels_EE_TBV[[5]],redmodel_EE_TBV[[5]])
partialR2_EE_TBV_V6 <- partialRsq(gamModels_EE_TBV[[6]],redmodel_EE_TBV[[6]])
partialR2_EE_TBV_V7 <- partialRsq(gamModels_EE_TBV[[7]],redmodel_EE_TBV[[7]])
partialR2_EE_TBV_V8 <- partialRsq(gamModels_EE_TBV[[8]],redmodel_EE_TBV[[8]])
partialR2_EE_TBV_V9 <- partialRsq(gamModels_EE_TBV[[9]],redmodel_EE_TBV[[9]])
partialR2_EE_TBV_V10 <- partialRsq(gamModels_EE_TBV[[10]],redmodel_EE_TBV[[10]])
partialR2_EE_TBV_V11 <- partialRsq(gamModels_EE_TBV[[11]],redmodel_EE_TBV[[11]])
partialR2_EE_TBV_V12 <- partialRsq(gamModels_EE_TBV[[12]],redmodel_EE_TBV[[12]])
partialR2_EE_TBV_V13 <- partialRsq(gamModels_EE_TBV[[13]],redmodel_EE_TBV[[13]])
partialR2_EE_TBV_V14 <- partialRsq(gamModels_EE_TBV[[14]],redmodel_EE_TBV[[14]])

#Merge Partial R2 values with F-stats and p-values
partialR2_EE_TBV<-as.data.frame(cbind(partialR2_EE_TBV_V1[[1]],partialR2_EE_TBV_V2[[1]],partialR2_EE_TBV_V3[[1]],partialR2_EE_TBV_V4[[1]],partialR2_EE_TBV_V5[[1]],
                                   partialR2_EE_TBV_V6[[1]],partialR2_EE_TBV_V7[[1]],partialR2_EE_TBV_V8[[1]],partialR2_EE_TBV_V9[[1]],partialR2_EE_TBV_V10[[1]],
                                   partialR2_EE_TBV_V11[[1]],partialR2_EE_TBV_V12[[1]],partialR2_EE_TBV_V13[[1]],partialR2_EE_TBV_V14[[1]]))
partialR2_EE_TBV<-as.data.frame(t(partialR2_EE_TBV))
partialR2_EE_TBV<-partialR2_EE_TBV %>%
  dplyr::rename(partialR2=V1)
cog_stats_TBV<-cbind(cog_stats_TBV,partialR2_EE_TBV)
cog_stats_TBV$partialR2<-as.numeric(as.character(cog_stats_TBV$partialR2))
cog_stats_TBV_order<-cog_stats_TBV[order(-cog_stats_TBV$partialR2),]


##Bar plot of EE partial R2 while controlling for TBV
cog_stats_TBV_order$bundles <- factor(cog_stats_TBV_order$bundles, levels = cog_stats_TBV_order$bundles)
figure3S<-ggplot(cog_stats_TBV_order, aes(x=bundles, y=partialR2, fill=bundles)) + geom_bar(stat="identity") + scale_fill_manual("Processing Method", values = c("CC (middle)" = "#F77F85", "SLF" = "#008B45", "Splenium" = "#EE3B3B", "Fornix" = "#4169E1", "Sup. CST" = "#9CB9F5", "Inf. CST" = "#6E91EB", "U-fibers" = "#9CCBB3", "Rostrum" = "#FFB6C1", "SLF (parietal)" = "#D0E0D8", "Uncinate" = "#4EAB7C", "Sup. Cerebellum" = "#F7EAAA", "Vermix" = "#F3D370", "Int. capsule" = "#CAE1FF", "Middle CP" = "#EEB422")) + labs(x ="Covariance Networks", y="Partial R2 for EE") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(face="bold",size=18, angle = 50, hjust=1), axis.text.y = element_text(face="bold",size=20), axis.title = element_text(size=28), axis.line= element_line(colour = 'black', size = 1.5), legend.position = 'none', legend.text = element_text(size=20), legend.title = element_text(size = 24), plot.title = element_text(face="bold",size = 20)) + ggtitle("Non-linear age effects on FDC") +scale_y_continuous(breaks=c(0,0.005,0.01,0.015)) +
  theme(plot.title = element_text(hjust = 0.5))
figure3S
ggsave(plot = figure3S,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure3S_bargraph_partialR2_EE_sensitivity.png",device = "png",width = 200,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


