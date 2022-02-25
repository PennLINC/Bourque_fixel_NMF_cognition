#This script is testing the development of white matter covariance networks in terms of Fiber density and cross-section (FDC) metric
#using Generelized Additive Models (GAMs). All GAMs include sex, mean DWI framewise displacement and number of DWI bad slices 
#(index of scan quality) as covariates. Included at the end of the script are sensitivity analyses in which we controlled 
#for Total Brain Volume. 

########################
#### LOAD LIBRARIES ####
########################
library(rlang)
library(mgcv)
library(stringr)
library(ggplot2)
library(fitdistrplus)
require(utils)
library(mgcViz)
library(cowplot)
library(dplyr)
library(olsrr)

######################
#### READ IN DATA ####
######################
#load components data
fdc <- read.csv("/cbica/projects/pnc_fixel_cs/GAMs/input_data/14comp_loadings_fdc_ltn.csv", header=F)
bblids<-read.csv("/cbica/projects/pnc_fixel_cs/nmf/mif_h5/cohort_files/ltn_FDC.csv",header=T)

#load in-scanner QC data
QC <- read.csv("/cbica/projects/pnc_fixel_cs/GAMs/input_data/QC_measures.csv", header=T, sep=",")

#load behavioral data
demo <- read.csv("/cbica/projects/GURLAB/dataFreezes/n1601_dataFreeze/demographics/n1601_demographics_go1_20161212.csv")

#load cognitive data
cnb_tm<-read.csv("/cbica/projects/GURLAB/dataFreezes/n1601_dataFreeze/cnb/n1601_cnb_factor_scores_tymoore_20151006.csv")

#load TBV data
TBV<-read.csv("/cbica/projects/pnc_fixel_cs/GAMs/input_data/n1601_ctVol20170412.csv")

#######################
#### ORGANIZE DATA ####
#######################
bblids<-bblids%>%
  dplyr::rename(bblid=subject)

comp_fdc<-cbind(bblids$bblid,fdc)
#Rename variable correctly
comp_fdc<-comp_fdc%>%
  dplyr::rename(bblid='bblids$bblid')

#Just keep necessary variables
TBV<-TBV %>%
  dplyr::select(bblid,mprage_antsCT_vol_TBV) %>%
  rename(TBV=mprage_antsCT_vol_TBV)

colnames(cnb_tm)
cog<-cnb_tm %>%
  dplyr::select(bblid,F1_Exec_Comp_Res_Accuracy,F3_Executive_Efficiency)

#Make sex an ordered variable and age in years
demo<-demo%>%
  mutate(oSex=sex) %>%
  mutate(Age=ageAtScan1/12)
demo$oSex<-ordered(demo$oSex)

#Join dataframes
df_fdc<-left_join(comp_fdc, demo, by="bblid")
df_fdc<-left_join(df_fdc, QC, by="bblid")
df_fdc<-left_join(df_fdc, cog, by="bblid")
df_fdc<-left_join(df_fdc, TBV, by="bblid")

#Get components numbering
comp_names<-df_fdc %>%
  dplyr::select(V1:V14)
Components <- names(comp_names)

createInteger <- function(f) {
  as.numeric(as.character(f))
}
comp_names <- as.data.frame(mapply(createInteger,comp_names))

#Check and remove missings for GAMs
sum(is.na(df_fdc$Age))
df_fdc<-df_fdc %>%
  filter(!is.na(Age)) %>%
  filter(!is.na(F1_Exec_Comp_Res_Accuracy)) %>%
  filter(!is.na(mean_fd)) %>%
  filter(!is.na(raw_num_bad_slices))

### FINAL SAMPLE FOR GAMs  IS N=939 (LOST N=2 DUE TO MISSING COGNITIVE DATA) ###
  
#List components' name
bundles<-c("Splenium","Fornix","Inf. CST","Int. capsule","SLF","CC (middle)","Rostrum","Sup. CST","Uncinate","SLF (parietal)","Middle CP",
           "Vermix","Sup. Cerebellum","U-fibers")
bundles<-as.data.frame(bundles)


########################
#### VISUALIZE DATA ####  ##Can skip this section##
########################
#Dependent variable (individuals' loadings on each of the 14 covariance networks)
foo <- function(x){
  require(ggplot2)
  ggplot(comp_names, aes(x = x)) + geom_histogram()
}
histograms<-lapply(comp_names,foo)
histogram<-plot_grid(plotlist=histograms,labels=c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10","V11","V12","V13","V14"),ncol=3, nrow=5)
print(histogram) 

#SCATTER PLOT
#Explore raw age effects for all covariance network
smooth_plot_age <- function(yvar){
  require(ggplot2)
  ggplot(df_fdc, aes_(x=~Age, y=as.name(yvar))) + geom_point() + geom_smooth(method = "gam", formula = y ~s(x))
}
comp_age_list<-lapply(names(df_fdc[c(2:15)]), smooth_plot_age)
plot_age_comp_smooth<-plot_grid(plotlist=comp_age_list)
plot_age_comp_smooth


########################################
#### BUILD GAM MODELS - AGE EFFECTS ####
########################################
#Full model
gamModels_age_fdc <- lapply(Components, function(x) {
  gam(substitute(i ~ s(Age) + oSex + raw_num_bad_slices + mean_fd, list(i = as.name(x))), method="REML", data = df_fdc)
})

#Look at model summaries
models_age_fdc <- lapply(gamModels_age_fdc, summary)

#Pull F-statistics
fstat_age_fdc <- sapply(gamModels_age_fdc, function(v) summary(v)$s.table[3])
fstat_age_fdc <- as.data.frame(fstat_age_fdc)

#Pull p-values - 
p_age_fdc <- sapply(gamModels_age_fdc, function(v) summary(v)$s.table[4])
p_age_fdc <- as.data.frame(p_age_fdc)
p_age_fdc <- round(p_age_fdc,3)

#FDR corrected p-values
p_age_fdc_fdr <- as.data.frame(p.adjust(p_age_fdc[,1], method="fdr"))
p_age_fdc_fdr <- round(p_age_fdc_fdr,3)

#Join F-stats and FDR corrected p-values for age effects
age_gam<-cbind(bundles,fstat_age_fdc,p_age_fdc_fdr)

#Get partial R2 (as a metric of effect size) based on Chenying script (https://github.com/PennLINC/ModelArray_paper/blob/enh/figures/notebooks/utils.R#L12)
#Reduced models
df_fdc_min<-df_fdc %>% #no missing data here on any variable, thus can use df_fdc_min data for reduced model testing
  dplyr::select(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,Age,oSex,mean_fd,raw_num_bad_slices)
count(df_fdc_min[rowSums(is.na(df_fdc_min))==0,])

redmodel_fdc <- lapply(Components, function(x) {
  gam(substitute(i ~ oSex + raw_num_bad_slices + mean_fd, list(i = as.name(x))), method="REML", data = df_fdc_min)
})

partialRsq <- function(fullmodel, redmodel) {
  # calculating SSE: used observed y (i.e. excluding observations with NA), and fitted values, directly from model object
  
  sse.full <- sum( (fullmodel$y - fullmodel$fitted.values)^2 )
  sse.red <- sum( (redmodel$y - redmodel$fitted.values)^2 )
  
  partialRsq <- (sse.red - sse.full) / sse.red
  
  toReturn <- list(partialRsq = partialRsq,
                   sse.full = sse.full,
                   sse.red = sse.red)
  return(toReturn)
}

partialR2_V1 <- partialRsq(gamModels_age_fdc[[1]],redmodel_fdc[[1]])
partialR2_V2 <- partialRsq(gamModels_age_fdc[[2]],redmodel_fdc[[2]])
partialR2_V3 <- partialRsq(gamModels_age_fdc[[3]],redmodel_fdc[[3]])
partialR2_V4 <- partialRsq(gamModels_age_fdc[[4]],redmodel_fdc[[4]])
partialR2_V5 <- partialRsq(gamModels_age_fdc[[5]],redmodel_fdc[[5]])
partialR2_V6 <- partialRsq(gamModels_age_fdc[[6]],redmodel_fdc[[6]])
partialR2_V7 <- partialRsq(gamModels_age_fdc[[7]],redmodel_fdc[[7]])
partialR2_V8 <- partialRsq(gamModels_age_fdc[[8]],redmodel_fdc[[8]])
partialR2_V9 <- partialRsq(gamModels_age_fdc[[9]],redmodel_fdc[[9]])
partialR2_V10 <- partialRsq(gamModels_age_fdc[[10]],redmodel_fdc[[10]])
partialR2_V11 <- partialRsq(gamModels_age_fdc[[11]],redmodel_fdc[[11]])
partialR2_V12 <- partialRsq(gamModels_age_fdc[[12]],redmodel_fdc[[12]])
partialR2_V13 <- partialRsq(gamModels_age_fdc[[13]],redmodel_fdc[[13]])
partialR2_V14 <- partialRsq(gamModels_age_fdc[[14]],redmodel_fdc[[14]])

#Merge Partial R2 values with F-stats and p-values
partialR2<-as.data.frame(cbind(partialR2_V1[[1]],partialR2_V2[[1]],partialR2_V3[[1]],partialR2_V4[[1]],partialR2_V5[[1]],
                               partialR2_V6[[1]],partialR2_V7[[1]],partialR2_V8[[1]],partialR2_V9[[1]],partialR2_V10[[1]],
                               partialR2_V11[[1]],partialR2_V12[[1]],partialR2_V13[[1]],partialR2_V14[[1]]))
partialR2<-as.data.frame(t(partialR2))
partialR2<-partialR2 %>%
  dplyr::rename(partial_R2=V1)
age_gam<-cbind(age_gam,partialR2)
age_gam$partial_R2<-as.numeric(as.character(age_gam$partial_R2))
age_gam<-age_gam[order(-age_gam$partial_R2),]


########################################################################
#### FIGURE 3B - BAR PLOT OF PARTIAL R2 FOR EACH COVARIANCE NETWORK ####
########################################################################
#First find colors similar to Figure 2 which depicts the 14 covariance networks on the brain
# Define the color ramp (returns a function object)
red <- colorRamp(c("brown2", "lightpink"))
# Define the ramp hex triplets
red.list <- rgb(red(seq(0, 1, length = 10)), max = 255)
print(red.list)
# Plot the color ramp
barplot(rep(1, 10), axes = FALSE, space = 0, col = red.list)

blue <- colorRamp(c("royalblue", "lightsteelblue1"))
blue.list<- rgb(blue(seq(0, 1, length = 10)), max = 255)
print(blue.list)
barplot(rep(1, 10), axes = FALSE, space = 0, col = blue.list)

green <-colorRamp(c("springgreen4", "grey92"))
green.list<- rgb(green(seq(0, 1, length = 10)), max = 255)
print(green.list)
barplot(rep(1, 10), axes = FALSE, space = 0, col = green.list)

yellow <-colorRamp(c("goldenrod2", "lightgoldenrodyellow"))
yellow.list<- rgb(yellow(seq(0, 1, length = 10)), max = 255)
print(yellow.list)
barplot(rep(1, 10), axes = FALSE, space = 0, col = yellow.list)

library(ggplot2)
age_gam$bundles <- factor(age_gam$bundles, levels = age_gam$bundles)
figure3b<-ggplot(age_gam, aes(x=bundles, y=partial_R2, fill=bundles)) + geom_bar(stat="identity") + scale_fill_manual("Processing Method", values = c("CC (middle)" = "#F77F85", "SLF" = "#008B45", "Splenium" = "#EE3B3B", "Fornix" = "#4169E1", "Sup. CST" = "#9CB9F5", "Inf. CST" = "#6E91EB", "U-fibers" = "#9CCBB3", "Rostrum" = "#FFB6C1", "SLF (parietal)" = "#D0E0D8", "Uncinate" = "#4EAB7C", "Sup. Cerebellum" = "#F7EAAA", "Vermix" = "#F3D370", "Int. capsule" = "#CAE1FF", "Middle CP" = "#EEB422")) + labs(x ="Covariance Networks", y="Partial R2") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(face="bold",size=18, angle = 50, hjust=1), axis.text.y = element_text(face="bold",size=20), axis.title = element_text(size=28), axis.line= element_line(colour = 'black', size = 1.5), legend.position = 'none', legend.text = element_text(size=20), legend.title = element_text(size = 24), plot.title = element_text(face="bold",size = 20)) + ggtitle("Non-linear age effects on FDC") +scale_y_continuous(breaks=c(0,0.05,0.10,0.15,0.20)) +
  theme(plot.title = element_text(hjust = 0.5))
figure3b
ggsave(plot = figure3b,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure3B_bargraph_partialR2.png",device = "png",width = 200,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


#########################################
#### COVARIANCE NETWORKS PREDICT AGE ####
#########################################
#Age ~ framewise_displacement + scan quality + 14NMF components, using linear modeling

Fullmodel_NMFpredictage<-lm(Age ~ raw_num_bad_slices + mean_fd + V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14, data=df_fdc)


NullModel_NMFpredictage<-lm(Age ~ raw_num_bad_slices + mean_fd ,data=df_fdc)


NMF_only_predict_age<-lm(Age ~ V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12+V13+V14, data=df_fdc)


#Make sure network components are not collinear
ols_vif_tol(NMF_only_predict_age)

Fullmodel_NMFpredictage_summary<-summary(Fullmodel_NMFpredictage)
NullModel_NMFpredictage_summary<-summary(NullModel_NMFpredictage)
NMF_only_predict_age_summary<-summary(NMF_only_predict_age)
Fullmodel_NMFpredictage_summary
NullModel_NMFpredictage_summary
NMF_only_predict_age_summary

#1.Compare Fullmodel and Nullmodel with an F test to test significant contribution of networks in predicting age
Ftest<-anova(Fullmodel_NMFpredictage,NullModel_NMFpredictage)
Ftest

#2.Correlation between actual age and predicted age by NMF_only_predict_age model
NMFCor <- cor.test(predict(NMF_only_predict_age), df_fdc$Age)
NMFCor #the r value will be added in Figure 3C

###########################################################################
#### FIGURE 3C - CORRELATION PLOT BETWEEN PREDICTED AGE AND ACTUAL AGE ####
###########################################################################
png(file='/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/Figures/Figure3C_NMF_predicted_vs_actual_age.png',width=2000,height=2000,res=300)
par(mar=c(5,5,2,2), cex.axis=2, bty="l")
plot(predict(NMF_only_predict_age), df_fdc$Age, ylab="Age", xlab="Predicted Age", cex.lab=2, pch=19, col="darkorange1")
abline(a=0,b=1, col = "darkorange3", lwd = 4)
dev.off()


#####################################################################
#### FIGURE 4 - DEVELOPMENTAL PLOTS FOR EACH COVARIANCE NETWORKS ####
#####################################################################
#This part of the code was taken from Adam Pines' 
#(https://github.com/PennBBL/multishell_diffusion/blob/master/PostProc/multishell_analyses.Rmd)

this_font_size=50
library(gratia)
library(dplyr)
library(svglite)
library(cowplot)

get_derivs_and_plot <- function(modobj,smooth_var,low_color=NULL,hi_color=NULL){
  #this_font_size = font_size*1.25
  if (is.null(low_color)){low_color = "white"}
  if (is.null(hi_color)){hi_color = "gray40"}
  derv<-derivatives(modobj,term=smooth_var)
  derv<- derv %>%
    mutate(sig = !(0 >lower & 0 < upper))
  derv$sig_deriv = derv$derivative*derv$sig
  cat(sprintf("\nSig change: %1.2f - %1.2f\n",min(derv$data[derv$sig==T]),max(derv$data[derv$sig==T])))
  d1<- ggplot(data=derv) + geom_tile(aes(x = data, y = .5, fill = sig_deriv))
  
  # will have to change this for negatively-sloped age-effects (swap low and high), color corresponds to metric family used
  d1 <- d1 + scale_fill_gradient(low = low_color, high = hi_color,limits = c(0,max(derv$derivative)))
  #    d1 <- d1 + scale_fill_gradient2(low = "#324c1d", midpoint = 0, mid = "white", high = "white",limits = c(min(derv$derivative),max(derv$derivative)))
  
  d1 <- d1 + 
    labs(x = smooth_var,fill = sprintf("")) + 
    theme(axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = this_font_size),
          axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          text = element_text(size=this_font_size),
          legend.text = element_text(size = this_font_size),
          axis.title = element_text(size = this_font_size),
          legend.key.width = unit(1,"cm"),
          legend.position = "right",
          plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    guides(fill = guide_colorbar(ticks=T, ticks.linewidth = 2, ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top")) +
    geom_rect(aes(ymin=0,ymax=1,xmin=min(data),xmax=max(data)),color="gray20",fill="white",alpha = 0)
  return(d1)
}

#Developmental plot for Covariance Network #1 (i.e.,V1)
p1<-ggplot(df_fdc,aes(x=Age,y=V1)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('#EE3B3B'), fill = "#EE3B3B", alpha = .8) + xlab("Age") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
p1<- p1+theme(text = element_text(size=this_font_size),
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

p1
# derivative calculation and plot object
d1 <- get_derivs_and_plot(modobj = gamModels_age_fdc[[1]],smooth_var = "s(Age)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 1.5), axis.line.y = element_line(colour = 'gray20', size = 1.5), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50)) + guides(fill = guide_colorbar(draw.ulim = F, frame.colour = "gray20", frame.linetype=1, frame.linewidth = 3 , ticks=T, ticks.linewidth = 2,  ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top"))  +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) #Top, left,Bottom, right
# stick plots together
scatter <- list(p1)
bar <- list(d1)
allplots <- c(scatter,bar)
pg<-plot_grid(rel_heights = c(16,2),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
final_plot <- pg
final_plot
# save
ggsave(plot = final_plot,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure4_derivative_V1.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)

#Developmental plot for Covariance Network #2 (i.e.,V2)
p2<-ggplot(df_fdc,aes(x=Age,y=V2)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('#4169E1'), fill = "#4169E1", alpha = .8) + xlab("Age") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
p2<- p2+theme(text = element_text(size=this_font_size),
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

p2
# derivative calculation and plot object
d2 <- get_derivs_and_plot(modobj = gamModels_age_fdc[[2]],smooth_var = "s(Age)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 1.5), axis.line.y = element_line(colour = 'gray20', size = 1.5), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50)) + guides(fill = guide_colorbar(draw.ulim = F, frame.colour = "gray20", frame.linetype=1, frame.linewidth = 3 , ticks=T, ticks.linewidth = 2,  ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top"))  +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) #Top, left,Bottom, right
# stick plots together
scatter <- list(p2)
bar <- list(d2)
allplots <- c(scatter,bar)
pg<-plot_grid(rel_heights = c(16,2),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
final_plot <- pg
final_plot
# save
ggsave(plot = final_plot,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure4_derivative_V2.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


#Developmental plot for Covariance Network #3 (i.e.,V3)
p3<-ggplot(df_fdc,aes(x=Age,y=V3)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('#6E91EB'), fill = "#6E91EB", alpha = .8) + xlab("Age") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
p3<- p3+theme(text = element_text(size=this_font_size),
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

p3
# derivative calculation and plot object
d3 <- get_derivs_and_plot(modobj = gamModels_age_fdc[[3]],smooth_var = "s(Age)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 1.5), axis.line.y = element_line(colour = 'gray20', size = 1.5), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50)) + guides(fill = guide_colorbar(draw.ulim = F, frame.colour = "gray20", frame.linetype=1, frame.linewidth = 3 , ticks=T, ticks.linewidth = 2,  ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top"))  +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) #Top, left,Bottom, right
# stick plots together
scatter <- list(p3)
bar <- list(d3)
allplots <- c(scatter,bar)
pg<-plot_grid(rel_heights = c(16,2),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
final_plot <- pg
final_plot
# save
ggsave(plot = final_plot,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure4_derivative_V3.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


#Developmental plot for Covariance Network #4 (i.e.,V4)
p4<-ggplot(df_fdc,aes(x=Age,y=V4)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('#CAE1FF'), fill = "#CAE1FF", alpha = .8) + xlab("Age") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
p4<- p4+theme(text = element_text(size=this_font_size),
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

p4
# derivative calculation and plot object
d4 <- get_derivs_and_plot(modobj = gamModels_age_fdc[[4]],smooth_var = "s(Age)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 1.5), axis.line.y = element_line(colour = 'gray20', size = 1.5), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50)) + guides(fill = guide_colorbar(draw.ulim = F, frame.colour = "gray20", frame.linetype=1, frame.linewidth = 3 , ticks=T, ticks.linewidth = 2,  ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top"))  +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) #Top, left,Bottom, right
# stick plots together
scatter <- list(p4)
bar <- list(d4)
allplots <- c(scatter,bar)
pg<-plot_grid(rel_heights = c(16,2),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
final_plot <- pg
final_plot
# save
ggsave(plot = final_plot,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure4_derivative_V4.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


#Developmental plot for Covariance Network #5 (i.e.,V5)
p5<-ggplot(df_fdc,aes(x=Age,y=V5)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('#008B45'), fill = "#008B45", alpha = .8) + xlab("Age") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
p5<- p5+theme(text = element_text(size=this_font_size),
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

p5
# derivative calculation and plot object
d5 <- get_derivs_and_plot(modobj = gamModels_age_fdc[[5]],smooth_var = "s(Age)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 1.5), axis.line.y = element_line(colour = 'gray20', size = 1.5), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50)) + guides(fill = guide_colorbar(draw.ulim = F, frame.colour = "gray20", frame.linetype=1, frame.linewidth = 3 , ticks=T, ticks.linewidth = 2,  ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top"))  +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) #Top, left,Bottom, right
# stick plots together
scatter <- list(p5)
bar <- list(d5)
allplots <- c(scatter,bar)
pg<-plot_grid(rel_heights = c(16,2),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
final_plot <- pg
final_plot
# save
ggsave(plot = final_plot,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure4_derivative_V5.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


#Developmental plot for Covariance Network #6 (i.e.,V6)
p6<-ggplot(df_fdc,aes(x=Age,y=V6)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('#F77F85'), fill = "#F77F85", alpha = .8) + xlab("Age") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
p6<- p6+theme(text = element_text(size=this_font_size),
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

p6
# derivative calculation and plot object
d6 <- get_derivs_and_plot(modobj = gamModels_age_fdc[[6]],smooth_var = "s(Age)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 1.5), axis.line.y = element_line(colour = 'gray20', size = 1.5), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50)) + guides(fill = guide_colorbar(draw.ulim = F, frame.colour = "gray20", frame.linetype=1, frame.linewidth = 3 , ticks=T, ticks.linewidth = 2,  ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top"))  +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) #Top, left,Bottom, right
# stick plots together
scatter <- list(p6)
bar <- list(d6)
allplots <- c(scatter,bar)
pg<-plot_grid(rel_heights = c(16,2),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
final_plot <- pg
final_plot
# save
ggsave(plot = final_plot,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure4_derivative_V6.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


#Developmental plot for Covariance Network #7 (i.e.,V7)
p7<-ggplot(df_fdc,aes(x=Age,y=V7)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('#FFB6C1'), fill = "#FFB6C1", alpha = .8) + xlab("Age") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
p7<- p7+theme(text = element_text(size=this_font_size),
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

p7
# derivative calculation and plot object
d7 <- get_derivs_and_plot(modobj = gamModels_age_fdc[[7]],smooth_var = "s(Age)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 1.5), axis.line.y = element_line(colour = 'gray20', size = 1.5), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50)) + guides(fill = guide_colorbar(draw.ulim = F, frame.colour = "gray20", frame.linetype=1, frame.linewidth = 3 , ticks=T, ticks.linewidth = 2,  ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top"))  +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) #Top, left,Bottom, right
# stick plots together
scatter <- list(p7)
bar <- list(d7)
allplots <- c(scatter,bar)
pg<-plot_grid(rel_heights = c(16,2),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
final_plot <- pg
final_plot
# save
ggsave(plot = final_plot,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure4_derivative_V7.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


#Developmental plot for Covariance Network #8 (i.e.,V8)
p8<-ggplot(df_fdc,aes(x=Age,y=V8)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('#9CB9F5'), fill = "#9CB9F5", alpha = .8) + xlab("Age") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
p8<- p8+theme(text = element_text(size=this_font_size),
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

p8
# derivative calculation and plot object
d8 <- get_derivs_and_plot(modobj = gamModels_age_fdc[[8]],smooth_var = "s(Age)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 1.5), axis.line.y = element_line(colour = 'gray20', size = 1.5), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50)) + guides(fill = guide_colorbar(draw.ulim = F, frame.colour = "gray20", frame.linetype=1, frame.linewidth = 3 , ticks=T, ticks.linewidth = 2,  ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top"))  +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) #Top, left,Bottom, right
# stick plots together
scatter <- list(p8)
bar <- list(d8)
allplots <- c(scatter,bar)
pg<-plot_grid(rel_heights = c(16,2),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
final_plot <- pg
final_plot
# save
ggsave(plot = final_plot,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure4_derivative_V8.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


#Developmental plot for Covariance Network #9 (i.e.,V9)
p9<-ggplot(df_fdc,aes(x=Age,y=V9)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('#4EAB7C'), fill = "#4EAB7C", alpha = .8) + xlab("Age") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
p9<- p9+theme(text = element_text(size=this_font_size),
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

p9
# derivative calculation and plot object
d9 <- get_derivs_and_plot(modobj = gamModels_age_fdc[[9]],smooth_var = "s(Age)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 1.5), axis.line.y = element_line(colour = 'gray20', size = 1.5), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50)) + guides(fill = guide_colorbar(draw.ulim = F, frame.colour = "gray20", frame.linetype=1, frame.linewidth = 3 , ticks=T, ticks.linewidth = 2,  ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top"))  +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) #Top, left,Bottom, right
# stick plots together
scatter <- list(p9)
bar <- list(d9)
allplots <- c(scatter,bar)
pg<-plot_grid(rel_heights = c(16,2),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
final_plot <- pg
final_plot
# save
ggsave(plot = final_plot,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure4_derivative_V9.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


#Developmental plot for Covariance Network #10 (i.e.,V10)
p10<-ggplot(df_fdc,aes(x=Age,y=V10)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('#D0E0D8'), fill = "#D0E0D8", alpha = .8) + xlab("Age") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
p10<- p10+theme(text = element_text(size=this_font_size),
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

p10
# derivative calculation and plot object
d10 <- get_derivs_and_plot(modobj = gamModels_age_fdc[[10]],smooth_var = "s(Age)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 1.5), axis.line.y = element_line(colour = 'gray20', size = 1.5), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50)) + guides(fill = guide_colorbar(draw.ulim = F, frame.colour = "gray20", frame.linetype=1, frame.linewidth = 3 , ticks=T, ticks.linewidth = 2,  ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top"))  +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) #Top, left,Bottom, right
# stick plots together
scatter <- list(p10)
bar <- list(d10)
allplots <- c(scatter,bar)
pg<-plot_grid(rel_heights = c(16,2),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
final_plot <- pg
final_plot
# save
ggsave(plot = final_plot,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure4_derivative_V10.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


#Developmental plot for Covariance Network #11 (i.e.,V11)
p11<-ggplot(df_fdc,aes(x=Age,y=V11)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('#EEB422'), fill = "#EEB422", alpha = .8) + xlab("Age") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
p11<- p11+theme(text = element_text(size=this_font_size),
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

p11
# derivative calculation and plot object
d11 <- get_derivs_and_plot(modobj = gamModels_age_fdc[[11]],smooth_var = "s(Age)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 1.5), axis.line.y = element_line(colour = 'gray20', size = 1.5), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50)) + guides(fill = guide_colorbar(draw.ulim = F, frame.colour = "gray20", frame.linetype=1, frame.linewidth = 3 , ticks=T, ticks.linewidth = 2,  ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top"))  +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) #Top, left,Bottom, right
# stick plots together
scatter <- list(p11)
bar <- list(d11)
allplots <- c(scatter,bar)
pg<-plot_grid(rel_heights = c(16,2),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
final_plot <- pg
final_plot
# save
ggsave(plot = final_plot,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure4_derivative_V11.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


#Developmental plot for Covariance Network #12 (i.e.,V12)
p12<-ggplot(df_fdc,aes(x=Age,y=V12)) + scale_y_continuous(breaks = seq(40, 100,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('#F3D370'), fill = "#F3D370", alpha = .8) + xlab("Age") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
p12<- p12+theme(text = element_text(size=this_font_size),
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

p12
# derivative calculation and plot object
d12 <- get_derivs_and_plot(modobj = gamModels_age_fdc[[12]],smooth_var = "s(Age)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 1.5), axis.line.y = element_line(colour = 'gray20', size = 1.5), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50)) + guides(fill = guide_colorbar(draw.ulim = F, frame.colour = "gray20", frame.linetype=1, frame.linewidth = 3 , ticks=T, ticks.linewidth = 2,  ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top"))  +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) #Top, left,Bottom, right
# stick plots together
scatter <- list(p12)
bar <- list(d12)
allplots <- c(scatter,bar)
pg<-plot_grid(rel_heights = c(16,2),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
final_plot <- pg
final_plot
# save
ggsave(plot = final_plot,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure4_derivative_V12.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


#Developmental plot for Covariance Network #13 (i.e.,V13)
p13<-ggplot(df_fdc,aes(x=Age,y=V13)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('#F7EAAA'), fill = "#F7EAAA", alpha = .8) + xlab("Age") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
p13<- p13+theme(text = element_text(size=this_font_size),
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

p13
# derivative calculation and plot object
d13 <- get_derivs_and_plot(modobj = gamModels_age_fdc[[13]],smooth_var = "s(Age)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 1.5), axis.line.y = element_line(colour = 'gray20', size = 1.5), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50)) + guides(fill = guide_colorbar(draw.ulim = F, frame.colour = "gray20", frame.linetype=1, frame.linewidth = 3 , ticks=T, ticks.linewidth = 2,  ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top"))  +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) #Top, left,Bottom, right
# stick plots together
scatter <- list(p13)
bar <- list(d13)
allplots <- c(scatter,bar)
pg<-plot_grid(rel_heights = c(16,2),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
final_plot <- pg
final_plot
# save
ggsave(plot = final_plot,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure4_derivative_V13.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


#Developmental plot for Covariance Network #14 (i.e.,V14)
p14<-ggplot(df_fdc,aes(x=Age,y=V14)) + scale_y_continuous(breaks = seq(80, 150,by=30)) + geom_point(size=2, colour="gray56") + geom_smooth(method = 'gam', formula = y~s(x,k=4), colour=('#9CCBB3'), fill = "#9CCBB3", alpha = .8) + xlab("Age") +ylab("FDC") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 2), axis.line.y = element_line(colour = 'gray20', size = 2),  axis.text = element_text(size=50), axis.title = element_text(size=50), axis.title.y = element_text(margin = margin(t = 0, r =55, b = 0, l = 0)))
p14<- p14+theme(text = element_text(size=this_font_size),
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

p14
# derivative calculation and plot object
d14 <- get_derivs_and_plot(modobj = gamModels_age_fdc[[14]],smooth_var = "s(Age)") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + theme(axis.line.x = element_line(colour = 'gray20', size = 1.5), axis.line.y = element_line(colour = 'gray20', size = 1.5), axis.ticks.length = unit(.25, "cm"), axis.text = element_text(size=50), axis.title = element_text(size=50)) + guides(fill = guide_colorbar(draw.ulim = F, frame.colour = "gray20", frame.linetype=1, frame.linewidth = 3 , ticks=T, ticks.linewidth = 2,  ticks.colour = "gray20", reverse = F,direction = "vertical",title.position = "top"))  +
  theme(plot.margin=unit(c(0,0,0,0),"cm")) #Top, left,Bottom, right
# stick plots together
scatter <- list(p14)
bar <- list(d14)
allplots <- c(scatter,bar)
pg<-plot_grid(rel_heights = c(16,2),plotlist = allplots,align = "v",axis = "lr",ncol = 1)
final_plot <- pg
final_plot
# save
ggsave(plot = final_plot,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure4_derivative_V14.png",device = "png",width = 320,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)



################################
#### SUPPLEMENTARY ANALYSES ####
################################

################################################################
#### 1.Controlling for Total Brain Volume (TBV) in the GAMs ####
################################################################
gamModels_age_fdc_TBV <- lapply(Components, function(x) {
  gam(substitute(i ~ s(Age) + oSex + TBV + raw_num_bad_slices + mean_fd, list(i = as.name(x))), method="REML", data = df_fdc)
})

#Look at model summaries
models_age_fdc_TBV <- lapply(gamModels_age_fdc_TBV, summary)

## MAIN EFFECT OF AGE WITH TBV ##
#Pull F-statistics
fstat_fdc_TBV <- sapply(gamModels_age_fdc_TBV, function(v) summary(v)$s.table[3])
fstat_fdc_TBV <- as.data.frame(fstat_fdc_TBV)

#Pull p-values - 
p_age_fdc_TBV <- sapply(gamModels_age_fdc_TBV, function(v) summary(v)$s.table[4])
p_age_fdc_TBV <- as.data.frame(p_age_fdc_TBV)
p_age_fdc_TBV <- round(p_age_fdc_TBV,3)

#FDR corrected p-values
p_age_fdc_fdr_TBV <- as.data.frame(p.adjust(p_age_fdc_TBV[,1], method="fdr"))
p_age_fdc_fdr_TBV <- round(p_age_fdc_fdr_TBV,3)

age_gam_TBV<-cbind(bundles,fstat_fdc_TBV,p_age_fdc_fdr_TBV)
age_gam_TBV<-age_gam_TBV%>%
  dplyr::rename(f_age=fstat_fdc_TBV)%>%
  dplyr::rename(p_fdr_age=`p.adjust(p_age_fdc_TBV[, 1], method = "fdr")`)

#Calculate Partial R2
df_fdc_TBV_min<-df_fdc %>% #no missing data here on any variable, thus can use df_fdc_min data for reduced model testing
  dplyr::select(V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13,V14,Age,oSex,mean_fd,raw_num_bad_slices,TBV)
count(df_fdc_TBV_min[rowSums(is.na(df_fdc_TBV_min))==0,])

redmodel_fdc_TBV <- lapply(Components, function(x) {
  gam(substitute(i ~ oSex + raw_num_bad_slices + mean_fd + TBV, list(i = as.name(x))), method="REML", data = df_fdc_TBV_min)
})

partialR2_TBV_V1 <- partialRsq(gamModels_age_fdc_TBV[[1]],redmodel_fdc_TBV[[1]])
partialR2_TBV_V2 <- partialRsq(gamModels_age_fdc_TBV[[2]],redmodel_fdc_TBV[[2]])
partialR2_TBV_V3 <- partialRsq(gamModels_age_fdc_TBV[[3]],redmodel_fdc_TBV[[3]])
partialR2_TBV_V4 <- partialRsq(gamModels_age_fdc_TBV[[4]],redmodel_fdc_TBV[[4]])
partialR2_TBV_V5 <- partialRsq(gamModels_age_fdc_TBV[[5]],redmodel_fdc_TBV[[5]])
partialR2_TBV_V6 <- partialRsq(gamModels_age_fdc_TBV[[6]],redmodel_fdc_TBV[[6]])
partialR2_TBV_V7 <- partialRsq(gamModels_age_fdc_TBV[[7]],redmodel_fdc_TBV[[7]])
partialR2_TBV_V8 <- partialRsq(gamModels_age_fdc_TBV[[8]],redmodel_fdc_TBV[[8]])
partialR2_TBV_V9 <- partialRsq(gamModels_age_fdc_TBV[[9]],redmodel_fdc_TBV[[9]])
partialR2_TBV_V10 <- partialRsq(gamModels_age_fdc_TBV[[10]],redmodel_fdc_TBV[[10]])
partialR2_TBV_V11 <- partialRsq(gamModels_age_fdc_TBV[[11]],redmodel_fdc_TBV[[11]])
partialR2_TBV_V12 <- partialRsq(gamModels_age_fdc_TBV[[12]],redmodel_fdc_TBV[[12]])
partialR2_TBV_V13 <- partialRsq(gamModels_age_fdc_TBV[[13]],redmodel_fdc_TBV[[13]])
partialR2_TBV_V14 <- partialRsq(gamModels_age_fdc_TBV[[14]],redmodel_fdc_TBV[[14]])

#Merge Partial R2 values with F-stats and p-values
partialR2_TBV<-as.data.frame(cbind(partialR2_TBV_V1[[1]],partialR2_TBV_V2[[1]],partialR2_TBV_V3[[1]],partialR2_TBV_V4[[1]],partialR2_TBV_V5[[1]],
                               partialR2_TBV_V6[[1]],partialR2_TBV_V7[[1]],partialR2_TBV_V8[[1]],partialR2_TBV_V9[[1]],partialR2_TBV_V10[[1]],
                               partialR2_TBV_V11[[1]],partialR2_TBV_V12[[1]],partialR2_TBV_V13[[1]],partialR2_TBV_V14[[1]]))
partialR2_TBV<-as.data.frame(t(partialR2_TBV))
partialR2_TBV<-partialR2_TBV %>%
  dplyr::rename(partialR2=V1)
age_gam_TBV<-cbind(age_gam_TBV,partialR2_TBV)
age_gam_TBV$partialR2<-as.numeric(as.character(age_gam_TBV$partialR2))
age_gam_TBV<-age_gam_TBV[order(-age_gam_TBV$partialR2),]

##Bar plot of partial R2 while controlling for TBV
age_gam_TBV$bundles <- factor(age_gam_TBV$bundles, levels = age_gam_TBV$bundles)
figure2S<-ggplot(age_gam_TBV, aes(x=bundles, y=partialR2, fill=bundles)) + geom_bar(stat="identity") + scale_fill_manual("Processing Method", values = c("CC (middle)" = "#F77F85", "SLF" = "#008B45", "Splenium" = "#EE3B3B", "Fornix" = "#4169E1", "Sup. CST" = "#9CB9F5", "Inf. CST" = "#6E91EB", "U-fibers" = "#9CCBB3", "Rostrum" = "#FFB6C1", "SLF (parietal)" = "#D0E0D8", "Uncinate" = "#4EAB7C", "Sup. Cerebellum" = "#F7EAAA", "Vermix" = "#F3D370", "Int. capsule" = "#CAE1FF", "Middle CP" = "#EEB422")) + labs(x ="Covariance Networks", y="Partial R2") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(face="bold",size=18, angle = 50, hjust=1), axis.text.y = element_text(face="bold",size=20), axis.title = element_text(size=28), axis.line= element_line(colour = 'black', size = 1.5), legend.position = 'none', legend.text = element_text(size=20), legend.title = element_text(size = 24), plot.title = element_text(face="bold",size = 20)) + ggtitle("Non-linear age effects on FDC") +scale_y_continuous(breaks=c(0,0.05,0.10,0.15,0.20)) +
  theme(plot.title = element_text(hjust = 0.5))
figure2S
ggsave(plot = figure2S,filename = "/Users/jbourque/UPENN/projects/PNC_fixel/gam_analysis/figures/Figure2S_bargraph_partialR2_sensitivity.png",device = "png",width = 200,height = 210,units = "mm")# Print the saved plot to the rmarkdown document (optional)


##############################################
#### 2.Main effect of sex (including TBV) ####
##############################################
models_age_fdc_TBV[[1]]
#Pull t-statistics
tstat_sex <- sapply(gamModels_age_fdc_TBV, function(v) summary(v)$p.table[2,3])
tstat_sex <- as.data.frame(tstat_sex)
#Positive t-stat is associated with increased TBV in females, relative to males

#Pull p-values - 
p_sex_fdc <- sapply(gamModels_age_fdc_TBV, function(v) summary(v)$p.table[2,4])
p_sex_fdc <- as.data.frame(p_sex_fdc)
p_sex_fdc <- round(p_sex_fdc,3)

#FDR corrected p-values
p_sex_fdc_fdr <- as.data.frame(p.adjust(p_sex_fdc[,1], method="fdr"))
p_sex_fdc_fdr <- round(p_sex_fdc_fdr,3)

#Join DFs
sex_gam<-cbind(bundles,tstat_sex,p_sex_fdc_fdr)
sex_gam<-sex_gam%>%
  dplyr::rename(t_sex=tstat_sex) %>%
  dplyr::rename(p_fdr_sex=`p.adjust(p_sex_fdc[, 1], method = "fdr")`)


####################################################
#### 3.Age by sex interactions? (including TBV) ####
####################################################
#Males are be the reference group
gamModels_fdc_int_TBV<-lapply(Components, function(x) {
  gam(substitute(i ~ s(Age) + oSex + s(Age,by=oSex) + mean_fd + raw_num_bad_slices + TBV, list(i=as.name(x))), method="REML", data=df_fdc)
})

#Look at models' summaries
models_int_TBV <- lapply(gamModels_fdc_int_TBV, summary)
models_int_TBV[[1]]

#extract p-value of age by sex interaction
p_sex_int_TBV <- sapply(gamModels_fdc_int_TBV, function(v) summary(v)$s.table[8])
p_sex_int_TBV <- as.data.frame(p_sex_int_TBV)
p_sex_int_fdr_TBV <- as.data.frame(p.adjust(p_sex_int_TBV[,1], method="fdr"))
p_sex_int_fdr_TBV <- round(p_sex_int_fdr_TBV,3)

#extract f-stat of age by sex interaction
f_sex_int_TBV <- sapply(gamModels_fdc_int_TBV, function(v) summary(v)$s.table[6])
f_sex_int_TBV <- as.data.frame(f_sex_int_TBV)

sex_int_fdr<-cbind(bundles,f_sex_int_TBV,p_sex_int_fdr_TBV)
colnames(sex_int_fdr)
sex_int_fdr<-sex_int_fdr%>%
  dplyr::rename(p_fdr_sex_int_TBV=`p.adjust(p_sex_int_TBV[, 1], method = "fdr")`)
#No age by sex interaction after FDR correction
