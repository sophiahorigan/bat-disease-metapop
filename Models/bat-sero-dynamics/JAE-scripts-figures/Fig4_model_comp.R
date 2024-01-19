
rm(list=ls())

library(epitools)
library(dplyr)
library(plyr)
library(cowplot)
library(deSolve)
library(mgcv)
library(lubridate)
library(matrixcalc)
library(Matrix)
library(ggplot2)
#library(IPMpack)
#library(splitstackshape)

#axdding in births distributed across biweeks means we need to make a giant matrix for both demographic and epidemic transitions

homewd="/Users/carabrook/Developer/bat-sero-dynamics"

setwd(paste0(homewd, "/JAE-scripts-figures/"))
load(paste0(homewd,"/JAE-scripts-figures/data/combine.fit.all.Rdata"))
load(paste0(homewd,"/JAE-scripts-figures/data/combine.mod.all.Rdata"))



get.seroprev.dat = function(data, vis_split, cutoff){
  
  visbin = seq(3,floor(max(data$age)), vis_split)
  #but breakdown the early years into more
  visbin = c(c(0,.5, 1, 2),   visbin)
  data$age_year <- NA
  
  for (i in 1:length(data$age_year)){
    tmp = visbin[data$age[i] > visbin]
    data$age_year[i] <- tmp[length(tmp)] 
  }
  
  if(cutoff=="mean"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev)/length(prev), seropos=sum(prev), count=length(prev))  
  }else if (cutoff=="uci"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev_uci)/length(prev_uci), seropos=sum(prev), count=length(prev_uci)) 
  }else if (cutoff=="lci"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev_lci)/length(prev_lci), seropos=sum(prev), count=length(prev_lci)) 
  }
  
  #you want the midpoint between either end of the age class now
  vect.age.yr = sort(unique(data$age_year))
  dat.sum2$age_plot = NA
  
  for (i in 1:(length(dat.sum2$age_year)-1)){
    dat.sum2$age_plot [i] = ((dat.sum2$age_year[i + 1] - dat.sum2$age_year[i])/2) + dat.sum2$age_year[i]
  }
  dat.sum2$age_plot[length(dat.sum2$age_plot)] =  (ceiling(max(data$age)) - dat.sum2$age_year[length(dat.sum2$age_year)])/2  + dat.sum2$age_year[length(dat.sum2$age_year)]
  
  names(dat.sum2)  <- c("real_age", "prevalence", "pos", "count", "age_year") #<- names(dat.sum2.tmp)
  
  # dat.sum2 <-  rbind(dat.sum2.tmp, dat.sum2)
  
  dat.sum3 = dat.sum2
  dat.sum3$class = "seropositive"
  rownames(dat.sum3) <- c()
  dat.sum3$species = unique(data$species)
  dat.sum3$type= unique(data$type)
  
  return(dat.sum3)
  
}
get.avg.age.sero <- function(dat){
  #round each age by year
  dat$age <- floor(dat$age)
  dat2 <- ddply(dat, .(model, age), summarize, mean_seroprev = mean(seroprevalence), lci_seroprev= mean(seroprev_lci), uci_seroprev= mean(seroprev_uci))
  dat2$species = unique(dat$species)
  dat2$type= unique(dat2$type)
  return(dat2)
}

load(paste0(homewd,"/JAE-scripts-figures/data/seasonal_prev.Rdata"))


dat.tot <- dplyr::select(dat.new, -(month), -(year), -(titer), -(titer_resid), -(titer_resid_lci), -(titer_resid_uci), -(doy))
head(dat.tot)

plot.annual.age.seroprev.all.mods.squiggle <- function(mod.dat, prev.dat, AIC.dat, type1, cutoff1, species1, do.save, filename){
  #first, choose
  prev.dat1 <- subset(prev.dat, species==species1 & type==type1)
  prev.dat1  <- subset(prev.dat1, !is.na(age))
  prev.dat1  <- subset(prev.dat1, fit==1)
  
  #get age.seroprev from data
  dat.seroprev = get.seroprev.dat(data=prev.dat1, vis_split=3, cutoff=cutoff1)
  
  mod.dat1 = subset(mod.dat, species==species1 & type==type1 & cutoff==cutoff1)
  
  AIC.dat1 = subset(AIC.dat, species==species1 & type==type1 & cutoff==cutoff1)
  #add an "age" value to fit.dat for plotting
  AIC.dat1$age <- NA
  if(species1=="Eidolon dupreanum"){
    AIC.dat1$age<- 18.2
  }else{
    AIC.dat1$age <- 13.5  
  }
  
  
  #and get deltaAIC relative to seroprev
  AIC.df <- subset(AIC.dat1, cutoff==cutoff1)
  min.AIC = min(AIC.df$AIC)
  AIC.df$deltaAIC =AIC.df$AIC-min.AIC
  
  
  max.a = ceiling(max(dat.seroprev$age_year))
  #if(cutoff1=="mean" | cutoff1=="lci"){
   # ymax =.5
  #}else if (cutoff1=="uci"){
   # ymax=
  #}
    
  #AIC.mult <- 1/(max(AIC.df$deltaAIC) + 1)
  #AIC.df$relAIC <- AIC.df$deltaAIC*AIC.mult
  AIC.mult <- ceiling(max(AIC.df$deltaAIC))
  AIC.df$relAIC <- AIC.df$deltaAIC/AIC.mult
  AIC.df$model <- factor(AIC.df$model, levels = c("MSIR", "MSRIR", "MSIRS", "MSIRN-matAB","MSIRN-matSus", "MSIRNR-matAB",  "MSIRNR-matSus"))
  AIC.df$outline <- NA
  AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
  AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
  AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
  
  
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 16, "EBOV" = 17)
  
  
  #then, reorder the models
  mod.dat1$model = factor(mod.dat1$model, levels = c("MSIR", "MSRIR", "MSIRS", "MSIRN-matAB","MSIRN-matSus", "MSIRNR-matAB",  "MSIRNR-matSus"))
  #with intra-annual squiggles:
  mod.dat2 <- ddply(mod.dat1, .(model, age), summarize, mean_seroprev = mean(seroprevalence), lci_seroprev= mean(seroprev_lci), uci_seroprev= mean(seroprev_uci), species=unique(species), type=unique(type)) 
  #mod.dat2 <- get.avg.age.sero(dat=mod.dat1) #no squiggles
  mod.dat2 = subset(mod.dat2, age<=15)
  
  #plot together
 # p1 <- ggplot(data=mod.dat2) + geom_line(aes(x=age, y=mean_seroprev), color="purple", size=1) + 
  #  facet_grid(~model) + geom_ribbon(aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), alpha=.3, fill="purple") +
   # geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence), linetype=2) +
    #geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), pch=1) +
    #coord_cartesian(xlim=c(0,15), ylim=c(0,1)) 
  
  #print(p1)
  label.dat = cbind.data.frame(c("MSIR", "MSRIR", "MSIRS", "MSIRN-matAB","MSIRN-matSus", "MSIRNR-matAB",  "MSIRNR-matSus"), c("A.", "B.", "C.", "D.", "E.", "F.", "G."))
  names(label.dat) = c("model", "label")
  #plot data
  if(species1=="Eidolon dupreanum"){
    
  p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1,  alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
    geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
    # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
    geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="purple", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
    geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "purple", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =18, size =4, color="navy", show.legend = FALSE) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 5, color="navy", size=5, show.legend = FALSE) +
    geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=5) +
    theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 16), 
          axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy", size=16), axis.text.y.right = element_text(colour = "navy", size=12),
          legend.text = element_blank(), plot.margin = unit(c(.1,.1,.1,.1), "cm")) + 
    ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=15), color="navy") + coord_cartesian(xlim = c(0,20), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
    scale_x_discrete(limits = c(0,10,15), labels = c("0", "10", "15")) + geom_segment(data=AIC.df, aes(x=0, xend=20, y=0, yend=0), linetype=2, color="navy")
  
  }else{
    mod.dat2 = subset(mod.dat2, age<=11)
    p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=2,  alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
      geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
      # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
      geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="lightgreen", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
      geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "lightgreen", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
      geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=5) +
      theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 16), 
            axis.title.x = element_text(size=16), axis.text = element_text(size=12), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy", size=16), axis.text.y.right = element_text(colour = "navy", size=12),
            legend.text = element_blank(), plot.margin = unit(c(.1,.1,.1,.1), "cm")) + 
      ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=11), color="navy") + coord_cartesian(xlim = c(0,15), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
      scale_x_discrete(limits = c(0,5,10), labels = c("0", "5", "10")) + geom_segment(data=AIC.df, aes(x=0, xend=15, y=0, yend=0), linetype=2, color="navy")
    
  }

  print(p2)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=140, 
           height=40, 
           scale=3, 
           dpi=300)
  }
  
  
}
plot.annual.age.seroprev.all.mods <- function(mod.dat, prev.dat, AIC.dat, type1, cutoff1, species1, do.save, filename){
  #first, choose
  prev.dat1 <- subset(prev.dat, species==species1 & type==type1)
  prev.dat1  <- subset(prev.dat1, !is.na(age))
  prev.dat1  <- subset(prev.dat1, fit==1)
  
  #get age.seroprev from data
  dat.seroprev = get.seroprev.dat(data=prev.dat1, vis_split=3, cutoff=cutoff1)
  
  mod.dat1 = subset(mod.dat, species==species1 & type==type1 & cutoff==cutoff1)
  
  AIC.dat1 = subset(AIC.dat, species==species1 & type==type1 & cutoff==cutoff1)
  #add an "age" value to fit.dat for plotting
  AIC.dat1$age <- NA
  if(species1=="Eidolon dupreanum"){
    AIC.dat1$age<- 18.2
  }else{
    AIC.dat1$age <- 13.5  
  }
  
  
  #and get deltaAIC relative to seroprev
  AIC.df <- subset(AIC.dat1, cutoff==cutoff1)
  min.AIC = min(AIC.df$AIC)
  AIC.df$deltaAIC =AIC.df$AIC-min.AIC
  
  
  max.a = ceiling(max(dat.seroprev$age_year))
  #if(cutoff1=="mean" | cutoff1=="lci"){
  # ymax =.5
  #}else if (cutoff1=="uci"){
  # ymax=
  #}
  
  #AIC.mult <- 1/(max(AIC.df$deltaAIC) + 1)
  #AIC.df$relAIC <- AIC.df$deltaAIC*AIC.mult
  AIC.mult <- ceiling(max(AIC.df$deltaAIC))
  AIC.df$relAIC <- AIC.df$deltaAIC/AIC.mult
  AIC.df$model <- factor(AIC.df$model, levels = c("MSIR", "MSRIR", "MSIRS", "MSIRN-matAB", "MSIRNR-matAB", "MSIRN-matSus", "MSIRNR-matSus"))
  AIC.df$outline <- NA
  AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
  AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
  AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
  
  
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 16, "EBOV" = 17)
  
  
  #then, reorder the models
  mod.dat1$model = factor(mod.dat1$model, levels = c("MSIR", "MSRIR", "MSIRS", "MSIRN-matAB", "MSIRNR-matAB", "MSIRN-matSus", "MSIRNR-matSus"))
  #with intra-annual squiggles:
  mod.dat2 <- ddply(mod.dat1, .(model, age), summarize, mean_seroprev = mean(seroprevalence), lci_seroprev= mean(seroprev_lci), uci_seroprev= mean(seroprev_uci), species=unique(species), type=unique(type)) 
  mod.dat2 <- get.avg.age.sero(dat=mod.dat1) #no squiggles
  mod.dat2 = subset(mod.dat2, age<=15)
  
  #plot together
  # p1 <- ggplot(data=mod.dat2) + geom_line(aes(x=age, y=mean_seroprev), color="purple", size=1) + 
  #  facet_grid(~model) + geom_ribbon(aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), alpha=.3, fill="purple") +
  # geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence), linetype=2) +
  #geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), pch=1) +
  #coord_cartesian(xlim=c(0,15), ylim=c(0,1)) 
  
  #print(p1)
  label.dat = cbind.data.frame(c("MSIR", "MSRIR", "MSIRS", "MSIRN-matAB", "MSIRNR-matAB", "MSIRN-matSus", "MSIRNR-matSus"), c("A.", "B.", "C.", "D.", "E.", "F.", "G."))
  names(label.dat) = c("model", "label")
  #plot data
  if (species1=="Eidolon dupreanum"){
  p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1,  alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
    geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
    # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
    geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="purple", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
    geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "purple", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =18, size =4, color="navy", show.legend = FALSE) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 5, color="navy", size=5, show.legend = FALSE) +
    geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=4.5) +
    theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 12), 
          axis.title.x = element_text(size=11), axis.text = element_text(size=10), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy"), axis.text.y.right = element_text(colour = "navy"),
          legend.text = element_blank(), plot.margin = unit(c(.2,.2,.2,1.3), "cm")) + 
    ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=15), color="navy") + coord_cartesian(xlim = c(0,20), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
    scale_x_discrete(limits = c(0,10,15), labels = c("0", "10", "15")) + geom_segment(data=AIC.df, aes(x=0, xend=20, y=0, yend=0), linetype=3, color="navy")
  
}else{
  mod.dat2 = subset(mod.dat2, age<=11)
  p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=2,  alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
    geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
    # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
    geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="lightgreen", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
    geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "lightgreen", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
    geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=4.5) +
    theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 12), 
          axis.title.x = element_text(size=11), axis.text = element_text(size=10), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy"), axis.text.y.right = element_text(colour = "navy"),
          legend.text = element_blank(), plot.margin = unit(c(.2,.2,.2,1.3), "cm")) + 
    ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=11), color="navy") + coord_cartesian(xlim = c(0,15), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
    scale_x_discrete(limits = c(0,5,10), labels = c("0", "5", "10")) + geom_segment(data=AIC.df, aes(x=0, xend=15, y=0, yend=0), linetype=2, color="navy")
  
}
  
  print(p2)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=100, 
           height=65, 
           scale=3, 
           dpi=300)
  }
  
  
}
plot.annual.age.seroprev.matSus.squiggle <- function(mod.dat, prev.dat, AIC.dat, type1, cutoff1, species1, do.save, filename){
  #first, choose
  prev.dat1 <- subset(prev.dat, species==species1 & type==type1)
  prev.dat1  <- subset(prev.dat1, !is.na(age))
  prev.dat1  <- subset(prev.dat1, fit==1)
  
  
  #get age.seroprev from data
  dat.seroprev = get.seroprev.dat(data=prev.dat1, vis_split=3, cutoff=cutoff1)
  
  mod.dat1 = subset(mod.dat, species==species1 & type==type1 & cutoff==cutoff1)
  mod.dat1$model[mod.dat1$model=="MSIRN-matSus"] = "MSIRN"
  mod.dat1$model[mod.dat1$model=="MSIRNR-matSus"] = "MSIRNR"
  mod.dat1 = subset(mod.dat1, model!="MSIRN-matAB")
  mod.dat1 = subset(mod.dat1, model!="MSIRNR-matAB")
  
  AIC.dat1 = subset(AIC.dat, species==species1 & type==type1 & cutoff==cutoff1)
  AIC.dat1$model[AIC.dat1$model=="MSIRN-matSus"] = "MSIRN"
  AIC.dat1$model[AIC.dat1$model=="MSIRNR-matSus"] = "MSIRNR"
  AIC.dat1 = subset(AIC.dat1, model!="MSIRN-matAB")
  AIC.dat1 = subset(AIC.dat1, model!="MSIRNR-matAB")
  #add an "age" value to fit.dat for plotting
  AIC.dat1$age <- NA
  if(species1=="Eidolon dupreanum"){
    AIC.dat1$age<- 18.2
  }else{
    AIC.dat1$age <- 13.5  
  }
  
  
  #and get deltaAIC relative to seroprev
  AIC.df <- subset(AIC.dat1, cutoff==cutoff1)
  min.AIC = min(AIC.df$AIC)
  AIC.df$deltaAIC =AIC.df$AIC-min.AIC
  
  
  max.a = ceiling(max(dat.seroprev$age_year))
  #if(cutoff1=="mean" | cutoff1=="lci"){
  # ymax =.5
  #}else if (cutoff1=="uci"){
  # ymax=
  #}
  
  #AIC.mult <- 1/(max(AIC.df$deltaAIC) + 1)
  #AIC.df$relAIC <- AIC.df$deltaAIC*AIC.mult
  AIC.mult <- ceiling(max(AIC.df$deltaAIC))
  AIC.df$relAIC <- AIC.df$deltaAIC/AIC.mult
  AIC.df$model <- factor(AIC.df$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  AIC.df$outline <- NA
  AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
  AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
  AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
  
  
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 16, "EBOV" = 17)
  
  
  #then, reorder the models
  mod.dat1$model = factor(mod.dat1$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  #with intra-annual squiggles:
  mod.dat2 <- ddply(mod.dat1, .(model, age), summarize, mean_seroprev = mean(seroprevalence), lci_seroprev= mean(seroprev_lci), uci_seroprev= mean(seroprev_uci), species=unique(species), type=unique(type)) 
  #mod.dat2 <- get.avg.age.sero(dat=mod.dat1) #no squiggles
  mod.dat2 = subset(mod.dat2, age<=15)
  #plot together
  # p1 <- ggplot(data=mod.dat2) + geom_line(aes(x=age, y=mean_seroprev), color="purple", size=1) + 
  #  facet_grid(~model) + geom_ribbon(aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), alpha=.3, fill="purple") +
  # geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence), linetype=2) +
  #geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), pch=1) +
  #coord_cartesian(xlim=c(0,15), ylim=c(0,1)) 
  
  #print(p1)
  label.dat = cbind.data.frame(c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"), c("A.", "B.", "C.", "D.", "E."))
  names(label.dat) = c("model", "label")
  #plot data
  if (species1=="Eidolon dupreanum"){
    
    p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1,  alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
      geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
      # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
      geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="purple", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
      geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "purple", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =18, size =4, color="navy", show.legend = FALSE) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 5, color="navy", size=5, show.legend = FALSE) +
      geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=4.5) +
      theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 12), 
            axis.title.x = element_text(size=11), axis.text = element_text(size=10), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy"), axis.text.y.right = element_text(colour = "navy"),
            legend.text = element_blank(), plot.margin = unit(c(.2,.2,.2,1.3), "cm")) + 
      ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=15), color="navy") + coord_cartesian(xlim = c(0,20), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
      scale_x_discrete(limits = c(0,10,15), labels = c("0", "10", "15")) + geom_segment(data=AIC.df, aes(x=0, xend=20, y=0, yend=0), linetype=2, color="navy")
    
  }else{
    mod.dat2 = subset(mod.dat2, age<=11)
    p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=2,  alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
      geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
      # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
      geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="lightgreen", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
      geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "lightgreen", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
      geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=4.5) +
      theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 12), 
            axis.title.x = element_text(size=11), axis.text = element_text(size=10), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy"), axis.text.y.right = element_text(colour = "navy"),
            legend.text = element_blank(), plot.margin = unit(c(.2,.2,.2,1.3), "cm")) + 
      ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=11), color="navy") + coord_cartesian(xlim = c(0,15), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
      scale_x_discrete(limits = c(0,5,10), labels = c("0", "5", "10")) + geom_segment(data=AIC.df, aes(x=0, xend=15, y=0, yend=0), linetype=2, color="navy")
    
  }
  
  
  print(p2)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=100, 
           height=65, 
           scale=3, 
           dpi=300)
  }
  
  
}
plot.annual.age.seroprev.matSus <- function(mod.dat, prev.dat, AIC.dat, type1, cutoff1, species1, do.save, filename){
  #first, choose
  prev.dat1 <- subset(prev.dat, species==species1 & type==type1)
  prev.dat1  <- subset(prev.dat1, !is.na(age))
  prev.dat1  <- subset(prev.dat1, fit==1)
  
  
  #get age.seroprev from data
  dat.seroprev = get.seroprev.dat(data=prev.dat1, vis_split=3, cutoff=cutoff1)
  
  mod.dat1 = subset(mod.dat, species==species1 & type==type1 & cutoff==cutoff1)
  mod.dat1$model = as.character(mod.dat1$model)
  mod.dat1$model[mod.dat1$model=="MSIRN-matSus"] = "MSIRN"
  mod.dat1$model[mod.dat1$model=="MSIRNR-matSus"] = "MSIRNR"
  mod.dat1 = subset(mod.dat1, model!="MSIRN-matAB")
  mod.dat1 = subset(mod.dat1, model!="MSIRNR-matAB")
  
  AIC.dat1 = subset(AIC.dat, species==species1 & type==type1 & cutoff==cutoff1)
  AIC.dat1$model = as.character(AIC.dat1$model)
  AIC.dat1$model[AIC.dat1$model=="MSIRN-matSus"] = "MSIRN"
  AIC.dat1$model[AIC.dat1$model=="MSIRNR-matSus"] = "MSIRNR"
  AIC.dat1 = subset(AIC.dat1, model!="MSIRN-matAB")
  AIC.dat1 = subset(AIC.dat1, model!="MSIRNR-matAB")
  #add an "age" value to fit.dat for plotting
  AIC.dat1$age <- NA
  if(species1=="Eidolon dupreanum"){
    AIC.dat1$age<- 18.2
  }else{
    AIC.dat1$age <- 13.5  
  }
  
  
  #and get deltaAIC relative to seroprev
  AIC.df <- subset(AIC.dat1, cutoff==cutoff1)
  min.AIC = min(AIC.df$AIC)
  AIC.df$deltaAIC =AIC.df$AIC-min.AIC
  
  
  max.a = ceiling(max(dat.seroprev$age_year))
  #if(cutoff1=="mean" | cutoff1=="lci"){
  # ymax =.5
  #}else if (cutoff1=="uci"){
  # ymax=
  #}
  
  #AIC.mult <- 1/(max(AIC.df$deltaAIC) + 1)
  #AIC.df$relAIC <- AIC.df$deltaAIC*AIC.mult
  AIC.mult <- ceiling(max(AIC.df$deltaAIC))
  AIC.df$relAIC <- AIC.df$deltaAIC/AIC.mult
  AIC.df$model <- factor(AIC.df$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  AIC.df$outline <- NA
  AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
  AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
  AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
  
  
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 16, "EBOV" = 17)
  
  label.dat = cbind.data.frame(c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"), c("A.", "B.", "C.", "D.", "E."))
  names(label.dat) = c("model", "label")
  #then, reorder the models
  mod.dat1$model = factor(mod.dat1$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  #with intra-annual squiggles:
  mod.dat2 <- ddply(mod.dat1, .(model, age), summarize, mean_seroprev = mean(seroprevalence), lci_seroprev= mean(seroprev_lci), uci_seroprev= mean(seroprev_uci), species=unique(species), type=unique(type)) 
  mod.dat2 <- get.avg.age.sero(dat=mod.dat1) #no squiggles
  mod.dat2 = subset(mod.dat2, age<=15)
  #plot together
  # p1 <- ggplot(data=mod.dat2) + geom_line(aes(x=age, y=mean_seroprev), color="purple", size=1) + 
  #  facet_grid(~model) + geom_ribbon(aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), alpha=.3, fill="purple") +
  # geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence), linetype=2) +
  #geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), pch=1) +
  #coord_cartesian(xlim=c(0,15), ylim=c(0,1)) 
  
  #print(p1)
  
  #plot data
  if (species1=="Eidolon dupreanum"){
    
    
    
    p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1,  alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
      geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
      # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
      geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="purple", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
      geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "purple", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
      geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=4.5) +
      theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 12), 
            axis.title.x = element_text(size=11), axis.text = element_text(size=10), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy"), axis.text.y.right = element_text(colour = "navy"),
            legend.text = element_blank(), plot.margin = unit(c(.2,.2,.2,1.3), "cm")) + 
      ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=15), color="navy") + coord_cartesian(xlim = c(0,20), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
      scale_x_discrete(limits = c(0,10,15), labels = c("0", "10", "15")) + geom_segment(data=AIC.df, aes(x=0, xend=20, y=0, yend=0), linetype=2, color="navy")
    
    
  }else{
    mod.dat2 = subset(mod.dat2, age<=11)
    p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=2,  alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
      geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
      # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
      geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="lightgreen", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
      geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "lightgreen", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
      geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=4.5) +
      theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 12), 
            axis.title.x = element_text(size=11), axis.text = element_text(size=10), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy"), axis.text.y.right = element_text(colour = "navy"),
            legend.text = element_blank(), plot.margin = unit(c(.2,.2,.2,1.3), "cm")) + 
      ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=11), color="navy") + coord_cartesian(xlim = c(0,15), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
      scale_x_discrete(limits = c(0,5,10), labels = c("0", "5", "10")) + geom_segment(data=AIC.df, aes(x=0, xend=15, y=0, yend=0), linetype=2, color="navy")
    
  }
  print(p2)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=100, 
           height=65, 
           scale=3, 
           dpi=300)
  }
  
  
}
plot.annual.age.seroprev.matAB <- function(mod.dat, prev.dat, AIC.dat, type1, cutoff1, species1, do.save, filename){
  #first, choose
  prev.dat1 <- subset(prev.dat, species==species1 & type==type1)
  prev.dat1  <- subset(prev.dat1, !is.na(age))
  prev.dat1  <- subset(prev.dat1, fit==1)
  
  
  #get age.seroprev from data
  dat.seroprev = get.seroprev.dat(data=prev.dat1, vis_split=3, cutoff=cutoff1)
  
  mod.dat1 = subset(mod.dat, species==species1 & type==type1 & cutoff==cutoff1)
  mod.dat1$model = as.character(mod.dat1$model)
  mod.dat1$model[mod.dat1$model=="MSIRN-matAB"] = "MSIRN"
  mod.dat1$model[mod.dat1$model=="MSIRNR-matAB"] = "MSIRNR"
  mod.dat1 = subset(mod.dat1, model!="MSIRN-matSus")
  mod.dat1 = subset(mod.dat1, model!="MSIRNR-matSus")
  
  AIC.dat1 = subset(AIC.dat, species==species1 & type==type1 & cutoff==cutoff1)
  AIC.dat1$model = as.character(AIC.dat1$model)
  AIC.dat1$model[AIC.dat1$model=="MSIRN-matAB"] = "MSIRN"
  AIC.dat1$model[AIC.dat1$model=="MSIRNR-matAB"] = "MSIRNR"
  AIC.dat1 = subset(AIC.dat1, model!="MSIRN-matSus")
  AIC.dat1 = subset(AIC.dat1, model!="MSIRNR-matSus")
  #add an "age" value to fit.dat for plotting
  AIC.dat1$age <- NA
  if(species1=="Eidolon dupreanum"){
    AIC.dat1$age<- 18.2
  }else{
    AIC.dat1$age <- 13.5  
  }
  
  
  #and get deltaAIC relative to seroprev
  AIC.df <- subset(AIC.dat1, cutoff==cutoff1)
  min.AIC = min(AIC.df$AIC)
  AIC.df$deltaAIC =AIC.df$AIC-min.AIC
  
  
  max.a = ceiling(max(dat.seroprev$age_year))
  #if(cutoff1=="mean" | cutoff1=="lci"){
  # ymax =.5
  #}else if (cutoff1=="uci"){
  # ymax=
  #}
  
  #AIC.mult <- 1/(max(AIC.df$deltaAIC) + 1)
  #AIC.df$relAIC <- AIC.df$deltaAIC*AIC.mult
  AIC.mult <- ceiling(max(AIC.df$deltaAIC))
  AIC.df$relAIC <- AIC.df$deltaAIC/AIC.mult
  AIC.df$model <- factor(AIC.df$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  AIC.df$outline <- NA
  AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
  AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
  AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
  
  
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 16, "EBOV" = 17)
  
  label.dat = cbind.data.frame(c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"), c("A.", "B.", "C.", "D.", "E."))
  names(label.dat) = c("model", "label")
  #then, reorder the models
  mod.dat1$model = factor(mod.dat1$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  #with intra-annual squiggles:
  mod.dat2 <- ddply(mod.dat1, .(model, age), summarize, mean_seroprev = mean(seroprevalence), lci_seroprev= mean(seroprev_lci), uci_seroprev= mean(seroprev_uci), species=unique(species), type=unique(type)) 
  mod.dat2 <- get.avg.age.sero(dat=mod.dat1) #no squiggles
  mod.dat2 = subset(mod.dat2, age<=15)
  #plot together
  # p1 <- ggplot(data=mod.dat2) + geom_line(aes(x=age, y=mean_seroprev), color="purple", size=1) + 
  #  facet_grid(~model) + geom_ribbon(aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), alpha=.3, fill="purple") +
  # geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence), linetype=2) +
  #geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), pch=1) +
  #coord_cartesian(xlim=c(0,15), ylim=c(0,1)) 
  
  #print(p1)
  
  #plot data
  if (species1=="Eidolon dupreanum"){
    
  
  
  p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1,  alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
    geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
    # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
    geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="purple", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
    geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "purple", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
    geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=4.5) +
    theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 12), 
          axis.title.x = element_text(size=11), axis.text = element_text(size=10), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy"), axis.text.y.right = element_text(colour = "navy"),
          legend.text = element_blank(), plot.margin = unit(c(.2,.2,.2,1.3), "cm")) + 
    ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=15), color="navy") + coord_cartesian(xlim = c(0,20), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
    scale_x_discrete(limits = c(0,10,15), labels = c("0", "10", "15")) + geom_segment(data=AIC.df, aes(x=0, xend=20, y=0, yend=0), linetype=2, color="navy")
  
  
  }else{
    mod.dat2 = subset(mod.dat2, age<=11)
    p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=2,  alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
      geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
      # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
      geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="lightgreen", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
      geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "lightgreen", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
      geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=4.5) +
      theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 12), 
            axis.title.x = element_text(size=11), axis.text = element_text(size=10), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy"), axis.text.y.right = element_text(colour = "navy"),
            legend.text = element_blank(), plot.margin = unit(c(.2,.2,.2,1.3), "cm")) + 
      ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=11), color="navy") + coord_cartesian(xlim = c(0,15), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
      scale_x_discrete(limits = c(0,5,10), labels = c("0", "5", "10")) + geom_segment(data=AIC.df, aes(x=0, xend=15, y=0, yend=0), linetype=2, color="navy")
    
  }
  print(p2)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=100, 
           height=65, 
           scale=3, 
           dpi=300)
  }
  
  
}
plot.annual.age.seroprev.matAB.squiggle <- function(mod.dat, prev.dat, AIC.dat, type1, cutoff1, species1, do.save, filename){
  #first, choose
  prev.dat1 <- subset(prev.dat, species==species1 & type==type1)
  prev.dat1  <- subset(prev.dat1, !is.na(age))
  prev.dat1  <- subset(prev.dat1, fit==1)
  
  
  #get age.seroprev from data
  dat.seroprev = get.seroprev.dat(data=prev.dat1, vis_split=3, cutoff=cutoff1)
  
  mod.dat1 = subset(mod.dat, species==species1 & type==type1 & cutoff==cutoff1)
  mod.dat1$model[mod.dat1$model=="MSIRN-matAB"] = "MSIRN"
  mod.dat1$model[mod.dat1$model=="MSIRNR-matAB"] = "MSIRNR"
  mod.dat1 = subset(mod.dat1, model!="MSIRN-matSus")
  mod.dat1 = subset(mod.dat1, model!="MSIRNR-matSus")
  
  AIC.dat1 = subset(AIC.dat, species==species1 & type==type1 & cutoff==cutoff1)
  AIC.dat1$model[AIC.dat1$model=="MSIRN-matAB"] = "MSIRN"
  AIC.dat1$model[AIC.dat1$model=="MSIRNR-matAB"] = "MSIRNR"
  AIC.dat1 = subset(AIC.dat1, model!="MSIRN-matSus")
  AIC.dat1 = subset(AIC.dat1, model!="MSIRNR-matSus")
  #add an "age" value to fit.dat for plotting
  AIC.dat1$age <- NA
  if(species1=="Eidolon dupreanum"){
    AIC.dat1$age<- 18.2
  }else{
    AIC.dat1$age <- 13.5  
  }
  
  
  #and get deltaAIC relative to seroprev
  AIC.df <- subset(AIC.dat1, cutoff==cutoff1)
  min.AIC = min(AIC.df$AIC)
  AIC.df$deltaAIC =AIC.df$AIC-min.AIC
  
  
  max.a = ceiling(max(dat.seroprev$age_year))
  #if(cutoff1=="mean" | cutoff1=="lci"){
  # ymax =.5
  #}else if (cutoff1=="uci"){
  # ymax=
  #}
  
  #AIC.mult <- 1/(max(AIC.df$deltaAIC) + 1)
  #AIC.df$relAIC <- AIC.df$deltaAIC*AIC.mult
  AIC.mult <- ceiling(max(AIC.df$deltaAIC))
  AIC.df$relAIC <- AIC.df$deltaAIC/AIC.mult
  AIC.df$model <- factor(AIC.df$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  AIC.df$outline <- NA
  AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
  AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
  AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
  
  
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 16, "EBOV" = 17)
  
  
  #then, reorder the models
  mod.dat1$model = factor(mod.dat1$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  #with intra-annual squiggles:
  mod.dat2 <- ddply(mod.dat1, .(model, age), summarize, mean_seroprev = mean(seroprevalence), lci_seroprev= mean(seroprev_lci), uci_seroprev= mean(seroprev_uci), species=unique(species), type=unique(type)) 
  #mod.dat2 <- get.avg.age.sero(dat=mod.dat1) #no squiggles
  mod.dat2 = subset(mod.dat2, age<=15)
  mod.dat2$model = factor(mod.dat2$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  #plot together
  # p1 <- ggplot(data=mod.dat2) + geom_line(aes(x=age, y=mean_seroprev), color="purple", size=1) + 
  #  facet_grid(~model) + geom_ribbon(aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), alpha=.3, fill="purple") +
  # geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence), linetype=2) +
  #geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), pch=1) +
  #coord_cartesian(xlim=c(0,15), ylim=c(0,1)) 
  
  #print(p1)
  label.dat = cbind.data.frame(c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"), c("(a)", "(b)", "(c)", "(d)", "(e)"))
  names(label.dat) = c("model", "label")
  label.dat$model <- factor(label.dat$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  #plot data
  if (species1=="Eidolon dupreanum"){
    
    p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1,  alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
      geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
      # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
      geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="purple", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
      geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "purple", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =18, size =4, color="navy", show.legend = FALSE) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 5, color="navy", size=5, show.legend = FALSE) +
      geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=5) +
      theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 16), 
            axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy", size=16), axis.text.y.right = element_text(colour = "navy", size=12),
            legend.text = element_blank(), plot.margin = unit(c(.1,.1,.1,.1), "cm")) + 
      ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=15), color="navy") + coord_cartesian(xlim = c(0,20), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
      scale_x_discrete(limits = c(0,10,15), labels = c("0", "10", "15")) + geom_segment(data=AIC.df, aes(x=0, xend=20, y=0, yend=0), linetype=2, color="navy")
    
  }else{
    mod.dat2 = subset(mod.dat2, age<=11)
    p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=2,  alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
      geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
      # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
      geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="lightgreen", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
      geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "lightgreen", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
      geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=4.5) +
      theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 16), 
            axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy", size=16), axis.text.y.right = element_text(colour = "navy", size=12),
            legend.text = element_blank(), plot.margin = unit(c(.1,.1,.1,.1), "cm")) + 
      ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=11), color="navy") + coord_cartesian(xlim = c(0,15), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
      scale_x_discrete(limits = c(0,5,10), labels = c("0", "5", "10")) + geom_segment(data=AIC.df, aes(x=0, xend=15, y=0, yend=0), linetype=2, color="navy")
    
  }
  
  
  print(p2)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=100, 
           height=40, 
           scale=3, 
           dpi=300)
  }
  
  
}



#for the main text:
plot.annual.age.seroprev.all.mods(mod.dat=combine.mod.all,
                                  prev.dat=dat.tot, 
                                  AIC.dat=combine.fit.all,
                                  cutoff1 = "mean",
                                  type1="NIVG",
                                  species1="Eidolon dupreanum",
                                  do.save=F,
                                  filename=NA)



plot.annual.age.seroprev.all.mods.squiggle(mod.dat=combine.mod.all,
                                  prev.dat=dat.tot, 
                                  AIC.dat=combine.fit.all,
                                  cutoff1 = "mean",
                                  type1="NIVG",
                                  species1="Eidolon dupreanum",
                                  do.save=F,
                                  filename=NA)



plot.annual.age.seroprev.matSus(mod.dat=combine.mod.all,
                                  prev.dat=dat.tot, 
                                  AIC.dat=combine.fit.all,
                                  type1="NIVG",
                                  cutoff1 = "mean",
                                  species1="Eidolon dupreanum",
                                  do.save=F,
                                  filename=NA)

plot.annual.age.seroprev.matAB.squiggle(mod.dat=combine.mod.all,
                                prev.dat=dat.tot, 
                                AIC.dat=combine.fit.all,
                                type1="NIVG",
                                cutoff1 = "mean",
                                species1="Eidolon dupreanum",
                                do.save=T,
                                filename=paste0(homewd, "/JAE-scripts-figures/figures/Fig4.pdf"))

plot.annual.age.seroprev.matAB.squiggle(mod.dat=combine.mod.all,
                                        prev.dat=dat.tot, 
                                        AIC.dat=combine.fit.all,
                                        type1="NIVG",
                                        cutoff1 = "mean",
                                        species1="Eidolon dupreanum",
                                        do.save=T,
                                        filename=paste0(homewd, "/JAE-scripts-figures/figures/Fig4.tiff"))

# 
# plot.annual.age.seroprev.all.mods.squiggle(mod.dat=combine.mod.all,
#                                            prev.dat=dat.tot, 
#                                            AIC.dat=combine.fit.all,
#                                            cutoff1 = "mean",
#                                            type1="NIVG",
#                                            species1="Eidolon dupreanum",
#                                            do.save=T,
#                                            filename=paste0(homewd, "/JAE-scripts-figures/figures/FigS10.pdf"))
# 
# plot.annual.age.seroprev.matAB.squiggle(mod.dat=combine.mod.all,
#                                         prev.dat=dat.tot, 
#                                         AIC.dat=combine.fit.all,
#                                         type1="EBOV",
#                                         cutoff1 = "mean",
#                                         species1="Pteropus rufus",
#                                         do.save=T,
#                                         filename=paste0(homewd, "/JAE-scripts-figures/figures/FigS11.pdf"))
# 





#and the combination at lci and uci
plot.annual.age.seroprev.matSus.both <- function(mod.dat, prev.dat, AIC.dat, cutoff1, do.save, filename){
  #first, choose
  #prev.dat1 <- subset(prev.dat, species==species1 & type==type1)
  
  prev.dat1  <- subset(prev.dat, !is.na(age))
  prev.dat1  <- subset(prev.dat1, fit==1)
  prev.dat.Eid = subset(prev.dat1, species=="Eidolon dupreanum" & type=="NIVG")
  prev.dat.Pter = subset(prev.dat1, species=="Pteropus rufus" & type=="EBOV")
  
  #get age.seroprev from data
  dat.seroprev.Eid = get.seroprev.dat(data=prev.dat.Eid, vis_split=3, cutoff=cutoff1)
  dat.seroprev.Pter= get.seroprev.dat(data=prev.dat.Pter, vis_split=3, cutoff=cutoff1)
  dat.seroprev = rbind(dat.seroprev.Eid, dat.seroprev.Pter)
  
  mod.dat1 = subset(mod.dat,  cutoff==cutoff1)
  mod.dat1$model = as.character(mod.dat1$model)
  mod.dat1$model[mod.dat1$model=="MSIRN-matSus"] = "MSIRN"
  mod.dat1$model[mod.dat1$model=="MSIRNR-matSus"] = "MSIRNR"
  mod.dat1 = subset(mod.dat1, model!="MSIRN-matAB")
  mod.dat1 = subset(mod.dat1, model!="MSIRNR-matAB")
  #mod.dat.Eid = subset(mod.dat1, species=="Eidolon dupreanum")
  #mod.dat.Pter = subset(mod.dat1, species!="Eidolon dupreanum")
  
  AIC.dat1 = subset(AIC.dat, cutoff==cutoff1)
  AIC.dat1$model = as.character(AIC.dat1$model)
  AIC.dat1$model[AIC.dat1$model=="MSIRN-matSus"] = "MSIRN"
  AIC.dat1$model[AIC.dat1$model=="MSIRNR-matSus"] = "MSIRNR"
  AIC.dat1 = subset(AIC.dat1, model!="MSIRN-matAB")
  AIC.dat1 = subset(AIC.dat1, model!="MSIRNR-matAB")
  #add an "age" value to fit.dat for plotting

  AIC.dat1$age<- 18.2
  
  
  
  #and get deltaAIC relative to seroprev
  AIC.df <- subset(AIC.dat1, cutoff==cutoff1)
  min.AIC.Eid = min(AIC.df$AIC[AIC.df$species=="Eidolon dupreanum"])
  min.AIC.Pter = min(AIC.df$AIC[AIC.df$species=="Pteropus rufus"])
  AIC.df$deltaAIC =NA
  AIC.df$deltaAIC[AIC.df$species=="Eidolon dupreanum"] = AIC.df$AIC[AIC.df$species=="Eidolon dupreanum"]-min.AIC.Eid
  AIC.df$deltaAIC[AIC.df$species=="Pteropus rufus"] = AIC.df$AIC[AIC.df$species=="Pteropus rufus"]-min.AIC.Pter
  
  
  max.a = ceiling(max(dat.seroprev$age_year))
  #if(cutoff1=="mean" | cutoff1=="lci"){
  # ymax =.5
  #}else if (cutoff1=="uci"){
  # ymax=
  #}
  
  AIC.mult <- ceiling(max(AIC.df$deltaAIC))
  AIC.df$relAIC <- AIC.df$deltaAIC/AIC.mult
  AIC.df$model <- factor(AIC.df$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  AIC.df$outline <- NA
  AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
  AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
  AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
  
  
 
  
  label.dat = cbind.data.frame(rep(c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"), 2), c("A.", "B.", "C.", "D.", "E.", "F.", "G.", "H.", "I.", "J."), rep(c("Eidolon dupreanum", "Pteropus rufus"), each=5))
  names(label.dat) = c("model", "label", "species")
  #then, reorder the models
  mod.dat1$model = factor(mod.dat1$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  #with intra-annual squiggles:
  mod.dat2 <- ddply(mod.dat1, .(species, model, age), summarize, mean_seroprev = mean(seroprevalence), lci_seroprev= mean(seroprev_lci), uci_seroprev= mean(seroprev_uci), species=unique(species), type=unique(type)) 
  mod.dat2.Eid <- get.avg.age.sero(dat=subset(mod.dat1, species=="Eidolon dupreanum")) #no squiggles
  mod.dat2.Pter <- get.avg.age.sero(dat=subset(mod.dat1, species=="Pteropus rufus"))
  mod.dat2 = rbind(mod.dat2.Eid, mod.dat2.Pter)
  mod.dat2 = subset(mod.dat2, age<=15)
  #plot together
  # p1 <- ggplot(data=mod.dat2) + geom_line(aes(x=age, y=mean_seroprev), color="purple", size=1) + 
  #  facet_grid(~model) + geom_ribbon(aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), alpha=.3, fill="purple") +
  # geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence), linetype=2) +
  #geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), pch=1) +
  #coord_cartesian(xlim=c(0,15), ylim=c(0,1)) 
  
  #print(p1)
  
  #plot data
 
  color.comp = c("Pteropus rufus" = "lightgreen", "Eidolon dupreanum" = "purple")
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 1, "EBOV" = 2)
    
    p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count, shape=type), color ="black", alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
      geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
      # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
      geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev, color=species), size=1.5, show.legend = FALSE) + scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) + scale_fill_manual(values=color.comp) +
      geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev, fill=species),  alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
      geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
      geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=4.5) +
      theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 12), 
            axis.title.x = element_text(size=11), axis.text = element_text(size=10), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy"), axis.text.y.right = element_text(colour = "navy"),
            legend.text = element_blank(), plot.margin = unit(c(.2,.2,.2,1.3), "cm")) + 
      ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=15), color="navy") + coord_cartesian(xlim = c(0,20), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
      scale_x_discrete(limits = c(0,10,15), labels = c("0", "10", "15")) + geom_segment(data=AIC.df, aes(x=0, xend=20, y=0, yend=0), linetype=2, color="navy")
    
  print(p2)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=100, 
           height=65, 
           scale=3, 
           dpi=300)
  }
  
  
}
plot.annual.age.seroprev.matSus.both.squiggle <- function(mod.dat, prev.dat, AIC.dat, cutoff1, do.save, filename){
  #first, choose
  #prev.dat1 <- subset(prev.dat, species==species1 & type==type1)
  
  prev.dat1  <- subset(prev.dat, !is.na(age))
  prev.dat1  <- subset(prev.dat1, fit==1)
  prev.dat.Eid = subset(prev.dat1, species=="Eidolon dupreanum" & type=="NIVG")
  prev.dat.Pter = subset(prev.dat1, species=="Pteropus rufus" & type=="EBOV")
  
  #get age.seroprev from data
  dat.seroprev.Eid = get.seroprev.dat(data=prev.dat.Eid, vis_split=3, cutoff=cutoff1)
  dat.seroprev.Pter= get.seroprev.dat(data=prev.dat.Pter, vis_split=3, cutoff=cutoff1)
  dat.seroprev = rbind(dat.seroprev.Eid, dat.seroprev.Pter)
  
  mod.dat1 = subset(mod.dat,  cutoff==cutoff1)
  mod.dat1$model = as.character(mod.dat1$model)
  mod.dat1$model[mod.dat1$model=="MSIRN-matSus"] = "MSIRN"
  mod.dat1$model[mod.dat1$model=="MSIRNR-matSus"] = "MSIRNR"
  mod.dat1 = subset(mod.dat1, model!="MSIRN-matAB")
  mod.dat1 = subset(mod.dat1, model!="MSIRNR-matAB")
  #mod.dat.Eid = subset(mod.dat1, species=="Eidolon dupreanum")
  #mod.dat.Pter = subset(mod.dat1, species!="Eidolon dupreanum")
  
  AIC.dat1 = subset(AIC.dat, cutoff==cutoff1)
  AIC.dat1$model = as.character(AIC.dat1$model)
  AIC.dat1$model[AIC.dat1$model=="MSIRN-matSus"] = "MSIRN"
  AIC.dat1$model[AIC.dat1$model=="MSIRNR-matSus"] = "MSIRNR"
  AIC.dat1 = subset(AIC.dat1, model!="MSIRN-matAB")
  AIC.dat1 = subset(AIC.dat1, model!="MSIRNR-matAB")
  #add an "age" value to fit.dat for plotting
  
  AIC.dat1$age<- 18.2
  
  
  
  #and get deltaAIC relative to seroprev
  AIC.df <- subset(AIC.dat1, cutoff==cutoff1)
  min.AIC.Eid = min(AIC.df$AIC[AIC.df$species=="Eidolon dupreanum"])
  min.AIC.Pter = min(AIC.df$AIC[AIC.df$species=="Pteropus rufus"])
  AIC.df$deltaAIC =NA
  AIC.df$deltaAIC[AIC.df$species=="Eidolon dupreanum"] = AIC.df$AIC[AIC.df$species=="Eidolon dupreanum"]-min.AIC.Eid
  AIC.df$deltaAIC[AIC.df$species=="Pteropus rufus"] = AIC.df$AIC[AIC.df$species=="Pteropus rufus"]-min.AIC.Pter
  
  
  max.a = ceiling(max(dat.seroprev$age_year))
  #if(cutoff1=="mean" | cutoff1=="lci"){
  # ymax =.5
  #}else if (cutoff1=="uci"){
  # ymax=
  #}
  
  AIC.mult <- ceiling(max(AIC.df$deltaAIC))
  AIC.df$relAIC <- AIC.df$deltaAIC/AIC.mult
  AIC.df$model <- factor(AIC.df$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  AIC.df$outline <- NA
  AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
  AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
  AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
  
  
  
  
  label.dat = cbind.data.frame(rep(c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"), 2), c("A.", "B.", "C.", "D.", "E.", "F.", "G.", "H.", "I.", "J."), rep(c("Eidolon dupreanum", "Pteropus rufus"), each=5))
  names(label.dat) = c("model", "label", "species")
  #then, reorder the models
  mod.dat1$model = factor(mod.dat1$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  #with intra-annual squiggles:
  mod.dat2 <- ddply(mod.dat1, .(species, model, age), summarize, mean_seroprev = mean(seroprevalence), lci_seroprev= mean(seroprev_lci), uci_seroprev= mean(seroprev_uci), species=unique(species), type=unique(type)) 
  #mod.dat2.Eid <- get.avg.age.sero(dat=subset(mod.dat1, species=="Eidolon dupreanum")) #no squiggles
  #mod.dat2.Pter <- get.avg.age.sero(dat=subset(mod.dat1, species=="Pteropus rufus"))
  #mod.dat2 = rbind(mod.dat2.Eid, mod.dat2.Pter)
  mod.dat2 = subset(mod.dat2, age<=15)
  #plot together
  # p1 <- ggplot(data=mod.dat2) + geom_line(aes(x=age, y=mean_seroprev), color="purple", size=1) + 
  #  facet_grid(~model) + geom_ribbon(aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), alpha=.3, fill="purple") +
  # geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence), linetype=2) +
  #geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), pch=1) +
  #coord_cartesian(xlim=c(0,15), ylim=c(0,1)) 
  
  #print(p1)
  
  #plot data
  
  color.comp = c("Pteropus rufus" = "lightgreen", "Eidolon dupreanum" = "purple")
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 1, "EBOV" = 2)
  
  p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count, shape=type), color ="black", alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
    geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
    # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
    geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev, color=species), size=1.5, show.legend = FALSE) + scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) + scale_fill_manual(values=color.comp) +
    geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev, fill=species),  alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
    geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=4.5) +
    theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 12), 
          axis.title.x = element_text(size=11), axis.text = element_text(size=10), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy"), axis.text.y.right = element_text(colour = "navy"),
          legend.text = element_blank(), plot.margin = unit(c(.2,.2,.2,1.3), "cm")) + 
    ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=15), color="navy") + coord_cartesian(xlim = c(0,20), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
    scale_x_discrete(limits = c(0,10,15), labels = c("0", "10", "15")) + geom_segment(data=AIC.df, aes(x=0, xend=20, y=0, yend=0), linetype=2, color="navy")
  
  print(p2)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=100, 
           height=65, 
           scale=3, 
           dpi=300)
  }
  
  
}
plot.annual.age.seroprev.matAB.both <- function(mod.dat, prev.dat, AIC.dat, cutoff1, do.save, filename){
  #first, choose
  #prev.dat1 <- subset(prev.dat, species==species1 & type==type1)
  
  prev.dat1  <- subset(prev.dat, !is.na(age))
  prev.dat1  <- subset(prev.dat1, fit==1)
  prev.dat.Eid = subset(prev.dat1, species=="Eidolon dupreanum" & type=="NIVG")
  prev.dat.Pter = subset(prev.dat1, species=="Pteropus rufus" & type=="EBOV")
  
  #get age.seroprev from data
  dat.seroprev.Eid = get.seroprev.dat(data=prev.dat.Eid, vis_split=3, cutoff=cutoff1)
  dat.seroprev.Pter= get.seroprev.dat(data=prev.dat.Pter, vis_split=3, cutoff=cutoff1)
  dat.seroprev = rbind(dat.seroprev.Eid, dat.seroprev.Pter)
  
  mod.dat1 = subset(mod.dat,  cutoff==cutoff1)
  mod.dat1$model = as.character(mod.dat1$model)
  mod.dat1$model[mod.dat1$model=="MSIRN-matAB"] = "MSIRN"
  mod.dat1$model[mod.dat1$model=="MSIRNR-matAB"] = "MSIRNR"
  mod.dat1 = subset(mod.dat1, model!="MSIRN-matSus")
  mod.dat1 = subset(mod.dat1, model!="MSIRNR-matSus")
  #mod.dat.Eid = subset(mod.dat1, species=="Eidolon dupreanum")
  #mod.dat.Pter = subset(mod.dat1, species!="Eidolon dupreanum")
  
  AIC.dat1 = subset(AIC.dat, cutoff==cutoff1)
  AIC.dat1$model = as.character(AIC.dat1$model)
  AIC.dat1$model[AIC.dat1$model=="MSIRN-matAB"] = "MSIRN"
  AIC.dat1$model[AIC.dat1$model=="MSIRNR-matAB"] = "MSIRNR"
  AIC.dat1 = subset(AIC.dat1, model!="MSIRN-matSus")
  AIC.dat1 = subset(AIC.dat1, model!="MSIRNR-matSus")
  #add an "age" value to fit.dat for plotting
  
  AIC.dat1$age<- 18.2
  
  
  
  #and get deltaAIC relative to seroprev
  AIC.df <- subset(AIC.dat1, cutoff==cutoff1)
  min.AIC.Eid = min(AIC.df$AIC[AIC.df$species=="Eidolon dupreanum"])
  min.AIC.Pter = min(AIC.df$AIC[AIC.df$species=="Pteropus rufus"])
  AIC.df$deltaAIC =NA
  AIC.df$deltaAIC[AIC.df$species=="Eidolon dupreanum"] = AIC.df$AIC[AIC.df$species=="Eidolon dupreanum"]-min.AIC.Eid
  AIC.df$deltaAIC[AIC.df$species=="Pteropus rufus"] = AIC.df$AIC[AIC.df$species=="Pteropus rufus"]-min.AIC.Pter
  
  
  max.a = ceiling(max(dat.seroprev$age_year))
  #if(cutoff1=="mean" | cutoff1=="lci"){
  # ymax =.5
  #}else if (cutoff1=="uci"){
  # ymax=
  #}
  
  AIC.mult <- ceiling(max(AIC.df$deltaAIC))
  AIC.df$relAIC <- AIC.df$deltaAIC/AIC.mult
  AIC.df$model <- factor(AIC.df$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  AIC.df$outline <- NA
  AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
  AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
  AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
  
  
  
  
  label.dat = cbind.data.frame(rep(c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"), 2), c("A.", "B.", "C.", "D.", "E.", "F.", "G.", "H.", "I.", "J."), rep(c("Eidolon dupreanum", "Pteropus rufus"), each=5))
  names(label.dat) = c("model", "label", "species")
  #then, reorder the models
  mod.dat1$model = factor(mod.dat1$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  #with intra-annual squiggles:
  mod.dat2 <- ddply(mod.dat1, .(species, model, age), summarize, mean_seroprev = mean(seroprevalence), lci_seroprev= mean(seroprev_lci), uci_seroprev= mean(seroprev_uci), species=unique(species), type=unique(type)) 
  mod.dat2.Eid <- get.avg.age.sero(dat=subset(mod.dat1, species=="Eidolon dupreanum")) #no squiggles
  mod.dat2.Pter <- get.avg.age.sero(dat=subset(mod.dat1, species=="Pteropus rufus"))
  mod.dat2 = rbind(mod.dat2.Eid, mod.dat2.Pter)
  mod.dat2 = subset(mod.dat2, age<=15)
  #plot together
  # p1 <- ggplot(data=mod.dat2) + geom_line(aes(x=age, y=mean_seroprev), color="purple", size=1) + 
  #  facet_grid(~model) + geom_ribbon(aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), alpha=.3, fill="purple") +
  # geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence), linetype=2) +
  #geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), pch=1) +
  #coord_cartesian(xlim=c(0,15), ylim=c(0,1)) 
  
  #print(p1)
  
  #plot data
  
  color.comp = c("Pteropus rufus" = "lightgreen", "Eidolon dupreanum" = "purple")
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 1, "EBOV" = 2)
  
  p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count, shape=type), color ="black", alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
    geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
    # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
    geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev, color=species), size=1.5, show.legend = FALSE) + scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) + scale_fill_manual(values=color.comp) +
    geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev, fill=species),  alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
    geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=5) +
    theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 16), 
          axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy", size=16), axis.text.y.right = element_text(colour = "navy", size=12),
          legend.text = element_blank(), plot.margin =unit(c(.1,.1,.1,.1), "cm")) + 
    ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=15), color="navy") + coord_cartesian(xlim = c(0,20), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
    scale_x_discrete(limits = c(0,10,15), labels = c("0", "10", "15")) + geom_segment(data=AIC.df, aes(x=0, xend=20, y=0, yend=0), linetype=2, color="navy")
  
  print(p2)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=100, 
           height=65, 
           scale=3, 
           dpi=300)
  }
  
  
}
plot.annual.age.seroprev.matAB.both.squiggle <- function(mod.dat, prev.dat, AIC.dat, cutoff1, do.save, filename){
  #first, choose
  #prev.dat1 <- subset(prev.dat, species==species1 & type==type1)
  
  prev.dat1  <- subset(prev.dat, !is.na(age))
  prev.dat1  <- subset(prev.dat1, fit==1)
  prev.dat.Eid = subset(prev.dat1, species=="Eidolon dupreanum" & type=="NIVG")
  prev.dat.Pter = subset(prev.dat1, species=="Pteropus rufus" & type=="EBOV")
  
  #get age.seroprev from data
  dat.seroprev.Eid = get.seroprev.dat(data=prev.dat.Eid, vis_split=3, cutoff=cutoff1)
  dat.seroprev.Pter= get.seroprev.dat(data=prev.dat.Pter, vis_split=3, cutoff=cutoff1)
  dat.seroprev = rbind(dat.seroprev.Eid, dat.seroprev.Pter)
  
  mod.dat1 = subset(mod.dat,  cutoff==cutoff1)
  mod.dat1$model = as.character(mod.dat1$model)
  mod.dat1$model[mod.dat1$model=="MSIRN-matAB"] = "MSIRN"
  mod.dat1$model[mod.dat1$model=="MSIRNR-matAB"] = "MSIRNR"
  mod.dat1 = subset(mod.dat1, model!="MSIRN-matSus")
  mod.dat1 = subset(mod.dat1, model!="MSIRNR-matSus")
  #mod.dat.Eid = subset(mod.dat1, species=="Eidolon dupreanum")
  #mod.dat.Pter = subset(mod.dat1, species!="Eidolon dupreanum")
  
  AIC.dat1 = subset(AIC.dat, cutoff==cutoff1)
  AIC.dat1$model = as.character(AIC.dat1$model)
  AIC.dat1$model[AIC.dat1$model=="MSIRN-matAB"] = "MSIRN"
  AIC.dat1$model[AIC.dat1$model=="MSIRNR-matAB"] = "MSIRNR"
  AIC.dat1 = subset(AIC.dat1, model!="MSIRN-matSus")
  AIC.dat1 = subset(AIC.dat1, model!="MSIRNR-matSus")
  #add an "age" value to fit.dat for plotting
  
  AIC.dat1$age<- 18.2
  
  
  
  #and get deltaAIC relative to seroprev
  AIC.df <- subset(AIC.dat1, cutoff==cutoff1)
  min.AIC.Eid = min(AIC.df$AIC[AIC.df$species=="Eidolon dupreanum"])
  min.AIC.Pter = min(AIC.df$AIC[AIC.df$species=="Pteropus rufus"])
  AIC.df$deltaAIC =NA
  AIC.df$deltaAIC[AIC.df$species=="Eidolon dupreanum"] = AIC.df$AIC[AIC.df$species=="Eidolon dupreanum"]-min.AIC.Eid
  AIC.df$deltaAIC[AIC.df$species=="Pteropus rufus"] = AIC.df$AIC[AIC.df$species=="Pteropus rufus"]-min.AIC.Pter
  
  
  max.a = ceiling(max(dat.seroprev$age_year))
  #if(cutoff1=="mean" | cutoff1=="lci"){
  # ymax =.5
  #}else if (cutoff1=="uci"){
  # ymax=
  #}
  
  AIC.mult <- ceiling(max(AIC.df$deltaAIC))
  AIC.df$relAIC <- AIC.df$deltaAIC/AIC.mult
  AIC.df$model <- factor(AIC.df$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  AIC.df$outline <- NA
  AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
  AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
  AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
  
  
  
  
  label.dat = cbind.data.frame(rep(c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"), 2), c("A.", "B.", "C.", "D.", "E.", "F.", "G.", "H.", "I.", "J."), rep(c("Eidolon dupreanum", "Pteropus rufus"), each=5))
  names(label.dat) = c("model", "label", "species")
  #then, reorder the models
  mod.dat1$model = factor(mod.dat1$model, levels = c("MSIR", "MSRIR", "MSIRS",  "MSIRN", "MSIRNR"))
  #with intra-annual squiggles:
  mod.dat2 <- ddply(mod.dat1, .(species, model, age), summarize, mean_seroprev = mean(seroprevalence), lci_seroprev= mean(seroprev_lci), uci_seroprev= mean(seroprev_uci), species=unique(species), type=unique(type)) 
  #mod.dat2.Eid <- get.avg.age.sero(dat=subset(mod.dat1, species=="Eidolon dupreanum")) #no squiggles
  #mod.dat2.Pter <- get.avg.age.sero(dat=subset(mod.dat1, species=="Pteropus rufus"))
  #mod.dat2 = rbind(mod.dat2.Eid, mod.dat2.Pter)
  mod.dat2 = subset(mod.dat2, age<=15)
  #plot together
  # p1 <- ggplot(data=mod.dat2) + geom_line(aes(x=age, y=mean_seroprev), color="purple", size=1) + 
  #  facet_grid(~model) + geom_ribbon(aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), alpha=.3, fill="purple") +
  # geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence), linetype=2) +
  #geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), pch=1) +
  #coord_cartesian(xlim=c(0,15), ylim=c(0,1)) 
  
  #print(p1)
  
  #plot data
  
  color.comp = c("Pteropus rufus" = "lightgreen", "Eidolon dupreanum" = "purple")
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 1, "EBOV" = 2)
  
  p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count, shape=type), color ="black", alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
    geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
    # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
    geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev, color=species), size=1.5, show.legend = FALSE) + scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) + scale_fill_manual(values=color.comp) +
    geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev, fill=species),  alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
    geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=5) +
    theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 16), 
          axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy", size=16), axis.text.y.right = element_text(colour = "navy", size=12),
          legend.text = element_blank(), plot.margin = unit(c(.1,.1,.1,.1), "cm")) + 
    ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=15), color="navy") + coord_cartesian(xlim = c(0,20), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
    scale_x_discrete(limits = c(0,10,15), labels = c("0", "10", "15")) + geom_segment(data=AIC.df, aes(x=0, xend=20, y=0, yend=0), linetype=2, color="navy")
  
  print(p2)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=100, 
           height=65, 
           scale=3, 
           dpi=300)
  }
  
  
}

plot.annual.age.seroprev.all.mods.both.squiggle <- function(mod.dat, prev.dat, AIC.dat, cutoff1, do.save, filename){
  #first, choose
  #prev.dat1 <- subset(prev.dat, species==species1 & type==type1)
  
  prev.dat1  <- subset(prev.dat, !is.na(age))
  prev.dat1  <- subset(prev.dat1, fit==1)
  prev.dat.Eid = subset(prev.dat1, species=="Eidolon dupreanum" & type=="NIVG")
  prev.dat.Pter = subset(prev.dat1, species=="Pteropus rufus" & type=="EBOV")
  
  #get age.seroprev from data
  dat.seroprev.Eid = get.seroprev.dat(data=prev.dat.Eid, vis_split=3, cutoff=cutoff1)
  dat.seroprev.Pter= get.seroprev.dat(data=prev.dat.Pter, vis_split=3, cutoff=cutoff1)
  dat.seroprev = rbind(dat.seroprev.Eid, dat.seroprev.Pter)
  
  mod.dat1 = subset(mod.dat,  cutoff==cutoff1)
  mod.dat1$model = as.character(mod.dat1$model)
  #mod.dat1$model[mod.dat1$model=="MSIRN-matAB"] = "MSIRN"
  #mod.dat1$model[mod.dat1$model=="MSIRNR-matAB"] = "MSIRNR"
  #mod.dat1 = subset(mod.dat1, model!="MSIRN-matSus")
  #mod.dat1 = subset(mod.dat1, model!="MSIRNR-matSus")
  #mod.dat.Eid = subset(mod.dat1, species=="Eidolon dupreanum")
  #mod.dat.Pter = subset(mod.dat1, species!="Eidolon dupreanum")
  
  AIC.dat1 = subset(AIC.dat, cutoff==cutoff1)
  AIC.dat1$model = as.character(AIC.dat1$model)
  #AIC.dat1$model[AIC.dat1$model=="MSIRN-matAB"] = "MSIRN"
  #AIC.dat1$model[AIC.dat1$model=="MSIRNR-matAB"] = "MSIRNR"
  #AIC.dat1 = subset(AIC.dat1, model!="MSIRN-matSus")
  #AIC.dat1 = subset(AIC.dat1, model!="MSIRNR-matSus")
  #add an "age" value to fit.dat for plotting
  
  AIC.dat1$age<- 18.2
  
  
  
  #and get deltaAIC relative to seroprev
  AIC.df <- subset(AIC.dat1, cutoff==cutoff1)
  min.AIC.Eid = min(AIC.df$AIC[AIC.df$species=="Eidolon dupreanum"])
  min.AIC.Pter = min(AIC.df$AIC[AIC.df$species=="Pteropus rufus"])
  AIC.df$deltaAIC =NA
  AIC.df$deltaAIC[AIC.df$species=="Eidolon dupreanum"] = AIC.df$AIC[AIC.df$species=="Eidolon dupreanum"]-min.AIC.Eid
  AIC.df$deltaAIC[AIC.df$species=="Pteropus rufus"] = AIC.df$AIC[AIC.df$species=="Pteropus rufus"]-min.AIC.Pter
  
  
  max.a = ceiling(max(dat.seroprev$age_year))
  #if(cutoff1=="mean" | cutoff1=="lci"){
  # ymax =.5
  #}else if (cutoff1=="uci"){
  # ymax=
  #}
  
  AIC.mult <- ceiling(max(AIC.df$deltaAIC))
  AIC.df$relAIC <- AIC.df$deltaAIC/AIC.mult
  AIC.df$model <- factor(AIC.df$model, levels = c("MSIR", "MSRIR", "MSIRS", "MSIRN-matAB","MSIRN-matSus", "MSIRNR-matAB",  "MSIRNR-matSus"))
  AIC.df$outline <- NA
  AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
  AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
  AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
  
  
  
  
  label.dat = cbind.data.frame(rep(c("MSIR", "MSRIR", "MSIRS", "MSIRN-matAB","MSIRN-matSus", "MSIRNR-matAB",  "MSIRNR-matSus"), 2), c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)", "(g)"), rep(c("Eidolon dupreanum", "Pteropus rufus"), each=7))
  names(label.dat) = c("model", "label", "species")
  
  label.dat$model = factor(label.dat$model, levels = c("MSIR", "MSRIR", "MSIRS", "MSIRN-matAB","MSIRN-matSus", "MSIRNR-matAB",  "MSIRNR-matSus"))
  
  #then, reorder the models
  mod.dat1$model = factor(mod.dat1$model, levels = c("MSIR", "MSRIR", "MSIRS", "MSIRN-matAB","MSIRN-matSus", "MSIRNR-matAB",  "MSIRNR-matSus"))
  #with intra-annual squiggles:
  mod.dat2 <- ddply(mod.dat1, .(species, model, age), summarize, mean_seroprev = mean(seroprevalence), lci_seroprev= mean(seroprev_lci), uci_seroprev= mean(seroprev_uci), species=unique(species), type=unique(type)) 
  #mod.dat2.Eid <- get.avg.age.sero(dat=subset(mod.dat1, species=="Eidolon dupreanum")) #no squiggles
  #mod.dat2.Pter <- get.avg.age.sero(dat=subset(mod.dat1, species=="Pteropus rufus"))
  #mod.dat2 = rbind(mod.dat2.Eid, mod.dat2.Pter)
  mod.dat2 = subset(mod.dat2, age<=15)
  
  mod.dat2$model = factor(mod.dat2$model, levels = c("MSIR", "MSRIR", "MSIRS", "MSIRN-matAB","MSIRN-matSus", "MSIRNR-matAB",  "MSIRNR-matSus"))
  #plot together
  # p1 <- ggplot(data=mod.dat2) + geom_line(aes(x=age, y=mean_seroprev), color="purple", size=1) + 
  #  facet_grid(~model) + geom_ribbon(aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), alpha=.3, fill="purple") +
  # geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence), linetype=2) +
  #geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), pch=1) +
  #coord_cartesian(xlim=c(0,15), ylim=c(0,1)) 
  
  #print(p1)
  
  #plot data
  
  color.comp = c("Pteropus rufus" = "lightgreen", "Eidolon dupreanum" = "purple")
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 1, "EBOV" = 2)
  
  dat.seroprev1 = dat.seroprev
  dat.seroprev1$type[dat.seroprev1$type=="NIVG"] = "NiV-G"
  dat.seroprev1$type[dat.seroprev1$type=="EBOV"] = "EBOV-Gp"
  dat.seroprev1$type = factor(dat.seroprev1$type, levels = c("NiV-G", "EBOV-Gp"))
  leg.grab <- ggplot() + geom_point(data=dat.seroprev1, aes(x=age_year, y=prevalence,  shape=type), color ="black", alpha=.9) + theme_bw() + facet_grid(species~model) +
    geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
    # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
    geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev, color=species), size=1.5) + scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) + scale_fill_manual(values=color.comp) +
    geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev, fill=species),  alpha=.3) + #scale_fill_manual(values=color.comp) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
    geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=5) +
    theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 16), 
          axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy", size=16), axis.text.y.right = element_text(colour = "navy", size=12),
          legend.direction = "horizontal", legend.position = "top", plot.margin = unit(c(0,0,0,0), "cm"), legend.text = element_text(face = "italic")) + 
    ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=15), color="navy") + coord_cartesian(xlim = c(0,20), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
    scale_x_discrete(limits = c(0,10,15), labels = c("0", "10", "15")) + geom_segment(data=AIC.df, aes(x=0, xend=20, y=0, yend=0), linetype=2, color="navy")
  
  leg1 = get_legend(leg.grab)
  
  p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count, shape=type), color ="black", alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
    geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
    # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
    geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev, color=species), size=1.5, show.legend = FALSE) + scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) + scale_fill_manual(values=color.comp) +
    geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev, fill=species),  alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
    geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
    geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=5) +
    theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 16), 
          axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy", size=16), axis.text.y.right = element_text(colour = "navy", size=12),
          legend.text = element_blank(), plot.margin = unit(c(.1,.1,.1,.1), "cm")) + 
    ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=15), color="navy") + coord_cartesian(xlim = c(0,20), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult,name = expression(Delta~AIC))) +
    scale_x_discrete(limits = c(0,10,15), labels = c("0", "10", "15")) + geom_segment(data=AIC.df, aes(x=0, xend=20, y=0, yend=0), linetype=2, color="navy")
  
  cowplot::plot_grid(leg1, p2, nrow=2, ncol=1, rel_heights = c(.05,1))
  
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=140, 
           height=70, 
           scale=3, 
           dpi=300)
  }
  
  
}


plot.annual.age.seroprev.all.mods.both.squiggle(mod.dat=combine.mod.all,
                                    prev.dat=dat.tot, 
                                    AIC.dat=combine.fit.all,
                                    cutoff1 = "mean",
                                    do.save=T,
                                    filename=paste0(homewd, "/JAE-scripts-figures/figures/FigS10.pdf"))

plot.annual.age.seroprev.all.mods.both.squiggle(mod.dat=combine.mod.all,
                                                prev.dat=dat.tot, 
                                                AIC.dat=combine.fit.all,
                                                cutoff1 = "lci",
                                                do.save=T,
                                                filename=paste0(homewd, "/JAE-scripts-figures/figures/FigS11.pdf"))

plot.annual.age.seroprev.all.mods.both.squiggle(mod.dat=combine.mod.all,
                                                prev.dat=dat.tot, 
                                                AIC.dat=combine.fit.all,
                                                cutoff1 = "uci",
                                                do.save=T,
                                                filename=paste0(homewd, "/JAE-scripts-figures/figures/FigS12.pdf"))

plot.annual.age.seroprev.matAB.both(mod.dat=combine.mod.all,
                                         prev.dat=dat.tot, 
                                         AIC.dat=combine.fit.all,
                                         cutoff1 = "uci",
                                         do.save=F,
                                         filename=NA)

plot.annual.age.seroprev.matAB.both(mod.dat=combine.mod.all,
                                    prev.dat=dat.tot, 
                                    AIC.dat=combine.fit.all,
                                    cutoff1 = "lci",
                                    do.save=F,
                                    filename=NA)


plot.annual.age.seroprev.matAB.both.squiggle(mod.dat=combine.mod.all,
                                    prev.dat=dat.tot, 
                                    AIC.dat=combine.fit.all,
                                    cutoff1 = "uci",
                                    do.save=T,
                                    filename=paste0(homewd, "/JAE-scripts-figures/figures/FigS12.pdf"))

plot.annual.age.seroprev.matAB.both.squiggle(mod.dat=combine.mod.all,
                                    prev.dat=dat.tot, 
                                    AIC.dat=combine.fit.all,
                                    cutoff1 = "lci",
                                    do.save=T,
                                    filename=paste0(homewd, "/JAE-scripts-figures/figures/FigS11.pdf"))
