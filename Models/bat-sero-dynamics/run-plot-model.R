# fit the model to these new 

rm(list=ls())
#set wd
#.libPaths("~/R/x86_64-pc-linux-gnu-library/3.4")

library(epitools)
library(dplyr)
library(plyr)
library(cowplot)
library(deSolve)
library(mgcv)
library(lubridate)
library(matrixcalc)
library(Matrix)
library(dplyr)
library(plyr)
library(ggplot2)
library(mvtnorm)

#load data
homewd = "/Users/carabrook/Developer/bat-sero-dynamics"

load(paste0(homewd, "/data/MSILI.Eid.mean.GhV.nobiwk.Rdata"))
load(paste0(homewd, "/data/MSIRN.Eid.mean.NiV.Rdata"))

mod.dat = MSIRN.EID.mean.NiV[[2]]
AIC.dat = MSIRN.EID.mean.NiV[[1]]
unique(mod.dat$biwk)

p1 <- ggplot(data=mod.dat) + geom_line(aes(x=age, y=seroprevalence, color=biwk)) +
      facet_wrap(~biwk) + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) +
      geom_ribbon(aes(x=age, ymin=seroprev_lci, ymax=seroprev_uci, fill=biwk), alpha=.3)

p1
#and take the mean and plot with the data



dat <- read.csv(file = paste0(homewd, "/data/EidNiVfit.csv"), header = T, stringsAsFactors = F)
names(dat)[names(dat)=="bat_species"] <- "species"
names(dat)[names(dat)=="antigen"] <- "type"
names(dat)[names(dat)=="seropos"] <- "prev"
names(dat)[names(dat)=="sampleid"] <- "sampleID"

#get data in age-seroprev form
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
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev)/length(prev), count=length(prev))  
  }else if (cutoff=="uci"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev_uci)/length(prev_uci), count=length(prev_uci)) 
  }else if (cutoff=="lci"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev_lci)/length(prev_lci), count=length(prev_lci)) 
  }
  
  #you want the midpoint between either end of the age class now
  vect.age.yr = sort(unique(data$age_year))
  dat.sum2$age_plot = NA
  
  for (i in 1:(length(dat.sum2$age_year)-1)){
    dat.sum2$age_plot [i] = ((dat.sum2$age_year[i + 1] - dat.sum2$age_year[i])/2) + dat.sum2$age_year[i]
  }
  dat.sum2$age_plot[length(dat.sum2$age_plot)] =  (ceiling(max(data$age)) - dat.sum2$age_year[length(dat.sum2$age_year)])/2  + dat.sum2$age_year[length(dat.sum2$age_year)]
  
  names(dat.sum2)  <- c("real_age", "prevalence", "count", "age_year") #<- names(dat.sum2.tmp)
  
  # dat.sum2 <-  rbind(dat.sum2.tmp, dat.sum2)
  
  dat.sum3 = dat.sum2
  dat.sum3$class = "seropositive"
  rownames(dat.sum3) <- c()
  
  return(dat.sum3)
  
}

dat.sum = get.seroprev.dat(data=dat[!is.na(dat$age),], vis_split = 3, cutoff = "mean")

  #plot with data
  max.a = ceiling(max(dat.sum$age_year))
  p2 <- ggplot() + geom_point(data=dat.sum, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1, alpha=.9) + theme_bw() +
    geom_line(data=mod.dat, aes(x=age, y=seroprevalence), color = "royalblue") + 
    geom_ribbon(data=mod.dat, aes(x=age, ymin=seroprev_lci, ymax=seroprev_uci), fill = "royalblue", alpha =.3) +
    theme(legend.position = c(.27,.9),  legend.title=element_blank(), legend.spacing = unit(.03,"cm"),  legend.text=element_text(size=10), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, max.a)) + ylab("seroprevalence") + xlab("age(yrs)")  + guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=11), axis.text = element_text(size=10), panel.grid = element_blank())
  print(p2)
  
  #and also plot a version where the model gets binned over .5-year increments
  mod.dat$age_year <- floor(mod.dat$age)
  new.mod <- ddply(mod.dat, .(age_year), summarize, seroprevalence= mean(seroprevalence), seroprev_lci=mean(seroprev_lci), seroprev_uci=mean(seroprev_uci))#, seroprev_lci2=mean(seroprev_lci2), seroprev_uci2=mean(seroprev_uci2))
  
  p3 <- ggplot() + geom_point(data=dat.sum, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1, alpha=.9) + theme_bw() +
    geom_line(data=new.mod, aes(x=age_year, y=seroprevalence), color = "royalblue") + 
    geom_ribbon(data=new.mod, aes(x=age_year, ymin=seroprev_lci, ymax=seroprev_uci), fill = "royalblue", alpha =.3) +
    theme(legend.position = c(.27,.9),  legend.title=element_blank(), legend.spacing = unit(.03,"cm"),  legend.text=element_text(size=10), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
    coord_cartesian(ylim=c(0,1), xlim=c(0, max.a)) + ylab("seroprevalence") + xlab("age(yrs)")  + guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
    theme(axis.title = element_text(size=11), axis.text = element_text(size=10), panel.grid = element_blank())
  print(p3)

  
  
  
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
plot.comp.squiggle <- function(mod.dat, prev.dat, AIC.dat, type1, cutoff1, species1, do.save, filename){
    #first, choose
    prev.dat1 <- subset(prev.dat, species==species1 & type==type1)
    prev.dat1  <- subset(prev.dat1, !is.na(age))
    #prev.dat1  <- subset(prev.dat1, fit==1)
    
    
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
    
    AIC.mult <- 1/(max(AIC.df$deltaAIC) + 1)
    AIC.df$relAIC <- AIC.df$deltaAIC*AIC.mult
    AIC.df <- subset(AIC.df,  model=="MSIRS" | model=="MSIRN" | model=="MSILI")
    AIC.df$model <- factor(AIC.df$model, levels = c(  "MSIRN","MSIRS", "MSILI"))
    AIC.df$outline <- NA
    AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
    AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
    AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
    
    
    shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 16, "EBOV" = 17)
    
    
    #then, reorder the models
    mod.dat1 <- subset(mod.dat1,   model=="MSIRS" | model=="MSIRN" | model=="MSILI")
    mod.dat1$model = factor(mod.dat1$model, levels = c( "MSIRN", "MSIRS", "MSILI"))
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
    label.dat = cbind.data.frame(c( "MSIRN", "MSIRS", "MSILI"), c("A.", "B.", "C."))
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
        geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=5) +
        theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 16), 
              axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy", size=16), axis.text.y.right = element_text(colour = "navy", size=12),
              legend.text = element_blank(), plot.margin = unit(c(.1,.1,.1,.1), "cm")) + 
        ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=15), color="navy") + coord_cartesian(xlim = c(0,20), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult*1000,name = expression(Delta~AIC))) +
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
        ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=11), color="navy") + coord_cartesian(xlim = c(0,15), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult*1000,name = expression(Delta~AIC))) +
        scale_x_discrete(limits = c(0,5,10), labels = c("0", "5", "10")) + geom_segment(data=AIC.df, aes(x=0, xend=15, y=0, yend=0), linetype=2, color="navy")
      
    }
    
    
    print(p2)
    
    if(do.save==TRUE){
      ggsave(file = filename,
             units="mm",  
             width=50, 
             height=30, 
             scale=3, 
             dpi=300)
    }
    
    
  }
plot.comp.squiggle.four <- function(mod.dat, prev.dat, AIC.dat, type1, cutoff1, species1, do.save, filename){
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
  
  AIC.mult <- 1/(max(AIC.df$deltaAIC) + 1)
  AIC.df$relAIC <- AIC.df$deltaAIC/3#AIC.mult
  AIC.df <- subset(AIC.df,  model=="MSIRS" | model=="MSIRN" | model=="MSILI" | model=="MSIR")
  AIC.df$model <- factor(AIC.df$model, levels = c("MSIR", "MSIRS","MSIRN", "MSILI"))
  AIC.df$outline <- NA
  AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
  AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
  AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
  
  
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 16, "EBOV" = 17)
  
  
  #then, reorder the models
  mod.dat1 <- subset(mod.dat1,   model=="MSIRS" | model=="MSIRN" | model=="MSILI"| model=="MSIR")
  mod.dat1$model = factor(mod.dat1$model, levels = c("MSIR",  "MSIRS", "MSIRN", "MSILI"))
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
  label.dat = cbind.data.frame(c("MSIR",  "MSIRS", "MSIRN", "MSILI"), c("A.", "B.", "C.", "D."))
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
      geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=5) +
      theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 16), 
            axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy", size=16), axis.text.y.right = element_text(colour = "navy", size=12),
            legend.text = element_blank(), plot.margin = unit(c(.1,.1,.1,.1), "cm")) + 
      ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=15), color="navy") + coord_cartesian(xlim = c(0,20), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*3,name = expression(Delta~AIC))) +
      scale_x_discrete(limits = c(0,5,10,15), labels = c("0","5", "10", "15")) + geom_segment(data=AIC.df, aes(x=0, xend=20, y=0, yend=0), linetype=2, color="navy")
    
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
      ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=11), color="navy") + coord_cartesian(xlim = c(0,15), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult*1000,name = expression(Delta~AIC))) +
      scale_x_discrete(limits = c(0,5,10), labels = c("0", "5", "10")) + geom_segment(data=AIC.df, aes(x=0, xend=15, y=0, yend=0), linetype=2, color="navy")
    
  }
  
  
  print(p2)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=70, 
           height=30, 
           scale=3, 
           dpi=300)
  }
  
  
}
plot.comp.squiggle.one <- function(mod.dat, prev.dat, AIC.dat, type1, cutoff1, species1, do.save, filename){
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
  
  AIC.mult <- 1/(max(AIC.df$deltaAIC) + 1)
  AIC.df$relAIC <- AIC.df$deltaAIC/3#AIC.mult
  AIC.df <- subset(AIC.df,   model=="MSIRN" )
  AIC.df$model <- factor(AIC.df$model, levels = c("MSIRN"))
  AIC.df$outline <- NA
  AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
  AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
  AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
  
  
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 16, "EBOV" = 17)
  
  
  #then, reorder the models
  mod.dat1 <- subset(mod.dat1,  model=="MSIRN" )
  mod.dat1$model = factor(mod.dat1$model, levels = c("MSIRN"))
  #with intra-annual squiggles:
  mod.dat2 <- ddply(mod.dat1, .(model, age), summarize, mean_seroprev = mean(seroprevalence), lci_seroprev= mean(seroprev_lci), uci_seroprev= mean(seroprev_uci), species=unique(species), type=unique(type)) 
  #mod.dat2 <- get.avg.age.sero(dat=mod.dat1) #no squiggles
  mod.dat2 = subset(mod.dat2, age<=16)
  #plot together
  # p1 <- ggplot(data=mod.dat2) + geom_line(aes(x=age, y=mean_seroprev), color="purple", size=1) + 
  #  facet_grid(~model) + geom_ribbon(aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), alpha=.3, fill="purple") +
  # geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence), linetype=2) +
  #geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), pch=1) +
  #coord_cartesian(xlim=c(0,15), ylim=c(0,1)) 
  
  #print(p1)
  #label.dat = cbind.data.frame(c("MSIR",  "MSIRS", "MSIRN", "MSILI"), c("A.", "B.", "C.", "D."))
  #names(label.dat) = c("model", "label")
  #plot data
  if (species1=="Eidolon dupreanum"){
    
    p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1,  alpha=.9) + theme_bw() + facet_grid(species~model) +
      geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence), linetype=2, color ="black",  size=.2) +
      # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
      geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="purple", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
      geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "purple", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
      #geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =18, size =4, color="navy", show.legend = FALSE) +
      #geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 5, color="navy", size=5, show.legend = FALSE) +
      #geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=5) +
      theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 16), 
            axis.title = element_text(size=16), 
            axis.text = element_text(size=12), 
            legend.title = element_blank(), legend.position = c(.9,.85),
            axis.title.y.right = element_text(colour = "navy", size=16), axis.text.y.right = element_text(colour = "navy", size=12),
            plot.margin = unit(c(.1,.1,.1,.1), "cm")) + 
      ylab("seroprevalence") + xlab("age(yrs)")   + coord_cartesian(xlim = c(0,16), ylim=c(0,1)) +
      scale_x_discrete(limits = c(0,10,15), labels = c("0", "10", "15")) 
    
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
      ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=11), color="navy") + coord_cartesian(xlim = c(0,15), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult*1000,name = expression(Delta~AIC))) +
      scale_x_discrete(limits = c(0,5,10), labels = c("0", "5", "10")) + geom_segment(data=AIC.df, aes(x=0, xend=15, y=0, yend=0), linetype=2, color="navy")
    
  }
  
  
  print(p2)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=50, 
           height=42, 
           scale=3, 
           dpi=300)
  }
  
  
}
#and just data
data.one <- function( prev.dat, type1, cutoff1, species1, do.save, filename){
  #first, choose
  prev.dat1 <- subset(prev.dat, species==species1 & type==type1)
  prev.dat1  <- subset(prev.dat1, !is.na(age))
  prev.dat1  <- subset(prev.dat1, fit==1)
  
  
  #get age.seroprev from data
  dat.seroprev = get.seroprev.dat(data=prev.dat1, vis_split=3, cutoff=cutoff1)
  
  #and get CIs
  dat.seroprev$lci = binom.exact(dat.seroprev$pos, dat.seroprev$count, conf.level = .95)$lower
  dat.seroprev$uci = binom.exact(dat.seroprev$pos, dat.seroprev$count, conf.level = .95)$upper
  
  
  
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 16, "EBOV" = 17)
  
  
  #facet_grid(~model) + geom_ribbon(aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), alpha=.3, fill="purple") +
  # geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence), linetype=2) +
  #geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), pch=1) +
  #coord_cartesian(xlim=c(0,15), ylim=c(0,1)) 
  
  #print(p1)
  #label.dat = cbind.data.frame(c("MSIR",  "MSIRS", "MSIRN", "MSILI"), c("A.", "B.", "C.", "D."))
  
  
  p2 <- ggplot() + 
    geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1) + theme_bw() + #facet_grid(species~model) +
    geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=1.1, show.legend = FALSE) +
    geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=lci, ymax=uci), linetype=3, color ="black",  size=.7, show.legend = FALSE) +
    theme(panel.grid = element_blank(),   strip.background = element_blank(), 
          strip.text.y = element_blank(), strip.text.x = element_text(size = 16), 
          axis.title = element_text(size=16), axis.text = element_text(size=14), legend.title = element_blank(), 
          axis.title.y.right = element_text(colour = "navy", size=16), axis.text.y.right = element_text(colour = "navy", size=12),
          legend.position = c(.9,.85), plot.margin = unit(c(.1,.1,.1,.1), "cm")) + 
    ylab("seroprevalence") + xlab("age(yrs)") + coord_cartesian(xlim = c(0,16), ylim=c(0,1)) +
    scale_x_discrete(limits = c(0,10,15), labels = c("0", "10", "15")) 
  
  
  
  print(p2)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=50, 
           height=42, 
           scale=3, 
           dpi=300)
  }
  
  
}

dat$fit <- 1


plot.comp.squiggle.one(mod.dat=mod.dat,
                   prev.dat=dat, 
                   AIC.dat=AIC.dat,
                   type1="NiV",
                   cutoff1 = "mean",
                   species1="Eidolon dupreanum",
                   do.save=T,
                   filename=paste0(homewd,"/figures/Eid-NiV-age-seroprev.png"))


data.one(prev.dat=dat, 
         type1="NiV",
         cutoff1 = "mean",
         species1="Eidolon dupreanum",
         do.save=T,
         filename=paste0(homewd,"/figures/Eid-NiV-age-seroprev-data.png"))

