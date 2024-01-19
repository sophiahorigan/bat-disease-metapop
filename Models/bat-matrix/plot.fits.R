rm(list=ls())

setwd("/Users/caraebrook/Documents/R/R_repositories/bat-matrix")

library(plyr)
library(dplyr)
library(ggplot2)
library(epitools)

load("seasonal_prev.Rdata")
load("refit.out.MSIRN.Eid.one.matAB.Rdata")
load("fit.out.MSIRS.Eid.one.Rdata") #rerunning because it was fit to wrong datase
load("fit.Eid.mean.MSILI.Ifit.Rdata")

#load("fit.Eid.lci.MSILI.Ifit.Rdata")
#load("fit.out.MSIRS.Eid.lci.Rdata")
#load()

# 
 mod.dat <- rbind(refit.out.MSIRN.Eid.one.matAB[[2]],
                  fit.Eid.mean.MSILI.Ifit[[2]],
                  fit.out.MSIRS.Eid.one [[2]])
# 
# mod.lci <- rbind(refit.out.MSIRN.Eid.one.matAB[[2]],#need to change
#                  fit.Eid.lci.MSILI.Ifit[[2]],
#                  fit.out.MSIRS.Eid.lci[[2]]) 
# mod.dat=mod.lci
#prev.dat = dat.new


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
plot.comp.squiggle.multi <- function(mod.dat, prev.dat, type1, 
                                   cutoff1, species1, do.save, filename){
  #first, choose
  prev.dat1 <- subset(prev.dat, species==species1 & type==type1)
  prev.dat1  <- subset(prev.dat1, !is.na(age))
  #prev.dat1  <- subset(prev.dat1, fit==1)
  
  
  #get age.seroprev from data
  dat.seroprev = get.seroprev.dat(data=prev.dat1, vis_split=3, cutoff=cutoff1)
  dat.seroprev$lci = binom.exact(x=dat.seroprev$pos, n=dat.seroprev$count, conf.level = .95)$lower
  dat.seroprev$uci = binom.exact(x=dat.seroprev$pos, n=dat.seroprev$count, conf.level = .95)$upper
  
  mod.dat1 = subset(mod.dat, species==species1 & type==type1)# & cutoff==cutoff1)
  mod.dat1$model[mod.dat1$model=="MSIRN-matAB"] = "MSIRN"
  mod.dat1$model[mod.dat1$model=="MSILI"] = "MSIRNI"
  mod.dat1$model[mod.dat1$model=="MSIRNR-matAB"] = "MSIRNR"
  mod.dat1 = subset(mod.dat1, model!="MSIRN-matSus")
  mod.dat1 = subset(mod.dat1, model!="MSIRNR-matSus")
  
  # AIC.dat1 = subset(AIC.dat, species==species1 & type==type1 & cutoff==cutoff1)
  # AIC.dat1$model[AIC.dat1$model=="MSIRN-matAB"] = "MSIRN"
  # AIC.dat1$model[AIC.dat1$model=="MSIRNR-matAB"] = "MSIRNR"
  # AIC.dat1 = subset(AIC.dat1, model!="MSIRN-matSus")
  # AIC.dat1 = subset(AIC.dat1, model!="MSIRNR-matSus")
  # #add an "age" value to fit.dat for plotting
  # AIC.dat1$age <- NA
  # if(species1=="Eidolon dupreanum"){
  #   AIC.dat1$age<- 18.2
  # }else{
  #   AIC.dat1$age <- 13.5  
  # }
  # 
  # 
  # #and get deltaAIC relative to seroprev
  # AIC.df <- subset(AIC.dat1, cutoff==cutoff1)
  # min.AIC = min(AIC.df$AIC)
  # AIC.df$deltaAIC =AIC.df$AIC-min.AIC
  # 
  
  max.a = ceiling(max(dat.seroprev$age_year))
  # #if(cutoff1=="mean" | cutoff1=="lci"){
  # # ymax =.5
  # #}else if (cutoff1=="uci"){
  # # ymax=
  # #}
  # 
  # AIC.mult <- 1/(max(AIC.df$deltaAIC) + 1)
  # AIC.df$relAIC <- AIC.df$deltaAIC/3#AIC.mult
  # AIC.df <- subset(AIC.df,   model=="MSIRN" )
  # AIC.df$model <- factor(AIC.df$model, levels = c("MSIRN"))
  # AIC.df$outline <- NA
  # AIC.df$outline[AIC.df$type == "NIVG"] <- "NiV-G"
  # AIC.df$outline[AIC.df$type == "EBOV"] <- "EBOV-Gp"
  # AIC.df$outline <- factor(AIC.df$outline, levels=c("NiV-G", "EBOV-Gp"))
  # 
  # 
  shapz <- c("NiV-G" = 1, "EBOV-Gp" =2, "NIVG" = 16, "EBOV" = 17)
  
  
  #then, reorder the models
  #mod.dat1 <- subset(mod.dat1,  model=="MSIRN" )
  mod.dat1$model = factor(mod.dat1$model, levels = c("MSIRN", "MSIRS", "MSIRNI"))
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
  #label.dat = cbind.data.frame(c("MSIR",  "MSIRS", "MSIRN", "MSILI"), c("A.", "B.", "C.", "D."))
  #names(label.dat) = c("model", "label")
  #plot data
  if (species1=="Eidolon dupreanum"){
    
    p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1,  alpha=.9) + theme_bw() + facet_grid(species~model) +
      geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
      geom_linerange(data = dat.seroprev, aes(x=age_year, ymin=lci, ymax=uci), linetype=3, size=.1) +
      # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
      geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="purple", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
      geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "purple", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
      #geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =18, size =4, color="navy", show.legend = FALSE) +
      #geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 5, color="navy", size=5, show.legend = FALSE) +
      #geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=5) +
      theme(panel.grid = element_blank(),   strip.background = element_rect(fill="white"), strip.text.y = element_blank(), strip.text.x = element_text(size=16), 
            axis.title = element_text(size=16), axis.text = element_text(size=14), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy", size=16), axis.text.y.right = element_text(colour = "navy", size=12),
            legend.position = c(.93,.86), plot.margin = unit(c(.1,.1,.1,.1), "cm")) + 
      ylab("seroprevalence") + xlab("age(yrs)")   + coord_cartesian(xlim = c(0,15), ylim=c(0,1)) +
      scale_x_discrete(limits = c(0,5,10,15), labels = c("0","5", "10", "15")) 
    
    
    
  }else{
    mod.dat2 = subset(mod.dat2, age<=11)
    p2 <- ggplot() + geom_point(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), color ="black", pch=2,  alpha=.9, show.legend = FALSE) + theme_bw() + facet_grid(species~model) +
      geom_line(data=dat.seroprev, aes(x=age_year, y=prevalence, size=count), linetype=2, color ="black",  size=.2, show.legend = FALSE) +
      # geom_linerange(data=dat.seroprev, aes(x=age_year, ymin=prev_lci, ymax=prev_uci), linetype=3, size=.1, color ="black",  size=.2, show.legend = FALSE) +
      geom_line(data=mod.dat2, aes(x=age, y=mean_seroprev), color="lightgreen", size=1.5, show.legend = FALSE) +# scale_color_manual(values=color.comp) + scale_shape_manual(values=shapz) +
      geom_ribbon(data=mod.dat2, aes(x=age, ymin=lci_seroprev, ymax=uci_seroprev), fill = "lightgreen", alpha=.3, show.legend = FALSE) + #scale_fill_manual(values=color.comp) +
      #geom_point(data=AIC.df, aes(x=age, y=relAIC),  shape =5, size =4, color="navy", show.legend = FALSE) +
      #geom_point(data=AIC.df, aes(x=age, y=relAIC), shape = 18, color="navy", size=5, show.legend = FALSE) +
      geom_label(data=label.dat, aes(x=2, y=1, label=label), label.size = NA, size=4.5) +
      theme(panel.grid = element_blank(),   strip.background = element_blank(), strip.text.y = element_blank(), strip.text.x = element_text(size = 16), 
            axis.title = element_text(size=16), axis.text = element_text(size=12), legend.title = element_blank(), axis.title.y.right = element_text(colour = "navy", size=16), axis.text.y.right = element_text(colour = "navy", size=12),
            legend.text = element_blank(), plot.margin = unit(c(.1,.1,.1,.1), "cm")) + 
      ylab("seroprevalence") + xlab("age(yrs)")  + geom_vline(data=AIC.df, aes(xintercept=11), color="navy") + coord_cartesian(xlim = c(0,15), ylim=c(0,1)) + scale_y_continuous(sec.axis = sec_axis(trans= ~.*AIC.mult*1000,name = expression(Delta~AIC))) +
      scale_x_discrete(limits = c(0,5,10), labels = c("0", "5", "10")) + geom_segment(data=AIC.df, aes(x=0, xend=15, y=0, yend=0), linetype=2, color="navy")
    
  }
  
  
  
  #print(p2)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           plot=p2,
           units="mm",  
           width=60, 
           height=30, 
           scale=3, 
           dpi=300)
    
    
  }
  
  
}



plot.comp.squiggle.multi(mod.dat= mod.dat, 
                         prev.dat=dat.new, 
                         type1="NIVG",
                         cutoff1 = "mean",
                         species1="Eidolon dupreanum", 
                         do.save = T,
                         filename = "NIH-DP2-Fig.png")
