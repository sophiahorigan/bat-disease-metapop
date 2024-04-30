#now fit age-seroprevalence curves
#first, load data and plot
rm(list=ls())
library(epitools)
library(dplyr)
library(plyr)
library(cowplot)
library(MASS)
library(popbio)
library(ggplot2)

#first get the age-freq distribution, binned by year (0-1, 1-2, and on up to 15 for each species...)
#survivorship = proportion of population that survives to a given age (l_x)
#l_x = exp(-mu*x) 


#set this to your home directory
#homewd="/Users/carabrook/Developer/bat-PVA"
homewd= "/Users/theresalaverty/Documents/R/R_repositories/bat-PVA"
setwd(homewd)


#load from the catching data - you can clone this from github to your own computer
#madabatwd = "/Users/carabrook/Developer/Madagascar-Bat-Data"
madabatwd = "/Users/theresalaverty/Documents/R/R_repositories/Madagascar-Bat-Data"
dat.new <- read.csv(file = paste0(madabatwd, "/Catching-Data/catching_data.csv"), header = T, stringsAsFactors = F)
head(dat.new)

#this paper focused just on data prior to 2018, so for consistency, we will subset it here--
#but it could be updated!
dat.new$collection_date <- as.Date(dat.new$collection_date, format = "%m/%d/%y")
#dat.new = subset(dat.new, collection_date < "2018-01-01")
unique(dat.new$bat_age_class)

#link these with ages
age.dat <- read.csv(file = paste0(madabatwd, "/Age-Data/aged_bats.csv"), header = T, stringsAsFactors = F)
age.dat$collection_date
head(age.dat)
length(unique(age.dat$sampleid))

age.dat[duplicated(age.dat$sampleid),] 
#these are these duplicate id bats from Anecia and Kim's first sampling. 
#they are two collection dates, so they are independent samples, 
#despite having duplicate ids


#now, to give ages to the catching data, link by both sampleid and collection date
age.link <- dplyr::select(age.dat, sampleid, collection_date, age)

#and merge
dat.new <- merge(dat.new, age.link, by = c("sampleid", "collection_date"), all.x = T) 

head(dat.new)

#now get the pregnancy status
#get female Eid and Pteropus caught in Sept/Oct and tabulate pregnancy rates from these data

dat.preg = subset(dat.new, bat_species == "Eidolon dupreanum" | bat_species == "Pteropus rufus")

#select relevant columns
dat.preg = dplyr::select(dat.preg, sampleid, roost_site, collection_date, bat_species, bat_sex, bat_age_class, age)
dat.preg = subset(dat.preg, bat_age_class!="J" & !is.na(bat_age_class) & bat_sex == "female")
library(lubridate)
dat.preg$month = month(dat.preg$collection_date)


dat.preg = subset(dat.preg, age>1 | is.na(age)) #189 female adult bats with ages that could or could not be pregnant
#now pick out months when you would expect to observe pregnancy or lactation
dat.preg = subset(dat.preg,   month ==11 | month ==12 ) #month == 9 (changes very little when Sept included...goes down a little when December included--probably due to bat death)

#give 0: 1 for reproductive vs not
dat.preg$bat_age_class[dat.preg$bat_age_class =="L" | dat.preg$bat_age_class =="P"] = 1
dat.preg$bat_age_class[dat.preg$bat_age_class!=1] = 0
dat.preg$bat_age_class = as.numeric(dat.preg$bat_age_class)

#sep by species and get percents...
dat.preg.Eid = subset(dat.preg, bat_species=="Eidolon dupreanum")
dat.preg.Pter= subset(dat.preg, bat_species=="Pteropus rufus")


sum(dat.preg.Eid$bat_age_class)/length(dat.preg.Eid$bat_age_class) 
#these percents below are from the older data (prior to 2018) - could be redone, but around 50% total
#9-12: .4222; 10-12: .422; 9-11: .46875; 10-11: .46875
#percents from newer data (through 2020)
#9-12: 0.3529412; 10-12: 0.4761905; 9-11: 0.4594595; 10-11: 0.4871795

sum(dat.preg.Pter$bat_age_class)/length(dat.preg.Pter$bat_age_class) 
#these percents below are from the older data (prior to 2018) - could be redone, but around 95% total
#9-12: .746; 10-12: .7878; 9-11: .8571; 10-11:.9;
#percents from newer data (through 2020)
#9-12: 0.5833333; 10-12: 0.7727273; 9-11: 0.6896552; 10-11: 0.8518519

#there are female bats between 1-2 years for both species that have been identified as
#pregnant or lactating, so we assume these animals birth starting in the second year of life 
#(that 1 is probably a 2 actually)
#and remove the irrelevant columns
#we're just interested in demography here so we can really cull this
dat.tot <- dplyr::select(dat.new, sampleid, roost_site, collection_date, bat_species, bat_sex, bat_age_class, age)
head(dat.tot)

#make an age subset
dat.age <- dat.tot[!is.na(dat.tot$age),]
head(dat.age)
dat.age <- dat.age[!duplicated(dat.age),]

#and make sure it's just your two big species
dat.age <- subset(dat.age, bat_species=="Eidolon dupreanum" | bat_species =="Pteropus rufus")

#plot prevalence by season for the other
#check we have the right numbers:
nrow(subset(dat.age, bat_species=="Eidolon dupreanum"))#125
nrow(subset(dat.age, bat_species=="Pteropus rufus"))#145

#JAE reports 109 and 142 respectively -- a few were added here from the juvenile aging
#remember this could be re-run using the newer data

#bin by age.
age.bins = 0:15
dat.age$age_bin <- NA
for (i in 1:length(dat.age$age)){
  dat.age$age_bin[i] <- tail(age.bins[dat.age$age[i] > age.bins],n=1)
}

#now plot as a histogram.

#and the better, ggplot version with error bars...
plot.hist <- function(tot.dat, do.save, filename){
  
  #get age bin count for histogram
  hist.count <- ddply(tot.dat, .(bat_species, age_bin), summarize, count = length(sampleid))
  
  hist.count$hist_year <- hist.count$age_bin+.5
  
  hist.sum.Eid <- sum(hist.count$count[hist.count$bat_species=="Eidolon dupreanum"])
  hist.sum.Pter <- sum(hist.count$count[hist.count$bat_species=="Pteropus rufus"])
  
  hist.count$lci <- hist.count$uci <- NA
  hist.count$lci[hist.count$bat_species=="Eidolon dupreanum"] <- binom.exact(x=hist.count$count[hist.count$bat_species=="Eidolon dupreanum"], n=hist.sum.Eid)$lower*hist.sum.Eid
  hist.count$uci[hist.count$bat_species=="Eidolon dupreanum"] <- binom.exact(x=hist.count$count[hist.count$bat_species=="Eidolon dupreanum"], n=hist.sum.Eid)$upper*hist.sum.Eid
  
  hist.count$lci[hist.count$bat_species=="Pteropus rufus"] <- binom.exact(x=hist.count$count[hist.count$bat_species=="Pteropus rufus"], n=hist.sum.Pter)$lower*hist.sum.Pter
  hist.count$uci[hist.count$bat_species=="Pteropus rufus"] <- binom.exact(x=hist.count$count[hist.count$bat_species=="Pteropus rufus"], n=hist.sum.Pter)$upper*hist.sum.Pter
  
  
  
  #and plot by species
  p1 <- ggplot() + geom_bar(data=hist.count, aes(x=hist_year, y = count), stat="identity", fill="white", col="black") + 
    geom_linerange(data=hist.count, aes(x=hist_year, ymin=lci, ymax=uci), size=.3, linetype=3) + facet_grid(bat_species~.) +
    xlab("age(yrs)") + ylab("frequency") + theme_bw() + coord_cartesian(xlim=c(0, 15), ylim = c(0,160)) + 
    theme(axis.title = element_text(size=30, face = "plain"), panel.grid = element_blank(), 
          axis.text = element_text(size=25), title = element_text(size = 30, face="italic"),
          strip.text = element_text(face = "italic", size = 20), strip.background = element_blank())   
  print(p1)
  
  
  
  if(do.save==TRUE){

      ggsave(file = filename,
             units="mm",  
             width=60, 
             height=80, 
             scale=3, 
             dpi=200)
      
  }
  
}

plot.hist(tot.dat = dat.age, do.save = FALSE, filename = NA)

#also plot the inferred age distribution
plot.age.dist <- function(tot.dat, do.save, filename){
  
  #get age bin count for histogram
  hist.count <- ddply(tot.dat, .(species, age_bin), summarize, count = length(sampleid))
  
  
  hist.count$hist_year <- hist.count$age_bin+.5
  
  hist.sum.Eid <- sum(hist.count$count[hist.count$bat_species=="Eidolon dupreanum"])
  hist.sum.Pter <- sum(hist.count$count[hist.count$bat_species=="Pteropus rufus"])
  
  hist.count$prop <- NA
  hist.count$prop[hist.count$bat_species=="Eidolon dupreanum"] <- hist.count$count[hist.count$bat_species=="Eidolon dupreanum"]/hist.sum.Eid
  hist.count$prop[hist.count$bat_species=="Pteropus rufus"] <-hist.count$count[hist.count$bat_species=="Pteropus rufus"]/hist.sum.Pter
  
  hist.count$lci <- hist.count$uci <- NA
  hist.count$lci[hist.count$bat_species=="Eidolon dupreanum"] <- binom.exact(x=hist.count$count[hist.count$bat_species=="Eidolon dupreanum"], n=hist.sum.Eid)$lower
  hist.count$uci[hist.count$bat_species=="Eidolon dupreanum"] <- binom.exact(x=hist.count$count[hist.count$bat_species=="Eidolon dupreanum"], n=hist.sum.Eid)$upper
  
  hist.count$lci[hist.count$bat_species=="Pteropus rufus"] <- binom.exact(x=hist.count$count[hist.count$bat_species=="Pteropus rufus"], n=hist.sum.Pter)$lower
  hist.count$uci[hist.count$bat_species=="Pteropus rufus"] <- binom.exact(x=hist.count$count[hist.count$bat_species=="Pteropus rufus"], n=hist.sum.Pter)$upper
  
  
  
  #and plot by species
  p1 <- ggplot() + geom_point(data=hist.count, aes(x=hist_year, y = prop, size=count), pch=1, col="black") + 
    geom_linerange(data=hist.count, aes(x=hist_year, ymin=lci, ymax=uci), size=.3, linetype=3) + facet_grid(bat_species~.) +
    xlab("age(yrs)") + ylab("survivorship") + theme_bw() + coord_cartesian(xlim=c(0, 15), ylim = c(0,.5)) + 
    theme(axis.title = element_text(size=30, face = "plain"), panel.grid = element_blank(), 
          axis.text = element_text(size=25), title = element_text(size = 30, face="italic"),
          strip.text = element_text(face = "italic", size = 20), strip.background = element_blank())   
  print(p1)
  
  
  
  if(do.save==TRUE){
    
    ggsave(file = filename,
           units="mm",  
           width=60, 
           height=80, 
           scale=3, 
           dpi=200)
    
  }
  
}
loglik=function(par, data, method_mod){
  #print(par)
  tmp = (length(par[par<0]) + length(par[par>100]))
  if(tmp>0){
    return(1000000)
  }else{
    
    if (method_mod=="constant"){
      lx = list()
      for(a in 1:length(data$prop)){ #for-loops over every individual in the dataset
        lx[[a]] = exp(-1*par*data$age_bin[a])#probability of surviving to a given age bin with the fitted slope (survival rate)
      }
      lx = c(unlist(lx))
    }else if (method_mod=="maturation"){
      lx1 = list()
      lx2 = list()
      for(a in 1:length(data$prop)){ #for-loops over every individual in the dataset
        lx1[[a]] = exp((-1*par[2]/par[3])*(1-exp(-1*par[3]*data$age_bin[a])))
        lx2[[a]] = exp(-1*par[1]*data$age_bin[a])#probability of surviving to a given age bin with the fitted slope (survival rate)
      }
      lx1=c(unlist(lx1))
      #then, we cull it at only the first age class (because maturation)
      lx1[2:length(lx1)] = 1
      lx2=c(unlist(lx2))
      lx=lx1*lx2
    }else if (method_mod=="senescence"){
      lx2 = list()
      lx3 = list()
      for(a in 1:length(data$prop)){ #for-loops over every individual in the dataset
        lx3[[a]] = exp((-1*par[2]/par[3])*(1-exp(-1*par[3]*data$age_bin[a])))
        lx2[[a]] = exp(-1*par[1]*data$age_bin[a])#probability of surviving to a given age bin with the fitted slope (survival rate)
      }
      lx2=c(unlist(lx2))
      lx3=c(unlist(lx3))
      #then, we cull it at only individuals over 5 (because senescence)
      lx3[1:5] = 1
      lx=lx2*lx3
    }else if (method_mod =="both"){
      lx1 = list()
      lx2 = list()
      lx3 = list()
      for(a in 1:length(data$prop)){ #for-loops over every individual in the dataset
        lx1[[a]] = exp((-1*par[2]/par[3])*(1-exp(-1*par[3]*data$age_bin[a])))
        lx2[[a]] = exp(-1*par[1]*data$age_bin[a])#probability of surviving to a given age bin with the fitted slope (survival rate)
        lx3[[a]] = exp((-1*par[4]/par[5])*(1-exp(-1*par[5]*data$age_bin[a])))
      }
      lx1=c(unlist(lx1))
      #then, we cull it at only the first age class (because maturation)
      lx1[2:length(lx1)] = 1
      lx2=c(unlist(lx2))
      lx3=c(unlist(lx3))
      #then, we cull it at only individuals over 5 (because senescence)
      lx3[1:5] = 1
      lx=lx1*lx2*lx3
    }
    
    Px = lx/ sum(lx)
    #print(Px)
    ll = dmultinom(x=data$count, prob = Px, log=T)
    return(-ll)
  }
}
compare.models = function(data){
  out = optim(par =.1, fn=loglik, data=data, method_mod = "constant", method="BFGS", hessian=T)
  out_sens = optim(par = c(.1, .1,.1), fn=loglik, data=data, method_mod = "senescence", method="BFGS", hessian=T)
  out_mat = optim(par = c(.1, .1,.1), fn=loglik, data=data, method_mod = "maturation", method="BFGS", hessian=T)
  out_both = optim(par = c(.1, .1,.1, .1, .1), fn=loglik, data=data, method_mod = "both", method="BFGS", hessian=T)
  
  mod.type = c("constant", "senescence", "maturation", "both")  
  loglik.dat = c(out$value, out_sens$value, out_mat$value, out_both$value)  
  k = c(1,3,3,5)
  cmp.df = cbind.data.frame( mod.type , -1*loglik.dat, k)
  names(cmp.df) = c("mod_type", "llik", "k")
  #and get AIC
  cmp.df$AIC =  2*cmp.df$k-2*cmp.df$llik
  return(cmp.df)
}
surv.siler.comp <- function(species1, tot.dat){
  
  #first, get data into proportion at each age bin to compare with model
  sub.dat <- subset(tot.dat, bat_species==species1)
  hist.count2 <- ddply(sub.dat, .(age_bin), summarize, count = length(sampleid))
  age.bins <- 0:14
  
  hist.count <- cbind.data.frame(age.bins, rep(0, length(age.bins)))#get the whole span of ages
  names(hist.count) <- c("age_bin", "count")
  
  #fill in data
  for (i in 1:length(hist.count$count)){
    hist.count$count[hist.count$age_bin==hist.count2$age_bin[i]] <- hist.count2$count[i]  
  }
  
  hist.sum <- sum(hist.count$count)
  hist.count$prop <- NA
  hist.count$prop <- hist.count$count/hist.sum
  hist.count$lci <- binom.exact(x=hist.count$count, n=hist.sum)$lower
  hist.count$uci <- binom.exact(x=hist.count$count, n=hist.sum)$upper
  hist.count$hist_year <- hist.count$age_bin +.5
  
  #then, look at just bats aged over 1 year
  adult.hist <- subset(hist.count, age_bin>0)
  
  #then, fit a simple exponenttial to infer survivorship. must remove 0s
  adult.hist.mod <- adult.hist[adult.hist$count>0,]
  
  #but have to remake proportion from 100
  adult.hist.mod$N = sum(adult.hist.mod$count)
  adult.hist.mod$prop = adult.hist.mod$count/adult.hist.mod$N
  adult.hist.mod$lci <- binom.exact(x=adult.hist.mod$count, n=adult.hist.mod$N)$lower
  adult.hist.mod$uci <- binom.exact(x=adult.hist.mod$count, n=adult.hist.mod$N)$upper
  
  
  #then, run and  compare
  
  
  compare.df = compare.models(data=adult.hist.mod)
  return(compare.df)
}
getJuvSurvAD <- function(adult_fecundity, adult.surv){
  i = (1-adult.surv)/(adult_fecundity*(adult.surv^2))
  juv.surv <- i*adult.surv
  return(juv.surv)
}
get.stab.age <- function(surv.par, adult_fecundity){
  
  #then, build your matrix
  bat.mat <- matrix(0,15,15) #we make a matrix of 15 age classes - assume individuals much older than that are VERY unlikely
  bat.mat[1,2:15] <- adult_fecundity #hardwire fertility - female pups born annually to each bat, starting at age 3...
  bat.mat[2,1] <- surv.par[1] #infant survivorship (first year of life 0 to 1)
  #bat.mat[cbind(3,2)] <-(surv.par[2]) #juvenile survivorship (year 1 to 2)
  bat.mat[cbind(3:15,2:14)] <- (surv.par[length(surv.par)]) #adult survivorship = everything else
  bat.mat[15,15] <- (surv.par[length(surv.par)]) #adult survivorship = everything else
  
  stab.age <- Re(eigen(bat.mat)$vector[,1])
  stab.age <- stab.age/sum(stab.age)
  return(stab.age )
}
R0Calc <- function (Pmatrix, Fmatrix) {
  ns <- nrow(Pmatrix)
  Imatrix <- matrix(0, ns,ns)
  diag(Imatrix) <- 1
  Nmatrix <- ginv(Imatrix - Pmatrix)
  Rmatrix <- Fmatrix %*% Nmatrix
  ave.R0 <- Re(eigen(Rmatrix)$values[1])
  return(ave.R0)
}
getR0 <- function(juv.surv, adult.surv, adult_fecundity){
  
  ## Pmatrix is transitions due to survival and growth, Fmatrix transitions due to reproduction
  ## so the latter might just have values in the top right hand corner. 
  Pmat <- matrix(0,15,15)
  Pmat[2,1] <- juv.surv #your guess which gets optimized
  Pmat[cbind(3:15,2:14)] <- adult.surv
  
  
  Fmat <- matrix(0,15,15)
  Fmat[1,2:15]<- adult_fecundity
  
  print(paste0("R0=",R0Calc(Pmat,Fmat)))
  
  R0.out <- R0Calc(Pmat,Fmat)
  
  return(R0.out)
}
constant.exp.fit <- function(species1, tot.dat, adult_fecundity, do.plot){
  
  #first, get data into proportion at each age bin to compare with model
  sub.dat <- subset(tot.dat, bat_species==species1)
  hist.count2 <- ddply(sub.dat, .(age_bin), summarize, count = length(sampleid))
  age.bins <- 0:14
  
  hist.count <- cbind.data.frame(age.bins, rep(0, length(age.bins)))#get the whole span of ages
  names(hist.count) <- c("age_bin", "count")
  
  #fill in data
  for (i in 1:length(hist.count$count)){
    hist.count$count[hist.count$age_bin==hist.count2$age_bin[i]] <- hist.count2$count[i]  
  }
  
  hist.sum <- sum(hist.count$count)
  hist.count$prop <- NA
  hist.count$prop <- hist.count$count/hist.sum
  hist.count$lci <- binom.exact(x=hist.count$count, n=hist.sum)$lower
  hist.count$uci <- binom.exact(x=hist.count$count, n=hist.sum)$upper
  hist.count$hist_year <- hist.count$age_bin +.5
  
  #then, look at just bats aged over 1 year
  adult.hist <- subset(hist.count, age_bin>0)
  
  #then, fit a simple exponenttial to infer survivorship. must remove 0s
  adult.hist.mod <- adult.hist[adult.hist$count>0,]
  
  #but have to remake proportion from 100
  adult.hist.mod$N = sum(adult.hist.mod$count)
  adult.hist.mod$prop = adult.hist.mod$count/adult.hist.mod$N
  adult.hist.mod$lci <- binom.exact(x=adult.hist.mod$count, n=adult.hist.mod$N)$lower
  adult.hist.mod$uci <- binom.exact(x=adult.hist.mod$count, n=adult.hist.mod$N)$upper
  
  
  #then, we use the constant model
  out = optim(par =.1, fn=loglik, data=  adult.hist.mod, method_mod = "constant", method="BFGS", hessian=T)
  
  m=out$par  #then, run each individual as a count
  llik= -1*out$value
  AIC = 2*1-2*llik
  adult.surv = 1-m
  
  
  #then, we fix fertility and optimize juvenile survivorship to make pop growth=1.
  inf.surv <- getJuvSurvAD(adult_fecundity=adult_fecundity, adult.surv=adult.surv)
  
  if(inf.surv>1){
    warning(paste0("Infant growth rate cannot be high enough to protect this species. Infant survival estimated at = "),round(inf.surv,3), ". Rate instead fixed at E. duprenanum assumption.")
   inf.surv = .68 
  }
  
  #can you ask for stable age fit?
  stab.struc.mod <- get.stab.age(surv.par=c(inf.surv, adult.surv), adult_fecundity=adult_fecundity)
  
  
  #sm_sq <- sum((hist.count$prop-stab.struc.mod)^2) 
  convergence <- out$convergence
  #then, plot stable age distribution with this, with data
  
  #then, print R0 (for a different matrix structure)
  R0 <- round(getR0(juv.surv = inf.surv, adult.surv = adult.surv, adult_fecundity = adult_fecundity), 2)
  
  est.surv <- c(inf.surv, adult.surv)
  names(est.surv) <- c()
  
  if(do.plot==TRUE){
    mod.stab.age <- get.stab.age(est.surv, adult_fecundity)
    mod.stab.age <- cbind.data.frame(hist.count$hist_year, mod.stab.age)
    names(mod.stab.age) <- c("hist_year", "prop")
    
    p1 <- ggplot() +  geom_point(data=hist.count, aes(x=hist_year, y = prop, size=count), pch=1, col="black") + 
      geom_linerange(data=hist.count, aes(x=hist_year, ymin=lci, ymax=uci), size=.3, linetype=3) + 
      geom_line(data = mod.stab.age, aes(x=hist_year, y = prop), col="red") +
      xlab("age(yrs)") + ylab("proportion of population") + theme_bw() + coord_cartesian(xlim=c(0, 15), ylim = c(0,.5)) + 
      theme(axis.title = element_text(size=30, face = "plain"), panel.grid = element_blank(), 
            axis.text = element_text(size=25), title = element_text(size = 30, face="italic"),
            strip.text = element_text(face = "italic", size = 20), strip.background = element_blank())   
    print(p1)
  }
  
  #and return the population growth rate
  
  new.dat <- cbind.data.frame(species1, adult_fecundity, est.surv[1], est.surv[2], R0, llik, AIC, convergence)
  names(new.dat) <- c("species", "adult_fecundity", "juv_survival", "adult_survival", "pop_growth_rate", "llik", "AIC", "convergence")
  
  return(new.dat)
}


out.comp.Eid  = surv.siler.comp(species1 = "Eidolon dupreanum", tot.dat = dat.age)
out.comp.Pter  = surv.siler.comp(species1 = "Pteropus rufus", tot.dat = dat.age)
#constant fit best for all


constant.exp.fit(species1 = "Eidolon dupreanum", tot.dat = dat.age, adult_fecundity = (  0.2205), do.plot = T) #if adult fecundity really is half of P. rufus, these guys could be in serious trouble too...but this is likely not the case
constant.exp.fit(species1 = "Eidolon dupreanum", tot.dat = dat.age, adult_fecundity = (.48), do.plot = T) 
constant.exp.fit(species1 = "Pteropus rufus", tot.dat = dat.age, adult_fecundity = .48, do.plot = T)
constant.exp.fit(species1 = "Pteropus rufus", tot.dat = dat.age, adult_fecundity = 0.424, do.plot = T)
 
#and double check these rates against the simple exponential model
#these rates are save in bat demo. now load and plot with data

load(paste0(homewd,"/prior-scripts/BiolConsBat/data/bat.demo.Oct.25.2018.Rdata"))
bat.demo
names(bat.demo)[names(bat.demo)== "infant_IBI_survival"] <-  "infant_annual_survival" 
names(bat.demo)[names(bat.demo)== "adult_IBI_survival"] <-  "adult_annual_survival" 
bat.demo$IBI_fecundity <- c(.48, .2205, .48, .424)
bat.demo$infant_annual_survival <- c(.5866705, 1.277,4.519, 5.115) #this is added in to recover a stable age population
bat.demo$adult_annual_survival <- c(.7802735,.7802735,.3155653,.3155653)


#and the better, ggplot version with error bars...
plot.hist.model <- function(tot.dat, mod.dat, do.save, filename){
  
  #just take "best"
  mod.dat = subset(mod.dat, scenario == "best")
  #write over juv estimates for Pteropus (over 1)
  mod.dat$infant_annual_survival[mod.dat$infant_annual_survival>1] <- mod.dat$infant_annual_survival[mod.dat$infant_annual_survival<1] 
  #get age bin count for histogram
  
  hist.count2 <- ddply(tot.dat, .(bat_species, age_bin), summarize, count = length(sampleid))
  age.bins <- 0:14
  
  hist.count <- cbind.data.frame(age.bins, rep(0, length(age.bins)))#get the whole span of ages
  names(hist.count) <- c("age_bin", "count")
  hist.count <- rbind(  hist.count,   hist.count)
  hist.count$species <- rep(c("Eidolon dupreanum", "Pteropus rufus"), each = length(age.bins))
  
  #fill in data
  for (i in 1:length(hist.count$count)){
    hist.count$count[hist.count$age_bin==hist.count2$age_bin[i] & hist.count$species==hist.count2$bat_species[i]] <- hist.count2$count[i]  
  }
  
  
  hist.count$hist_year <- hist.count$age_bin+.5
  
  hist.sum.Eid <- sum(hist.count$count[hist.count$species=="Eidolon dupreanum"])
  hist.sum.Pter <- sum(hist.count$count[hist.count$species=="Pteropus rufus"])
  
  hist.count$lci <- hist.count$uci <- NA
  hist.count$lci[hist.count$species=="Eidolon dupreanum"] <- binom.exact(x=hist.count$count[hist.count$species=="Eidolon dupreanum"], n=hist.sum.Eid)$lower*hist.sum.Eid
  hist.count$uci[hist.count$species=="Eidolon dupreanum"] <- binom.exact(x=hist.count$count[hist.count$species=="Eidolon dupreanum"], n=hist.sum.Eid)$upper*hist.sum.Eid
  
  hist.count$lci[hist.count$species=="Pteropus rufus"] <- binom.exact(x=hist.count$count[hist.count$species=="Pteropus rufus"], n=hist.sum.Pter)$lower*hist.sum.Pter
  hist.count$uci[hist.count$species=="Pteropus rufus"] <- binom.exact(x=hist.count$count[hist.count$species=="Pteropus rufus"], n=hist.sum.Pter)$upper*hist.sum.Pter
  
  #and build model dataset
  
  mod.stab.age.Eid <- get.stab.age(unique(c(mod.dat$infant_annual_survival[mod.dat$species=="Eidolon dupreanum"], mod.dat$adult_annual_survival[mod.dat$species=="Eidolon dupreanum"])), adult_fecundity= unique(mod.dat$IBI_fecundity[mod.dat$species=="Eidolon dupreanum"]))
  mod.stab.age.Pter <- get.stab.age(unique(c(mod.dat$infant_annual_survival[mod.dat$species=="Pteropus rufus"], mod.dat$adult_annual_survival[mod.dat$species=="Pteropus rufus"])), adult_fecundity= unique(mod.dat$IBI_fecundity[mod.dat$species=="Pteropus rufus"]))
  mod.stab.age <- cbind.data.frame(hist.count$hist_year, c(mod.stab.age.Eid, mod.stab.age.Pter))
  names(mod.stab.age) <- c("hist_year", "prop")
  mod.stab.age$species <- rep(c("Eidolon dupreanum", "Pteropus rufus"), each = length(age.bins))
  mod.stab.age$count <- NA
  mod.stab.age$count[mod.stab.age$species=="Eidolon dupreanum"] <- round(mod.stab.age$prop[mod.stab.age$species=="Eidolon dupreanum"]*hist.sum.Eid ,0)
  mod.stab.age$count[mod.stab.age$species=="Pteropus rufus"] <- round(mod.stab.age$prop[mod.stab.age$species=="Pteropus rufus"]*hist.sum.Pter,0)
   
  #and lci uci for model
  mod.stab.age$prop_lci<-  mod.stab.age$prop_uci<- mod.stab.age$count_lci <- mod.stab.age$count_uci <- NA
  
  mod.stab.age$prop_lci[mod.stab.age$species=="Eidolon dupreanum"] <- binom.exact(mod.stab.age$count[mod.stab.age$species=="Eidolon dupreanum"], hist.sum.Eid)$lower
  mod.stab.age$prop_uci[mod.stab.age$species=="Eidolon dupreanum"] <- binom.exact(mod.stab.age$count[mod.stab.age$species=="Eidolon dupreanum"], hist.sum.Eid)$upper
  
  mod.stab.age$prop_lci[mod.stab.age$species=="Pteropus rufus"] <- binom.exact(mod.stab.age$count[mod.stab.age$species=="Pteropus rufus"], hist.sum.Pter)$lower
  mod.stab.age$prop_uci[mod.stab.age$species=="Pteropus rufus"] <- binom.exact(mod.stab.age$count[mod.stab.age$species=="Pteropus rufus"], hist.sum.Pter)$upper
  
  
  mod.stab.age$count_lci[mod.stab.age$species=="Eidolon dupreanum"] <- mod.stab.age$prop_lci[mod.stab.age$species=="Eidolon dupreanum"]*hist.sum.Eid 
  mod.stab.age$count_uci[mod.stab.age$species=="Eidolon dupreanum"] <- mod.stab.age$prop_uci[mod.stab.age$species=="Eidolon dupreanum"]*hist.sum.Eid 
  
  mod.stab.age$count_lci[mod.stab.age$species=="Pteropus rufus"] <- mod.stab.age$prop_lci[mod.stab.age$species=="Pteropus rufus"]*hist.sum.Pter
  mod.stab.age$count_uci[mod.stab.age$species=="Pteropus rufus"] <- mod.stab.age$prop_uci[mod.stab.age$species=="Pteropus rufus"]*hist.sum.Pter
  
  #and include the relevant parametric data as text
  mod.dat$Fa = mod.dat$IBI_fecundity
  mod.dat$Fa <- as.character(as.expression(paste0("F[A]==", mod.dat$Fa)))
  mod.dat$sa = mod.dat$adult_annual_survival
  mod.dat$sa <- as.character(as.expression(paste0("s[A]==", round(mod.dat$sa,3))))
  mod.dat$sj = mod.dat$infant_annual_survival
  mod.dat$sj <- as.character(as.expression(paste0("s[J]==", round(mod.dat$sj,3))))
  #and plot by species
  
  p1 <- ggplot() + geom_bar(data=hist.count, aes(x=hist_year, y = count), stat="identity", fill="white", col="black") + 
    geom_linerange(data=hist.count, aes(x=hist_year, ymin=lci, ymax=uci), size=.3, linetype=3) + facet_grid(species~.) +
    geom_line(data=mod.stab.age, aes(x=hist_year, y=count), col="red", size=1)+
    geom_label(data=mod.dat, aes(x=14, y=120, label =sj), parse=T, label.size = NA, col="red")+
    geom_label(data=mod.dat, aes(x=14, y=135, label =sa), parse=T, label.size = NA, col="red")+
    geom_label(data=mod.dat, aes(x=14, y=150, label =Fa), parse=T, label.size = NA, col="red")+
    geom_ribbon(data=mod.stab.age, aes(x=hist_year, ymax=count_uci, ymin=count_lci), fill="red", alpha=.3)+
    xlab("age(yrs)") + ylab("frequency") + theme_bw() + coord_cartesian(xlim=c(0, 15), ylim = c(0,160)) + 
    theme(axis.title = element_text(size=30, face = "plain"), panel.grid = element_blank(), 
          axis.text = element_text(size=25), title = element_text(size = 30, face="italic"),
          strip.text = element_text(face = "italic", size = 20), strip.background = element_blank())   
  print(p1)
  
  
  
  if(do.save==TRUE){
    
    ggsave(file = filename,
           units="mm",  
           width=60, 
           height=80, 
           scale=3, 
           dpi=200)
    
  }
  
}
plot.hist.model(tot.dat = dat.age, mod.dat=bat.demo, do.save = FALSE, filename = NA)


plot.hist.model(tot.dat = dat.age, mod.dat=bat.demo, do.save = TRUE, filename = paste0(homewd, "/prior-scripts/BiolConsBat/figures/Fig2_age_freq_mod_dat.pdf"))

##########################################################################################
##########################################################################################

#now compare against subfossil data
#now get age at death from samples ages and compare:

# Simulate age at death from our Leslie survival parameters (above)
#if we use the high birth rate for Eid:

get.sample = function(juv_surv, adult_surv){

  surv.par <- c(juv_surv, adult_surv) #juvenile, adult
  nlive <- nstart <- 1000
  ndead <- rep(NA,1000)
  
  for (j in 1:1000) { 
    nlive <- rbinom(1,nstart,surv.par[min(j,length(surv.par))])
    ndead[j] <- nstart-nlive
    nstart <- nlive
  }
  
  ndead = ndead[ndead>0]
  data.f = cbind.data.frame(1:length(ndead), ndead)
  names(data.f) = c("age", "num_dead")
  
  age.sample = list()
  for (i in 1:length(data.f$age)){
    #print(i)
    tmp =  rep(data.f$age[i], data.f$num_dead[i])  
    age.sample[[i]] = tmp
  }
  age.sample = c(unlist(age.sample))
  
  return(age.sample)
}
#age.at.death.mod.Eid = get.sample(juv_surv = .68, adult_surv =  0.753853)
age.at.death.mod.Eid = get.sample(juv_surv = .68, adult_surv =  0.753853)


#then, pull data from the subfossil
#madabatwd = "/Users/carabrook/Developer/Madagascar-Bat-Data"
madabatwd = "/Users/theresalaverty/Documents/R/R_repositories/Madagascar-Bat-Data"

df <- read.csv( paste0(madabatwd, "/Age-Data/Raw-Matsons-Data/subfossil_age_data.csv"), header = T)
head(df)
unique(df$species)
dat = subset(df, sample_type == "subfossil" & species=="Eidolon dupreanum") #34
dat = dplyr::select(dat, species, tooth_type, sampleid, age, cc, age_range)
length(dat$sampleid) #34
length(unique(dat$sampleid)) #34 - there are not multiple teeth from the same bat- but we can't be sure
dat$age = as.numeric(dat$age)
dat = subset(dat, !is.na(age)) #19
#age should be whatever is specified - look at age-freq distribution, by single year

pdf(paste0(homewd, "/prior-scripts/BiolConsBat/figures/age_at_death.pdf"))
par(cex = (1.5))
plot(density(dat$age), ylab="Density", xlab=c("Age at death"), xlim=c(0,15), ylim=c(0,.5), col="blue", main = NA,lwd=1.5)
lines(density(age.at.death.mod.Eid), col="red", lwd=1.5)
legend(x=8.5,y=.5, legend = c("modern (modeled)", "subfossil"), col=c("red", "blue"), lwd=1.5, cex = .7)
dev.off()

#and compare statistically
ks.test(x=age.at.death.mod.Eid, y= dat$age)
#D = 0.165, p-value = 0.6904
#alternative hypothesis: two-sided
#not significant difference

#and add in the Pteropus distribution:
#age.at.death.mod.Pter = get.sample(juv_surv = .68, adult_surv =  0.3849284)
age.at.death.mod.Pter = get.sample(juv_surv = .68, adult_surv =  0.3849284)

par(cex = (1.5))
plot(density(dat$age), ylab="Density", xlab=c("Age at death"), xlim=c(0,15), ylim=c(0,.5), col="blue", main = NA,lwd=1.5)
lines(density(age.at.death.mod.Eid), col="red", lwd=1.5)
lines(density(age.at.death.mod.Pter), col="purple", lwd=1.5)
legend(x=8.5,y=.5, legend = c("modern (modeled)", "subfossil"), col=c("red", "blue"), lwd=1.5, cex = .7)
