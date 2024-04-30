#first, make some phase plane plots for your bats, based on literature-derived demographic parameters. 
#then, plot mortalities on top of them

rm(list=ls())

#first, load bat ages

library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
library(stringr)

#set this to your home directory
homewd="/Users/carabrook/Developer/bat-PVA"
setwd(homewd)

#load demographic info
load(paste0(homewd, "/prior-scripts/BiolConsBat/data/bat.demo.Oct.25.2018.Rdata"))

#goal is to print all combinations for each species and lambda with a range
#based on each par
calc.IBI.fec.iso <- function(df){
  
  #IBI fecundity is pretty much always going to be 1 (for species that twin) or .5 (for species that don't twin)
  new.rows <- data.frame(matrix(NA, nrow=2, ncol=ncol(df)))
  names(new.rows) <- names(df)
  new.rows$parameter_type <- "IBI_fecundity"
  new.rows$stage <- "adult"
  new.rows$species <- as.character(unique(df$species))
  new.rows$scenario <- c("best", "worst")
  new.rows$value[new.rows$scenario=="best"] <- as.numeric(df$value[df$parameter_type=="litter_size" & df$scenario=="best"])
  new.rows$value[new.rows$scenario=="worst"] <- as.numeric(df$value[df$parameter_type=="litter_size" & df$scenario=="worst"])
  
  new.rows$iso_sim <- NA
  #now add the letters
  if(length(unique(new.rows$value)) > 1){
    new.rows$iso_sim[new.rows$scenario=="best"] <- "A"
    new.rows$iso_sim[new.rows$scenario=="worst"] <- "B"
  }
  if(length(unique(new.rows$value)) == 1){
    new.rows$iso_sim <- "C"  
  }
  
  #let's make it female fecundity and divide by 2!!!
  new.rows$value <- new.rows$value/2
  
  #and bind to old
  df <- rbind(df, new.rows)
  return(df)  
}
choose.par <- function(df, type){
  #split df by parameter type
  if(type=="isocline"){
    #we don't need mortality data
    df <- df[df$parameter_type!="annual_mortality",] 
    #then, we jut split by type
    par.list <- dlply(df, .(parameter_type))  
  }
  if(type!="isocline"){#otherwise, we need to keep the stage info too
    par.list <- dlply(df, .(parameter_type, stage))  
  }
  par.list <- lapply(par.list, cull.par, type=type)
  #now par.list as dataframe
  df2 <- do.call("rbind", par.list)
  return(df2)
}
cull.par <- function(df, type){
  if(type=="worst"){
    df1 <- df[1,]
    if (unique(df$parameter_type)=="IBI"){
      df1$value <- max(df$value)  
    }
    if (unique(df$parameter_type)=="age_1st_rep"){
      df1$value <- max(df$value)  
    }
    else{
      df1$value <- min(df$value) 
    }
    
    df <- df1
  }
  
  if(type=="mean"){
    df1 <- df[1,]
    df1$value <- mean(df$value)
    df <- df1
  }
  if(type=="best"){
    df1 <- df[1,]
    if (unique(df$parameter_type)=="IBI"){
      df1$value <- min(df$value)  
    }
    if (unique(df$parameter_type)=="age_1st_rep"){
      df1$value <- min(df$value)  
    }
    else{
      df1$value <- max(df$value) 
    }
    df <- df1
  }
  if(type=="isocline"){
    #make two places this time
    #and take either A and B or C
    df1 <- rbind(df[1,], df[1,])
    df1$scenario <- c("best", "worst")
    if (unique(df$parameter_type)=="IBI"){
      #write over iso_sim and value categories
      df1$value[1] <- df$value[df$iso_sim=="A" | df$iso_sim=="C" ]  #best
      df1$iso_sim[1] <- df$iso_sim[df$iso_sim=="A" | df$iso_sim=="C" ] 
      df1$value[2] <- df$value[df$iso_sim=="B" | df$iso_sim=="C"]  #worst
      df1$iso_sim[2] <- df$iso_sim[df$iso_sim=="B" | df$iso_sim=="C"]  #worst
    }
    if (unique(df$parameter_type)=="age_1st_rep"){
      df1$value[1] <- min(df$value)  #best
      df1$value[2] <- max(df$value)  #worst
    }
    if (unique(df$parameter_type)=="litter_size"){
      df1$value[1] <- max(df$value)  #best
      df1$value[2] <- min(df$value)  #worst
    }
    if (unique(df$parameter_type)=="avg_lifespan"){
      df1$value[1] <- max(df$value)  #best
      df1$value[2] <- min(df$value)  #worst
    }
    df <- df1
  }
  return(df)
}
wrap.choose.par <- function(df, type){
  df.list <- dlply(df, .(species))
  df.list <- lapply(df.list, choose.par, type=type)
  df <- do.call("rbind", df.list)
  rownames(df) <- c()
  return(df)
}
reshape.data <- function(df, type){
  #now select the columns you need for leslie matrix and rename to be consistent with pop_model
  if(type=="isocline"){
    df2 <- dcast(df, scenario + species ~ parameter_type, value.var="value")
    df3 <- select(df2, c(species, scenario, IBI_fecundity, age_1st_rep, avg_lifespan, IBI, litter_size))
  }
  if(type!="isocline"){
    df2 <- dcast(df, stage + species ~ parameter_type, value.var="value")
    df3 <- select(df2, c(species, stage, annual_mortality, annual_fecundity, age_1st_rep, avg_lifespan, dispersal_prob, IBI, litter_size))
  }
  #and return to send to our model
  return(df3)
}
wrap.pop.data <- function(df, type){
  #select only the columns we need for the mat
  df <- select(df, c(parameter_type, stage, species, iso_sim, value))
  head(df)
  
  #then, if type=isocline, remove the NAs
  if(type=="isocline"){
    df <- df[complete.cases(df),]  
  }
  
  #first need to pick just one (worst, best, or mean) of every parameter
  #split by species
  df2 <- wrap.choose.par(df, type=type)
  
  #then, calculate annual fecundity (for those species that lack it)
  #from IBI and litter size. Then, reshape for population model.
  df.list <- dlply(df2, .(species))
  
  if (type=="isocline"){
    df.list <- lapply(df.list, calc.IBI.fec.iso)  
  }
  if (type!="isocline"){
    df.list <- lapply(df.list, calc.annual.fec)  
  }
  
  #bind back
  df3 <- do.call("rbind", df.list)
  rownames(df3) <- c()
  
  
  #now reshape table into the one that you need
  #need the following columns: species, stage, stage_survivorship, annual_fecundity, age_1st_rep, 
  #avg_lifespan, dispersal_prob (by stage)
  #also, just to have it for the lefkovitch matrices, calculate the fecundity
  #by the interbirth interval - should just be the littersize (divided by 2 for females)
  #also, attach the IBI to each species because we may need it later for calculating 
  #Lefkovitch survivorship
  
  df4 <- reshape.data(df3, type=type)
  
  #when not making isoclines, itgenerates NAs because it makes a column
  #for every parameter for all stages. We can fill these in by repeating the adult values
  if (type!="isocline"){
    df4.list <-  dlply(df4, .(species))
    df4.list <- lapply(df4.list, rep.adult.par)
    
    #bind and return
    pop.data <- do.call("rbind", df4.list)
    
    
    #now overwrite names for litter size and divide by 2 because we are only modeling females
    names(pop.data)[length(names(pop.data))] <- "IBI_fecundity"
    pop.data$IBI_fecundity <- as.numeric(pop.data$IBI_fecundity)/2
    
    #then make sure that everything that all the parameters are numeric
    pop.data$annual_mortality <- as.numeric(pop.data$annual_mortality)
    pop.data$annual_fecundity <- as.numeric(pop.data$annual_fecundity)
    # Note: We here list fecundity for infants of both sexes produced.
    # Must halve fecundity if we are just modeling females.
    # Do this right away.
    pop.data$annual_fecundity <- pop.data$annual_fecundity/2
    pop.data$age_1st_rep <- as.numeric(pop.data$age_1st_rep)
    pop.data$avg_lifespan <- as.numeric(pop.data$avg_lifespan)
    pop.data$dispersal_prob <- as.numeric(pop.data$dispersal_prob)
    rownames(pop.data) <- c()
  }
  
  #isocline data is already ready to go 
  if (type=="isocline"){
    pop.data <- df4
  }
  #and return
  return(pop.data)
  
}
isocline.wrap <- function(df, p_sequence){
  #new version essentially plots isocline wrap based on the age at first rep for each species
  
  zero_growth_line.list <- list()
  for(i in 1:2){
    zero_growth_line.list[[i]] <- lapply(p_sequence, get_isoclines, IBI_fecundity=as.numeric(as.character(df$IBI_fecundity[i])), IBI=as.numeric(as.character(df$IBI[i])), age_1st_rep=as.numeric(as.character(df$age_1st_rep[i])))           
  }
  zero_growth_line <- matrix(unlist(zero_growth_line.list[[1]]), ncol=2, byrow=T)
  zero_growth_line <- data.frame(zero_growth_line)
  names(zero_growth_line) <- c("IBI_infant_mortality_rate", "IBI_adult_mortality_rate")
  zero_growth_line$scenario <- "best"
  zero_growth_line_worst <- matrix(unlist(zero_growth_line.list[[2]]), ncol=2, byrow=T)
  zero_growth_line_worst <- data.frame(zero_growth_line_worst)  
  names(zero_growth_line_worst) <- c("IBI_infant_mortality_rate", "IBI_adult_mortality_rate")
  zero_growth_line_worst$scenario <- "worst"
  
  #bind to other
  zero_growth_line <- rbind(zero_growth_line, zero_growth_line_worst)
  
  #now attach some identifier information to this
  zero_growth_line$species <- unique(df$species)
  zero_growth_line$IBI[zero_growth_line$scenario=="best"] <- df$IBI[df$scenario=="best"]
  zero_growth_line$IBI[zero_growth_line$scenario=="worst"] <- df$IBI[df$scenario=="worst"]
  
  zero_growth_line$IBI_fecundity[zero_growth_line$scenario=="best"] <- df$IBI_fecundity[df$scenario=="best"]
  zero_growth_line$IBI_fecundity[zero_growth_line$scenario=="worst"] <- df$IBI_fecundity[df$scenario=="worst"]
  
  zero_growth_line$age_1st_rep[zero_growth_line$scenario=="best"] <- df$age_1st_rep[df$scenario=="best"]
  zero_growth_line$age_1st_rep[zero_growth_line$scenario=="worst"] <- df$age_1st_rep[df$scenario=="worst"]
  
  #zero_growth_line$avg_lifespan[zero_growth_line$scenario=="best"] <- df$avg_lifespan[df$scenario=="best"]
  #zero_growth_line$avg_lifespan[zero_growth_line$scenario=="worst"] <- df$avg_lifespan[df$scenario=="worst"]
  
  zero_growth_line$litter_size[zero_growth_line$scenario=="best"] <- df$litter_size[df$scenario=="best"]
  zero_growth_line$litter_size[zero_growth_line$scenario=="worst"] <- df$litter_size[df$scenario=="worst"]
  
  #and for good measure, include the annual rate too
  #zero_growth_line$lit_infant_annual_mortality_rate <- (df$stage_mortality_rate[df$stage=="infant"])/IBI
  #zero_growth_line$lit_juvenile_annual_mortality_rate <- (df$stage_mortality_rate[df$stage=="juvenile"])/IBI
  #zero_growth_line$lit_adult_annual_mortality_rate <- (df$stage_mortality_rate[df$stage=="adult"])/IBI
  
  
  return(zero_growth_line)
}
isocline.survival.wrap <- function(df, p_sequence){
  #new version essentially plots isocline wrap based on the age at first rep for each species
  
  zero_growth_line.list <- list()
  for(i in 1:2){
    zero_growth_line.list[[i]] <- lapply(p_sequence, get_isoclines_survival, IBI_fecundity=as.numeric(as.character(df$IBI_fecundity[i])), IBI=as.numeric(as.character(df$IBI[i])), age_1st_rep=as.numeric(as.character(df$age_1st_rep[i])))           
  }
  zero_growth_line <- matrix(unlist(zero_growth_line.list[[1]]), ncol=2, byrow=T)
  zero_growth_line <- data.frame(zero_growth_line)
  names(zero_growth_line) <- c("IBI_infant_survival_rate", "IBI_adult_survival_rate")
  zero_growth_line$scenario <- "best"
  zero_growth_line_worst <- matrix(unlist(zero_growth_line.list[[2]]), ncol=2, byrow=T)
  zero_growth_line_worst <- data.frame(zero_growth_line_worst)  
  names(zero_growth_line_worst) <- c("IBI_infant_survival_rate", "IBI_adult_survival_rate")
  zero_growth_line_worst$scenario <- "worst"
  
  #bind to other
  zero_growth_line <- rbind(zero_growth_line, zero_growth_line_worst)
  
  #now attach some identifier information to this
  zero_growth_line$species <- unique(df$species)
  zero_growth_line$IBI[zero_growth_line$scenario=="best"] <- df$IBI[df$scenario=="best"]
  zero_growth_line$IBI[zero_growth_line$scenario=="worst"] <- df$IBI[df$scenario=="worst"]
  
  zero_growth_line$IBI_fecundity[zero_growth_line$scenario=="best"] <- df$IBI_fecundity[df$scenario=="best"]
  zero_growth_line$IBI_fecundity[zero_growth_line$scenario=="worst"] <- df$IBI_fecundity[df$scenario=="worst"]
  
  zero_growth_line$age_1st_rep[zero_growth_line$scenario=="best"] <- df$age_1st_rep[df$scenario=="best"]
  zero_growth_line$age_1st_rep[zero_growth_line$scenario=="worst"] <- df$age_1st_rep[df$scenario=="worst"]
  
  #zero_growth_line$avg_lifespan[zero_growth_line$scenario=="best"] <- df$avg_lifespan[df$scenario=="best"]
  #zero_growth_line$avg_lifespan[zero_growth_line$scenario=="worst"] <- df$avg_lifespan[df$scenario=="worst"]
  
  zero_growth_line$litter_size[zero_growth_line$scenario=="best"] <- df$litter_size[df$scenario=="best"]
  zero_growth_line$litter_size[zero_growth_line$scenario=="worst"] <- df$litter_size[df$scenario=="worst"]
  
  #and for good measure, include the annual rate too
  #zero_growth_line$lit_infant_annual_mortality_rate <- (df$stage_mortality_rate[df$stage=="infant"])/IBI
  #zero_growth_line$lit_juvenile_annual_mortality_rate <- (df$stage_mortality_rate[df$stage=="juvenile"])/IBI
  #zero_growth_line$lit_adult_annual_mortality_rate <- (df$stage_mortality_rate[df$stage=="adult"])/IBI
  
  
  return(zero_growth_line)
}
get_isoclines <- function(p, IBI_fecundity, IBI, age_1st_rep){
  
  #first, choose 
  #first record timesteps of adulthood
  #dur_rep <- ((avg_lifespan/IBI)-(age_1st_rep/IBI)) #duration of reproduction in IBI timesteps
  #dur_mat <- ((avg_lifespan/IBI)-1) # durating of juvenile or adult class in IBI timesteps
  
  #and here is where we actually choose our proper IBI fecundity
  
  #and get all our variables in order
  #s was fed in through the function
  #then, adult stage reproductive output (must be multiplied by the probability of escaping this class)
  
  Fa = IBI_fecundity
  
  #our system somehwat breaks down in the case where age at 1st rep = 2 because the first birth cohort must st
  #*(dur_rep)
  #then, age at first rep - this has been carefully chosen to be an even number divisible by the IBI
  a = age_1st_rep/IBI
  
  #L = avg_lifespan/IBI
  
  #then, we solve for i:
  
  #i <- (((-p^(L))+sqrt(p^(2*L)+4*Fa*p^(a+1)))/(2*(Fa*p^(a+1))))
  i <- (1-p)/(Fa*(p^a))
  
  #however, if age at 1st rep=1, then we have no juvenile class and we have to make amends
  #in this case our equation is given by:
  #maybe not. I think our equation still holds. but the leslie matrix will change for sure
  #if(a==1){
   # i <- -(1/Fa)
  #}
  
  
  #adult_surv =  p #IBI adult mortality rate
  #infant_surv = (p*i)
  
  adult_mort =  1-p
  infant_mort = 1-(p*i)
  #and then return them!
  #return(c(infant_surv,adult_surv))
  return(c(infant_mort,adult_mort))
  #return(c(i,p)) #if you would rather plot annual mortaluty...
}
get_isoclines_survival <- function(p, IBI_fecundity, IBI, age_1st_rep){
  
  #first, choose 
  #first record timesteps of adulthood
  #dur_rep <- ((avg_lifespan/IBI)-(age_1st_rep/IBI)) #duration of reproduction in IBI timesteps
  #dur_mat <- ((avg_lifespan/IBI)-1) # durating of juvenile or adult class in IBI timesteps
  
  #and here is where we actually choose our proper IBI fecundity
  
  #and get all our variables in order
  #s was fed in through the function
  #then, adult stage reproductive output (must be multiplied by the probability of escaping this class)
  
  Fa = IBI_fecundity
  
  #our system somehwat breaks down in the case where age at 1st rep = 2 because the first birth cohort must st
  #*(dur_rep)
  #then, age at first rep - this has been carefully chosen to be an even number divisible by the IBI
  a = age_1st_rep/IBI
  
  #L = avg_lifespan/IBI
  
  #then, we solve for i:
  
  #i <- (((-p^(L))+sqrt(p^(2*L)+4*Fa*p^(a+1)))/(2*(Fa*p^(a+1))))
  i <- (1-p)/(Fa*(p^a))
  
  #however, if age at 1st rep=1, then we have no juvenile class and we have to make amends
  #in this case our equation is given by:
  #maybe not. I think our equation still holds. but the leslie matrix will change for sure
  #if(a==1){
  # i <- -(1/Fa)
  #}
  
  
  #adult_surv =  p #IBI adult mortality rate
  #infant_surv = (p*i)
  
  adult_surv =  p
  infant_surv = (p*i)
  #and then return them!
  #return(c(infant_surv,adult_surv))
  return(c(infant_surv,adult_surv))
  #return(c(i,p)) #if you would rather plot annual mortaluty...
}
wrap.isocline.by.multiplier <- function(dat1, type,  p_sequence, do.plot, do.save, filename){
  
  dat2 <- wrap.pop.data(dat1, type=type) 
  
  #now split to list
  dat.list <- dlply(dat2, .(species))
  
  #now make list of isoclines per species
  #start with
  isocline.list <- lapply(dat.list, isocline.wrap, p_sequence=p_sequence)
  #and bind them all together
  isocline.df <- do.call("rbind", isocline.list)
  rownames(isocline.df) <- c()
  
  #add a grouping ID for each line
  isocline.df$groupID <- paste(isocline.df$species, isocline.df$scenario, sep= "-")
  
  #remove any NAs
  isocline.df <- isocline.df[complete.cases(isocline.df),]
  
  #and set any rates below zero to zero:
  isocline.df$IBI_infant_mortality_rate[isocline.df$IBI_infant_mortality_rate < 0] <- 0
  isocline.df$IBI_infant_mortality_rate[isocline.df$IBI_adult_mortality_rate<0] <- 0
  
  #add a column for multiplier
  isocline.df$multiplier <- as.factor(isocline.df$age_1st_rep/isocline.df$IBI)
  isocline.df$IBI <- factor(isocline.df$IBI)
  isocline.df$age_1st_rep <- factor(isocline.df$age_1st_rep)
  #and we plot
  if(do.plot==TRUE){
    
    #and we'll shade above line and below, so let's add a mock variable here for each
    plot.upper.df <- subset(isocline.df, scenario=="best")
    plot.upper.df$ymax=.8
    #and the next box
    max.df <- ddply(plot.upper.df, .(species), summarize, x=c(max(IBI_adult_mortality_rate[IBI_infant_mortality_rate>0]),max(IBI_adult_mortality_rate[IBI_infant_mortality_rate>0]),.35,.35),  y=c(0,.8,.8,0))
    
    #and the lower 
    plot.lower.df <- subset(isocline.df, scenario=="worst")
    plot.lower.df$ymin = 0
    
    #and let's fill the blocks in our graph
    species <- rep(unique(isocline.df$species), each=4)
    x <- c(.83,.83,1,1,0,0,0,0,0,0,0,0,.93,.93,1,1,0,0,0,0,0,0,0,0,.53,.53,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    y <- rep(c(0,1,1,0), 11)
    group <- rep(seq(1,11,1),each=4)
    dat6 <- cbind.data.frame(species,x,y,group)
    
    #want to fill between each lines
    p <-  ggplot(data=isocline.df) + 
      geom_line(aes(x=IBI_adult_mortality_rate,y=IBI_infant_mortality_rate, group=groupID, linetype=multiplier, colour=age_1st_rep)) + facet_grid(species~.) +
      theme_classic() + ylim(c(0,1)) +xlim(c(0,1))+
      geom_segment(aes(x=0,xend=0,y=0, yend=1), colour="black", size=.1) +
      geom_segment(aes(x=0,xend=1,y=0, yend=0), colour="black", size=.1) +
      xlab("adult mortality rate per IBI") + ylab("infant mortality rate per IBI") +
      annotate("text", x=.8, y=.8, label="population extinction", colour="black", size=2.5) + annotate("text", x=.1, y=.2, label="population growth", size=2.5) +
      ggtitle("Zero-Growth Isoclines with Species-Specific Parameters")  
    
    print(p) 
    
    
    
  }
  return(list(isocline.df, p))
}
wrap.isocline.bats <- function(dat1, type, p_sequence, do.plot, do.save, filename){
  
  
  dat1$adult_annual_mortality <- 1- dat1$adult_annual_survival 
  dat1$infant_annual_mortality <- 1- dat1$infant_annual_survival 
  dat1$infant_annual_mortality[dat1$infant_annual_mortality<0 & dat1$scenario=="best"] <- 0
  dat1$infant_annual_mortality[dat1$infant_annual_mortality<0 & dat1$scenario=="worst"] <- dat1$infant_annual_mortality[dat1$species=="Eidolon dupreanum" & dat1$scenario=="worst"]
  
  #now split to list
  dat.list <- dlply(dat1, .(species))
  
  #now make list of isoclines per species
  #start with
  isocline.list <- lapply(dat.list, isocline.wrap, p_sequence=p_sequence)
  #and bind them all together
  isocline.df <- do.call("rbind", isocline.list)
  rownames(isocline.df) <- c()
  
  #add a grouping ID for each line
  isocline.df$groupID <- paste(isocline.df$species, isocline.df$scenario, sep= "-")
  
  #remove any NAs
  isocline.df <- isocline.df[complete.cases(isocline.df),]
  
  #and set any rates below zero to zero:
  isocline.df$IBI_infant_mortality_rate[isocline.df$IBI_infant_mortality_rate < 0] <- 0
  isocline.df$IBI_infant_mortality_rate[isocline.df$IBI_adult_mortality_rate<0] <- 0
  
  #add a column for multiplier
  isocline.df$multiplier <- as.factor(isocline.df$age_1st_rep/isocline.df$IBI)
  isocline.df$IBI <- factor(isocline.df$IBI)
  isocline.df$age_1st_rep <- factor(isocline.df$age_1st_rep)
  
  
  #and we plot
  if(do.plot==TRUE){
    
    #factor species by body size
    isocline.df$species <- factor(isocline.df$species, levels = c("Eidolon dupreanum", "Pteropus rufus"))
    
    #and we'll shade above line and below, so let's add a mock variable here for each
    plot.upper.df <- subset(isocline.df, scenario=="best")
    plot.upper.df$ymax=.8
    #and the next box
    max.df <- ddply(plot.upper.df, .(species), summarize, x=c(max(IBI_adult_mortality_rate[IBI_infant_mortality_rate>0]),max(IBI_adult_mortality_rate[IBI_infant_mortality_rate>0]),.35,.35),  y=c(0,.8,.8,0))
    
    #and the lower 
    plot.lower.df <- subset(isocline.df, scenario=="worst")
    plot.lower.df$ymin = 0
    
    #and let's fill the blocks in our graph
    species <- rep(unique(isocline.df$species), each=4)
    x <- c(.99,.99,1,1,.99,.99,1,1)
    y <- rep(c(0,1,1,0), 2)
    group <- rep(seq(1,2,1),each=4)
    dat6 <- cbind.data.frame(species,x,y,group)
    #dat6$y[dat6$species=="Microcebus sp." & dat6$x==.99 & dat6$y == 0] <- .34
    #dat6$y[dat6$species=="Microcebus sp." & dat6$x==1 & dat6$y == 0] <- .34
    
    #want to fill between each lines
    
    
    #redo the levels
    isocline.df$species <- factor(isocline.df$species, levels= c("Eidolon dupreanum", "Pteropus rufus"))
    
    plot.upper.df$species <- factor(plot.upper.df$species, levels= c("Eidolon dupreanum", "Pteropus rufus"))
    
    dat6$species <- factor(dat6$species, levels= c("Eidolon dupreanum", "Pteropus rufus"))
    
    dat1$adult_annual_mortality = 1- dat1$adult_annual_survival
    dat1$infant_annual_mortality = 1- dat1$infant_annual_survival
    
    
    quartz()
    
    p <-  ggplot(data=isocline.df) + 
      geom_line(aes(IBI_adult_mortality_rate,IBI_infant_mortality_rate, 
                    group=groupID), colour="white", alpha=0) +
      geom_ribbon(data =  plot.upper.df, aes(x= IBI_adult_mortality_rate, ymin=IBI_infant_mortality_rate, ymax=1), colour="black", fill="black") + #color above line
      geom_polygon(data=dat6, aes(x=x,y=y), color="black", fill="black") +
      geom_ribbon(data =  plot.lower.df, aes(x= IBI_adult_mortality_rate, ymin=ymin, ymax=IBI_infant_mortality_rate), fill="gray85", colour="black") + #color above line
      theme_classic() + ylim(c(0,1)) +xlim(c(0,1)) + 
      theme(strip.text = element_text(size=20, face = "italic"), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), 
            axis.text.x=element_text(size=17), axis.text.y=element_text(size=17)) + 
      geom_segment(aes(x=0,xend=0,y=0, yend=1), colour="black", size=.1) +
      geom_segment(aes(x=0,xend=1,y=0, yend=0), colour="black", size=.1) +
      xlab("adult annual mortality") + ylab("juvenile annual mortality") +
      annotate("text", x=.85, y=.7, label="population\nextinction", colour="white", size=5, fontface="bold") + annotate("text", x=.11, y=.25, label="population\ngrowth", size=5, fontface="bold") +
      facet_grid(species~.) +
      geom_segment(data=dat1, aes(x=adult_annual_mortality, y=0, xend=adult_annual_mortality, yend=1), linetype=2, size=1, color="magenta") +
      geom_point(data=dat1, aes(x=adult_annual_mortality, y=infant_annual_mortality), shape=18, size=6, color="magenta")
    
    print(p)  
    
    if(do.save==TRUE){
      
      ggsave(file = filename,
             units="mm",  
             width=60, 
             height=60, 
             scale=3, 
             dpi=300)
    }
    
  }
  return(list(isocline.df, p))
}
wrap.isocline.bats.raw.mort <- function(dat1, type, p_sequence, do.plot, do.save, filename){
  
  
  dat1$adult_annual_mortality <- 1- dat1$adult_annual_survival 
  dat1$infant_annual_mortality <- 1- dat1$infant_annual_survival 
  dat1$infant_annual_mortality[dat1$infant_annual_mortality<0 & dat1$scenario=="best"] <- 0
  dat1$infant_annual_mortality[dat1$infant_annual_mortality<0 & dat1$scenario=="worst"] <- dat1$infant_annual_mortality[dat1$species=="Eidolon dupreanum" & dat1$scenario=="worst"]
  
  #now split to list
  dat.list <- dlply(dat1, .(species))
  
  #now make list of isoclines per species
  #start with
  isocline.list <- lapply(dat.list, isocline.wrap, p_sequence=p_sequence)
  #and bind them all together
  isocline.df <- do.call("rbind", isocline.list)
  rownames(isocline.df) <- c()
  
  #add a grouping ID for each line
  isocline.df$groupID <- paste(isocline.df$species, isocline.df$scenario, sep= "-")
  
  #remove any NAs
  isocline.df <- isocline.df[complete.cases(isocline.df),]
  
  #and set any rates below zero to zero:
  isocline.df$IBI_infant_mortality_rate[isocline.df$IBI_infant_mortality_rate < 0] <- 0
  isocline.df$IBI_infant_mortality_rate[isocline.df$IBI_adult_mortality_rate<0] <- 0
  
  #add a column for multiplier
  isocline.df$multiplier <- as.factor(isocline.df$age_1st_rep/isocline.df$IBI)
  isocline.df$IBI <- factor(isocline.df$IBI)
  isocline.df$age_1st_rep <- factor(isocline.df$age_1st_rep)
  
  
  #and we plot
  if(do.plot==TRUE){
    
    #factor species by body size
    isocline.df$species <- factor(isocline.df$species, levels = c("Eidolon dupreanum", "Pteropus rufus"))
    
    #and we'll shade above line and below, so let's add a mock variable here for each
    plot.upper.df <- subset(isocline.df, scenario=="best")
    plot.upper.df$ymax=.8
    #and the next box
    max.df <- ddply(plot.upper.df, .(species), summarize, x=c(max(IBI_adult_mortality_rate[IBI_infant_mortality_rate>0]),max(IBI_adult_mortality_rate[IBI_infant_mortality_rate>0]),.35,.35),  y=c(0,.8,.8,0))
    
    #and the lower 
    plot.lower.df <- subset(isocline.df, scenario=="worst")
    plot.lower.df$ymin = 0
    
    #and let's fill the blocks in our graph
    species <- rep(unique(isocline.df$species), each=4)
    x <- c(.99,.99,1,1,.99,.99,1,1)
    y <- rep(c(0,1,1,0), 2)
    group <- rep(seq(1,2,1),each=4)
    dat6 <- cbind.data.frame(species,x,y,group)
    #dat6$y[dat6$species=="Microcebus sp." & dat6$x==.99 & dat6$y == 0] <- .34
    #dat6$y[dat6$species=="Microcebus sp." & dat6$x==1 & dat6$y == 0] <- .34
    
    #want to fill between each lines
    
    
    #redo the levels
    isocline.df$species <- factor(isocline.df$species, levels= c("Eidolon dupreanum", "Pteropus rufus"))
    
    plot.upper.df$species <- factor(plot.upper.df$species, levels= c("Eidolon dupreanum", "Pteropus rufus"))
    
    dat6$species <- factor(dat6$species, levels= c("Eidolon dupreanum", "Pteropus rufus"))
    
    
    
    quartz()
    
    p <-  ggplot(data=isocline.df) + 
      geom_line(aes(IBI_adult_mortality_rate,IBI_infant_mortality_rate, 
                    group=groupID), colour="white", alpha=0) +
      geom_ribbon(data =  plot.upper.df, aes(x= IBI_adult_mortality_rate, ymin=IBI_infant_mortality_rate, ymax=1), colour="black", fill="black") + #color above line
      geom_polygon(data=dat6, aes(x=x,y=y), color="black", fill="black") +
      geom_ribbon(data =  plot.lower.df, aes(x= IBI_adult_mortality_rate, ymin=ymin, ymax=IBI_infant_mortality_rate), fill="gray85", colour="black") + #color above line
      theme_classic() + ylim(c(0,1)) +xlim(c(0,1)) + 
      theme(strip.text = element_text(size=20, face = "italic"), axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), 
            axis.text.x=element_text(size=17), axis.text.y=element_text(size=17)) + 
      geom_segment(aes(x=0,xend=0,y=0, yend=1), colour="black", size=.1) +
      geom_segment(aes(x=0,xend=1,y=0, yend=0), colour="black", size=.1) +
      xlab("adult annual mortality") + ylab("juvenile annual mortality") +
      annotate("text", x=.85, y=.7, label="population\nextinction", colour="white", size=5, fontface="bold") + annotate("text", x=.11, y=.25, label="population\ngrowth", size=5, fontface="bold") +
      facet_grid(species~.) #+
      #geom_segment(data=dat1, aes(x=adult_annual_mortality, y=0, xend=adult_annual_mortality, yend=1), linetype=2, size=1, color="magenta") +
      #geom_point(data=dat1, aes(x=adult_annual_mortality, y=infant_annual_mortality), shape=18, size=6, color="magenta")
    
    print(p)  
    
    if(do.save==TRUE){
      
      ggsave(file = filename,
             units="mm",  
             width=55, 
             height=60, 
             scale=3, 
             dpi=300)
    }
    
  }
  return(list(isocline.df, p))
}
wrap.isocline.bats.survival <- function(dat1, type, p_sequence, do.plot, do.save, filename){
  
  
  dat1$adult_IBI_mortality <- 1- dat1$adult_annual_survival 
  dat1$infant_IBI_mortality <- 1- dat1$infant_annual_survival 
  dat1$infant_IBI_mortality[dat1$infant_IBI_mortality<0 & dat1$scenario=="best"] <- 0
  dat1$infant_IBI_mortality[dat1$infant_IBI_mortality<0 & dat1$scenario=="worst"] <- dat1$infant_IBI_mortality[dat1$species=="Eidolon dupreanum" & dat1$scenario=="worst"]
  
  #now split to list
  dat.list <- dlply(dat1, .(species))
  
  #now make list of isoclines per species
  #start with
  isocline.list <- lapply(dat.list, isocline.survival.wrap, p_sequence=p_sequence)
  #and bind them all together
  isocline.df <- do.call("rbind", isocline.list)
  rownames(isocline.df) <- c()
  
  #add a grouping ID for each line
  isocline.df$groupID <- paste(isocline.df$species, isocline.df$scenario, sep= "-")
  
  #remove any NAs
  isocline.df <- isocline.df[complete.cases(isocline.df),]
  
  #and set any rates below zero to zero:
  isocline.df$IBI_infant_survival_rate[isocline.df$IBI_infant_surival_rate < 0] <- 0
  isocline.df$IBI_adult_survival_rate[isocline.df$IBI_adult_survival_rate<0] <- 0
  
  #add a column for multiplier
  isocline.df$multiplier <- as.factor(isocline.df$age_1st_rep/isocline.df$IBI)
  isocline.df$IBI <- factor(isocline.df$IBI)
  isocline.df$age_1st_rep <- factor(isocline.df$age_1st_rep)
  
  
  #and we plot
  if(do.plot==TRUE){
    
    #factor species by body size
    isocline.df$species <- factor(isocline.df$species, levels = c("Eidolon dupreanum", "Pteropus rufus"))
    
    #cull those that are over 1
    #isocline.df = subset(isocline.df, IBI_infant_survival_rate<=1)
    
    
    #and we'll shade above line and below, so let's add a mock variable here for each
    #these get reverse from mortality
    plot.lower.df <- subset(isocline.df, scenario=="best") #lowest survival rates to allow for persistence
    #plot.lower.df$ymin=0
    #and the next box
    #max.df <- ddply(plot.lower.df, .(species), summarize, x=c(max(IBI_adult_survival_rate[IBI_infant_survival_rate>0]),max(IBI_adult_survival_rate[IBI_infant_survival_rate>0]),.35,.35),  y=c(0,.8,.8,0))
    
    #and the lower 
    plot.upper.df <- subset(isocline.df, scenario=="worst")
   # plot.upper.df$ymax = 0
    
    
    
    
    #want to fill between each lines
    
    
    #redo the levels
    isocline.df$species <- factor(isocline.df$species, levels= c("Eidolon dupreanum", "Pteropus rufus"))
    
    plot.upper.df$species <- factor(plot.upper.df$species, levels= c("Eidolon dupreanum", "Pteropus rufus"))
    
    
    
    #you want the adult survival rate where infant survival = 1 (or as close as possible)
    #adult_min_df = subset(isocline.df, IBI_infant_survival_rate<=1 & scenario == "worst")
    #infant_max_df = ddply(adult_min_df, .(species), summarize, IBI_infant_survival_rate = max(IBI_infant_survival_rate))#, IBI_adult_survival_rate = unique(IBI_adult_survival_rate))
     
    #find minimum survival required per species
    #adult.min.Eid = adult_min_df$IBI_adult_survival_rate[ adult_min_df$species == "Eidolon dupreanum" & adult_min_df$IBI_infant_survival_rate==infant_max_df$IBI_infant_survival_rate[infant_max_df$species=="Eidolon dupreanum"]]
    #adult.min.Pter = adult_min_df$IBI_adult_survival_rate[adult_min_df$species == "Pteropus rufus" & adult_min_df$IBI_infant_survival_rate==infant_max_df$IBI_infant_survival_rate[infant_max_df$species=="Pteropus rufus"]]
                                                         
    
 
    #and let's fill the blocks in our graph
    #species <- rep(unique(isocline.df$species), each=4)
    #x.Eid <- c(adult.min.Eid ,adult.min.Eid ,0,0)
    #x.Pter <- c(adult.min.Pter ,adult.min.Pter ,0,0)
    #x = c(x.Eid, x.Pter)
    #y <- rep(c(0,1,1,0), 2)
    #group <- rep(seq(1,2,1),each=4)
    #dat6 <- cbind.data.frame(species,x,y,group)
    #dat6$species <- factor(dat6$species, levels= c("Eidolon dupreanum", "Pteropus rufus"))
    #dat6$y[dat6$species=="Microcebus sp." & dat6$x==.99 & dat6$y == 0] <- .34
    #dat6$y[dat6$species=="Microcebus sp." & dat6$x==1 & dat6$y == 0] <- .34
    
    
    #for max, make it just above the infant mort
    dat1 = subset(dat1, scenario=="best")
    #then overwrite Pteropus to use Eidolon
    #dat1$infant_annual_survival[dat1$species=="Pteropus rufus"] = dat1$infant_annual_survival[dat1$species=="Eidolon dupreanum"]
    
    quartz()
    
    p <-  ggplot(data=isocline.df) + facet_grid(species~.) +
     geom_ribbon(data =  plot.upper.df, aes(x= IBI_adult_survival_rate, ymin=IBI_infant_survival_rate, ymax=1), colour="black", fill="gray85") + #color above line
      geom_ribbon(data =  plot.lower.df, aes(x= IBI_adult_survival_rate, ymin=0, ymax=IBI_infant_survival_rate), fill="black", colour="black") + #color below line
      geom_line(aes(IBI_adult_survival_rate,IBI_infant_survival_rate, 
                    group=groupID), colour="white", alpha=0) +
      theme_classic()  +coord_cartesian(xlim = c(.25,1), ylim = c(0,1), expand=FALSE) +
      theme(strip.text = element_text(size=20, face = "italic"), 
            axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), 
            axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), panel.spacing = unit(.5, "cm")) + 
      xlab("adult annual survival") + ylab("juvenile annual survival") +
      annotate("text", x=.93, y=.8, label="population\ngrowth", size=4.5, fontface="bold") + 
      annotate("text", x=.4, y=.5, label="population\nextinction",colour="white", size=4.5, fontface="bold") +
      geom_segment(data=dat1, aes(x=adult_annual_survival, y=0, xend=adult_annual_survival, yend=1), linetype=2, size=1, color="magenta") +
      geom_point(data=dat1, aes(x=adult_annual_survival, y=infant_annual_survival), shape=18, size=6, color="magenta")
    
    print(p)  
    
    if(do.save==TRUE){
      
      ggsave(file = filename,
             units="mm",  
             width=55, 
             height=60, 
             scale=3, 
             dpi=300)
    }
    
  }
  return(list(isocline.df, p))
}
wrap.isocline.bats.survival.BW <- function(dat1, type, p_sequence, do.plot, do.save, filename){
  
  
  dat1$adult_IBI_mortality <- 1- dat1$adult_annual_survival 
  dat1$infant_IBI_mortality <- 1- dat1$infant_annual_survival 
  dat1$infant_IBI_mortality[dat1$infant_IBI_mortality<0 & dat1$scenario=="best"] <- 0
  dat1$infant_IBI_mortality[dat1$infant_IBI_mortality<0 & dat1$scenario=="worst"] <- dat1$infant_IBI_mortality[dat1$species=="Eidolon dupreanum" & dat1$scenario=="worst"]
  
  #now split to list
  dat.list <- dlply(dat1, .(species))
  
  #now make list of isoclines per species
  #start with
  isocline.list <- lapply(dat.list, isocline.survival.wrap, p_sequence=p_sequence)
  #and bind them all together
  isocline.df <- do.call("rbind", isocline.list)
  rownames(isocline.df) <- c()
  
  #add a grouping ID for each line
  isocline.df$groupID <- paste(isocline.df$species, isocline.df$scenario, sep= "-")
  
  #remove any NAs
  isocline.df <- isocline.df[complete.cases(isocline.df),]
  
  #and set any rates below zero to zero:
  isocline.df$IBI_infant_survival_rate[isocline.df$IBI_infant_surival_rate < 0] <- 0
  isocline.df$IBI_adult_survival_rate[isocline.df$IBI_adult_survival_rate<0] <- 0
  
  #add a column for multiplier
  isocline.df$multiplier <- as.factor(isocline.df$age_1st_rep/isocline.df$IBI)
  isocline.df$IBI <- factor(isocline.df$IBI)
  isocline.df$age_1st_rep <- factor(isocline.df$age_1st_rep)
  
  
  #and we plot
  if(do.plot==TRUE){
    
    #factor species by body size
    isocline.df$species <- factor(isocline.df$species, levels = c("Eidolon dupreanum", "Pteropus rufus"))
    
    #cull those that are over 1
    #isocline.df = subset(isocline.df, IBI_infant_survival_rate<=1)
    
    
    #and we'll shade above line and below, so let's add a mock variable here for each
    #these get reverse from mortality
    plot.lower.df <- subset(isocline.df, scenario=="best") #lowest survival rates to allow for persistence
    #plot.lower.df$ymin=0
    #and the next box
    #max.df <- ddply(plot.lower.df, .(species), summarize, x=c(max(IBI_adult_survival_rate[IBI_infant_survival_rate>0]),max(IBI_adult_survival_rate[IBI_infant_survival_rate>0]),.35,.35),  y=c(0,.8,.8,0))
    
    #and the lower 
    plot.upper.df <- subset(isocline.df, scenario=="worst")
    # plot.upper.df$ymax = 0
    
    
    
    
    #want to fill between each lines
    
    
    #redo the levels
    isocline.df$species <- factor(isocline.df$species, levels= c("Eidolon dupreanum", "Pteropus rufus"))
    
    plot.upper.df$species <- factor(plot.upper.df$species, levels= c("Eidolon dupreanum", "Pteropus rufus"))
    
    
    
    #you want the adult survival rate where infant survival = 1 (or as close as possible)
    #adult_min_df = subset(isocline.df, IBI_infant_survival_rate<=1 & scenario == "worst")
    #infant_max_df = ddply(adult_min_df, .(species), summarize, IBI_infant_survival_rate = max(IBI_infant_survival_rate))#, IBI_adult_survival_rate = unique(IBI_adult_survival_rate))
    
    #find minimum survival required per species
    #adult.min.Eid = adult_min_df$IBI_adult_survival_rate[ adult_min_df$species == "Eidolon dupreanum" & adult_min_df$IBI_infant_survival_rate==infant_max_df$IBI_infant_survival_rate[infant_max_df$species=="Eidolon dupreanum"]]
    #adult.min.Pter = adult_min_df$IBI_adult_survival_rate[adult_min_df$species == "Pteropus rufus" & adult_min_df$IBI_infant_survival_rate==infant_max_df$IBI_infant_survival_rate[infant_max_df$species=="Pteropus rufus"]]
    
    
    
    #and let's fill the blocks in our graph
    #species <- rep(unique(isocline.df$species), each=4)
    #x.Eid <- c(adult.min.Eid ,adult.min.Eid ,0,0)
    #x.Pter <- c(adult.min.Pter ,adult.min.Pter ,0,0)
    #x = c(x.Eid, x.Pter)
    #y <- rep(c(0,1,1,0), 2)
    #group <- rep(seq(1,2,1),each=4)
    #dat6 <- cbind.data.frame(species,x,y,group)
    #dat6$species <- factor(dat6$species, levels= c("Eidolon dupreanum", "Pteropus rufus"))
    #dat6$y[dat6$species=="Microcebus sp." & dat6$x==.99 & dat6$y == 0] <- .34
    #dat6$y[dat6$species=="Microcebus sp." & dat6$x==1 & dat6$y == 0] <- .34
    
    
    #for max, make it just above the infant mort
    dat1 = subset(dat1, scenario=="best")
    #then overwrite Pteropus to use Eidolon
    #dat1$infant_annual_survival[dat1$species=="Pteropus rufus"] = dat1$infant_annual_survival[dat1$species=="Eidolon dupreanum"]
    
    quartz()
    
    p <-  ggplot(data=isocline.df) + facet_grid(species~.) +
      geom_ribbon(data =  plot.upper.df, aes(x= IBI_adult_survival_rate, ymin=IBI_infant_survival_rate, ymax=1), colour="black", fill="gray85") + #color above line
      geom_ribbon(data =  plot.lower.df, aes(x= IBI_adult_survival_rate, ymin=0, ymax=IBI_infant_survival_rate), fill="black", colour="black") + #color below line
      geom_line(aes(IBI_adult_survival_rate,IBI_infant_survival_rate, 
                    group=groupID), colour="white", alpha=0) +
      theme_classic()  +coord_cartesian(xlim = c(.25,1), ylim = c(0,1), expand=FALSE) +
      theme(strip.text = element_text(size=20, face = "italic"), 
            axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), 
            axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), panel.spacing = unit(.5, "cm")) + 
      xlab("adult annual survival") + ylab("juvenile annual survival") +
      annotate("text", x=.93, y=.8, label="population\ngrowth", size=4.5, fontface="bold") + 
      annotate("text", x=.4, y=.5, label="population\nextinction",colour="white", size=4.5, fontface="bold") +
      geom_segment(data=dat1, aes(x=adult_annual_survival, y=0, xend=adult_IBI_survival, yend=1), linetype=1, size=1.1, color="black") +
      geom_segment(data=dat1, aes(x=adult_annual_survival, y=0, xend=adult_IBI_survival, yend=1), linetype=2, size=1, color="gray85") +
      geom_point(data=dat1, aes(x=adult_annual_survival, y=infant_IBI_survival), shape=18, size=7, color="black") +
      geom_point(data=dat1, aes(x=adult_annual_survival, y=infant_IBI_survival), shape=18, size=6, color="gray85")
    
    
    print(p)  
    
    if(do.save==TRUE){
      
      ggsave(file = filename,
             units="mm",  
             width=55, 
             height=60, 
             scale=3, 
             dpi=300)
    }
    
  }
  return(list(isocline.df, p))
}

#rename since IBI for bats is always 1
names(bat.demo)[names(bat.demo)=="adult_IBI_survival"] <- "adult_annual_survival"
names(bat.demo)[names(bat.demo)=="infant_IBI_survival"] <- "infant_annual_survival"

wrap.isocline.bats.raw.mort(dat1=bat.demo, type="isocline", p_sequence=seq(0, 1, .01), 
                            do.plot=TRUE, do.save=TRUE, filename=paste0(homewd, "/prior-scripts/BiolConsBat/figures/raw_bat_isoclines_mort.pdf"))

wrap.isocline.bats(dat1=bat.demo, type="isocline", p_sequence=seq(0, 1, .01), do.plot=TRUE, 
                   do.save=TRUE, filename=paste0(homewd, "/prior-scripts/BiolConsBat/figures/bat_isoclines_mort.pdf"))

wrap.isocline.bats.survival(dat1=bat.demo, type="isocline", p_sequence=seq(0, 1, .01), do.plot=TRUE, do.save=FALSE, filename=NA)




wrap.isocline.bats.raw.mort(dat1=bat.demo, type="isocline", p_sequence=seq(0, 1, .01), 
                   do.plot=TRUE, do.save=TRUE, filename=paste0(homewd, "/prior-scripts/BiolConsBat/figures/raw_bat_isoclines_mort.pdf"))

wrap.isocline.bats.survival(dat1=bat.demo, type="isocline", p_sequence=seq(0, 1, .01), 
                   do.plot=TRUE, do.save=TRUE, filename=paste0(homewd, "/prior-scripts/BiolConsBat/figures/Fig2_bat_survival_isoclines.pdf"))


wrap.isocline.bats.survival(dat1=bat.demo, type="isocline", p_sequence=seq(0, 1, .01), 
                            do.plot=TRUE, do.save=TRUE, filename=paste0(homewd, "/prior-scripts/BiolConsBat/figures/Fig2_bat_survival_isoclines.pdf"))



wrap.isocline.bats.survival(dat1=bat.demo, type="isocline", p_sequence=seq(0, 1, .01), 
                            do.plot=TRUE, do.save=TRUE, filename=paste0(homewd, "/prior-scripts/BiolConsBat/figures/Fig2_bat_survival_isoclines.pdf"))


wrap.isocline.bats.survival.BW(dat1=bat.demo, type="isocline", p_sequence=seq(0, 1, .01), 
                            do.plot=TRUE, do.save=TRUE, filename=paste0(homewd, "/prior-scripts/BiolConsBat/figures/Fig2_bat_survival_isoclines_BW.pdf"))



#get lambda from population estimates
mat.Eid = matrix(c(0,(.48*.793),.544,.793 ), byrow=T, nrow=2)
mat.Pter = matrix(c(0,(.48*.511),.544,.511 ), byrow=T, nrow=2)
Re(eigen(mat.Eid)$values)[1] #1
Re(eigen(mat.Pter)$values)[1] #.70


#and for no range (perfect isocline)
bat.demo = subset(bat.demo, scenario=="best")
bat.demo2 = bat.demo
bat.demo2$scenario = "worst"

bat.demo = rbind(bat.demo, bat.demo2)

wrap.isocline.bats.raw.mort(dat1=bat.demo, type="isocline", p_sequence=seq(0, 1, .01), 
                            do.plot=TRUE, do.save=TRUE, filename=paste0(homewd, "/prior-scripts/BiolConsBat/figures/raw_bat_isoclines_mort_one_line.pdf"))

wrap.isocline.bats(dat1=bat.demo, type="isocline", p_sequence=seq(0, 1, .01), do.plot=TRUE, 
                   do.save=TRUE, filename=paste0(homewd, "/prior-scripts/BiolConsBat/figures/bat_isoclines_mort_one_line.pdf"))





  