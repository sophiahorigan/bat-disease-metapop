########################
## SIMULATE MSIRN ##
#######################

# Created by Sophia Horigan (shorigan@uchicago.edu)
# Started: 01-17-24
# Last updated: 01-19-24

# Goal: Get MSIRM model working 

# Project Overview: This project seeks to use existing demographic data and the current best-fit transmission model
# to explore persistence threshholds for hypothetical pathogens under data-based metapopulation structure of Pteropus
# rufus in Madagascar. 

## PART 1 : SIMULATE MSIRN MODEL
# set model parameters
# function : run model
# plotting : 
# statistics : AP (annual persistence: with prob > 50% infection will persist in population for 1 year)
#           : LP (long-term persisence: with prob > 50% infection will persist in population for 100 years)

#rm(list=ls())

# Accomplishments
 # 01-17-24 removed all model simulation except MSRIN, and removed extraneous functions 
 # 01-19-24 modified main simulation matrix to also produce a dataframe tracking age class dynamics, including changing variable names


# To Do
 # Figure out how to return both disease and population dataframes from main function (different dims)
 # move helper functions to a different script to get them out of the way?

library(dplyr)
library(plyr)
library(cowplot)
library(deSolve)
library(lubridate)
library(matrixcalc)
library(Matrix)
library(ggplot2)
library(reshape2)
library(beepr)


# helper functions

build.pop.mat = function(surv, surv_juv, s, adult_fec){ 
  pop.mat = matrix(0,  nrow=(s-1), ncol = (s-1))
  diag(pop.mat) = surv
  diag(pop.mat)[1] = surv_juv
  col_s = c(rep(0, s-2), surv)
  pop.mat = cbind(pop.mat, col_s)
  row1 = c(0, rep((adult_fec*surv), (s-1))) #bats reproduce for the first time at the end of the second year of life. good for E. dup and P. ruf
  pop.mat = rbind(row1,pop.mat)
  return(pop.mat)
  
}#for stable age distribution

get.age.struct.M = function(pop, c){ #function that collapses population vector in disease state form to age structure
  age.mat <- mat_split(M=matrix(pop, ncol=1), r=c, c=1)
  age.mat.list <- c()
  for (i in 1:dim(age.mat)[3]){
    age.mat.list[[i]] <- age.mat[,,i]
  }
  #zage.mat <- Reduce('+', age.mat.list)
  age.mat <- lapply(age.mat.list, sum)
  age.mat.dat <- c( c(1:length(age.mat)), unlist(age.mat))
  names(age.mat.dat) <- c("age", "pop")
  return(age.mat)
}

find.biweek = function(t, times){
  
  biwks <- sort(unique(round(revtrunc(times),4)))
  this.wk = round(revtrunc(times[t]), 4)
  this.biwk <- which(biwks==this.wk)
  
  
  return(this.biwk)
  
  
}

revtrunc = function(x){
  newx = x - floor(x)
  return(newx)
}

buildTMat_MSIRN_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, mu.sick, rho, beta, recov, sigma, wane, add.inf.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I comes #3
  #foi = 1-exp(-beta*sum(I_m))
  foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  
  
  mat1 <- matrix(0,5,5) # MSIRN
  
  Tmat <- matrix(0,5*nage,5*nage) #MSIRN for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[3,5] <- rho
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1-rho
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],5);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    if (add.inf.mort==TRUE){
      surv[3] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
      #surv[4] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*5+1):(j*5+5),((j-1)*5+1):(j*5)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}

buildFMatrix_MSIRN <- function(age.classes, adult_fec, surv.biwk, biwk, N_stat){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  
  if(biwk==1){ #peak births
    new.fec = adult_fec*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = adult_fec*surv.biwk*.25
  }else if (biwk==3|biwk==25){
    new.fec = adult_fec*surv.biwk*.1
    #}else if (biwk==4|biwk ==24){
    #new.fec = adult_fec*surv.biwk*.05
  }else{
    new.fec = 0
  }
  s <- nage <- length(age.classes)
  
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  
  
  #make matrix the same size as the transition
  Fmat <- matrix(0,5*nage,5*nage)
  
  for (j in 1:nage) { #no fertility in first age class (mat immune)
    if(N_stat=="matAB"){
      #try it assuming N moms produce maternally immune pups
      Fmat[1,((j-1)*5+1):(j*5)] <- c(0,0,fert.biwk[j],fert.biwk[j], fert.biwk[j]) #the mat immune (so mom = I and R but not N)
      Fmat[2,((j-1)*5+1):(j*5)]<-  c(fert.biwk[j],fert.biwk[j],0,0, 0) #the susceptible (mom didn't get sick, so mom= N and S)
      
    }else if (N_stat=="matSus"){
      Fmat[1,((j-1)*5+1):(j*5)] <- c(0,0,fert.biwk[j],fert.biwk[j], 0) #the mat immune (so mom = I and R but not N)
      Fmat[2,((j-1)*5+1):(j*5)]<- c(fert.biwk[j],fert.biwk[j],0,0, fert.biwk[j]) #the susceptible (mom didn't get sick, so mom= N and S)
    }
    
    
    
  }
  
  return(Fmat)
}

transform.vect <- function(vec, s, c){ # vector, rows, columns
  vec2 <- t(commutation.matrix(r=s,c=c))%*%vec
  return(vec2)
}

mat_split <- function(M, r, c){
  nr <- ceiling(nrow(M)/r)
  nc <- ceiling(ncol(M)/c)
  newM <- matrix(NA, nr*r, nc*c)
  newM[1:nrow(M), 1:ncol(M)] <- M
  
  div_k <- kronecker(matrix(seq_len(nr*nc), nr, byrow = TRUE), matrix(1, r, c))
  matlist <- split(newM, div_k)
  N <- length(matlist)
  mats <- unlist(matlist)
  dim(mats)<-c(r, c, N)
  return(mats)
}#for matrix slicing. give it the number of rows and columns you want in the resulting matrices

sim.met.MSIRN.age <- function(burnin, sim_pop, yrs, ntyr, s, beta, age.brk, recov, mort, mort_juv, adult_fec, wane, rho, slope.wane, sigma, mu.sick, add.inf.mort, N_stat, printpop){
  
  c=5 # num columns - in this case M, S, I, R, N
  
  #can take more timesteps than the data
  times <- seq(0, yrs, by = 1/ntyr) 

  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  N_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRN_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, rho=rho,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRN(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  ##################################
  ## Create vector for age classes
  #################################


  N_pop_pop <- transform.vect(vec=N_pop_ts, s=s, c=c)
  print(dim(N_pop_pop))
  
  
  N_split_2 = mat_split(N_pop_pop, r=5, c=ncol(N_pop_pop))
  print(dim(N_split_2))
  
  #transform array into list
  N_split_pop = list()
  print(dim(N_split_2)[3])
  for (i in 1:dim(N_split_2)[3]){
    N_split_pop[[i]] = N_split_2[,,i]
  }
  
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total_pop = lapply(X=N_split_pop, FUN=colSums)
  print(length(N_total_pop))
  # print(N_total)
  

  
  ###############################################
  ## Create vector for infection state variables
  ###############################################
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  print(dim(N_pop_ts))
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  print(dim(N_split_1))
  
  #transform array into list
  N_split = list()
  print(dim(N_split_1)[3])
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  print(length(N_total))
 # print(N_total)

  
  ######################################################
  ## Compile infection data into dataframe and return
  ######################################################
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S", "I", "R", "N")
  dat.tot$tot_pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #print(dat.tot$tot_pop)
  
  dat.long <- melt(dat.tot, measure.vars = c("M", "S", "I", "R", "N", "tot_pop"), variable.name = "state", value.name = "count")
  dat.sub = subset(dat.long, state!="tot_pop")
  dat.pop = subset(dat.long, state=="tot_pop")
  
  dat.pop <- dplyr::select(dat.pop, -(state))
  names(dat.pop)[names(dat.pop)=="count"] <- "tot_pop"
  
  dat.out <- merge(dat.sub, dat.pop, by ="time", all.x = T, sort = F)
  head(dat.out)
  dat.out$proportion <- dat.out$count/dat.out$tot_pop
  
  dat.out$state = factor(dat.out$state, levels=c("M", "S", "I", "R", "N"))
  
  dat.out = subset(dat.out, time >= burnin)
  
  
  
  ######################################################
  ## Compile population data into dataframe and return
  ######################################################
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot.pop = cbind.data.frame(times,N_total_pop)
  names(dat.tot.pop) = c("time", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")
  dat.tot.pop$tot_pop = rowSums((dat.tot.pop[,2:ncol(dat.tot.pop)]))
  #print(dat.tot.pop$tot_pop)
  
  dat.long.pop <- melt(dat.tot.pop, measure.vars = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "tot_pop"), variable.name = "age", value.name = "count")
  dat.sub.pop = subset(dat.long.pop, age!="tot_pop")
  dat.pop.pop = subset(dat.long.pop, age=="tot_pop")
  
  dat.pop.pop <- dplyr::select(dat.pop.pop, -(age))
  names(dat.pop.pop)[names(dat.pop.pop)=="count"] <- "tot_pop"
  
  dat.out.pop <- merge(dat.sub.pop, dat.pop.pop, by ="time", all.x = T, sort = F)
  head(dat.out.pop)
  dat.out.pop$proportion <- dat.out.pop$count/dat.out.pop$tot_pop
  
  dat.out.pop$age = factor(dat.out.pop$age, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20"))
  
  # assign adult and juvenile categorization
  for(i in 1:nrow(dat.out.pop)){
    if(dat.out.pop[i,"age"] == 1){dat.out.pop[i, "class"] = 'J'}
    if(dat.out.pop[i,"age"] != 1){dat.out.pop[i, "class"] = 'A'}
  }
  
  
  # drop burnin
  dat.out.pop = subset(dat.out.pop, time >= burnin)
  
  
  # CAN ONLY RETURN ONE ARRAY IN R SO DUMB
  # NEED TO FIGURE OUT HOW TO COMBINE THEM BAH
  if(printpop == TRUE){
    return(dat.out.pop)
  }
  else{
    return(dat.out)
  }

}






run.sim.MSIRN <- sim.met.MSIRN.age(burnin=0,                  # years to ditch at start of sim
                                   sim_pop=1000,               # population size
                                   yrs=20,                     # years for simulation
                                   ntyr=26,                   # model timestep - greatly changes dynamics!
                                   s=20,                       # number of age classes
                                   beta = 2.203968,            # transmission
                                   age.brk = 20,               # how many distinct age classes there are
                                   recov=1,                    # recovery rate
                                   mort=.207,                  # adult natural mortality
                                   mort_juv=0.456,             # mortality juvenile
                                   adult_fec=.48,              # adult fecundity
                                   wane=0.0850142487956748,    # rate of waning maternal immunity - fixed
                                   slope.wane=1,               # slope for waning maternal immunity - logistic regression
                                   sigma=0.00748049787194744,  # waning adult humoral immunity 
                                   rho=0,                      # this is new from the paper, allows for N-class individuals to wane back to I (not S) so represents a persistent infection. in the future, could have this occur only seasonally
                                   mu.sick = 1,                # increased mortality rate for infected individuals
                                   add.inf.mort = FALSE,       # increased mortatlity rate for infected juveniles also?
                                   N_stat = "matAB",             # N class mothers produce M class young
                                   printpop = TRUE)            # print age class dynamics instead of infection dynamics
beep()

#######################
### PLOTTING INFECTION
#######################

colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N' = "navy")

p1 <- ggplot(data=run.sim.MSIRN) + geom_line(aes(x=time, y=proportion, color=state)) + 
  theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
  scale_color_manual(values=colz) +
  xlab("Year")
p1 


#######################
## PLOTTING POPULATION
######################

# by age
p2 <- ggplot(data=run.sim.MSIRN) + geom_line(aes(x=time, y=count, color=age)) + 
  theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
  #scale_color_manual(values=colz) +
  xlab("Year")
p2 + geom_line(aes(x=time, y=tot_pop))

# by class
agg_adjuv <- aggregate(run.sim.MSIRN$count, by=list(run.sim.MSIRN$time, run.sim.MSIRN$class), FUN=sum)
colnames(agg_adjuv) = c("time", "class", "count")

p3 <- ggplot(data=agg_adjuv) + geom_line(aes(x=time, y=count, color=class)) + 
  theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
  xlab("Year")

p3


#ggsave(file = "MSIRN.png",
 #      plot=p1,
  #     units="mm",  
   #    width=30, 
    #   height=60, 
     #  scale=3, 
      # dpi=300)

