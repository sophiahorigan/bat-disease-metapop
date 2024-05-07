####################################
## SIMULATE BAT DISEASE - PART 2 ##
###################################

# Created by Sophia Horigan (shorigan@uchicago.edu)
# Started: 01-30-24
# Last updated: 02-06-24

# Project Overview: This project seeks to use existing demographic data and the candidate transmission models
# to explore persistence threshholds for hypothetical pathogens under data-based metapopulation structure of Pteropus
# rufus in Madagascar. 

## PART 2 : ADD STOCHASTICITY
# for all pathogen and demographic params, add optional stochasticity
# this means instead of using a fixed value, I take an input value and use it as the
# mean in a distribution and then draw a random parameter value from
# that distribution

# Accomplishments
# 1-30-24 Added statements where stochasticity should go based on two new arguments: stoch_disease, stoch_demo
# 2-05-24 Started adding stochasticity - functions take in mean and sd, draw univariate random normal for each param
# 2-06-24 Finished adding stochasticity - modified code to take parameter array. in future, can populate parameter array in a loop
# 2-06-24 Added simulation loop

# To Do
# create plotting for stochastic simulations - mean and sd

rm(list=ls())

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


##########################
## Helper Functions
##########################
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

####################
## Matrix Functions
###################
# population matrix
build.pop.mat = function(surv, surv_juv, s, fecundity){ 
  pop.mat = matrix(0,  nrow=(s-1), ncol = (s-1))
  diag(pop.mat) = surv
  diag(pop.mat)[1] = surv_juv
  col_s = c(rep(0, s-2), surv)
  pop.mat = cbind(pop.mat, col_s)
  row1 = c(0, rep((fecundity*surv), (s-1))) #bats reproduce for the first time at the end of the second year of life. good for E. dup and P. ruf
  pop.mat = rbind(row1,pop.mat)
  return(pop.mat)
  
}#for stable age distribution

############
## MSIR
##########
# transmission matrix 
buildTMat_age <- function(model, stoch_disease, c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, param.array, add.inf.mort, age.rate){
  
  s <- nage <- length(age.classes)
  
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  
  # SET DISEASE PARAMETERS  
  if(stoch_disease == TRUE){
    beta <- rnorm(nage, param.array$beta_mean[model], param.array$beta_sd[model])
    sigma <- rnorm(nage, param.array$sigma_mean[model], param.array$sigma_sd[model])
    wane <- rnorm(nage, param.array$wane_mean[model], param.array$wane_sd[model])
    mu.sick <- rnorm(nage, param.array$mu.sick_mean[model], param.array$mu.sick_sd[model])
    recov <- rnorm(1, param.array$recov_mean[model], param.array$recov_sd[model])
  }
  else if(stoch_disease == FALSE){
    beta <- rep(param.array$beta_mean[model], nage) 
    sigma <- rep(param.array$sigma_mean[model], nage)
    mu.sick <- rep(param.array$mu.sick_mean[model], nage)
    wane <- rep(param.array$wane_mean[model], nage)
    recov <- param.array$recov_mean[model]
  }
  waning.maternal = wane
  age.rate <- rep(age.rate, nage)
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I comes #3
  foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi))) #vector

  mat1 <- matrix(0,4,4) # MSIR
  
  Tmat <- matrix(0,4*nage,4*nage) #MSIR for every age cohort
  for (j in 1:nage) {
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]#-psi
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    #mat1[4,2] <- psi
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[2,4] <- sigma[j]
    
    surv <- rep(1-mort.vect[j],4);
    
    if (add.inf.mort==TRUE){
      surv[3] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
    } 

    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*4+1):(j*4+4),((j-1)*4+1):(j*4)] <- mat1*surv*age.rate[j]
      
      #here, we ammend the epidemic transition matrix accordingly
      mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up
      Tmat[((j-1)*4+1):(j*4),((j-1)*4+1):(j*4)] <- mat2*surv*(1-age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*4+1):(j*4),((j-1)*4+1):(j*4)] <- mat1*surv
    }
  }
  return(Tmat)
}

# fecundity matrix
buildFMatrix <- function(age.classes, fecundity, surv.biwk, biwk){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  if(biwk==1){ #peak births
    new.fec = fecundity*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = fecundity*surv.biwk*.2
  }else if (biwk==3|biwk==25){
    new.fec = fecundity*surv.biwk*.1
  }else if (biwk==4|biwk ==24){
    new.fec = fecundity*surv.biwk*.05
  }else{
    new.fec = 0
  }
  s <- nage <- length(age.classes)
  
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  #make matrix the same size as the transition
  Fmat <- matrix(0,4*nage,4*nage)
  
  for (j in 1:nage) {
    Fmat[1,((j-1)*4+1):(j*4)] <- c(0,0,fert.biwk[j],fert.biwk[j]) #the mat immune
    Fmat[2,((j-1)*4+1):(j*4)]<- c(fert.biwk[j],fert.biwk[j],0,0) #the susceptible (mom didn't get sick)
  }
  
  return(Fmat)
}

# simulate model
sim.met.MSIR.age <- function(model, burnin, sim_pop, yrs, ntyr, age.brk, s, param.array, add.inf.mort, printpop, stoch_disease, stoch_demo){
  #################################
  ## Generate stable age structure
  ################################
  # num columns
  c=4 # M, S, I, R
  
  # generate time seq
  times <- seq(0, yrs, by = 1/ntyr) 
  
  # SET DEMOGRAPHIC PARAMETERS
  if(stoch_demo == TRUE){
    mort_ad <- rnorm(1, param.array$mort_ad_mean[model], param.array$mort_ad_sd[model])
    mort_juv <- rnorm(1, param.array$mort_juv_mean[model], param.array$mort_juv_sd[model])
    fecundity <- rnorm(1, param.array$fecundity_mean[model], param.array$fecundity_sd[model])
  }
  else if(stoch_demo == FALSE){
    mort_ad <- param.array$mort_ad_mean[model]
    mort_juv <- param.array$mort_juv_mean[model]
    fecundity <- param.array$fecundity_mean[model]
  }
  
  # Take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort_ad), surv_juv=(1-mort_juv), s=(s), fecundity = fecundity)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  # plot(stab.struct, xlab="Age", ylab="Proportion")
  
  ####################################################################
  ## Use stable age structure to generate bat population for time seq
  ####################################################################
  
  # check lambda
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  
  # gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  ######################################################
  ## Initiate infection and run dynamics to equilibrium
  ######################################################
  
  # introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 
  M_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  N_tot = cbind(M_init, S_init, I_init, R_init)
  
  N_pop = vec(N_tot) # by disease status
  M_pop = vec(t(N_tot)) # by age class
  
  # make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  ##################################
  ## Update stable age structure?? why?
  ##################################
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion")
  
  ##################################
  ## Iterate through timesteps
  ##################################
  
  for (i in 1:(length(times)-1)){
    # build matrices anew each time because foi depends on # infected
    # calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    
    # generate SIR matrix
    Tmat <- buildTMat_age(model = model, stoch_disease = stoch_disease, param.array = param.array, c = c, Npop = N_pop_ts[,i], age.classes = 1:s, age.brk = age.brk, surv.biwk = (1-mort_ad)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), add.inf.mort = add.inf.mort, age.rate = 1/ntyr)
    
    # generate fecundity matrix
    Fmat <- buildFMatrix(age.classes = 1:s, fecundity = fecundity, surv.biwk = (1-mort_ad)^(1/ntyr), biwk = biwk1)
    
    # make transition mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
  }
  
  ##############################
  ## Update stable structure
  #############################
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion")  
  
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  ######################
  ## Sum age classes
  #####################
  # split matrix into age classes
  N_split_pop = mat_split(N_pop_ts, r=c, c=ncol(N_pop_ts))
  
  # transform array into list and populate
  pop.list = list()
  for (i in 1:dim(N_split_pop)[3]){
    pop.list[[i]] = N_split_pop[,,i]
  }
  
  # then take column sums of each 
  pop.sum = lapply(X=pop.list, FUN = colSums)
  
  ##################################
  ## Sum disease state variables
  #################################
  # split matrix into state variables
  N_split_disease = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  # transform array into list and populate
  disease.list = list()
  for (i in 1:dim(N_split_disease)[3]){
    disease.list[[i]] = N_split_disease[,,i]
  }
  
  # then take column sums of each 
  disease.sum = lapply(X=disease.list, FUN = colSums)
  
  ################################
  ## Clean population dataframe
  ###############################
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot.pop = cbind.data.frame(times, pop.sum)
  names(dat.tot.pop) = c("time", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20")
  dat.tot.pop$tot_pop = rowSums((dat.tot.pop[,2:ncol(dat.tot.pop)]))
  
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
  
  dat.out.pop = subset(dat.out.pop, time >= burnin)

  ################################
  ## Clean disease dataframe
  ###############################
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot.disease = cbind.data.frame(times, disease.sum)
  names(dat.tot.disease) = c("time", "M", "S", "I", "R")
  dat.tot.disease$tot_pop = rowSums((dat.tot.disease[,2:ncol(dat.tot.disease)]))
  #print(dat.tot$tot_pop)
  
  dat.long.disease <- melt(dat.tot.disease, measure.vars = c("M", "S", "I", "R", "tot_pop"), variable.name = "state", value.name = "count")
  dat.sub.disease = subset(dat.long.disease, state!="tot_pop")
  dat.pop.disease = subset(dat.long.disease, state=="tot_pop")
  
  dat.pop.disease <- dplyr::select(dat.pop.disease, -(state))
  names(dat.pop.disease)[names(dat.pop.disease)=="count"] <- "tot_pop"
  
  dat.out.disease <- merge(dat.sub.disease, dat.pop.disease, by ="time", all.x = T, sort = F)
  #head(dat.out.disease)
  dat.out.disease$proportion <- dat.out.disease$count/dat.out.disease$tot_pop
  
  dat.out.disease$state = factor(dat.out.disease$state, levels=c("M", "S", "I", "R"))
  
  dat.out.disease = subset(dat.out.disease, time >= burnin)

  # CAN ONLY RETURN ONE ARRAY IN R SO DUMB
  # NEED TO FIGURE OUT HOW TO COMBINE THEM BAH
  if(printpop == TRUE){
    return(dat.out.pop)
  }
  if(printpop == FALSE){
    return(dat.out.disease)
  }
}

##########################
## MSIRN / MSIRNI
##########################
# age structure
get.age.struct.M = function(pop, c){ # function that collapses population vector in disease state form to age structure
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

# transmission matrix
buildTMat_MSIRN_age <- function(model, c, stoch_disease, Npop, age.classes, surv.biwk, age.brk, param.array, surv.juv.biwk, add.inf.mort, age.rate){
  
  ##########################
  ## Set up Parameters
  ##########################
  # num age classes
  s <- nage <- length(age.classes) 
  
  # mortality
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  
  # SET DISEASE PARAMETERS  
  if(stoch_disease == TRUE){
    beta <- rnorm(nage, param.array$beta_mean[model], param.array$beta_sd[model])
    sigma <- rnorm(nage, param.array$sigma_mean[model], param.array$sigma_sd[model])
    wane <- rnorm(nage, param.array$wane_mean[model], param.array$wane_sd[model])
    mu.sick <- rnorm(nage, param.array$mu.sick_mean[model], param.array$mu.sick_sd[model])
    recov <- rnorm(1, param.array$recov_mean[model], param.array$recov_sd[model])
    rho <- rnorm(1, param.array$rho_mean[model], param.array$rho_sd)
  }
  else if(stoch_disease == FALSE){
    beta <- rep(param.array$beta_mean[model], nage) 
    sigma <- rep(param.array$sigma_mean[model], nage)
    mu.sick <- rep(param.array$mu.sick_mean[model], nage)
    wane <- rep(param.array$wane_mean[model], nage)
    recov <- param.array$recov_mean[model]
    rho <- param.array$rho_mean[model]
  }
  waning.maternal = wane
  age.rate <- rep(age.rate, nage)
  
  # density dependent transmission
  # get population size
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  
  # calculate infecteds
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] # I is in position 3
  
  # vector of force of infection 
  foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi)))
  
  ##############################
  ## Populate Infection Matrix
  ##############################
  mat1 <- matrix(0,5,5) # MSIRN
  
  Tmat <- matrix(0,5*nage,5*nage) # MSIRN for every age cohort
  
  for (j in 1:nage) {     # fill in epi matrix for each age class. this gives probability of transmission per biweek for the model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[3,5] <- rho
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j]
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1-rho
    
    # survival 
    surv <- rep(1-mort.vect[j],5); # equal across all age groups
    
    # optional : add infection-induced mortality
    if (add.inf.mort==TRUE){
      surv[3] <- pmax(1-(mu.sick[j]*mort.vect[j]),0) # 3 is infected position
    } 
    
    ##################################
    ## Update State Variable Matrix
    #################################
    if (j!=nage) {       # if you are not in the last age class, you go into next age class
      
      # infection matrix, survival, and aging
      Tmat[(j*5+1):(j*5+5),((j-1)*5+1):(j*5)] <- mat1*surv*age.rate[j]
      
      # infection matrix, survival, NO aging
      Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(1-age.rate[j]) # not sure how this is decided
    } 
    else {
      # stay in this age class if you are at the peak age
      Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv
    }
  }
  return(Tmat)
}

# fecundity matrix
buildFMatrix_MSIRN <- function(age.classes, fecundity, surv.biwk, biwk, N_stat){ 	
  
  ################################
  ## Set timing of birth rates
  ###############################
  if(biwk==1){ 
    new.fec = fecundity*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = fecundity*surv.biwk*.25
  }else if (biwk==3|biwk==25){
    new.fec = fecundity*surv.biwk*.1
    #}else if (biwk==4|biwk ==24){
    #new.fec = fecundity*surv.biwk*.05
  }else{
    new.fec = 0
  }
  
  ####################
  ## Set up matrix
  ###################
  s <- nage <- length(age.classes)
  
  # compile into fecundity for each age class
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  Fmat <- matrix(0,5*nage,5*nage)
  
  ############################
  ## Update fertility matrix
  ###########################
  
  for (j in 1:nage) { # no fertility in first age class
    # try it assuming N moms produce maternally immune pups
    Fmat[1,((j-1)*5+1):(j*5)] <- c(0,0,fert.biwk[j],fert.biwk[j], fert.biwk[j]) # the mat immune (so mom = I and R but not N)
    Fmat[2,((j-1)*5+1):(j*5)]<-  c(fert.biwk[j],fert.biwk[j],0,0, 0) # the susceptible (mom didn't get sick, so mom= N and S)
  }
  return(Fmat)
}

# simulate model
sim.met.MSIRN.age <- function(model, burnin, sim_pop, yrs, ntyr, s, age.brk, param.array, add.inf.mort, N_stat, printpop, stoch_demo, stoch_disease){
  
  #################################
  ## Generate stable age structure
  ################################
  # num columns - in this case M, S, I, R, N
  c=5
  
  # generate time seq
  times <- seq(0, yrs, by = 1/ntyr) 
  
  # SET DEMOGRAPHIC PARAMETERS
  if(stoch_demo == TRUE){
    mort_ad <- rnorm(1, param.array$mort_ad_mean[model], param.array$mort_ad_sd[model])
    mort_juv <- rnorm(1, param.array$mort_juv_mean[model], param.array$mort_juv_sd[model])
    fecundity <- rnorm(1, param.array$fecundity_mean[model], param.array$fecundity_sd[model])
  }
  else if(stoch_demo == FALSE){
    mort_ad <- param.array$mort_ad_mean[model]
    mort_juv <- param.array$mort_juv_mean[model]
    fecundity <- param.array$fecundity_mean[model]
  }
  
  # Take juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort_ad), surv_juv=(1-mort_juv), s=(s), fecundity = fecundity)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  ####################################################################
  ## Use stable age structure to generate bat population for time seq
  ####################################################################
  # check lambda
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  #print(paste0("lambda = ", lambda))
  
  # gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  ######################################################
  ## Initiate infection and run dynamics to equilibrium
  ######################################################
  # introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 # does this initial value change dynamics much?
  M_init = rep(0, s)
  N_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init)
  
  N_pop = vec(N_tot) # by disease status
  M_pop = vec(t(N_tot)) # by age class
  
  # make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  ##################################
  ## Update stable age structure?? why?
  ##################################
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  ##################################
  ## Iterate through timesteps
  ##################################
  for (i in 1:(length(times)-1)){
    # build matrices anew each time because foi depends on # infected
    # calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times) # what is a biweek? different than selecting 26 for ntyr?
    
    # generate SIR matrix
    Tmat <- buildTMat_MSIRN_age(model = model, stoch_disease = stoch_disease, param.array = param.array, c = c, Npop = N_pop_ts[,i], age.classes = 1:s, age.brk = age.brk, surv.biwk = (1-mort_ad)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), add.inf.mort = add.inf.mort, age.rate = 1/ntyr)
    
    # generate fecundity matrix
    Fmat <- buildFMatrix_MSIRN(age.classes = 1:s, fecundity = fecundity, surv.biwk = (1-mort_ad)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    # make trans mat
    transMat <- Tmat + Fmat 
    
    # move forward in time
    nt1 <- (transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
  }
  
  ##############################
  ## Update stable structure
  #############################
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  ######################
  ## Sum age classes
  #####################
  # split matrix into age classes
  N_split_pop = mat_split(N_pop_ts, r=5, c=ncol(N_pop_ts))
  
  # transform matrix into list and populate
  pop.list = list()
  for (i in 1:dim(N_split_pop)[3]){
    pop.list[[i]] = N_split_pop[,,i]
  }
  
  # then take column sums of each and plot the class totals over time
  pop.sum = lapply(X=pop.list, FUN=colSums)
  
  ##################################
  ## Sum disease state variables
  #################################
  # split matrix into state variables
  N_split_disease = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  disease.list = list()
  for (i in 1:dim(N_split_disease)[3]){
    disease.list[[i]] = N_split_disease[,,i]
  }
  
  # then take column sums of each and plot the state variable totals over time
  disease.sum = lapply(X=disease.list, FUN=colSums)
  
  ################################
  ## Clean population dataframe
  ###############################
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot.pop = cbind.data.frame(times, pop.sum)
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
  
  dat.out.pop = subset(dat.out.pop, time >= burnin)

  ################################
  ## Clean disease dataframe
  ###############################
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot.disease = cbind.data.frame(times, disease.sum)
  names(dat.tot.disease) = c("time", "M", "S", "I", "R", "N")
  dat.tot.disease$tot_pop = rowSums((dat.tot.disease[,2:ncol(dat.tot.disease)]))
  #print(dat.tot$tot_pop)
  
  dat.long.disease <- melt(dat.tot.disease, measure.vars = c("M", "S", "I", "R", "N", "tot_pop"), variable.name = "state", value.name = "count")
  dat.sub.disease = subset(dat.long.disease, state!="tot_pop")
  dat.pop.disease = subset(dat.long.disease, state=="tot_pop")
  
  dat.pop.disease <- dplyr::select(dat.pop.disease, -(state))
  names(dat.pop.disease)[names(dat.pop.disease)=="count"] <- "tot_pop"
  
  dat.out.disease <- merge(dat.sub.disease, dat.pop.disease, by ="time", all.x = T, sort = F)
  #head(dat.out.disease)
  dat.out.disease$proportion <- dat.out.disease$count/dat.out.disease$tot_pop
  
  dat.out.disease$state = factor(dat.out.disease$state, levels=c("M", "S", "I", "R", "N"))
  
  dat.out.disease = subset(dat.out.disease, time >= burnin)

  # CAN ONLY RETURN ONE ARRAY IN R SO DUMB
  # NEED TO FIGURE OUT HOW TO COMBINE THEM BAH
  if(printpop == TRUE){
    return(dat.out.pop)
  }
  if(printpop == FALSE){
    return(dat.out.disease)
  }
}

##################
## RUN ONE SIMULATION
#################

run.sim.MSIRS <- sim.met.MSIR.age(model = 1,
                                  burnin = 0, 
                                  sim_pop = 10000, 
                                  yrs = 1, 
                                  ntyr = 26, 
                                  s = 20, 
                                  age.brk = 20, 
                                  param.array = param.array,
                                  add.inf.mort = FALSE,
                                  printpop = FALSE,
                                  stoch_demo = TRUE,
                                  stoch_disease = TRUE)

run.sim.MSIRN <- sim.met.MSIRN.age(model = 2,
                                   burnin=0,                  # years to ditch at start of sim
                                   sim_pop=1000,               # population size
                                   yrs=1,                     # years for simulation
                                   ntyr=26,                    # model timestep - greatly changes dynamics!
                                   s=20,                       # number of age classes
                                   age.brk = 20,               # how many distinct age classes there are
                                   param.array = param.array,
                                   add.inf.mort = FALSE,       # increased mortality rate for infected juveniles also?
                                   N_stat = "matAB",           # N class mothers produce M class young
                                   printpop = FALSE,        # print age class dynamics instead of infection dynamics
                                   stoch_demo = TRUE,         # stochasticity for demographic parameters
                                   stoch_disease = TRUE)      # stochasticity for disease parameters


run.sim.MSIRNI <- sim.met.MSIRN.age(model = 3,
                                    burnin = 0,
                                    sim_pop = 10000,
                                    yrs = 1,
                                    ntyr = 26,
                                    s = 20,
                                    age.brk = 20,
                                    param.array = param.array,
                                    add.inf.mort = FALSE, 
                                    N_stat = "matAB",
                                    printpop = FALSE,
                                    stoch_demo = TRUE,
                                    stoch_disease = TRUE)

beep()

###################
## PARAM ARRAY
##################

setwd("/Users/sophiahorigan/Documents/GitHub/Bat-disease-metapop/bat-disease-metapop/Models/Working")
param.array <- read.csv("params_brook2019.csv")

#######################
### PLOTTING INFECTION
#######################
## Seperate
## MSIRN
if (printpop == FALSE){
  
  colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N' = "navy")
  
  p1_msirn <- ggplot(data=run.sim.MSIRN) + geom_line(aes(x=time, y=proportion, color=state)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    scale_color_manual(values=colz) +
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRN$time), 1)) + ggtitle("MSIRN") +
    xlab("Year")
  p1_msirn 
}

## MSIRNI
if (printpop == FALSE){
  
  colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N' = "navy")
  
  p1_msirni <- ggplot(data=run.sim.MSIRNI) + geom_line(aes(x=time, y=proportion, color=state)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    scale_color_manual(values=colz) +
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRNI$time), 1)) + ggtitle("MSIRNI")
  xlab("Year")
  p1_msirni 
  
}

## MSIRS
if (printpop == FALSE){
  
  colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue")
  
  p1_msirs <- ggplot(data=run.sim.MSIRS) + geom_line(aes(x=time, y=proportion, color=state)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    scale_color_manual(values=colz) +
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRS$time), 1)) + ggtitle("MSIR") +
    xlab("Year")
  p1_msirs 
  
}

## Together



#######################
## PLOTTING POPULATION
######################
# MSIRS
if (printpop == TRUE){
  # by age
  p2_msirs <- ggplot(data=run.sim.MSIRS) + geom_line(aes(x=time, y=count, color=age)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    #scale_color_manual(values=colz) +
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRS$time), 1)) + ggtitle("MSIRS") +
    xlab("Year")
  p2_msirs + geom_line(aes(x=time, y=tot_pop))
  p2_msirs
}
if (printpop == TRUE){
  # by class
  agg_adjuv <- aggregate(run.sim.MSIRS$count, by=list(run.sim.MSIRS$time, run.sim.MSIRS$class), FUN=sum)
  colnames(agg_adjuv) = c("time", "class", "count")
  
  p3_msirs <- ggplot(data=agg_adjuv) + geom_line(aes(x=time, y=count, color=class)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRS$time), 1)) + ggtitle("MSIRS") + 
    xlab("Year")
  
  p3_msirs
}


# MSIRN
if (printpop == TRUE){
  # by age
  p2_msirn <- ggplot(data=run.sim.MSIRN) + geom_line(aes(x=time, y=count, color=age)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    #scale_color_manual(values=colz) +
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRN$time), 1)) + ggtitle("MSIRN") +
    xlab("Year")
  p2_msirn + geom_line(aes(x=time, y=tot_pop))
  p2_msirn
  
  # by class
  agg_adjuv <- aggregate(run.sim.MSIRN$count, by=list(run.sim.MSIRN$time, run.sim.MSIRN$class), FUN=sum)
  colnames(agg_adjuv) = c("time", "class", "count")
  
  p3_msirn <- ggplot(data=agg_adjuv) + geom_line(aes(x=time, y=count, color=class)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRN$time), 1)) + ggtitle("MSIRN") + 
    xlab("Year")
  
  p3_msirn
}

# MSIRNI
if (printpop == TRUE){
  # by age
  p2_msirni <- ggplot(data=run.sim.MSIRNI) + geom_line(aes(x=time, y=count, color=age)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    #scale_color_manual(values=colz) +
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRNI$time), 1)) + ggtitle("MSIRNI") +
    xlab("Year")
  p2_msirni + geom_line(aes(x=time, y=tot_pop))
  p2_msirni
  
  # by class
  agg_adjuv <- aggregate(run.sim.MSIRNI$count, by=list(run.sim.MSIRNI$time, run.sim.MSIRNI$class), FUN=sum)
  colnames(agg_adjuv) = c("time", "class", "count")
  
  p3_msirni <- ggplot(data=agg_adjuv) + geom_line(aes(x=time, y=count, color=class)) + 
    theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    scale_x_continuous(breaks=seq(0, max(run.sim.MSIRNI$time), 1)) + ggtitle("MSIRNI") +
    xlab("Year")
  
  p3_msirni
}


#####################
## SIMULATION WRAPPER
####################
## WIP

simulations <- function(model, numsims, param.array, burnin, sim_pop, yrs){
  
  # generate param array
  param.array = param.array # could have it change year to year 
  
  # generate output dataframe
  datalist = vector("list", length = numsims)

  # run sims, saving output
  if(model == 1){
    for(i in 1:numsims){
      output <- sim.met.MSIR.age(model = model, 
                       burnin = burnin,
                       sim_pop = sim_pop,
                       yrs = yrs,
                       ntyr = 26,
                       s = 20,
                       age.brk = 20,
                       param.array = param.array,
                       add.inf.mort = FALSE, 
                       printpop = FALSE,
                       stoch_demo = TRUE,
                       stoch_disease = TRUE)
      output$sim <- i
      datalist[[i]] <- output
    }
  }
  else if(model == 2 || model == 3){
    for(i in 1:numsims){
      output <- sim.met.MSIRN.age(model = model,
                        burnin = burnin,
                        sim_pop = sim_pop,
                        yrs = yrs,
                        ntyr = 26,
                        s = 20,
                        age.brk = 20,
                        param.array = param.array,
                        add.inf.mort = FALSE, 
                        N_stat = "matAB",
                        printpop = FALSE,
                        stoch_demo = TRUE,
                        stoch_disease = TRUE)
      output$sim <- i
      datalist[[i]] <- output
    }
  }
  
  sims <<- do.call(rbind, datalist)
}
test <- simulations(model = 1, numsims = 10, param.array = param.array, burnin = 0, sim_pop = 10000, yrs = 100)
beep()


colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue")
colz2 = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N' = "navy")

p1_msirni <- ggplot(data=test) + geom_point(aes(x=time, y=proportion, color=state)) + 
  theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
  scale_color_manual(values=colz2) +
  scale_x_continuous(breaks=seq(0, max(test$time), 1)) + ggtitle("MSIRNI") +
  xlab("Year")
p1_msirni 

ggsave(file = "MSIRNI_stochexample.png",
      plot=p1_msirni,
     units="mm",  
    width=70, 
   height=60, 
  scale=3, 
 dpi=300)

