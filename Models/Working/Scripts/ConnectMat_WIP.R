##############################
## CONNECTIVITY - PART 4 ##
############################

# Created by Sophia Horigan (shorigan@uchicago.edu)
# Started: 01-30-24
# Last updated: 04-20-24

# Project Overview: This project seeks to use existing demographic data and the candidate transmission models
# to explore persistence threshholds for hypothetical pathogens under data-based metapopulation structure of Pteropus
# rufus in Madagascar. 

## PART 4 : ADD IN CONNECTIVITY STRUCTURE
# Two mechanisms of connectivity
# 1. dispersal i.e. bat moves to different subpop
# 2. interminging i.e. interaction (read: potential transmission) at communal foraging site
# ASSUMMPTIONS:
  # Probability of dispersal and interaction is same regardless of disease state
  # Only adults disperse

# Accomplishments
# 02-12-24 Set up dispersal.df - movement of individuals
# 02-13-24 Realized I was making it harder than it needed to be - can keep connect.df the same if assuming equal prob of dispersal 
# 04-17-24 Improved saving and plotting features
# 04-20-24 Distinguished dispersal vs intermingling in code language
# 04-20-24 Started to implement intermingling
  # need to make all matrices multidimensional by num_subpops...
  # maybe should fix the number of subpops as a grid. if I do that, should make all variables an array that reflects position in the grid
# 04-20-24 Brainstormed grid structure - next evolution of model

# TO DO
# figure out how to list all init pop sizes on figure and output file - or print an associated params/inputs file for each figure??
# because the results depend so much on initial conditions, will need to simulate under a range of initpop, pop growth, and disease
  # heatmap of persistence, then mark where our data suggests the pop is 
# Develop BuildDMat and BuildIMat 
# Make Npop and other subpop tracking states into arrays instead of lists to accommodate grid locations (EXPAND TO GRID STRUCTURE)
  # input: grid dimensions
  # input: grid size (? or just say one patch holds one population. but need general size to overlay on map of Mada)
  # input: percentage of occupied patches starting
  # input: patch structure (randomly occupied, clustered, evenly spaced, etc)
# figure out how to combine buildTmat and BuildTmatMSIRN - same for fecundity matrix
# does R have something like a Struct in C that makes it easier to pass everything through functions?

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

setwd("/Users/sophiahorigan/Documents/GitHub/Bat-disease-metapop/bat-disease-metapop/Models/Working/")


##########################
## Helper Functions
##########################
TransformVect <- function(vec, s, c){ # vector, rows, columns
  vec2 <- t(commutation.matrix(r=s,c=c))%*%vec
  return(vec2)
}

MatSplit <- function(M, r, c){
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
} # for matrix slicing. give it the number of rows and columns you want in the resulting matrices

FindBiweek = function(t, times){
  
  biwks <- sort(unique(round(RevTrunc(times),4)))
  this.wk = round(RevTrunc(times[t]), 4)
  this.biwk <- which(biwks==this.wk)
  
  
  return(this.biwk)
  
  
}

RevTrunc = function(x){
  newx = x - floor(x)
  return(newx)
}

####################
## Matrix Functions
###################
# population matrix
BuildPopMat = function(surv, surv_juv, s, fecundity){ 
  pop.mat = matrix(0,  nrow=(s-1), ncol = (s-1))
  diag(pop.mat) = surv
  diag(pop.mat)[1] = surv_juv
  col_s = c(rep(0, s-2), surv)
  pop.mat = cbind(pop.mat, col_s)
  row1 = c(0, rep((fecundity*surv), (s-1))) #bats reproduce for the first time at the end of the second year of life. good for E. dup and P. ruf
  pop.mat = rbind(row1,pop.mat)
  return(pop.mat)
  
}#for stable age distribution

# transmission matrix # generate for all subpops at the same time
BuildTMatMSIR <- function(model, stoch_disease, c, Npop_epi, Ipop, sum_intermingle, age.classes, surv.biwk, age.brk, surv.juv.biwk, param.array, add.inf.mort, age.rate){
  
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
  
  # calculate force of infection
  foi = 1-exp(-beta*((sum(Ipop)/sum(Npop_epi)) + sum_intermingle))
  
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
BuildFMatMSIR <- function(age.classes, fecundity, surv.biwk, biwk){ 	#one for every age class
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

# age structure
GetAgeStructM = function(pop, c){ # function that collapses population vector in disease state form to age structure
  age.mat <- MatSplit(M=matrix(pop, ncol=1), r=c, c=1)
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
BuildTMatMSIRN <- function(model, c, stoch_disease, Npop_epi, Ipop, sum_intermingle, age.classes, surv.biwk, age.brk, param.array, surv.juv.biwk, add.inf.mort, age.rate){
  
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
  
  # calculate force of infection 
  foi = 1-exp(-beta*((sum(Ipop)/sum(Npop_epi)) + sum_intermingle))
  
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
BuildFMatMSIRN <- function(age.classes, fecundity, surv.biwk, biwk){ 	
  
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

# dispersal matrix # for larger simulations only
BuildDMat <- function(nage, nclass, dispersal.df, num_subpop, biweek){
  # maybe just fix connect.df in here
  ####################
  ## Set up matrix
  ###################
  if(nrow(connect.df) != num_subpop){
    print("Wrong dimension for dispersal matrix! Fix.")
  }
  
  connect.df
  # can modify base values here, or generate new ones
  # add stochasticity in values here
  # stochasticity allows for open population -- can 'gain' or 'lose' indls by having dispersal prob not total 1
  
  ###########################################
  ## EXPAND TO ALL AGE AND INFECTION CLASSES
  ##########################################

  #that blank matrix:
  n.dim2 <- nclass*nage*num_subpop
  
  connect.matrix2 <-matrix(NA,nrow=n.dim2, ncol=n.dim2) 
  
  #take each element of the matrix and make it the diagonal of its own matrix that is 
  #the same length as the number of age classes
  tot.mat.list <- list()
  for (i in 1:num_subpop){
    #take site "i" and make a list of dispersal probabilities to the other sites
    tmp <- connect.df[i,]
    
    #expand each probability by the number of age and infection classes
    matrix.list <- c()
    for (j in 1:length(tmp)){
      matrix.list[[j]] <- matrix(0, nrow=nage*nclass, ncol=nage*nclass)
      diag(matrix.list[[j]]) <- tmp[j]
    }
    
    #bind them together end-to-end in a long matrix
    tmp2 <- matrix(unlist(matrix.list), nrow=nage*nclass, ncol=nage*nclass*num_subpop)
    
    #store them in a list and then bind them later by rows
    tot.mat.list[[i]] <- tmp2
  }
  
  #now write over the columns with infants such that they will not disperse at all
  inf.columns <- c(1:nclass)
  
  for (m in 1:num_subpop){
    for (n in 1:length(inf.columns)){
      for (k in 1:length(inf.columns)){
        if (n == k){ # stay in same population
          tot.mat.list[[m]][n,k] = 1
        }
        else{
          tot.mat.list[[m]][n,k] = 0
        }
      }
    }
  }
  
  connect.matrix2 <- do.call("rbind", tot.mat.list)
  
  connect.matrix2 <- as.matrix(connect.matrix2)
  
  #and return the named matrix
  return(connect.matrix2)
}

# intermingling matrix # for larger simulations only
BuildIMat <- function(nage, nclass, intermingle.df, num_subpop, biweek){
  # maybe just fix connect.df in here
  ####################
  ## Set up matrix
  ###################
  if(nrow(connect.df) != num_subpop){
    print("Wrong dimension for dispersal matrix! Fix.")
  }
  
  connect.df
  # can modify base values here, or generate new ones
  # add stochasticity in values here
  # stochasticity allows for open population -- can 'gain' or 'lose' indls by having dispersal prob not total 1
  
  ###########################################
  ## EXPAND TO ALL AGE AND INFECTION CLASSES
  ##########################################
  
  #that blank matrix:
  n.dim2 <- nclass*nage*num_subpop
  
  connect.matrix2 <-matrix(NA,nrow=n.dim2, ncol=n.dim2) 
  
  #take each element of the matrix and make it the diagonal of its own matrix that is 
  #the same length as the number of age classes
  tot.mat.list <- list()
  for (i in 1:num_subpop){
    #take site "i" and make a list of dispersal probabilities to the other sites
    tmp <- connect.df[i,]
    
    #expand each probability by the number of age and infection classes
    matrix.list <- c()
    for (j in 1:length(tmp)){
      matrix.list[[j]] <- matrix(0, nrow=nage*nclass, ncol=nage*nclass)
      diag(matrix.list[[j]]) <- tmp[j]
    }
    
    #bind them together end-to-end in a long matrix
    tmp2 <- matrix(unlist(matrix.list), nrow=nage*nclass, ncol=nage*nclass*num_subpop)
    
    #store them in a list and then bind them later by rows
    tot.mat.list[[i]] <- tmp2
  }
  
  #now write over the columns with infants such that they will not disperse at all
  inf.columns <- c(1:nclass)
  
  for (m in 1:num_subpop){
    for (n in 1:length(inf.columns)){
      for (k in 1:length(inf.columns)){
        if (n == k){ # stay in same population
          tot.mat.list[[m]][n,k] = 1
        }
        else{
          tot.mat.list[[m]][n,k] = 0
        }
      }
    }
  }
  
  connect.matrix2 <- do.call("rbind", tot.mat.list)
  
  connect.matrix2 <- as.matrix(connect.matrix2)
  
  #and return the named matrix
  return(connect.matrix2)
}

################
## LOAD PARAMS
###############

# param array
param.array <- read.csv("Input/params_brook2019.csv")
param.array

# probability they switch to a different roost (also from telem? from other data?)
dispersal.df <- as.matrix(read.csv("Input/dispersal_static_low.csv", header = FALSE))
dispersal.df <- as.matrix(read.csv("Input/dispersal_static_high.csv", header = FALSE))
dispersal.df <- as.matrix(read.csv("Input/dispersal_static_NONE.csv", header = FALSE))
dispersal.df <- as.matrix(read.csv("Input/dispersal_static_extralow.csv", header = FALSE))
dispersal.df

# probability they will encounter an individual (home range overlap)
intermingle.df <- as.data.frame(read.csv("Input/intermingle_static_low.csv", header = FALSE))
intermingle.df

##########################
## SIMULATE MODEL
##########################

SimOneModel <- function(model, num_subpop, burnin, sim_pop, yrs, ntyr, age.brk, s, param.array, add.inf.mort, printpop, stoch_disease, stoch_demo, intermingling, dispersal){
  #################################
  ## Generate stable age structure
  ################################
  
  if(model==1){
    nclass=4 # M, S, I, R
  }
  else if(model==2 || model==3){
    nclass=5 # M, S, I, R, N
  }
  
  # generate time seq
  times <- seq(0, yrs, by = 1/ntyr) 
  
  # define lists
  # TO DO: MAKE THESE INTO ARRAYS TO MATCH PATCH GRID LOCATION
  # size - grid dimensions
  mort_ad <- list()
  mort_juv <- list()
  fecundity <- list()
  N_pop_ts <- list()
  out.pop <- list()
  out.disease <- list()
  Npop_epi <- list()
  Ipop <- list()
  Tmat <- list()
  
  ## START SUBPOP LOOP
  for(i in 1:num_subpop){
    
    # SET DEMOGRAPHIC PARAMETERS
    if(stoch_demo == TRUE){
      mort_ad[i] <- as.numeric(rnorm(1, param.array$mort_ad_mean[model], param.array$mort_ad_sd[model]))
      mort_juv[i] <- as.numeric(rnorm(1, param.array$mort_juv_mean[model], param.array$mort_juv_sd[model]))
      fecundity[i] <- as.numeric(rnorm(1, param.array$fecundity_mean[model], param.array$fecundity_sd[model]))
    }
    else if(stoch_demo == FALSE){
      mort_ad[i] <- as.numeric(param.array$mort_ad_mean[model])
      mort_juv[i] <- as.numeric(param.array$mort_juv_mean[model])
      fecundity[i] <- as.numeric(param.array$fecundity_mean[model])
    }
    
    # convert list to numeric
    mort_ad <- as.numeric(mort_ad)
    mort_juv <- as.numeric(mort_juv)
    fecundity <- as.numeric(fecundity)
    
    if (length(sim_pop) != num_subpop){
      print("Did not pass an init pop size for each subpop!")
      stop()
    }
    
    # Take our juvenile and adult survival rates and use them to get the stable age structure using subpop specific demographic parameters
    mat1 = BuildPopMat(surv=(1-mort_ad[i]), surv_juv=(1-mort_juv[i]), s=(s), fecundity = fecundity[i])
    
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
    bat.mat = stab.struct*sim_pop[i]
    
    ######################################################
    ## Initiate infection and run dynamics to equilibrium
    ######################################################
    
    # introduce a few infecteds and run it out to equilibrium before you grab the data
    
    if(model==1){ # MSIRS
      R_init = rep(0, s)
      I_init = rep(0, s); I_init[3] = 5 
      M_init = rep(0, s)
      S_init = bat.mat - I_init 
      
      N_tot = cbind(M_init, S_init, I_init, R_init)
    }
    else if(model==2 || model==3){ # MSIRN / MSIRNI
      R_init = rep(0, s)
      I_init = rep(0, s); I_init[3] = 5 # does this initial value change dynamics much?
      M_init = rep(0, s)
      N_init = rep(0, s)
      S_init = bat.mat - I_init 
      
      N_tot = cbind(M_init, S_init, I_init, R_init, N_init)
    }
    
    N_pop = vec(N_tot) # by disease status
    M_pop = vec(t(N_tot)) # by age class
    
    # make list of matricies to store your population as you go that holds all subpops
    tmp <- matrix(NA, ncol = length(times), nrow(N_pop))
    tmp[,1] <- M_pop # by age class
    N_pop_ts[[i]] <- tmp # fill into large array
    
    ##################################
    ## Update stable age structure?? why?
    ##################################
    
    stab.struct <- GetAgeStructM(M_pop, c=nclass)
    stab.struct <- c(unlist(stab.struct))
    stab.struct <- stab.struct/(sum(stab.struct))
    #plot(stab.struct, xlab="Age", ylab="Proportion")
  }
  
  ##################################
  ## Iterate through timesteps
  ##################################
  
  for (j in 1:(length(times)-1)){
    
    # generate dispersal matrix
    
    # generate intermingling matrix
    
    ####################
    ## INTERMINGLING
    ###################
    # calculate num infected and total pop for each subpop
      for(i in 1:num_subpop){
        
        tmp = TransformVect(vec=N_pop_ts[[i]][,j], s=s, c=nclass)
        Npop_epi[[i]] <- tmp
        
        tmp = MatSplit(matrix(Npop_epi[[i]], ncol=1), r=s, c=1)[,,3] # I is in position 3
        Ipop[[i]] <- tmp
      }
    
    # calculate intermingling factor
      sum_intermingle <- list(rep(0, as.numeric(num_subpop)))
      
      for (i in 1:num_subpop){  
        for (j in 1:num_subpop){
          if (i != j){
            if (intermingling == TRUE){
              sum_intermingle[[i]] = sum_intermingle[[i]] + (sum(Ipop[[j]])/sum(Npop_epi[[j]]) * intermingling.df[i][j]) 
            }
            else if (intermingling == FALSE){
              sum_intermingle[[i]] = 0
            }
          }
        }
      }
    
    # compute matrices and iterate one timestep
    for(i in 1:num_subpop){

      biwk1 <- FindBiweek(t=j, times=times)
      
      if(model==1){ # MSIRS
        # generate SIR matrix
         tmp = BuildTMatMSIR(model = model, Npop_epi = Npop_epi[[i]], Ipop = Ipop[[i]], sum_intermingle = sum_intermingle[[i]], stoch_disease = stoch_disease, param.array = param.array, c = nclass, age.classes = 1:s, age.brk = age.brk, surv.biwk = (1-mort_ad[i])^(1/ntyr), surv.juv.biwk = (1-mort_juv[i])^(1/ntyr), add.inf.mort = add.inf.mort, age.rate = 1/ntyr)
         Tmat[[i]] <- tmp
        
        # generate fecundity matrix
        Fmat <- BuildFMatMSIR(age.classes = 1:s, fecundity = fecundity[i], surv.biwk = (1-mort_ad[i])^(1/ntyr), biwk = biwk1)

      }
      else if(model==2 || model==3){ # MSIRN, MSIRNI
        # generate SIR matrix
        tmp <- BuildTMatMSIRN(model = model, Npop_epi = Npop_epi[[i]], Ipop = Ipop[[i]], sum_intermingle = sum_intermingle[[i]], stoch_disease = stoch_disease, param.array = param.array, c = nclass, age.classes = 1:s, age.brk = age.brk, surv.biwk = (1-mort_ad[i])^(1/ntyr), surv.juv.biwk = (1-mort_juv[i])^(1/ntyr), add.inf.mort = add.inf.mort, age.rate = 1/ntyr)
        Tmat[[i]] <- tmp
        
        # generate fecundity matrix
        Fmat <- BuildFMatMSIRN(age.classes = 1:s, fecundity = fecundity[i], surv.biwk = (1-mort_ad[i])^(1/ntyr), biwk = biwk1)
        
      }
      
      # make transition matrix
      transMat <- Tmat[[i]] + Fmat
      
      #move forward in time
      N_pop_ts[[i]][,j+1] <- transMat %*% N_pop_ts[[i]][,j]
      #print(paste0("Pre dispersal subpop", i, "= ",  N_pop_ts[[i]][,j+1]))
      #print(i)
      #print(N_pop_ts[[i]][,j+1])
      #print(class(N_pop_ts[[i]][,j+1]))
    }
    
    # stitch current pops into matrix
    tmp <- list()
    for (i in 1:num_subpop){
      tmp[[i]] <- N_pop_ts[[i]][,j+1]
    }
    N_pop <- do.call(rbind, tmp)
    
    ##############################
    ## DISPERSAL
    #############################
    if (dispersal == TRUE){
      N_pop <- dispersal.df %*% N_pop
    }
    
    # split resulting matrix and fill in pop list
    for (i in 1:num_subpop){
      N_pop_ts[[i]][,j+1] <- N_pop[i,]
    }
  }
  
  for(i in 1:num_subpop){
    ##############################
    ## Update stable structure
    #############################
    stab.struct <- GetAgeStructM(pop=N_pop_ts[[i]][,ncol(N_pop_ts[[i]])], c=nclass)
    stab.struct <- c(unlist(stab.struct))
    stab.struct <- stab.struct/(sum(stab.struct))
    #plot(stab.struct, xlab="Age", ylab="Proportion")  
    
    tmp <- TransformVect(vec=N_pop_ts[[i]], s=s, c=nclass)
    
    ######################
    ## Sum age classes
    #####################
    # split matrix into age classes
    N_split_pop = MatSplit(tmp, r=nclass, c=ncol(tmp))
    
    # transform array into list and populate
    pop.list = list()
    for (j in 1:dim(N_split_pop)[3]){
      pop.list[[j]] = N_split_pop[,,j]
    }
    
    # then take column sums of each 
    pop.sum = lapply(X=pop.list, FUN = colSums)
    
    ##################################
    ## Sum disease state variables
    #################################
    # split matrix into state variables
    N_split_disease = MatSplit(tmp, r=s, c=ncol(tmp))
    
    # transform array into list and populate
    disease.list = list()
    for (j in 1:dim(N_split_disease)[3]){
      disease.list[[j]] = N_split_disease[,,j]
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
    for(j in 1:nrow(dat.out.pop)){
      if(dat.out.pop[j,"age"] == 1){dat.out.pop[j, "class"] = 'J'}
      if(dat.out.pop[j,"age"] != 1){dat.out.pop[j, "class"] = 'A'}
    }
    
    dat.out.pop$subpop = i
    
    out.pop[[i]] = subset(dat.out.pop, time >= burnin)
    
    ################################
    ## Clean disease dataframe
    ###############################
    #times=seq(0,yrs,by =1/ntyr)
    
    if(model==1){ # MSIRS
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
    }
    else if(model==2 || model==3){ # MSIRN, MSIRNI
      dat.tot.disease = cbind.data.frame(times, disease.sum)
      names(dat.tot.disease) = c("time", "M", "S", "I", "R","N")
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
    }
    
    dat.out.disease$subpop = i
    
    out.disease[[i]] = subset(dat.out.disease, time >= burnin)
  }
  
  final.disease <- do.call(rbind, out.disease)
  final.pop <- do.call(rbind, out.pop)
  
  # CAN ONLY RETURN ONE ARRAY IN R SO DUMB
  # NEED TO FIGURE OUT HOW TO COMBINE THEM BAH
  if(printpop == FALSE){
    return(final.disease)
  }
  if(printpop == TRUE){
    return(final.pop)
  }
}

SimAllModels <- function(yrs, numsims, do.save){
  for (i in 1:3){
    
    SimOneModel(model = i,
                numsims = numsims,
                num_subpop = 3,
                sim_pop = c(2000,110,200), 
                param.array = param.array, 
                burnin = 0, 
                yrs = yrs, 
                do.plot = TRUE)
    beep(i)
    
  }
  p4 <- plot_grid(p1_1, p2_1, p3_1)
  print(p4)
  
  p5 <- plot_grid(p1_2, p2_2, p3_2)
  print(p4)
  
  if(do.save==TRUE){
    ggsave(file = paste0("AllModels_NumSims", numsims, "_NumYrs", yrs, "_", Sys.time(), ".png"),
           plot=p4,
           units="mm",  
           width=70, 
           height=60, 
           scale=3, 
           dpi=300)
    
    ggsave(file = paste0("CountAllModels_NumSims", numsims, "_NumYrs", yrs, "_", Sys.time(), ".png"),
           plot=p5,
           units="mm",  
           width=70, 
           height=60, 
           scale=3, 
           dpi=300)
  }
}

#####################
## SIMULATION WRAPPER
####################

SimWrap <- function(model, sim_pop, numsims, num_subpop, param.array, burnin, yrs, do.plot, do.save){
  
  # generate param array
  param.array = param.array # could have it change year to year 
  
  # can loop through sim_pops
  
  # generate output dataframe
  datalist = vector("list", length = numsims)
  
  # run sims, saving output
  for(i in 1:numsims){
    output <- SimOneModel(model = model,
                      num_subpop = num_subpop,
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
                      stoch_disease = TRUE,
                      intermingling = FALSE,
                      dispersal = FALSE)
    output$model <- model
    output$sim <- i
    datalist[[i]] <- output
  }
  
  sims <<- do.call(rbind, datalist)
  
  if(do.plot==TRUE){
    # get mean and CI for stochastic realizations
    sim.summary <- sims %>%
      group_by(time, state, subpop) %>%
      dplyr::summarise(mean.proportion = mean(proportion, na.rm = TRUE),
                       sd.proportion = sd(proportion, na.rm = TRUE),
                       n.proportion = n(),
                       mean.count = mean(count, na.rm = TRUE),
                       sd.count = sd(count, na.rm = TRUE),
                       n.count = n(),
                       mean.totpop = mean(tot_pop, na.rm = TRUE),
                       sd.totpop = sd(tot_pop, na.rm = TRUE),
                       n.totpop = n()) %>%
      mutate(se.proportion = sd.proportion / sqrt(n.proportion),
             lower.ci.proportion = mean.proportion - qt(1 - (0.05 / 2), n.proportion - 1) * se.proportion,
             upper.ci.proportion = mean.proportion + qt(1 - (0.05 / 2), n.proportion - 1) * se.proportion,
             se.count = sd.count / sqrt(n.count),
             lower.ci.count = mean.count - qt(1 - (0.05 / 2), n.count - 1) * se.count,
             upper.ci.count = mean.count + qt(1 - (0.05 / 2), n.count - 1) * se.count,
             se.totpop = sd.totpop / sqrt(n.totpop),
             lower.ci.totpop = mean.totpop - qt(1 - (0.05 / 2), n.totpop - 1) * se.totpop,
             upper.ci.totpop = mean.totpop + qt(1 - (0.05 / 2), n.totpop - 1) * se.totpop)
    
    # and plot by model type
    if(model == 1){
      colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue")
      p1_1 <<- ggplot(data=sim.summary) + geom_line(aes(x=time, y=mean.proportion, color=state)) + facet_grid(~subpop) +
        theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
        geom_ribbon(aes(x=time, ymin = lower.ci.proportion, ymax = upper.ci.proportion, fill = state), alpha = 0.2) +
        scale_color_manual(values=colz) +
        scale_x_continuous(breaks=seq(0, max(sim.summary$time), 1)) + ggtitle(paste0("MSIRS. initpop = ", sim_pop, " numsims = ", numsims)) +
        xlab("Year")
      print(p1_1)
      
      p1_2 <<- ggplot(data=sim.summary) + geom_line(aes(x=time, y=mean.totpop)) + facet_grid(~subpop) +
        geom_line(aes(x=time, y=mean.count, color=state)) +
        theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
        geom_ribbon(aes(x=time, ymin = lower.ci.count, ymax = upper.ci.count, fill = state), alpha = 0.2) +
        geom_ribbon(aes(x=time, ymin = lower.ci.totpop, ymax = upper.ci.totpop), alpha = 0.2) +
        scale_color_manual(values=colz) +
        scale_x_continuous(breaks=seq(0, max(sim.summary$time), 1)) + ggtitle(paste0("MSIRS COUNT. initpop = ", sim_pop, " numsims = ", numsims)) +
        xlab("Year")
      print(p1_2)
      
      if(do.save==TRUE){
        ggsave(path = "Output", file = paste0("Frac_Model_", model, "_InitPop_", sim_pop[model], "_NumSims_", numsims, "_Yrs_", yrs, "_", Sys.time(), ".png"),
               plot=p1_1,
               units="mm",  
               width=70, 
               height=60, 
               scale=3, 
               dpi=300)
        ggsave(path = "Output", file = paste0("Count_Model_", model, "_InitPop_", sim_pop[model], "_NumSims_", numsims, "_Yrs_", yrs, "_", Sys.time(), ".png"),
               plot=p1_2,
               units="mm",  
               width=70, 
               height=60, 
               scale=3, 
               dpi=300)
      }
    }
    if(model == 2){
      colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N' = "navy")
      p2_1 <<- ggplot(data=sim.summary) + geom_line(aes(x=time, y=mean.proportion, color=state)) + facet_grid(~subpop) +
        theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
        geom_ribbon(aes(x=time, ymin = lower.ci.proportion, ymax = upper.ci.proportion, fill = state), alpha = 0.2) +
        scale_color_manual(values=colz) +
        scale_x_continuous(breaks=seq(0, max(sim.summary$time), 1)) + ggtitle(paste0("MSIRN. initpop = ", sim_pop, " numsims = ", numsims)) +
        xlab("Year")
      print(p2_1)
      
      p2_2 <<- ggplot(data=sim.summary) + geom_line(aes(x=time, y=mean.totpop)) + facet_grid(~subpop) +
        geom_line(aes(x=time, y=mean.count, color=state)) +
        theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
        geom_ribbon(aes(x=time, ymin = lower.ci.count, ymax = upper.ci.count, fill = state), alpha = 0.2) +
        geom_ribbon(aes(x=time, ymin = lower.ci.totpop, ymax = upper.ci.totpop), alpha = 0.2) +
        scale_color_manual(values=colz) +
        scale_x_continuous(breaks=seq(0, max(sim.summary$time), 1)) + ggtitle(paste0("MSIRN COUNT. initpop = ", sim_pop, " numsims = ", numsims)) +
        xlab("Year")
      print(p2_2)
      
      if(do.save==TRUE){
        ggsave(path = "Output", file = paste0("Frac_Model_", model, "_InitPop_", sim_pop[model], "_NumSims_", numsims, "_Yrs_", yrs, "_", Sys.time(), ".png"),
               plot=p2_1,
               units="mm",  
               width=70, 
               height=60, 
               scale=3, 
               dpi=300)
        ggsave(path = "Output", file = paste0("Count_Model_", model, "_InitPop_", sim_pop[model], "_NumSims_", numsims, "_Yrs_", yrs, "_", Sys.time(), ".png"),
               plot=p2_2,
               units="mm",  
               width=70, 
               height=60, 
               scale=3, 
               dpi=300)
      }
    }
    if( model == 3){
      colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N' = "navy")
      p3_1 <<- ggplot(data=sim.summary) + geom_line(aes(x=time, y=mean.proportion, color=state)) + facet_grid(~subpop) +
        theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
        geom_ribbon(aes(x=time, ymin = lower.ci.proportion, ymax = upper.ci.proportion, fill = state), alpha = 0.2) +
        scale_color_manual(values=colz) +
        scale_x_continuous(breaks=seq(0, max(sim.summary$time), 1)) + ggtitle(paste0("MSIRNI. initpop = ", sim_pop, " numsims = ", numsims)) +
        xlab("Year")
      print(p3_1)
      
      p3_2 <<- ggplot(data=sim.summary) + geom_line(aes(x=time, y=mean.totpop)) + facet_grid(~subpop) +
        geom_line(aes(x=time, y=mean.count, color=state)) +
        theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
        geom_ribbon(aes(x=time, ymin = lower.ci.count, ymax = upper.ci.count, fill = state), alpha = 0.2) +
        geom_ribbon(aes(x=time, ymin = lower.ci.totpop, ymax = upper.ci.totpop), alpha = 0.2) +
        scale_color_manual(values=colz) +
        scale_x_continuous(breaks=seq(0, max(sim.summary$time), 1)) + ggtitle(paste0("MSIRNI COUNT. initpop = ", sim_pop, " numsims = ", numsims)) +
        xlab("Year")
      print(p3_2)
      
      if(do.save==TRUE){
        ggsave(path = "Output", file = paste0("Frac_Model_", model, "_InitPop_", sim_pop[model], "_NumSims_", numsims, "_Yrs_", yrs, "_", Sys.time(), ".png"),
               plot=p3_1,
               units="mm",  
               width=70, 
               height=60, 
               scale=3, 
               dpi=300)
        
        ggsave(path = "Output", file = paste0("Count_Model_", model, "_InitPop_", sim_pop[model], "_NumSims_", numsims, "_Yrs_", yrs, "_", Sys.time(), ".png"),
               plot=p3_2,
               units="mm",  
               width=70, 
               height=60, 
               scale=3, 
               dpi=300)
      }
    }
  }
}

#-------------------------------------------------------------------------------------------------------------------------#

###################
## RUN ONE MODEL
##################

test <- SimWrap(model = 1,
                  numsims = 5,
                  num_subpop = 2,
                  sim_pop = c(1000, 1000), 
                  param.array = param.array, 
                  burnin = 0, 
                  yrs = 1, 
                  do.plot = TRUE,
                  do.save = FALSE)


beep()

##############################
## SUPER LOOP - ALL 3 MODELS
#############################

SimAllModels(yrs = 1, numsims = 5, do.save=TRUE)  



