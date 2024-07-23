##############################
## CONNECTIVITY MATRICIES PT 2
############################

# Created by Sophia Horigan (shorigan@uchicago.edu)
# Started: 05-22-24
# Last updated: 07-22-24

# Project Overview: This project seeks to use existing demographic data and the candidate transmission models
# to explore persistence threshholds for hypothetical pathogens under data-based metapopulation structure of Eidolon
# dupreanum in Madagascar. 

## PART 6: Connectivity Pt 2
# Now that I know my two connectivity mechanisms, dispersal and intermingling, and have implemented them,
# I need to develop the functions that take in parameter estimates and make flexible sized
# dispersal and intermingling probabilities

# Dispersal - sub x sub matrix with probability
# Intermingling - sub x sub matrix with probability 

# Accomplishments
# Built DistMat and IMat functions that take in probabilities and return matricies
# added placement and conditions to generate dispersal matricies in code
# updated param values based on brook 2019
# added plotting for just infecteds over time - starting to explore figures demonstrating persistence
# converted param value stochastic draws to lognormal to avoid zero values
# floored all bat pop counts to avoid partial or negative bats
# made GetDistanceMat - pairwise distance between each subpop based on fixed grid size
# set (temporary) upper limits on param draws (uci) to avoid population explosions
# added super loop for varying simulation initial conditions
# added ScaleDistMat which will allow for dispersal prob to be based on distance between subpops

# TO DO
# figure out intermingling. I don't think it makes any sense.
# modify code to use mapply and Rcpp

rm(list=ls())

set.seed(123)

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
library(gridExtra)
library(Rcpp)
library(RcppArmadillo)
library(purrr)

# laptop
#setwd("/Users/sophiahorigan/Documents/GitHub/Bat-disease-metapop/bat-disease-metapop/Models/WIP/")
# desktop
setwd("/Users/shorigan/Documents/GitHub/Bat-disease-metapop/bat-disease-metapop/Models/WIP/")
# midway
#setwd("/home/shorigan/WIP/")

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
BuildTMatMSIRS <- function(model, stoch_disease, c, Npop_epi, Ipop, sum_intermingle, age.classes, surv.biwk, age.brk, surv.juv.biwk, param.array, add.inf.mort, age.rate){
  
  recov = 1 # fixed at 1 biweek for all models
  mu.sick = 1 # fixed for all models
  
  s <- nage <- length(age.classes)
  
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  
  # create param vectors
  beta <- as.vector(rep(100, nage))
  psi <- as.vector(rep(100, nage))
  omega <- as.vector(rep(100, nage))
  
  if(stoch_disease == TRUE){
    for (m in 1:nage){
      while (beta[m] > param.array$beta_max[model]){
        beta[m] <- rlnorm(1, param.array$beta_logmean[model], param.array$beta_logvar[model])
      }
      while (psi[m] > param.array$psi_max[model]){
        psi[m] <- rlnorm(1, param.array$psi_logmean[model], param.array$psi_logvar[model])
      }
      while (omega[m] > param.array$omega_max[model]){
        omega[m] <- rlnorm(1, param.array$omega_logmean[model], param.array$omega_logvar[model])
      }
    }
  }
  
  else if(stoch_disease == FALSE){
    beta <- rep(param.array$beta_mean[model], nage) 
    psi <- rep(param.array$psi_mean[model], nage)
    omega <- rep(param.array$omega_mean[model], nage)
  }
  
  age.rate <- rep(age.rate, nage)
  
  # calculate force of infection
  foi = 1-exp(-beta*((sum(Ipop)/sum(Npop_epi)) + sum_intermingle))
  
  mat1 <- matrix(0,4,4) # MSIRS
  
  Tmat <- matrix(0,4*nage,4*nage) #MSIRS for every age cohort
  for (j in 1:nage) {
    mat1[] <- 0
    mat1[1,1] <- 1- omega[j]
    mat1[2,1] <- omega[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-psi[j]
    mat1[2,4] <- psi[j]
    
    surv <- rep(1-mort.vect[j], 4); # 4 classes
    
    if (add.inf.mort==TRUE){ # additional death from infection (3rd slot in vect)
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
  #print(Tmat)
  #readline(prompt="Press [enter] to continue")
  return(Tmat)
}

# fecundity matrix
BuildFMatMSIRS <- function(age.classes, fecundity, surv.biwk, biwk){ 	#one for every age class
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
  
  recov = 1 # fixed at 1 biweek for all models
  mu.sick = 1 # fixed for all models
  
  # num age classes
  s <- nage <- length(age.classes) 
  
  # create param vectors
  beta <- as.vector(rep(100, nage))
  sigma <- as.vector(rep(100, nage))
  omega <- as.vector(rep(100, nage))
  
  if (model == 2){ # no rho
    rho <- as.vector(rep(0, nage))
  }
  if (model == 3){
    rho <- as.vector(rep(100, nage))
  }
  
  # mortality
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  
  # SET DISEASE PARAMETERS  
  if(stoch_disease == TRUE){
    for (m in 1:nage){
      while (beta[m] > param.array$beta_max[model]){
        beta[m] <- rlnorm(1, param.array$beta_logmean[model], param.array$beta_logvar[model])
      }
      while (sigma[m] > param.array$sigma_max[model]){
        sigma[m] <- rlnorm(1, param.array$sigma_logmean[model], param.array$sigma_logvar[model])
      }
      while (omega[m] > param.array$omega_max[model]){
        omega[m] <- rlnorm(1, param.array$omega_logmean[model], param.array$omega_logvar[model])
      }
      while (rho[m] > param.array$rho_max[model]){
        rho[m] <- rlnorm(1, param.array$rho_logmean[model], param.array$rho_logvar[model])
      }
    }
  }
  
  else if(stoch_disease == FALSE){
    beta <- rep(param.array$beta_mean[model], nage) 
    sigma <- rep(param.array$sigma_mean[model], nage)
    omega <- rep(param.array$omega_mean[model], nage)
    rho <- rep(param.array$rho_mean[model], nage)
  }
  
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
    mat1[1,1] <- 1-omega[j]
    mat1[2,1] <- omega[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[3,5] <- rho[j]
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] # problem
    mat1[5,4] <- sigma[j] # problem
    mat1[5,5] <- 1-rho[j]
    
    # survival 
    surv <- rep(1-mort.vect[j],5); # equal across all classes
    
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
    }else if (biwk==4|biwk ==24){
    new.fec = fecundity*surv.biwk*.05
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

# dispersal matrix 
BuildDMat <- function(num_patches, d.mode, biweek){
  
  if (d.mode == 'seasonal'){
    # seasonal changes in dispersal probability # for example
    if (biweek > 1 && biweek < 13){
      disp.vec = rlnorm(num_patches, param.array$prob_disp_mean_DRY[model], param.array$prob_disp_sd_DRY[model])
    }
    if (biweek > 13 && biweek < 26){
      disp.vec = rlnorm(num_patches, param.array$prob_disp_mean_WET[model], param.array$prob_disp_sd_WET[model])
    }
  }
  
  else if (d.mode == 'annual'){
    #disp.vec = rlnorm(num_patches, param.array$prob_disp_mean[model], param.array$prob_disp_sd[model])
    disp.vec = 0.2772918 # prob of dispersing within a biweek
  }

  # fill in matrix with dispersal out
  for (i in 1:num_patches){
    for (j in 1:num_patches){
      if (i == j){
        DMat[i,j] <<- 0 
      }
      else{
        DMat[i,j] <<- disp.vec * ScaleDistMat[i,j] # dispersal probability to each subpop adjusted by distance 
      }
    }
  }
  
  # calculate probability staying in
  for (i in 1:num_patches){
    for (j in 1:num_patches){
      if (i == j){
        DMat[i,j] <<- 1 - sum(DMat[i,]) # make less than 1 if you want loss of individuals (i.e. 0.001 change of emigrating)
      }
    }
  }

  return(DMat)
}

# intermingling matrix # for larger simulations only
BuildIMat <- function(model, num_patches, DistMat, biweek){
 
  # create tmp matricies
  OMat <- matrix(NA, ncol = num_patches, nrow = num_patches) # overlap
  
  for(i in 1:num_patches){
    for(j in 1:num_patches){
      if (DistMat[i,j] > param.array$int_rad_10_0[model]^2){ # no overlap
        OMat[i,j] = 0
        printf("no overlap")
        print(i)
        print(j)
      }
      else if (DistMat[i,j] < param.array$int_rad_10_0[model]^2){
        OMat[i,j] = param.array$int_rad_10_0[model]^2 - DistMat[i,j]
      }
    }
  }
  for(i in 1:num_patches){
    for(j in 1:num_patches){
      if (OMat[i,j] > 0){
        if (OMat[i,j] < param.array$int_rad_90_85[model]){ # 100 - 0 %
          IMat[i,j] <<- runif(1, min = 0.0, max = 1)
        }
        if (OMat[i,j] > param.array$int_rad_90_85[model] && OMat[i,j] < param.array$int_rad_85_80[model]){ # 89 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.89)
        }
        if (OMat[i,j] > param.array$int_rad_85_80[model] && OMat[i,j] < param.array$int_rad_80_75[model]){ # 84 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.84)
        }
        if (OMat[i,j] > param.array$int_rad_80_75[model] && OMat[i,j] < param.array$int_rad_75_70[model]){ # 79 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.79)
        }
        if (OMat[i,j] > param.array$int_rad_75_70[model] && OMat[i,j] < param.array$int_rad_70_65[model]){ # 74 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.74)
        }
        if (OMat[i,j] > param.array$int_rad_70_65[model] && OMat[i,j] < param.array$int_rad_65_60[model]){ # 69 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.69)
        }
        if (OMat[i,j] > param.array$int_rad_65_60[model] && OMat[i,j] < param.array$int_rad_60_55[model]){ # 64 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.64)
        }
        if (OMat[i,j] > param.array$int_rad_60_55[model] && OMat[i,j] < param.array$int_rad_55_50[model]){ # 59 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.59)
        } 
        if (OMat[i,j] > param.array$int_rad_55_50[model] && OMat[i,j] < param.array$int_rad_50_45[model]){ # 54 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.54)
        }
        if (OMat[i,j] > param.array$int_rad_50_45[model] && OMat[i,j] < param.array$int_rad_45_40[model]){ # 49 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.49)
        }
        if (OMat[i,j] > param.array$int_rad_45_40[model] && OMat[i,j] < param.array$int_rad_40_35[model]){ # 44 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.44)
        }
        if (OMat[i,j] > param.array$int_rad_40_35[model] && OMat[i,j] < param.array$int_rad_35_30[model]){ # 39 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.39)
        }
        if (OMat[i,j] > param.array$int_rad_35_30[model] && OMat[i,j] < param.array$int_rad_30_25[model]){ # 34 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.30)
        }
        if (OMat[i,j] > param.array$int_rad_30_25[model] && OMat[i,j] < param.array$int_rad_25_20[model]){ # 29 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.29)
        }
        if (OMat[i,j] > param.array$int_rad_25_20[model] && OMat[i,j] < param.array$int_rad_20_15[model]){ # 24 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.24)
        }
        if (OMat[i,j] > param.array$int_rad_20_15[model] && OMat[i,j] < param.array$int_rad_15_10[model]){ # 19 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.19)
        }
        if (OMat[i,j] > param.array$int_rad_15_10[model] && OMat[i,j] < param.array$int_rad_10_0[model]){ # 14 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.14)
        }
        if (OMat[i,j] > param.array$int_rad_10_0[model]){ # 9 - 0 %
          IMat[i,j] <<- runif(1, min = 0, max = 0.09)
        }
      }
      else if (OMat[i,j] == 0){
        IMat[i,j] <<- 0
      }
    }
  }

  return(IMat)
}

# initiate and randomly populate
GeneratePops <- function(num_patches, prop_occupied_patches, grid_mode, pop_mode, pop_low, pop_high){
  
  n.sims = length(pop_low)
  
  # check for even number of subpops for grid structure (num_subpop/2 x num_subpop/2 grid)
  if (sqrt(num_patches) %% 1 != 0){
    print("ERROR: need a square divisible number of grid patches.")
    stop()
  }
  
  # make pop model array to hold 1 for occupied, 0 for not occupied
  # same across all sims
  occupation_vector <<- c(rep(0, num_patches))
  
  # determine number of initial occupied patches
  # same across all sims
  num_occupied_patches <<- round(num_patches * prop_occupied_patches, digits = 0)
  
  # check mode for determining which patches are occupied
  if (grid_mode == 'stochastic'){
    # generate 'num_occupied_grids' random numbers from 1 - num_subpops
    random_pop_list = sample(1:num_patches, num_occupied_patches)
    init_pop_list <- c(random_pop_list)
    
    for (i in 1:num_patches){
      if (i %in% init_pop_list){
        occupation_vector[i] = 1
      }
    }
    # cat(num_occupied_patches, " initial occupied grids at locations: ", occupation_vector, "\n")
    
  }
  else if (grid_mode == 'static'){
    # generate 'num_occupied_grids' evenly spaced numbers from 1 - num_subpops
    # determine stepwise interval for occupied patches
    interval = round(num_patches / num_occupied_patches, digits = 0)
    
    # in a list, patches with an initial pop get a 1, without get a 0 
    for (i in 1:num_patches){
      if (i == 1){
        occupation_vector[i] = 1
        ticker = i
      }
      if (i == (ticker + interval) && sum(occupation_vector) < num_occupied_patches){
        
        occupation_vector[i] = 1
        ticker = i
      }
    }
   # cat(num_occupied_patches, " initial occupied grids at locations ", occupation_vector, "\n")
  }
  else {
    print("ERROR: options for grid_mode are 'static or 'stochastic'")
    stop()
  }
  
  # link each patch with a population size
  init_pop_vector <- c(rep(list(rep(0, num_patches)), n.sims))
  
  if (pop_mode == 'stochastic'){ ## THIS ISN'T WORKING YET
    for (i in 1:length(init_pop_vector)){ # patches and sims
      for (j in 1:num_patches){ # patches
        if (occupation_vector[j] == 1){ # patches
          init_pop_vector[[i]][[j]] = sample(1:(pop_low[[i]][[j]]*num_patches), 1) 
        }
      }
    }
     print(init_pop_vector)
  }
  else if (pop_mode == 'static'){
    init_pop_static = pop_low[[1]][[1]]
    for (i in 1:length(init_pop_vector)){
      for (j in 1:num_patches){
        if (occupation_vector[j] == 1){
          init_pop_vector[[i]][[j]] = init_pop_static 
        }
      }
    }
    print(init_pop_vector)
  }
  else {
    print("ERROR: options for grid_mode are 'static or 'stochastic'")
    stop()
  }
  
  return(init_pop_vector)
  
}

# generate pairwise distance matrix between all grids
GenerateDistMat <- function(patch_dim, grid_size){

  ID.key <- matrix(0, ncol = 2, nrow = patch_dim*patch_dim)
  
  # assign grid ID to coordinates
  for (i in 1:patch_dim){
    for (j in 1:patch_dim){
      grid_id = (i * j) + ((patch_dim - i) * (j - 1))
      ID.key[grid_id, 1] <- i
      ID.key[grid_id, 2] <- j
    }
  }
  
  # calculate distance between all points and fill matrix
  dist = grid_size # distance between center points in roost
  
  for (k in 1:(patch_dim*patch_dim)){
    for (l in 1:(patch_dim*patch_dim)){
      if (k == l){
        DistMat[k,l] <<- 0 # no distance between same roost
      }
      else{
        DistMat[k,l] <<- sqrt((ID.key[l, 1] - ID.key[k, 1])^2 + (ID.key[l, 2] - ID.key[k, 2])^2) * dist
      }
    }
  }
    return(DistMat)
}

# scale distance matrix so all values total 1
GenerateScaleDistMat <- function(num_patches){
  
  # rescale to add up to 1
  for(i in 1:num_patches){
    ScaleDistMat[i,] <<- DistMat[i,]/sum(DistMat[i,])
  }
  
  
  #now make everything a probability of dispersal by taking 1- the matrix
  ScaleDistMat <<- 1 - ScaleDistMat

  #and rescale again as a proportion of one:
  for (i in 1:nrow(ScaleDistMat)){
    tmp <- sum(ScaleDistMat[i,])
    past <- ScaleDistMat[i,] 
    ScaleDistMat[i,] <<- past/tmp
    if(round(sum(ScaleDistMat[i,]), digits=4) !=1){
      warning(paste("at i=", i, "rowSums probabilities do not add to 1", sep=" "))
    }
  }

  #and once more the other way
  for (i in 1:ncol(ScaleDistMat)){
    tmp <- sum(ScaleDistMat[,i])
    past <- ScaleDistMat[,i] 
    ScaleDistMat[,i] <<- past/tmp
    
    if(round(sum(ScaleDistMat[,i]), digits=4) !=1){
      warning(paste("at i=", i, "colSums probabilities do not add to 1", sep=" "))
    } 
  }
  
  return(ScaleDistMat)
}

################
## LOAD PARAMS
###############

# param array
param.array <- read.csv("Input/fit_parameters.csv")
param.array

##########################
## SIMULATE MODEL
##########################

SimOneModel <- function(model, num_patches, patch_dim, burnin, sim_pop, yrs, ntyr, age.brk, s, param.array, init_fracI, add.inf.mort, printpop, stoch_disease, stoch_demo, intermingling, dispersal, grid_size, prop_occupied_patches, sim){

  ################
  ## Grid prep
  ##############
  # create global matricies
  DMat <<- matrix(0, ncol = num_patches, nrow = num_patches) # dispersal
  IMat <<- matrix(0, ncol = num_patches, nrow = num_patches) # intermingling
  DistMat <<- matrix(0, ncol = num_patches, nrow = num_patches) # distance
  ScaleDistMat <<- matrix(0, ncol = num_patches, nrow = num_patches) # scaled distance

  # generate distance matrix for grid
  GenerateDistMat(patch_dim = patch_dim, grid_size = grid_size)
  # scale distance matrix
  GenerateScaleDistMat(num_patches = num_patches)
  
  #################################
  ## Generate stable age structure
  ################################
  
  if (model == 1){
    nclass = 4 # M, S, I, R
  }
  else if (model == 2 || model == 3){
    nclass = 5 # M, S, I, R, N
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
  for(i in 1:num_patches){
    
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
    
    if (length(sim_pop) != num_patches){
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
   # print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
    
    # gives counts of bats per age year
    bat.mat = stab.struct*sim_pop[i]
    
    ######################################################
    ## Initiate infection and run dynamics to equilibrium
    ######################################################
    
    if(model==1){ # MSIRS
      R_init = rep(0, s)
      I_init = rep(0, s)
      I_init[3] = floor(init_fracI[i] * sum(bat.mat))
      M_init = rep(0, s)
      S_init = bat.mat - I_init 
      
      N_tot = cbind(M_init, S_init, I_init, R_init)
    }
    else if(model==2 || model==3){ # MSIRN / MSIRNI
      R_init = rep(0, s)
      I_init = rep(0, s)
      I_init[3] = floor(init_fracI[i] * sum(bat.mat))
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
    
    # get biweek to have seasonal changes in dispersal
    biwk1 <- FindBiweek(t=j, times=times)
    
    # generate dispersal matrix
    if (dispersal == TRUE){
      BuildDMat(num_patches = num_patches, d.mode = 'annual', biweek = biwk1)
    }
    
    # generate intermingling matrix
    if (intermingling == TRUE){
      BuildIMat(model = model, num_patches = num_patches, DistMat = DistMat, biweek = biwk1)
    }
    
    ####################
    ## INTERMINGLING
    ###################
    # calculate num infected and total pop for each subpop
    for(i in 1:num_patches){
      
      tmp = TransformVect(vec=N_pop_ts[[i]][,j], s=s, c=nclass)
      Npop_epi[[i]] <- tmp
      
      tmp = MatSplit(matrix(Npop_epi[[i]], ncol=1), r=s, c=1)[,,3] # I is in position 3
      Ipop[[i]] <- tmp
    }
    
    sum_intermingle <- matrix(0, nrow = num_patches, ncol = num_patches)
    
    for (i in 1:num_patches){  
      for (ii in 1:num_patches){
        if (i != ii){
          if (intermingling == TRUE){
            sum_intermingle[i] = sum_intermingle[i] + (sum(Ipop[[ii]])/sum(Npop_epi[[ii]]) * IMat[i,ii]) 
          }
          else if (intermingling == FALSE){
            sum_intermingle[i] = 0
          }
        }
      }
    }
    
    # compute matrices and iterate one timestep
    for(i in 1:num_patches){
      
      if(model==1){ # MSIRS
        # generate SIR matrix
        tmp = BuildTMatMSIRS(model = model, Npop_epi = Npop_epi[[i]], Ipop = Ipop[[i]], sum_intermingle = sum_intermingle[i], stoch_disease = stoch_disease, param.array = param.array, c = nclass, age.classes = 1:s, age.brk = age.brk, surv.biwk = (1-mort_ad[i])^(1/ntyr), surv.juv.biwk = (1-mort_juv[i])^(1/ntyr), add.inf.mort = add.inf.mort, age.rate = 1/ntyr)
        Tmat[[i]] <- tmp
        #readline(prompt="Press [enter] to continue")
        
        # generate fecundity matrix
        Fmat <- BuildFMatMSIRS(age.classes = 1:s, fecundity = fecundity[i], surv.biwk = (1-mort_ad[i])^(1/ntyr), biwk = biwk1)
        
      }
      else if(model==2 || model==3){ # MSIRN, MSIRNI
        # generate SIR matrix
        tmp <- BuildTMatMSIRN(model = model, Npop_epi = Npop_epi[[i]], Ipop = Ipop[[i]], sum_intermingle = sum_intermingle[i], stoch_disease = stoch_disease, param.array = param.array, c = nclass, age.classes = 1:s, age.brk = age.brk, surv.biwk = (1-mort_ad[i])^(1/ntyr), surv.juv.biwk = (1-mort_juv[i])^(1/ntyr), add.inf.mort = add.inf.mort, age.rate = 1/ntyr)
        Tmat[[i]] <- tmp
        
        # generate fecundity matrix
        Fmat <- BuildFMatMSIRN(age.classes = 1:s, fecundity = fecundity[i], surv.biwk = (1-mort_ad[i])^(1/ntyr), biwk = biwk1)
        
      }
      
      # make transition matrix
      transMat <- Tmat[[i]] + Fmat
      
      #move forward in time
      N_pop_ts[[i]][,j+1] <- transMat %*% N_pop_ts[[i]][,j]
      
    }
    
    # stitch current pops into matrix
    tmp2 <- list()
    for (i in 1:num_patches){
      tmp <- N_pop_ts[[i]][,j+1]
      tmp2[[i]] <- tmp
    }
    N_pop <- do.call(rbind, tmp2)
    #print(N_pop)
    
    ##############################
    ## DISPERSAL
    #############################
    
    if (dispersal == TRUE){
      N_pop <- DMat %*% N_pop
    }
    
    # round to nearest individual and remove negative values
    N_pop <- pmax(floor(N_pop), 0)
    # convert nans to 0 
    N_pop[is.nan(N_pop)] = 0
    
    # split resulting matrix and fill in pop list
    for (i in 1:num_patches){
      N_pop_ts[[i]][,j+1] <- N_pop[i,]
      # print("N pop at t+1")
      #print(N_pop_ts[[i]][,j+1])
    }
  }
  
  for(i in 1:num_patches){
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
  
  ## add other important columns
  final.disease$grid_size = grid_size
  final.disease$disp = dispersal
  final.disease$inter = intermingling
  final.disease$num_patch = num_patches
  final.disease$init_patch_occ = prop_occupied_patches
  final.disease$model = model
  final.disease$sim = sim
  
  final.pop$grid_size = grid_size
  final.pop$disp = dispersal
  final.pop$inter = intermingling
  final.pop$num_patch = num_patches
  final.pop$init_patch_occ = prop_occupied_patches
  final.pop$model = model
  final.pop$sim = sim
  
  cat("sim ", sim, "\n")
  
  # CAN ONLY RETURN ONE ARRAY IN R SO DUMB
  # NEED TO FIGURE OUT HOW TO COMBINE THEM BAH
  
  
  if(printpop == FALSE){
   return(final.disease)
  }
  if(printpop == TRUE){
    return(final.pop)
  }
}

Plot <- function(sims, do.save.plot){
  
  # get mean and CI for stochastic realizations
  sims <- sims %>% mutate(proportion = ifelse(is.na(proportion), 0, proportion))
  
  sim.summary <<- sims %>%
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
  
  # plot by model type
  if(sims$model[1] == 1){
    colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue")
    
    # PROPORTION
    p1_1 <<- ggplot(data=sim.summary) + geom_line(aes(x=time, y=mean.proportion, color=state)) + facet_wrap(~subpop, ncol = sqrt(as.numeric(max(sim.summary$subpop)))) +
      theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
      geom_ribbon(aes(x=time, ymin = lower.ci.proportion, ymax = upper.ci.proportion, fill = state), alpha = 0.2) +
      scale_color_manual(values=colz) +
      scale_x_continuous(breaks=seq(0, max(sim.summary$time), 1)) + ggtitle(paste0("MSIRS PROP. numsims = ", max(sims$sim))) +
      xlab("Year")
    print(p1_1)
    
    # COUNT
    p1_2 <<- ggplot(data=sim.summary) + geom_line(aes(x=time, y=mean.totpop)) + facet_wrap(~subpop, ncol = sqrt(as.numeric(max(sim.summary$subpop)))) +
      geom_line(aes(x=time, y=mean.count, color=state)) +
      theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
      geom_ribbon(aes(x=time, ymin = lower.ci.count, ymax = upper.ci.count, fill = state), alpha = 0.2) +
      geom_ribbon(aes(x=time, ymin = lower.ci.totpop, ymax = upper.ci.totpop), alpha = 0.2) +
      scale_color_manual(values=colz) +
      scale_x_continuous(breaks=seq(0, max(sim.summary$time), 1)) + ggtitle(paste0("MSIRS COUNT")) +
      xlab("Year")
    print(p1_2)
    
    # INFECTED NOT AVERAGED
    sim.inf <<- sims %>% filter(state == 'I') 
    sim.inf.summary <<- sim.summary %>% filter(state == 'I') 
    
    #    p1_3 <<- ggplot(data = sim.inf, aes(x = time, y = proportion, color = as.factor(sim))) + geom_line() + facet_wrap(~subpop, ncol = patch_dim) +
    #     theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    #    theme(plot.title = element_text(hjust = 0.5)) +
    #   scale_x_continuous(breaks=seq(0, max(sims$time), 10)) + ggtitle(paste0("MSIRS")) +
    #  xlab("Year") + ylab("Proportion Infected")
    #  print(p1_3)
    
    # INFECTEDS AVERAGED
    p1_4 <<- ggplot(data = sim.inf.summary, aes(x = time, y = mean.proportion)) + geom_line() + facet_wrap(~subpop, ncol = sqrt(as.numeric(max(sim.summary$subpop)))) +
      geom_ribbon(aes(x=time, ymin = lower.ci.proportion, ymax = upper.ci.proportion, fill = state), alpha = 0.2) +
      theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(breaks=seq(0, max(sims$time), 10)) + ggtitle(paste0("MSIRS")) +
      xlab("Year") + ylab("Proportion Infected")
    print(p1_4)
    
    if(do.save.plot==TRUE){
      ggsave(path = "Output", file = paste0("Frac_Model_", model, "_Yrs_", yrs, "_", Sys.time(), ".png"),
             plot=p1_1,
             units="mm",  
             width=70, 
             height=60, 
             scale=3, 
             dpi=300)
      ggsave(path = "Output", file = paste0("Count_Model_", model, "_Yrs_", yrs, "_", Sys.time(), ".png"),
             plot=p1_2,
             units="mm",  
             width=70, 
             height=60, 
             scale=3, 
             dpi=300)
      ggsave(path = "Output", file = paste0("Inf_Model_Avg", model, "_Yrs_", yrs, "_", Sys.time(), ".png"),
             plot=p1_3,
             units="mm",  
             width=70, 
             height=60, 
             scale=3, 
             dpi=300)
      ggsave(path = "Output", file = paste0("Inf_Model_", model, "_Yrs_", yrs, "_", Sys.time(), ".png"),
             plot=p1_4,
             units="mm",  
             width=70, 
             height=60, 
             scale=3, 
             dpi=300)
    }
  }
  if(sims$model[1] == 2){
    colz <<- c('M'="violet", 'S' = "mediumseagreen", 'I' = "tomato", 'R' = "cornflowerblue", 'N' = "navy")
    # PROPORTION
    p2_1 <<- ggplot(data = sim.summary) + geom_line(aes(x=time, y=mean.proportion, color=state)) + facet_wrap(~subpop, ncol = sqrt(as.numeric(max(sim.summary$subpop)))) +
      theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
      geom_ribbon(aes(x=time, ymin = lower.ci.proportion, ymax = upper.ci.proportion, fill = state), alpha = 0.2) +
      scale_color_manual(values=colz) +
      scale_x_continuous(breaks=seq(0, max(sim.summary$time), 1)) + ggtitle(paste0("MSIRN PROP")) +
      xlab("Year")
    print(p2_1)
    
    # COUNT
    p2_2 <<- ggplot(data=sim.summary) + geom_line(aes(x=time, y=mean.totpop)) + facet_wrap(~subpop, ncol = sqrt(as.numeric(max(sim.summary$subpop)))) +
      geom_line(aes(x=time, y=mean.count, color=state)) +
      theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
      geom_ribbon(aes(x=time, ymin = lower.ci.count, ymax = upper.ci.count, fill = state), alpha = 0.2) +
      geom_ribbon(aes(x=time, ymin = lower.ci.totpop, ymax = upper.ci.totpop), alpha = 0.2) +
      scale_color_manual(values=colz) +
      scale_x_continuous(breaks=seq(0, max(sim.summary$time), 1)) + ggtitle(paste0("MSIRN COUNT")) +
      xlab("Year")
    print(p2_2)
    
    # INFECTED NOT AVERAGED
    sim.inf <<- sims %>% filter(state == 'I') 
    sim.inf.summary <<- sim.summary %>% filter(state == 'I') 
    
    #  p2_3 <<- ggplot(data = sim.inf, aes(x = time, y = proportion, color = as.factor(sim))) + geom_line() + facet_wrap(~subpop, ncol = patch_dim) +
    #     theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    #    theme(plot.title = element_text(hjust = 0.5)) +
    #   scale_x_continuous(breaks=seq(0, max(sims$time), 10)) + ggtitle(paste0("MSIRN")) +
    #  xlab("Year") + ylab("Proportion Infected")
    #print(p2_3)
    
    # INFECTEDS AVERAGED
    p2_4 <<- ggplot(data = sim.inf.summary, aes(x = time, y = mean.proportion)) + geom_line() + facet_wrap(~subpop, ncol = sqrt(as.numeric(max(sim.summary$subpop)))) +
      geom_ribbon(aes(x=time, ymin = lower.ci.proportion, ymax = upper.ci.proportion, fill = state), alpha = 0.2) +
      theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(breaks=seq(0, max(sims$time), 10)) + ggtitle(paste0("MSIRN")) +
      xlab("Year") + ylab("Proportion Infected")
    print(p2_4)
    
    
    if(do.save.plot==TRUE){
      ggsave(path = "Output", file = paste0("Frac_Model_", model, "_Yrs_", yrs, "_", Sys.time(), ".png"),
             plot=p2_1,
             units="mm",  
             width=70, 
             height=60, 
             scale=3, 
             dpi=300)
      ggsave(path = "Output", file = paste0("Count_Model_", model, "_Yrs_", yrs, "_", Sys.time(), ".png"),
             plot=p2_2,
             units="mm",  
             width=70, 
             height=60, 
             scale=3, 
             dpi=300)
      ggsave(path = "Output", file = paste0("Inf_Model_Avg", model, "_Yrs_", yrs, "_", Sys.time(), ".png"),
             plot=p2_3,
             units="mm",  
             width=70, 
             height=60, 
             scale=3, 
             dpi=300)
      ggsave(path = "Output", file = paste0("Inf_Model_", model, "_Yrs_", yrs, "_", Sys.time(), ".png"),
             plot=p2_4,
             units="mm",  
             width=70, 
             height=60, 
             scale=3, 
             dpi=300)
    }
  }
  if(sims$model[1] == 3){
    colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N' = "navy")
    
    # PROPORTION
    p3_1 <<- ggplot(data=sim.summary) + geom_line(aes(x=time, y=mean.proportion, color=state)) + facet_wrap(~subpop, ncol = sqrt(as.numeric(max(sim.summary$subpop)))) +
      theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
      geom_ribbon(aes(x=time, ymin = lower.ci.proportion, ymax = upper.ci.proportion, fill = state), alpha = 0.2) +
      scale_color_manual(values=colz) +
      scale_x_continuous(breaks=seq(0, max(sim.summary$time), 1)) + ggtitle("MSIRNI PROP") +
      xlab("Year")
    print(p3_1)
    
    # COUNT
    p3_2 <<- ggplot(data=sim.summary) + geom_line(aes(x=time, y=mean.totpop)) + facet_wrap(~subpop, ncol = sqrt(as.numeric(max(sim.summary$subpop)))) +
      geom_line(aes(x=time, y=mean.count, color=state)) +
      theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
      geom_ribbon(aes(x=time, ymin = lower.ci.count, ymax = upper.ci.count, fill = state), alpha = 0.2) +
      geom_ribbon(aes(x=time, ymin = lower.ci.totpop, ymax = upper.ci.totpop), alpha = 0.2) +
      scale_color_manual(values=colz) +
      scale_x_continuous(breaks=seq(0, max(sim.summary$time), 1)) + ggtitle("MSIRNI COUNT") +
      xlab("Year")
    print(p3_2)
    
    # INFECTED NOT AVERAGED
    sim.inf <<- sims %>% filter(state == 'I') 
    sim.inf.summary <<- sim.summary %>% filter(state == 'I') 
    
    #  p3_3 <<- ggplot(data = sim.inf, aes(x = time, y = proportion, color = as.factor(sim))) + geom_line() + facet_wrap(~subpop, ncol = patch_dim) +
    #   theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
    #  theme(plot.title = element_text(hjust = 0.5)) +
    # scale_x_continuous(breaks=seq(0, max(sims$time), 10)) + ggtitle(paste0("MSIRNI")) +
    #xlab("Year") + ylab("Proportion Infected")
    #  print(p3_3)
    
    # INFECTEDS AVERAGED
    p3_4 <<- ggplot(data = sim.inf.summary, aes(x = time, y = mean.proportion)) + geom_line() + facet_wrap(~subpop, ncol = sqrt(as.numeric(max(sim.summary$subpop)))) +
      geom_ribbon(aes(x=time, ymin = lower.ci.proportion, ymax = upper.ci.proportion, fill = state), alpha = 0.2) +
      theme_bw() + theme(panel.grid = element_blank(), strip.background = element_rect(fill="white")) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(breaks=seq(0, max(sims$time), 10)) + ggtitle(paste0("MSIRNI")) +
      xlab("Year") + ylab("Proportion Infected")
    print(p3_4)
    
    if(do.save.plot==TRUE){
      ggsave(path = "Output", file = paste0("Frac_Model_", model, "_Yrs_", yrs, "_", Sys.time(), ".png"),
             plot=p3_1,
             units="mm",  
             width=70, 
             height=60, 
             scale=3, 
             dpi=300)
      
      ggsave(path = "Output", file = paste0("Count_Model_", model, "_Yrs_", yrs, "_", Sys.time(), ".png"),
             plot=p3_2,
             units="mm",  
             width=70, 
             height=60, 
             scale=3, 
             dpi=300)
      ggsave(path = "Output", file = paste0("Inf_Model_Avg", model, "_Yrs_", yrs, "_", Sys.time(), ".png"),
             plot=p3_3,
             units="mm",  
             width=70, 
             height=60, 
             scale=3, 
             dpi=300)
      ggsave(path = "Output", file = paste0("Inf_Model_", model, "_Yrs_", yrs, "_", Sys.time(), ".png"),
             plot=p3_4,
             units="mm",  
             width=70, 
             height=60, 
             scale=3, 
             dpi=300)
    }
  }
}

#####################
## RUN MODEL
####################
sims <- data.frame()

# RCPP
## set number of sims
yrs = 2
numsims = 5
prop_occupied_patches = 1
do.save.data = FALSE
burnin = 0
stoch_demo = TRUE
stoch_disease = TRUE

model_vect <- as.vector(c(1))
grid_vect = as.vector(c(100))
patch_vect = as.vector(c(2))
pop_size = as.vector(c(10000))
intermingle_vect <- as.vector(c("TRUE", "FALSE"))
dispersal_vect <- as.vector(c("TRUE", "FALSE"))

 
## mapply
iter = 1
for (i in 1:length(model_vect)){ # models
  for (j in 1:length(grid_vect)){ # grid sizes
    for (k in 1:length(patch_vect)){ # patch dim
      for (m in 1:length(pop_size)){ # init population size
        for (n in 1:length(dispersal_vect)){
            
          cat("Model: ", i, " Grid size: ", grid_vect[j], " Patch dim: ", patch_vect[k], " Pop size: ", pop_size[m], " Dispersal: ", dispersal_vect[n], " Intermingling: ", intermingle_vect[n], "\n")
          
          
          model.list = c(rep(i, numsims))
          num_patches.list = c(rep(patch_vect[k]^2, numsims))
          patch_dim.list = c(rep(patch_vect[k], numsims))
          grid_size.list = c(rep(grid_vect[j], numsims))
          prop_occupied_patches.list = c(rep(prop_occupied_patches, numsims))
          yrs.list = c(rep(yrs, numsims))
          pop_low.list = rep(list(rep(pop_size[m]/patch_vect[k]^2, patch_vect[k]^2)), numsims)
          pop_high.list = rep(list(rep(pop_size[m]/patch_vect[k]^2, patch_vect[k]^2)), numsims)
          sim_pop.list = GeneratePops(num_patches = patch_vect[k]^2, prop_occupied_patches = prop_occupied_patches, grid_mode = "stochastic", pop_mode = "static", pop_low = pop_low.list, pop_high = pop_high.list)
          intermingling.list = c(rep(intermingle_vect[n], numsims))
          dispersal.list = c(rep(dispersal_vect[n], numsims))
          do.save.data.list = c(rep(do.save.data, numsims))
          sim.list = c(rep(1:numsims))
          burnin.list = c(rep(burnin, numsims))
          ntyr.list = c(rep(26, numsims))
          s.list = c(rep(20, numsims))
          age.brk.list = c(rep(20, numsims))
          init_fracI.list = rep(list(runif(patch_vect[k]^2, min = 0.1, max = 0.45)), numsims) 
          add.inf.mort.list = c(rep(FALSE, numsims))
          printpop.list = c(rep(FALSE, numsims))
          stoch_demo.list = c(rep(stoch_demo, numsims))
          stoch_disease.list = c(rep(stoch_disease, numsims))
          param.array.list = list(param.array)
          
          print("iteration")
          print(iter)
          
          if (iter == 1){
            output <<- mapply(SimOneModel, 
                              model = model.list,
                              num_patches = num_patches.list,
                              patch_dim = patch_dim.list,
                              grid_size = grid_size.list,
                              prop_occupied_patches = prop_occupied_patches.list,
                              yrs = yrs.list,
                              sim_pop = sim_pop.list,
                              intermingling = intermingling.list,
                              dispersal = dispersal.list,
                              burnin = burnin.list,
                              ntyr = ntyr.list,
                              s = s.list,
                              age.brk = age.brk.list,
                              param.array = param.array.list,
                              init_fracI = init_fracI.list,
                              add.inf.mort = add.inf.mort.list,
                              printpop = printpop.list,
                              stoch_demo = stoch_demo.list,
                              stoch_disease = stoch_disease.list,
                              sim = sim.list,
                              SIMPLIFY = FALSE)
            output <<- map_dfr(output, ~as.data.frame(.x))
          }
          else{
            tmp <<- mapply(SimOneModel, 
                           model = model.list,
                           num_patches = num_patches.list,
                           patch_dim = patch_dim.list,
                           grid_size = grid_size.list,
                           prop_occupied_patches = prop_occupied_patches.list,
                           yrs = yrs.list,
                           sim_pop = sim_pop.list,
                           intermingling = intermingling.list,
                           dispersal = dispersal.list,
                           burnin = burnin.list,
                           ntyr = ntyr.list,
                           s = s.list,
                           age.brk = age.brk.list,
                           param.array = param.array.list,
                           init_fracI = init_fracI.list,
                           add.inf.mort = add.inf.mort.list,
                           printpop = printpop.list,
                           stoch_demo = stoch_demo.list,
                           stoch_disease = stoch_disease.list,
                           sim = sim.list,
                           SIMPLIFY = FALSE)
            tmp <<- map_dfr(tmp, ~as.data.frame(.x))
            output <<- rbind(tmp, output)
          }
          iter = iter+1
        }
      }
    }
  }
}




  
# plot
Plot(sims, do.save.plot = FALSE)


