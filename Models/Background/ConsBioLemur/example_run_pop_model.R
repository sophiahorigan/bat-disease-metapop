#add in connectivity
#simulation plots

#first compare lambda, then simulation of an annual lefkovitch matrix vs. an explicit
#leslie matrix for a birth pulse population

#first build lambda for the two different structures and compare
#then simulate with density dependence, first with hunting-only mortality, then
#with elevated mortality.
# then see if connectivity can save it

#then apply hunt in three different seasons
#highlight which scenario applies best to which species (i.e. Indri and Aye-aye are different)

#Finally, take 1 scenario and ask if connectivity can save it

#first compare lambda
rm(list=ls())

homewd= "/Users/shorigan/Desktop/Mada-Bat-Metapopulation/Models"

#set wd
setwd(paste0(homewd, "/prior-scripts/ConsBioLemur/"))
#.libPaths("/home/caraeb/R/x86_64-redhat-linux-gnu-library/3.4/")

library(reshape2)
library(dplyr)
library(plyr)
library(ggplot2)
library(stringr)
library(plyr)
library(matrixStats)
library(RColorBrewer)
library(MASS)

#load data

load(file=paste0(homewd, "/prior-scripts/ConsBioLemur/data/local_lambdas.Rdata"))
head(lambda.df) 
#this loads a table with the population characteristics for each lemur species,
#as well as the lat/long of wach site.

#now load all the functions
noise <- function(param){ # no noise to zero params but shakes up 
  if(param == 0){             #the others (both fecundity and survivorship) quite a bit
    return(param)
  }else{
    return(rnorm(1, mean=param, sd = .1))
  }
} # noise function for stochastic simulations. jiters .01 std deviations in each direction
find.prop.zero <- function(Nnext2, a){
  if (a==1){
    site.list <- slice(Nnext2, 2)
    site.sum <- lapply(site.list, sum)
    site.sum <- c(unlist(site.sum))
    prop.zero=length(site.sum[site.sum!=0])/length(site.sum)
  }
  if (a>=2){
    site.list <- slice(Nnext2, 3)
    site.sum <- lapply(site.list, sum)
    site.sum <- c(unlist(site.sum))
    prop.zero = length(site.sum[site.sum!=0])/length(site.sum)
  }
  
  return(prop.zero)
}
collapse.adults <- function(Nmat.str, n.sites){
  
  infants <- list()
  juveniles <- list()
  adults <- list()
  
  Nmeta.pop.list <- list()
  for (i in 1:n.sites){
    infants[[i]] <- Nmat.str[1,,i]
    juveniles[[i]] <- Nmat.str[2,,i]
    adults[[i]] <- colSums(Nmat.str[3:4,,i])
    
    Nmeta.pop.list[[i]] <- rbind(infants[[i]], juveniles[[i]], adults[[i]])
  }
  return(Nmeta.pop.list)
  
}
wrap.noise <- function(pop.mat){
  matrix(sapply(pop.mat,noise), 
         byrow = F, ncol=ncol(pop.mat))
} #wrapper to apply over every matrix
wrap.lambda <- function(df, nat.mort, use.uci, leslie, seas){
  #first, sort by species, site, scenario
  df <- arrange(df, species, site, scenario)
  #build a separate population matrix for every site and scenario
  df.list <- dlply(df, .(species, site, scenario))
  
  #apply pop mat function. 
  #let it depend on whether you want a lefkovitch or a leslie matrix
  #default is lefkovitch. need to specify if you want seasonal or annual leslie matrix
  if(leslie=="annual"){
    pop.mat.list <- lapply(df.list, build.leslie.annual, nat.mort=nat.mort, use.uci=use.uci)  
  }
  if(leslie=="seasonal"){
    pop.mat.list <- lapply(df.list, build.leslie.seasonal, nat.mort=nat.mort, use.uci=use.uci, seas=seas)  
  }
  if(leslie==FALSE){
    pop.mat.list <- lapply(df.list, build.pop.mat, nat.mort=nat.mort, use.uci=use.uci)  
  }
  
  
  #apply get_lambda function
  lambda.list <- lapply(pop.mat.list, get_lambda)
  lambdas <- c(unlist(lambda.list))
  #now attach back to the site dataframe
  df$lambda <-  lambdas 
  return(df)
}
build.pop.mat <- function(df, nat.mort, use.uci, dens.dep, N, K){
  #add in the option to use the upper confidence limit of the mortality rate
  if(use.uci==TRUE){
    df$IBI_adult_mortality = df$IBI_upr
  }
  #add in the option to elevate mortality beyond that captured in the hunt
  df$IBI_adult_mortality =  df$IBI_adult_mortality + nat.mort
  #remember that we can't have mortality above 100%, so correct for that here
  df$IBI_adult_mortality[df$IBI_adult_mortality>=1] <- .9
  
  p = 1-as.numeric(as.character(df$IBI_adult_mortality))
  pi= 1-as.numeric(as.character(df$IBI_infant_mortality))
  a=  as.numeric(as.character(df$age_1st_rep))/as.numeric(as.character(df$IBI)) #age 1st rep in IBI units. all scenarios have been retsricted to instances in which a is a multiple of IBI
  Fap= as.numeric(as.character(df$IBI_fecundity))*p
  Fap1= (as.numeric(as.character(df$IBI_fecundity)))*p #first year breeders (same as others unles density-dependence comes into play)
  Faip=(as.numeric(as.character(df$IBI_fecundity)))*pi
  
  if(dens.dep==TRUE){
    Faip <- Faip*(1-N/K)
    Faip[Faip<0] <- 0 #no negative births allowed
    
    Fap1 <- Fap1*(1-N/K)
    Fap1[Fap1<0] <- 0
    
    #and all other reproduction is acted on less strongly - but we need a way to hold it steady
    Fap <- Fap*(1-(N/(K*2)))
    Fap[Fap<0] <- 0
  }
  
  #now build your pop matrix depending on your age at 1st rep
  pop.mat <- NA
  
  if(a==1){
    pop.mat <- matrix(c(Faip,Fap,
                        pi,p), nrow=2, ncol=2, byrow=TRUE)
  }
  if(a==2){
    pop.mat <- matrix(c(0,Fap1,Fap,
                        pi,0,0,
                        0,p, p), nrow=3, ncol=3, byrow=TRUE)
  }
  if(a>=3){
    pop.mat <- matrix(c(0,0,Fap1,Fap,
                        pi,0,0,0,
                        0, p,0,0,
                        0,0, p^(a-2),p), nrow=4, ncol=4, byrow=TRUE)
  }
  
  return(pop.mat)
}
build.leslie.annual <- function(df, nat.mort, use.uci, dens.dep, N, K){
  #add in the option to use the upper confidence limit of the mortality rate
  if(use.uci==TRUE){
    df$IBI_adult_mortality = df$IBI_upr
  }
  #add in the option to elevate mortality beyond that captured in the hunt
  df$IBI_adult_mortality =  df$IBI_adult_mortality + nat.mort
  #remember that we can't have mortality above 1===00%, so correct for that here
  df$IBI_adult_mortality[df$IBI_adult_mortality>=1] <- .9
  
  #now convert back to annual:
  adult_mort <- as.numeric(as.character(df$IBI_adult_mortality))/as.numeric(as.character(df$IBI))
  infant_mort <- as.numeric(as.character(df$IBI_infant_mortality))/as.numeric(as.character(df$IBI))
  
  p = 1-adult_mort
  pi= 1-infant_mort
  a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
  Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p #make into years again
  Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p #first year breeders
  Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
  L = as.numeric(as.character(df$avg_lifespan))
  
  if(dens.dep==TRUE){
    Faip <- Faip*(1-N/K) #only applicable if a==1
    Faip[Faip<0] <- 0 #no negative births allowed
    
    Fap1 <- Fap1*(1-N/K)
    Fap1[Fap1<0] <- 0
    
    #and all other reproduction is acted on less strongly - but we need a way to hold it steady
    Fap <- Fap*(1-(N/(K*2)))
    Fap[Fap<0] <- 0
  }
  
  #now build your pop matrix depending on your age at 1st rep
  #sketch the skeleton 
  pop.mat <- matrix(c(rep(0,length=(L)),rep(0,length=(L))), nrow=(L),ncol=(L))
  
  #add  the survivorship part
  pop.mat.s <- matrix(c(rep(0,length=(L-1)),rep(0,length=(L-1))), nrow=(L-1),ncol=(L-1))
  #we have to assume that infancy lasts for a full year in this case.
  #we fill the diagonal - first is always infant, then juvenile for 1, 2, or 3 slots (or more!), then the rest
  #since juveniles and adults have same survivorship, this bit is easy
  diag(pop.mat.s) <-c(pi, rep(p, length=(L-2)))
  #ends 2 short because you tack on the last adult later
  
  #make vectors for the top and side:
  #structure of top depends on a:
  if(a==1){
    top.row <- c(Faip, rep(Fap, length=(L-1)))
  }
  if(a==2){
    top.row <- c(0,Fap1, rep(Fap, length=(L-2)))
  }
  if(a>=3){
    top.row <- c(rep(0, length=(a-1)),Fap1,rep(Fap, length=(L-a)))
  }
  
  final.col <- c(Fap,rep(0, length=L-2), p)
  
  #and fill the skeleton:
  pop.mat[1,] <- top.row
  pop.mat[,ncol(pop.mat)] <- final.col
  pop.mat[2:nrow(pop.mat),1:(ncol(pop.mat)-1)] <- pop.mat.s
  
  #and return
  
  return(pop.mat)
}
build.leslie.seasonal <- function(df, nat.mort, use.uci, seas, dens.dep, N, K){
  #add in the option to use the upper confidence limit of the mortality rate
  if(use.uci==TRUE){
    df$IBI_adult_mortality = df$IBI_upr
  }
  #add in the option to elevate mortality beyond that captured in the hunt
  df$IBI_adult_mortality =  df$IBI_adult_mortality + nat.mort
  #remember that we can't have mortality above 100%, so correct for that here
  df$IBI_adult_mortality[df$IBI_adult_mortality>=1] <- .9
  
  #now convert back to annual and then make seasonal
  adult_mort <- (as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))/seas
  infant_mort <- (as.numeric(df$IBI_infant_mortality)/as.numeric(df$IBI))/seas
  
  p = 1-adult_mort
  pi= 1-infant_mort
  a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
  Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p #make into years again
  Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p #first year breeders
  Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
  L = as.numeric(as.character(df$avg_lifespan))
  
  if(dens.dep==TRUE){
    Faip <- Faip*(1-N/K) #only applicable if a==1
    Faip[Faip<0] <- 0 #no negative births allowed
    
    Fap1 <- Fap1*(1-N/K)
    Fap1[Fap1<0] <- 0
    
    #and all other reproduction is acted on less strongly - but we need a way to hold it steady
    Fap <- Fap*(1-(N/(K*2)))
    Fap[Fap<0] <- 0
  }
  
  #now build your pop matrix depending on your age at 1st rep
  #sketch the skeleton 
  pop.mat <- matrix(c(rep(0,length=(L*seas)),rep(0,length=(L*seas))), nrow=(L*seas),ncol=(L*seas))
  
  #add  the survivorship part
  pop.mat.s <- matrix(c(rep(0,length=(L*seas-1)),rep(0,length=(L*seas-1))), nrow=(L*seas-1),ncol=(L*seas-1))
  #we have to assume that infancy lasts for a full year in this case.
  #we fill the diagonal - first is always infant, then juvenile for 1, 2, or 3 slots (or more!), then the rest
  #since juveniles and adults have same survivorship, this bit is easy
  diag(pop.mat.s) <-c(rep(pi, times=seas), rep(p, length=(L*seas-seas-1)))
  #ends 2 short because you tack on the last adult later
  
  #make vectors for the top and side:
  #structure of top depends on a:
  if(a==1){
    top.row <- c(c(rep(0,length=(seas-1)), Faip), rep(c(rep(0,length=(seas-1)), Fap), times=(L-1)))
  }
  if(a==2){
    top.row <- c(c(rep(0,length=(seas*2-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-2)))
  }
  if(a>=3){
    top.row <- c(c(rep(0,length=(seas*a-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-a)))
  }
  
  final.col <- c(Fap,rep(0, length=L*seas-2), p)
  
  #and fill the skeleton:
  pop.mat[1,] <- top.row
  pop.mat[,ncol(pop.mat)] <- final.col
  pop.mat[2:nrow(pop.mat),1:(ncol(pop.mat)-1)] <- pop.mat.s
  
  #and return
  
  return(pop.mat)
}
build.leslie.seasonal.hunt <- function(df, nat.mort, use.uci, seas, dens.dep, N, K, hunt_sim){
  #add in the option to use the upper confidence limit of the mortality rate
  if(use.uci==TRUE){
    df$IBI_adult_mortality = df$IBI_upr
  }
  #add in the option to elevate mortality beyond that captured in the hunt
  df$IBI_adult_mortality =  df$IBI_adult_mortality + nat.mort
  #remember that we can't have mortality above 100%, so correct for that here
  df$IBI_adult_mortality[df$IBI_adult_mortality>=1] <- .9
  
  #now convert back to annual and then make seasonal
  adult_mort <- (as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))/seas
  infant_mort <- (as.numeric(df$IBI_infant_mortality)/as.numeric(df$IBI))/seas
  
  #now build matrix based on hunt_sim
  
  if (hunt_sim==1){
    adult_mort_hi <- (as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)/3
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    #fecundity also affected, depending on 
    Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #make into years again
    Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #first year breeders
    Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
    L = as.numeric(as.character(df$avg_lifespan))
    
    if(dens.dep==TRUE){
      Faip <- Faip*(1-N/K) #only applicable if a==1
      Faip[Faip<0] <- 0 #no negative births allowed
      
      Fap1 <- Fap1*(1-N/K)
      Fap1[Fap1<0] <- 0
      
      #and all other reproduction is acted on less strongly - but we need a way to hold it steady
      Fap <- Fap*(1-(N/(K*2)))
      Fap[Fap<0] <- 0
    }
    
    #now build your pop matrix depending on your age at 1st rep
    #sketch the skeleton 
    pop.mat <- matrix(c(rep(0,length=(L*seas)),rep(0,length=(L*seas))), nrow=(L*seas),ncol=(L*seas))
    
    #add  the survivorship part
    pop.mat.s <- matrix(c(rep(0,length=(L*seas-1)),rep(0,length=(L*seas-1))), nrow=(L*seas-1),ncol=(L*seas-1))
    #we have to assume that infancy lasts for a full year in this case.
    #we fill the diagonal - first is always infant, then juvenile for 1, 2, or 3 slots (or more!), then the rest
    
    #we have to be careful here since our adult mortality varies by sim. 
    #for hunt_sim=1, it is every 1 of 4
    diag(pop.mat.s) <-(c(rep(pi, times=seas), rep(c(p_hi_mort, p_low_mort,p_low_mort, p_low_mort), times=(L-1)))[-1])
    
    #ends 2 short because you tack on the last adult later
    
    #make vectors for the top and side:
    #structure of top depends on a:
    if(a==1){
      top.row <- c(c(rep(0,length=(seas-1)), Faip), rep(c(rep(0,length=(seas-1)), Fap), times=(L-1)))
    }
    if(a==2){
      top.row <- c(c(rep(0,length=(seas*2-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-2)))
    }
    if(a>=3){
      top.row <- c(c(rep(0,length=(seas*a-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-a)))
    }
    
    final.col <- c(Fap,rep(0, length=L*seas-2), p_low_mort)
  }
  if (hunt_sim==2){
    adult_mort_hi <- (as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)/3
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    #fecundity also affected, depending on 
    Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #make into years again
    Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #first year breeders
    Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
    L = as.numeric(as.character(df$avg_lifespan))
    
    if(dens.dep==TRUE){
      Faip <- Faip*(1-N/K) #only applicable if a==1
      Faip[Faip<0] <- 0 #no negative births allowed
      
      Fap1 <- Fap1*(1-N/K)
      Fap1[Fap1<0] <- 0
      
      #and all other reproduction is acted on less strongly - but we need a way to hold it steady
      Fap <- Fap*(1-(N/(K*2)))
      Fap[Fap<0] <- 0
    }
    
    #now build your pop matrix depending on your age at 1st rep
    #sketch the skeleton 
    pop.mat <- matrix(c(rep(0,length=(L*seas)),rep(0,length=(L*seas))), nrow=(L*seas),ncol=(L*seas))
    
    #add  the survivorship part
    pop.mat.s <- matrix(c(rep(0,length=(L*seas-1)),rep(0,length=(L*seas-1))), nrow=(L*seas-1),ncol=(L*seas-1))
    #we have to assume that infancy lasts for a full year in this case.
    #we fill the diagonal - first is always infant, then juvenile for 1, 2, or 3 slots (or more!), then the rest
    
    #we have to be careful here since our adult mortality varies by sim. 
    #for hunt_sim=2, it is every 2 of 4
    diag(pop.mat.s) <-(c(rep(pi, times=seas), rep(c(p_low_mort, p_hi_mort, p_low_mort, p_low_mort), times=(L-1)))[-1])
    
    #ends 2 short because you tack on the last adult later
    
    #make vectors for the top and side:
    #structure of top depends on a:
    if(a==1){
      top.row <- c(c(rep(0,length=(seas-1)), Faip), rep(c(rep(0,length=(seas-1)), Fap), times=(L-1)))
    }
    if(a==2){
      top.row <- c(c(rep(0,length=(seas*2-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-2)))
    }
    if(a>=3){
      top.row <- c(c(rep(0,length=(seas*a-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-a)))
    }
    
    final.col <- c(Fap,rep(0, length=L*seas-2), p_low_mort)
  }
  if (hunt_sim==3){
    adult_mort_hi <- (as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)/3
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    #fecundity also affected, depending on 
    Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #make into years again
    Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #first year breeders
    Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
    L = as.numeric(as.character(df$avg_lifespan))
    
    if(dens.dep==TRUE){
      Faip <- Faip*(1-N/K) #only applicable if a==1
      Faip[Faip<0] <- 0 #no negative births allowed
      
      Fap1 <- Fap1*(1-N/K)
      Fap1[Fap1<0] <- 0
      
      #and all other reproduction is acted on less strongly - but we need a way to hold it steady
      Fap <- Fap*(1-(N/(K*2)))
      Fap[Fap<0] <- 0
    }
    
    #now build your pop matrix depending on your age at 1st rep
    #sketch the skeleton 
    pop.mat <- matrix(c(rep(0,length=(L*seas)),rep(0,length=(L*seas))), nrow=(L*seas),ncol=(L*seas))
    
    #add  the survivorship part
    pop.mat.s <- matrix(c(rep(0,length=(L*seas-1)),rep(0,length=(L*seas-1))), nrow=(L*seas-1),ncol=(L*seas-1))
    #we have to assume that infancy lasts for a full year in this case.
    #we fill the diagonal - first is always infant, then juvenile for 1, 2, or 3 slots (or more!), then the rest
    
    #we have to be careful here since our adult mortality varies by sim. 
    #for hunt_sim=2, it is every 2 of 4
    diag(pop.mat.s) <-(c(rep(pi, times=seas), rep(c(p_low_mort, p_low_mort, p_hi_mort, p_low_mort), times=(L-1)))[-1])
    
    #ends 2 short because you tack on the last adult later
    
    #make vectors for the top and side:
    #structure of top depends on a:
    if(a==1){
      top.row <- c(c(rep(0,length=(seas-1)), Faip), rep(c(rep(0,length=(seas-1)), Fap), times=(L-1)))
    }
    if(a==2){
      top.row <- c(c(rep(0,length=(seas*2-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-2)))
    }
    if(a>=3){
      top.row <- c(c(rep(0,length=(seas*a-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-a)))
    }
    
    final.col <- c(Fap,rep(0, length=L*seas-2), p_low_mort)
  }
  if (hunt_sim==4){
    adult_mort_hi <- (as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)/3
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    #fecundity also affected, depending on 
    Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #make into years again
    Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #first year breeders
    Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
    L = as.numeric(as.character(df$avg_lifespan))
    
    if(dens.dep==TRUE){
      Faip <- Faip*(1-N/K) #only applicable if a==1
      Faip[Faip<0] <- 0 #no negative births allowed
      
      Fap1 <- Fap1*(1-N/K)
      Fap1[Fap1<0] <- 0
      
      #and all other reproduction is acted on less strongly - but we need a way to hold it steady
      Fap <- Fap*(1-(N/(K*2)))
      Fap[Fap<0] <- 0
    }
    
    #now build your pop matrix depending on your age at 1st rep
    #sketch the skeleton 
    pop.mat <- matrix(c(rep(0,length=(L*seas)),rep(0,length=(L*seas))), nrow=(L*seas),ncol=(L*seas))
    
    #add  the survivorship part
    pop.mat.s <- matrix(c(rep(0,length=(L*seas-1)),rep(0,length=(L*seas-1))), nrow=(L*seas-1),ncol=(L*seas-1))
    #we have to assume that infancy lasts for a full year in this case.
    #we fill the diagonal - first is always infant, then juvenile for 1, 2, or 3 slots (or more!), then the rest
    
    #we have to be careful here since our adult mortality varies by sim. 
    #for hunt_sim=2, it is every 2 of 4
    diag(pop.mat.s) <-(c(rep(pi, times=seas), rep(c(p_low_mort, p_low_mort, p_low_mort,  p_hi_mort), times=(L-1)))[-1])
    
    #ends 2 short because you tack on the last adult later
    
    #make vectors for the top and side:
    #structure of top depends on a:
    if(a==1){
      top.row <- c(c(rep(0,length=(seas-1)), Faip), rep(c(rep(0,length=(seas-1)), Fap), times=(L-1)))
    }
    if(a==2){
      top.row <- c(c(rep(0,length=(seas*2-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-2)))
    }
    if(a>=3){
      top.row <- c(c(rep(0,length=(seas*a-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-a)))
    }
    
    final.col <- c(Fap,rep(0, length=L*seas-2),  p_hi_mort)#this one has hi mort at end
  }
  if (hunt_sim==5){
    #this one puts hunt across two seasons
    adult_mort_hi <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9)/2
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)/2
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    #fecundity also affected, depending on 
    Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #make into years again
    Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #first year breeders
    Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
    L = as.numeric(as.character(df$avg_lifespan))
    
    if(dens.dep==TRUE){
      Faip <- Faip*(1-N/K) #only applicable if a==1
      Faip[Faip<0] <- 0 #no negative births allowed
      
      Fap1 <- Fap1*(1-N/K)
      Fap1[Fap1<0] <- 0
      
      #and all other reproduction is acted on less strongly - but we need a way to hold it steady
      Fap <- Fap*(1-(N/(K*2)))
      Fap[Fap<0] <- 0
    }
    
    #now build your pop matrix depending on your age at 1st rep
    #sketch the skeleton 
    pop.mat <- matrix(c(rep(0,length=(L*seas)),rep(0,length=(L*seas))), nrow=(L*seas),ncol=(L*seas))
    
    #add  the survivorship part
    pop.mat.s <- matrix(c(rep(0,length=(L*seas-1)),rep(0,length=(L*seas-1))), nrow=(L*seas-1),ncol=(L*seas-1))
    #we have to assume that infancy lasts for a full year in this case.
    #we fill the diagonal - first is always infant, then juvenile for 1, 2, or 3 slots (or more!), then the rest
    
    #we have to be careful here since our adult mortality varies by sim. 
    #for hunt_sim=2, it is every 2 of 4
    diag(pop.mat.s) <-(c(rep(pi, times=seas), rep(c(p_hi_mort, p_hi_mort, p_low_mort, p_low_mort), times=(L-1)))[-1])
    
    #ends 2 short because you tack on the last adult later
    
    #make vectors for the top and side:
    #structure of top depends on a:
    if(a==1){
      top.row <- c(c(rep(0,length=(seas-1)), Faip), rep(c(rep(0,length=(seas-1)), Fap), times=(L-1)))
    }
    if(a==2){
      top.row <- c(c(rep(0,length=(seas*2-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-2)))
    }
    if(a>=3){
      top.row <- c(c(rep(0,length=(seas*a-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-a)))
    }
    
    final.col <- c(Fap,rep(0, length=L*seas-2),  p_low_mort)
  }
  if (hunt_sim==6){
    #this one puts hunt across two seasons
    adult_mort_hi <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9)/2
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)/2
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    #fecundity also affected, depending on 
    Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #make into years again
    Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #first year breeders
    Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
    L = as.numeric(as.character(df$avg_lifespan))
    
    if(dens.dep==TRUE){
      Faip <- Faip*(1-N/K) #only applicable if a==1
      Faip[Faip<0] <- 0 #no negative births allowed
      
      Fap1 <- Fap1*(1-N/K)
      Fap1[Fap1<0] <- 0
      
      #and all other reproduction is acted on less strongly - but we need a way to hold it steady
      Fap <- Fap*(1-(N/(K*2)))
      Fap[Fap<0] <- 0
    }
    
    #now build your pop matrix depending on your age at 1st rep
    #sketch the skeleton 
    pop.mat <- matrix(c(rep(0,length=(L*seas)),rep(0,length=(L*seas))), nrow=(L*seas),ncol=(L*seas))
    
    #add  the survivorship part
    pop.mat.s <- matrix(c(rep(0,length=(L*seas-1)),rep(0,length=(L*seas-1))), nrow=(L*seas-1),ncol=(L*seas-1))
    #we have to assume that infancy lasts for a full year in this case.
    #we fill the diagonal - first is always infant, then juvenile for 1, 2, or 3 slots (or more!), then the rest
    
    #we have to be careful here since our adult mortality varies by sim. 
    #for hunt_sim=2, it is every 2 of 4
    diag(pop.mat.s) <-(c(rep(pi, times=seas), rep(c(p_low_mort, p_hi_mort, p_hi_mort, p_low_mort), times=(L-1)))[-1])
    
    #ends 2 short because you tack on the last adult later
    
    #make vectors for the top and side:
    #structure of top depends on a:
    if(a==1){
      top.row <- c(c(rep(0,length=(seas-1)), Faip), rep(c(rep(0,length=(seas-1)), Fap), times=(L-1)))
    }
    if(a==2){
      top.row <- c(c(rep(0,length=(seas*2-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-2)))
    }
    if(a>=3){
      top.row <- c(c(rep(0,length=(seas*a-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-a)))
    }
    
    final.col <- c(Fap,rep(0, length=L*seas-2),  p_low_mort)
  }
  if (hunt_sim==7){
    #this one puts hunt across two seasons
    adult_mort_hi <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9)/2
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)/2
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    #fecundity also affected, depending on 
    Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #make into years again
    Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #first year breeders
    Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
    L = as.numeric(as.character(df$avg_lifespan))
    
    if(dens.dep==TRUE){
      Faip <- Faip*(1-N/K) #only applicable if a==1
      Faip[Faip<0] <- 0 #no negative births allowed
      
      Fap1 <- Fap1*(1-N/K)
      Fap1[Fap1<0] <- 0
      
      #and all other reproduction is acted on less strongly - but we need a way to hold it steady
      Fap <- Fap*(1-(N/(K*2)))
      Fap[Fap<0] <- 0
    }
    
    #now build your pop matrix depending on your age at 1st rep
    #sketch the skeleton 
    pop.mat <- matrix(c(rep(0,length=(L*seas)),rep(0,length=(L*seas))), nrow=(L*seas),ncol=(L*seas))
    
    #add  the survivorship part
    pop.mat.s <- matrix(c(rep(0,length=(L*seas-1)),rep(0,length=(L*seas-1))), nrow=(L*seas-1),ncol=(L*seas-1))
    #we have to assume that infancy lasts for a full year in this case.
    #we fill the diagonal - first is always infant, then juvenile for 1, 2, or 3 slots (or more!), then the rest
    
    #we have to be careful here since our adult mortality varies by sim. 
    #for hunt_sim=2, it is every 2 of 4
    diag(pop.mat.s) <-(c(rep(pi, times=seas), rep(c(p_low_mort, p_low_mort, p_hi_mort, p_hi_mort), times=(L-1)))[-1])
    
    #ends 2 short because you tack on the last adult later
    
    #make vectors for the top and side:
    #structure of top depends on a:
    if(a==1){
      top.row <- c(c(rep(0,length=(seas-1)), Faip), rep(c(rep(0,length=(seas-1)), Fap), times=(L-1)))
    }
    if(a==2){
      top.row <- c(c(rep(0,length=(seas*2-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-2)))
    }
    if(a>=3){
      top.row <- c(c(rep(0,length=(seas*a-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-a)))
    }
    
    final.col <- c(Fap,rep(0, length=L*seas-2),  p_hi_mort)
  }
  if (hunt_sim==8){
    #this one puts hunt across two seasons
    adult_mort_hi <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9)/2
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)/2
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    #fecundity also affected, depending on 
    Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #make into years again
    Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #first year breeders
    Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
    L = as.numeric(as.character(df$avg_lifespan))
    
    if(dens.dep==TRUE){
      Faip <- Faip*(1-N/K) #only applicable if a==1
      Faip[Faip<0] <- 0 #no negative births allowed
      
      Fap1 <- Fap1*(1-N/K)
      Fap1[Fap1<0] <- 0
      
      #and all other reproduction is acted on less strongly - but we need a way to hold it steady
      Fap <- Fap*(1-(N/(K*2)))
      Fap[Fap<0] <- 0
    }
    
    #now build your pop matrix depending on your age at 1st rep
    #sketch the skeleton 
    pop.mat <- matrix(c(rep(0,length=(L*seas)),rep(0,length=(L*seas))), nrow=(L*seas),ncol=(L*seas))
    
    #add  the survivorship part
    pop.mat.s <- matrix(c(rep(0,length=(L*seas-1)),rep(0,length=(L*seas-1))), nrow=(L*seas-1),ncol=(L*seas-1))
    #we have to assume that infancy lasts for a full year in this case.
    #we fill the diagonal - first is always infant, then juvenile for 1, 2, or 3 slots (or more!), then the rest
    
    #we have to be careful here since our adult mortality varies by sim. 
    #for hunt_sim=2, it is every 2 of 4
    diag(pop.mat.s) <-(c(rep(pi, times=seas), rep(c(p_hi_mort, p_low_mort, p_low_mort, p_hi_mort), times=(L-1)))[-1])
    
    #ends 2 short because you tack on the last adult later
    
    #make vectors for the top and side:
    #structure of top depends on a:
    if(a==1){
      top.row <- c(c(rep(0,length=(seas-1)), Faip), rep(c(rep(0,length=(seas-1)), Fap), times=(L-1)))
    }
    if(a==2){
      top.row <- c(c(rep(0,length=(seas*2-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-2)))
    }
    if(a>=3){
      top.row <- c(c(rep(0,length=(seas*a-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-a)))
    }
    
    final.col <- c(Fap,rep(0, length=L*seas-2),  p_hi_mort)
  }
  if (hunt_sim==9){
    #this one puts hunt across three seasons
    adult_mort_hi <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9)/3
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    #fecundity also affected, depending on 
    Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #make into years again
    Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #first year breeders
    Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
    L = as.numeric(as.character(df$avg_lifespan))
    
    if(dens.dep==TRUE){
      Faip <- Faip*(1-N/K) #only applicable if a==1
      Faip[Faip<0] <- 0 #no negative births allowed
      
      Fap1 <- Fap1*(1-N/K)
      Fap1[Fap1<0] <- 0
      
      #and all other reproduction is acted on less strongly - but we need a way to hold it steady
      Fap <- Fap*(1-(N/(K*2)))
      Fap[Fap<0] <- 0
    }
    
    #now build your pop matrix depending on your age at 1st rep
    #sketch the skeleton 
    pop.mat <- matrix(c(rep(0,length=(L*seas)),rep(0,length=(L*seas))), nrow=(L*seas),ncol=(L*seas))
    
    #add  the survivorship part
    pop.mat.s <- matrix(c(rep(0,length=(L*seas-1)),rep(0,length=(L*seas-1))), nrow=(L*seas-1),ncol=(L*seas-1))
    #we have to assume that infancy lasts for a full year in this case.
    #we fill the diagonal - first is always infant, then juvenile for 1, 2, or 3 slots (or more!), then the rest
    
    #we have to be careful here since our adult mortality varies by sim. 
    #for hunt_sim=2, it is every 2 of 4
    diag(pop.mat.s) <-(c(rep(pi, times=seas), rep(c(p_hi_mort, p_hi_mort, p_hi_mort, p_low_mort), times=(L-1)))[-1])
    
    #ends 2 short because you tack on the last adult later
    
    #make vectors for the top and side:
    #structure of top depends on a:
    if(a==1){
      top.row <- c(c(rep(0,length=(seas-1)), Faip), rep(c(rep(0,length=(seas-1)), Fap), times=(L-1)))
    }
    if(a==2){
      top.row <- c(c(rep(0,length=(seas*2-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-2)))
    }
    if(a>=3){
      top.row <- c(c(rep(0,length=(seas*a-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-a)))
    }
    
    final.col <- c(Fap,rep(0, length=L*seas-2),  p_low_mort)
  }
  if (hunt_sim==10){
    #this one puts hunt across three seasons
    adult_mort_hi <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9)/3
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    #fecundity also affected, depending on 
    Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #make into years again
    Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #first year breeders
    Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
    L = as.numeric(as.character(df$avg_lifespan))
    
    if(dens.dep==TRUE){
      Faip <- Faip*(1-N/K) #only applicable if a==1
      Faip[Faip<0] <- 0 #no negative births allowed
      
      Fap1 <- Fap1*(1-N/K)
      Fap1[Fap1<0] <- 0
      
      #and all other reproduction is acted on less strongly - but we need a way to hold it steady
      Fap <- Fap*(1-(N/(K*2)))
      Fap[Fap<0] <- 0
    }
    
    #now build your pop matrix depending on your age at 1st rep
    #sketch the skeleton 
    pop.mat <- matrix(c(rep(0,length=(L*seas)),rep(0,length=(L*seas))), nrow=(L*seas),ncol=(L*seas))
    
    #add  the survivorship part
    pop.mat.s <- matrix(c(rep(0,length=(L*seas-1)),rep(0,length=(L*seas-1))), nrow=(L*seas-1),ncol=(L*seas-1))
    #we have to assume that infancy lasts for a full year in this case.
    #we fill the diagonal - first is always infant, then juvenile for 1, 2, or 3 slots (or more!), then the rest
    
    #we have to be careful here since our adult mortality varies by sim. 
    #for hunt_sim=2, it is every 2 of 4
    diag(pop.mat.s) <-(c(rep(pi, times=seas), rep(c(p_low_mort, p_hi_mort, p_hi_mort, p_hi_mort), times=(L-1)))[-1])
    
    #ends 2 short because you tack on the last adult later
    
    #make vectors for the top and side:
    #structure of top depends on a:
    if(a==1){
      top.row <- c(c(rep(0,length=(seas-1)), Faip), rep(c(rep(0,length=(seas-1)), Fap), times=(L-1)))
    }
    if(a==2){
      top.row <- c(c(rep(0,length=(seas*2-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-2)))
    }
    if(a>=3){
      top.row <- c(c(rep(0,length=(seas*a-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-a)))
    }
    
    final.col <- c(Fap,rep(0, length=L*seas-2),  p_hi_mort)
  }
  if (hunt_sim==11){
    #this one puts hunt across three seasons
    adult_mort_hi <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9)/3
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    #fecundity also affected, depending on 
    Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #make into years again
    Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #first year breeders
    Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
    L = as.numeric(as.character(df$avg_lifespan))
    
    if(dens.dep==TRUE){
      Faip <- Faip*(1-N/K) #only applicable if a==1
      Faip[Faip<0] <- 0 #no negative births allowed
      
      Fap1 <- Fap1*(1-N/K)
      Fap1[Fap1<0] <- 0
      
      #and all other reproduction is acted on less strongly - but we need a way to hold it steady
      Fap <- Fap*(1-(N/(K*2)))
      Fap[Fap<0] <- 0
    }
    
    #now build your pop matrix depending on your age at 1st rep
    #sketch the skeleton 
    pop.mat <- matrix(c(rep(0,length=(L*seas)),rep(0,length=(L*seas))), nrow=(L*seas),ncol=(L*seas))
    
    #add  the survivorship part
    pop.mat.s <- matrix(c(rep(0,length=(L*seas-1)),rep(0,length=(L*seas-1))), nrow=(L*seas-1),ncol=(L*seas-1))
    #we have to assume that infancy lasts for a full year in this case.
    #we fill the diagonal - first is always infant, then juvenile for 1, 2, or 3 slots (or more!), then the rest
    
    #we have to be careful here since our adult mortality varies by sim. 
    #for hunt_sim=2, it is every 2 of 4
    diag(pop.mat.s) <-(c(rep(pi, times=seas), rep(c(p_hi_mort, p_low_mort, p_hi_mort, p_hi_mort), times=(L-1)))[-1])
    
    #ends 2 short because you tack on the last adult later
    
    #make vectors for the top and side:
    #structure of top depends on a:
    if(a==1){
      top.row <- c(c(rep(0,length=(seas-1)), Faip), rep(c(rep(0,length=(seas-1)), Fap), times=(L-1)))
    }
    if(a==2){
      top.row <- c(c(rep(0,length=(seas*2-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-2)))
    }
    if(a>=3){
      top.row <- c(c(rep(0,length=(seas*a-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-a)))
    }
    
    final.col <- c(Fap,rep(0, length=L*seas-2),  p_hi_mort)
  }
  if (hunt_sim==12){
    #this one puts hunt across three seasons
    adult_mort_hi <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9)/3
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    #fecundity also affected, depending on 
    Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #make into years again
    Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #first year breeders
    Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
    L = as.numeric(as.character(df$avg_lifespan))
    
    if(dens.dep==TRUE){
      Faip <- Faip*(1-N/K) #only applicable if a==1
      Faip[Faip<0] <- 0 #no negative births allowed
      
      Fap1 <- Fap1*(1-N/K)
      Fap1[Fap1<0] <- 0
      
      #and all other reproduction is acted on less strongly - but we need a way to hold it steady
      Fap <- Fap*(1-(N/(K*2)))
      Fap[Fap<0] <- 0
    }
    
    #now build your pop matrix depending on your age at 1st rep
    #sketch the skeleton 
    pop.mat <- matrix(c(rep(0,length=(L*seas)),rep(0,length=(L*seas))), nrow=(L*seas),ncol=(L*seas))
    
    #add  the survivorship part
    pop.mat.s <- matrix(c(rep(0,length=(L*seas-1)),rep(0,length=(L*seas-1))), nrow=(L*seas-1),ncol=(L*seas-1))
    #we have to assume that infancy lasts for a full year in this case.
    #we fill the diagonal - first is always infant, then juvenile for 1, 2, or 3 slots (or more!), then the rest
    
    #we have to be careful here since our adult mortality varies by sim. 
    #for hunt_sim=2, it is every 2 of 4
    diag(pop.mat.s) <-(c(rep(pi, times=seas), rep(c(p_hi_mort, p_hi_mort, p_low_mort, p_hi_mort), times=(L-1)))[-1])
    
    #ends 2 short because you tack on the last adult later
    
    #make vectors for the top and side:
    #structure of top depends on a:
    if(a==1){
      top.row <- c(c(rep(0,length=(seas-1)), Faip), rep(c(rep(0,length=(seas-1)), Fap), times=(L-1)))
    }
    if(a==2){
      top.row <- c(c(rep(0,length=(seas*2-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-2)))
    }
    if(a>=3){
      top.row <- c(c(rep(0,length=(seas*a-1)), Fap1), rep(c(rep(0,length=(seas-1)), Fap), times=(L-a)))
    }
    
    final.col <- c(Fap,rep(0, length=L*seas-2),  p_hi_mort)
  }
  
  
  #and fill the skeleton:
  pop.mat[1,] <- top.row
  pop.mat[,ncol(pop.mat)] <- final.col
  pop.mat[2:nrow(pop.mat),1:(ncol(pop.mat)-1)] <- pop.mat.s
  
  #and return
  
  return(pop.mat)
}
build.leslie.periodic <- function(df, nat.mort, hunt.override, set_hunt, use.uci,  dens.dep, N, K){
  #add in the option to use the upper confidence limit of the mortality rate
  if(use.uci==TRUE){
    df$IBI_adult_mortality = df$IBI_upr
  }
  
  if(hunt.override==TRUE){
    df$IBI_adult_mortality = set_hunt
  }else{
    #add in the option to elevate mortality beyond that captured in the hunt
    df$IBI_adult_mortality =  df$IBI_adult_mortality + nat.mort
    #remember that we can't have mortality above 100%, so correct for that here
    df$IBI_adult_mortality[df$IBI_adult_mortality>=1] <- .9
    df$IBI_adult_mortality[df$IBI_adult_mortality<0] <- 0
  }
  #now convert back to annual and then make seasonal
  adult_mort <- (as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))/4
  infant_mort <- (as.numeric(df$IBI_infant_mortality)/as.numeric(df$IBI))/4
  
  #and build your respective matrices for each season
  
  p = 1-adult_mort
  pi= 1-infant_mort
  a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
  Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p #make into years again
  Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p #first year breeders
  Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
  L = as.numeric(as.character(df$avg_lifespan))
  
  if(dens.dep==TRUE){
    Faip <- Faip*(1-N/K) #only applicable if a==1
    Faip[Faip<0] <- 0 #no negative births allowed
    
    Fap1 <- Fap1*(1-N/K)
    Fap1[Fap1<0] <- 0
    
    #and all other reproduction is acted on less strongly - but we need a way to hold it steady
    Fap <- Fap*(1-(N/(K*2)))
    Fap[Fap<0] <- 0
  }
  
  #build seasonal survival matrices for each quarter of the year
  qtr1 = matrix(0,nrow=(L), ncol=(L))
  diag(qtr1) = c(pi, rep(p, times= (length(diag(qtr1))-1))) #add in sites
  #qtr1[nrow(qtr1), ncol(qtr1)] = p
  
  qtr2 = qtr3 = qtr1  #could alter later if you want to change the survival rate for infants halfway through
  
  #and qtr4 (with reproduction)
  qtr4 = matrix(0,nrow=(L-1), ncol=(L))
  diag(qtr4) = c(pi, rep(p, times= (length(diag(qtr4))-1))) #add in sites
  qtr4[nrow(qtr4), ncol(qtr4)] = p
  
  
  if(a==1){
    top.row <- c(Faip, rep(Fap, length=(L-1)))
  }
  if(a==2){
    top.row <- c(0,Fap1, rep(Fap, length=(L-2)))
  }
  if(a>=3){
    top.row <- c(rep(0, length=(a-1)),Fap1,rep(Fap, length=(L-a)))
  }
  
  
  
  qtr4 = rbind(top.row, qtr4) #survives 4th to next year
  
  #and get the total population level matrix as the product of these quarters to test for lambda
  
  
  #and finally the birth rate for the 4th quarter (where rep takes place)
  #all but those which do not reproduce have a fec. term
  
  #pop.mat =  qtr4 %*%  qtr1 %*%  qtr2 %*% qtr3
  
  pop.mat =   qtr1 %*%  qtr2 %*% qtr3 %*% qtr4 
  
  
  
  #and return
  
  return(pop.mat)
}
build.leslie.periodic.hunt <- function(df, hunt.override, set_hunt, nat.mort, use.uci,  dens.dep, N, K, hunt_sim){
  #add in the option to use the upper confidence limit of the mortality rate
  if(use.uci==TRUE){
    df$IBI_adult_mortality = df$IBI_upr
  }
  #add in the option to elevate mortality beyond that captured in the hunt
  if(hunt.override==TRUE){
    df$IBI_adult_mortality = set_hunt
  }else{
    
    
    df$IBI_adult_mortality =  df$IBI_adult_mortality + nat.mort
    #remember that we can't have mortality above 100%, so correct for that here
    df$IBI_adult_mortality[df$IBI_adult_mortality>=1] <- .9
    df$IBI_adult_mortality[df$IBI_adult_mortality<0] <- 0
  }
  
  if (hunt_sim==1 | hunt_sim==2 | hunt_sim==3 | hunt_sim==4){
    #and build your respective matrics for each season
    
    #these tabulations should hold for all 1 quarter hunts
    infant_mort <- (as.numeric(df$IBI_infant_mortality)/as.numeric(df$IBI))/4
    adult_mort_hi <- (as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)/3
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    
    if (hunt_sim==1){
      #fecundity also affected, depending on the season of hunting (when hunting is in the 4th quarter)
      Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #make into years again
      Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #first year breeders
      Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
      L = as.numeric(as.character(df$avg_lifespan))
      
      if(dens.dep==TRUE){
        Faip <- Faip*(1-N/K) #only applicable if a==1
        Faip[Faip<0] <- 0 #no negative births allowed
        
        Fap1 <- Fap1*(1-N/K)
        Fap1[Fap1<0] <- 0
        
        #and all other reproduction is acted on less strongly - but we need a way to hold it steady
        Fap <- Fap*(1-(N/(K*2)))
        Fap[Fap<0] <- 0
      }
      
      #build seasonal survival matrices for each quarter - this varies by hunt sim.
      
      qtr1 = matrix(0,nrow=L, ncol=L) 
      qtr2 = qtr1
      diag(qtr1) = c(pi, rep(p_hi_mort, times= (length(diag(qtr1))-1))) #hi hunt is concentrated in 1st quarter
      diag(qtr2) = c(pi, rep(p_low_mort, times= (length(diag(qtr2))-1))) #hi hunt is concentrated in 1st quarter
      #qtr1[nrow(qtr1), ncol(qtr1)] = p_hi_mort
      #qtr2[nrow(qtr2), ncol(qtr2)] = p_low_mort
      
      qtr3 = qtr2  #could alter later if you want to change the survival rate for infants halfway through
      
      #and add the 0s for fecundity to these seasons
      #qtr2 = rbind(rep(0, ncol(qtr2)), qtr2) #survives 3rd to 4th
      #qtr3 = rbind(rep(0, ncol(qtr3)), qtr3) #survives 2nd to 3rd
      #qtr1 = rbind(rep(0, ncol(qtr1)), qtr1) #survives 1st to 2nd
      
      
      if(a==1){
        top.row <- c(Faip, rep(Fap, length=(L-1)))
      }
      if(a==2){
        top.row <- c(0,Fap1, rep(Fap, length=(L-2)))
      }
      if(a>=3){
        top.row <- c(rep(0, length=(a-1)),Fap1,rep(Fap, length=(L-a)))
      }
      
      
      qtr4 = matrix(0,nrow=(L-1), ncol=(L))
      diag(qtr4) = c(pi, rep(p_low_mort, times= (length(diag(qtr4))-1))) #add in sites
      qtr4[nrow(qtr4), ncol(qtr4)] = p_low_mort
      
      
      qtr4 = rbind(top.row, qtr4) #survives 4th to next year
      
      #and get the total population level matrix as the product of these quarters to test for lambda
      
      
      #and finally the birth rate for the 4th quarter (where rep takes place)
      #all but those which do not reproduce have a fec. term
      
      
      pop.mat =   qtr1 %*%  qtr2 %*% qtr3 %*% qtr4
      
      
    }else if (hunt_sim==2){ # hi hunt in 2nd quarter
      #fecundity also affected, depending on the season of hunting (when hunting is in the 4th quarter)
      Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #make into years again
      Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #first year breeders
      Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
      L = as.numeric(as.character(df$avg_lifespan))
      
      if(dens.dep==TRUE){
        Faip <- Faip*(1-N/K) #only applicable if a==1
        Faip[Faip<0] <- 0 #no negative births allowed
        
        Fap1 <- Fap1*(1-N/K)
        Fap1[Fap1<0] <- 0
        
        #and all other reproduction is acted on less strongly - but we need a way to hold it steady
        Fap <- Fap*(1-(N/(K*2)))
        Fap[Fap<0] <- 0
      }
      
      #build seasonal survival matrices for each quarter - this varies by hunt sim.
      
      qtr1 = matrix(0,nrow=L, ncol=L) 
      qtr2 = qtr1
      diag(qtr1) = c(pi, rep(p_low_mort, times= (length(diag(qtr1))-1))) #hi hunt is concentrated in 1st quarter
      diag(qtr2) = c(pi, rep(p_hi_mort, times= (length(diag(qtr2))-1))) #hi hunt is concentrated in 1st quarter
      #qtr1[nrow(qtr1), ncol(qtr1)] = p_hi_mort
      #qtr2[nrow(qtr2), ncol(qtr2)] = p_low_mort
      
      qtr3 = qtr1  #could alter later if you want to change the survival rate for infants halfway through
      
      
      if(a==1){
        top.row <- c(Faip, rep(Fap, length=(L-1)))
      }
      if(a==2){
        top.row <- c(0,Fap1, rep(Fap, length=(L-2)))
      }
      if(a>=3){
        top.row <- c(rep(0, length=(a-1)),Fap1,rep(Fap, length=(L-a)))
      }
      
      
      qtr4 = matrix(0,nrow=(L-1), ncol=(L))
      diag(qtr4) = c(pi, rep(p_low_mort, times= (length(diag(qtr4))-1))) #add in sites
      qtr4[nrow(qtr4), ncol(qtr4)] = p_low_mort
      
      
      qtr4 = rbind(top.row, qtr4) #survives 4th to next year
      
      #and get the total population level matrix as the product of these quarters to test for lambda
      
      
      #and finally the birth rate for the 4th quarter (where rep takes place)
      #all but those which do not reproduce have a fec. term
      
      
      pop.mat =   qtr1 %*%  qtr2 %*% qtr3 %*% qtr4
      
      
    }else if (hunt_sim==3){
      
      #fecundity also affected, depending on the season of hunting (when hunting is in the 4th quarter)
      Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #make into years again
      Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #first year breeders
      Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
      L = as.numeric(as.character(df$avg_lifespan))
      
      if(dens.dep==TRUE){
        Faip <- Faip*(1-N/K) #only applicable if a==1
        Faip[Faip<0] <- 0 #no negative births allowed
        
        Fap1 <- Fap1*(1-N/K)
        Fap1[Fap1<0] <- 0
        
        #and all other reproduction is acted on less strongly - but we need a way to hold it steady
        Fap <- Fap*(1-(N/(K*2)))
        Fap[Fap<0] <- 0
      }
      
      #build seasonal survival matrices for each quarter - this varies by hunt sim.
      
      qtr1 = matrix(0,nrow=L, ncol=L) 
      qtr3 = qtr1
      diag(qtr1) = c(pi, rep(p_low_mort, times= (length(diag(qtr1))-1))) #hi hunt is concentrated in 1st quarter
      diag(qtr3) = c(pi, rep(p_hi_mort, times= (length(diag(qtr3))-1))) #hi hunt is concentrated in 1st quarter
      #qtr1[nrow(qtr1), ncol(qtr1)] = p_hi_mort
      #qtr2[nrow(qtr2), ncol(qtr2)] = p_low_mort
      
      qtr2 = qtr1  #could alter later if you want to change the survival rate for infants halfway through
      
      
      if(a==1){
        top.row <- c(Faip, rep(Fap, length=(L-1)))
      }
      if(a==2){
        top.row <- c(0,Fap1, rep(Fap, length=(L-2)))
      }
      if(a>=3){
        top.row <- c(rep(0, length=(a-1)),Fap1,rep(Fap, length=(L-a)))
      }
      
      
      qtr4 = matrix(0,nrow=(L-1), ncol=(L))
      diag(qtr4) = c(pi, rep(p_low_mort, times= (length(diag(qtr4))-1))) #add in sites
      qtr4[nrow(qtr4), ncol(qtr4)] = p_low_mort
      
      
      qtr4 = rbind(top.row, qtr4) #survives 4th to next year
      
      #and get the total population level matrix as the product of these quarters to test for lambda
      
      
      #and finally the birth rate for the 4th quarter (where rep takes place)
      #all but those which do not reproduce have a fec. term
      
      
      pop.mat =   qtr1 %*%  qtr2 %*% qtr3 %*% qtr4
      
      
      
      
    } else if (hunt_sim==4){
      
      #fecundity also affected, depending on the season of hunting (when hunting is in the 4th quarter)
      Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #make into years again
      Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #first year breeders
      Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
      L = as.numeric(as.character(df$avg_lifespan))
      
      if(dens.dep==TRUE){
        Faip <- Faip*(1-N/K) #only applicable if a==1
        Faip[Faip<0] <- 0 #no negative births allowed
        
        Fap1 <- Fap1*(1-N/K)
        Fap1[Fap1<0] <- 0
        
        #and all other reproduction is acted on less strongly - but we need a way to hold it steady
        Fap <- Fap*(1-(N/(K*2)))
        Fap[Fap<0] <- 0
      }
      
      #build seasonal survival matrices for each quarter - this varies by hunt sim.
      
      qtr1 = matrix(0,nrow=L, ncol=L) 
      #qtr2 = qtr1
      diag(qtr1) = c(pi, rep(p_low_mort, times= (length(diag(qtr1))-1))) #hi hunt is concentrated in 1st quarter
      #diag(qtr2) = c(pi, rep(p_hi_mort, times= (length(diag(qtr2))-1))) #hi hunt is concentrated in 1st quarter
      #qtr1[nrow(qtr1), ncol(qtr1)] = p_hi_mort
      #qtr2[nrow(qtr2), ncol(qtr2)] = p_low_mort
      
      qtr2 = qtr3 = qtr1  #could alter later if you want to change the survival rate for infants halfway through
      
      
      if(a==1){
        top.row <- c(Faip, rep(Fap, length=(L-1)))
      }
      if(a==2){
        top.row <- c(0,Fap1, rep(Fap, length=(L-2)))
      }
      if(a>=3){
        top.row <- c(rep(0, length=(a-1)),Fap1,rep(Fap, length=(L-a)))
      }
      
      
      qtr4 = matrix(0,nrow=(L-1), ncol=(L))
      diag(qtr4) = c(pi, rep(p_hi_mort, times= (length(diag(qtr4))-1))) #add in sites
      qtr4[nrow(qtr4), ncol(qtr4)] = p_hi_mort
      
      
      qtr4 = rbind(top.row, qtr4) #survives 4th to next year
      
      #and get the total population level matrix as the product of these quarters to test for lambda
      
      
      #and finally the birth rate for the 4th quarter (where rep takes place)
      #all but those which do not reproduce have a fec. term
      
      
      pop.mat =   qtr1 %*%  qtr2 %*% qtr3 %*% qtr4
      
      
    }
  }else if (hunt_sim==5 | hunt_sim==6 | hunt_sim==7 | hunt_sim==8){ # now try the scenarios where the hunting is spread over two seasons
    #this one puts hunt across two seasons
    infant_mort <- (as.numeric(df$IBI_infant_mortality)/as.numeric(df$IBI))/4
    adult_mort_hi <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9)/2
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)/2
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    
    
    if(hunt_sim==5){
      #fecundity also affected, depending on the season of hunting (when hunting is in the 4th quarter)
      Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #make into years again
      Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #first year breeders
      Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
      L = as.numeric(as.character(df$avg_lifespan))
      
      if(dens.dep==TRUE){
        Faip <- Faip*(1-N/K) #only applicable if a==1
        Faip[Faip<0] <- 0 #no negative births allowed
        
        Fap1 <- Fap1*(1-N/K)
        Fap1[Fap1<0] <- 0
        
        #and all other reproduction is acted on less strongly - but we need a way to hold it steady
        Fap <- Fap*(1-(N/(K*2)))
        Fap[Fap<0] <- 0
      }
      
      #build seasonal survival matrices for each quarter - this varies by hunt sim.
      
      qtr1 = matrix(0,nrow=L, ncol=L) 
      qtr3 = qtr1
      diag(qtr1) = c(pi, rep(p_hi_mort, times= (length(diag(qtr1))-1))) #hi hunt is concentrated in 1st quarter
      diag(qtr3) = c(pi, rep(p_low_mort, times= (length(diag(qtr3))-1))) #hi hunt is concentrated in 1st quarter
      #qtr1[nrow(qtr1), ncol(qtr1)] = p_hi_mort
      #qtr2[nrow(qtr2), ncol(qtr2)] = p_low_mort
      
      qtr2 = qtr1  #could alter later if you want to change the survival rate for infants halfway through
      
      
      if(a==1){
        top.row <- c(Faip, rep(Fap, length=(L-1)))
      }
      if(a==2){
        top.row <- c(0,Fap1, rep(Fap, length=(L-2)))
      }
      if(a>=3){
        top.row <- c(rep(0, length=(a-1)),Fap1,rep(Fap, length=(L-a)))
      }
      
      
      qtr4 = matrix(0,nrow=(L-1), ncol=(L))
      diag(qtr4) = c(pi, rep(p_low_mort, times= (length(diag(qtr4))-1))) #add in sites
      qtr4[nrow(qtr4), ncol(qtr4)] = p_low_mort
      
      
      qtr4 = rbind(top.row, qtr4) #survives 4th to next year
      
      
      pop.mat =   qtr1 %*%  qtr2 %*% qtr3 %*% qtr4
      
      
    } else if (hunt_sim==6){
      
      Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #make into years again
      Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #first year breeders
      Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
      L = as.numeric(as.character(df$avg_lifespan))
      
      if(dens.dep==TRUE){
        Faip <- Faip*(1-N/K) #only applicable if a==1
        Faip[Faip<0] <- 0 #no negative births allowed
        
        Fap1 <- Fap1*(1-N/K)
        Fap1[Fap1<0] <- 0
        
        #and all other reproduction is acted on less strongly - but we need a way to hold it steady
        Fap <- Fap*(1-(N/(K*2)))
        Fap[Fap<0] <- 0
      }
      
      #build seasonal survival matrices for each quarter - this varies by hunt sim.
      
      qtr1 = matrix(0,nrow=L, ncol=L) 
      qtr2 = qtr1
      diag(qtr2) = c(pi, rep(p_hi_mort, times= (length(diag(qtr2))-1))) #hi hunt is concentrated in 1st quarter
      diag(qtr1) = c(pi, rep(p_low_mort, times= (length(diag(qtr1))-1))) #hi hunt is concentrated in 1st quarter
      #qtr1[nrow(qtr1), ncol(qtr1)] = p_hi_mort
      #qtr2[nrow(qtr2), ncol(qtr2)] = p_low_mort
      
      qtr3 = qtr2  #could alter later if you want to change the survival rate for infants halfway through
      
      
      if(a==1){
        top.row <- c(Faip, rep(Fap, length=(L-1)))
      }
      if(a==2){
        top.row <- c(0,Fap1, rep(Fap, length=(L-2)))
      }
      if(a>=3){
        top.row <- c(rep(0, length=(a-1)),Fap1,rep(Fap, length=(L-a)))
      }
      
      
      qtr4 = matrix(0,nrow=(L-1), ncol=(L))
      diag(qtr4) = c(pi, rep(p_low_mort, times= (length(diag(qtr4))-1))) #add in sites
      qtr4[nrow(qtr4), ncol(qtr4)] = p_low_mort
      
      
      qtr4 = rbind(top.row, qtr4) #survives 4th to next year
      
      
      pop.mat =   qtr1 %*%  qtr2 %*% qtr3 %*% qtr4
      
    } else if (hunt_sim==7){
      Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #make into years again
      Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #first year breeders
      Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
      L = as.numeric(as.character(df$avg_lifespan))
      
      if(dens.dep==TRUE){
        Faip <- Faip*(1-N/K) #only applicable if a==1
        Faip[Faip<0] <- 0 #no negative births allowed
        
        Fap1 <- Fap1*(1-N/K)
        Fap1[Fap1<0] <- 0
        
        #and all other reproduction is acted on less strongly - but we need a way to hold it steady
        Fap <- Fap*(1-(N/(K*2)))
        Fap[Fap<0] <- 0
      }
      
      #build seasonal survival matrices for each quarter - this varies by hunt sim.
      
      qtr1 = matrix(0,nrow=L, ncol=L) 
      qtr3 = qtr1
      diag(qtr1) = c(pi, rep(p_low_mort, times= (length(diag(qtr1))-1))) #hi hunt is concentrated in 1st quarter
      diag(qtr3) = c(pi, rep(p_hi_mort, times= (length(diag(qtr3))-1))) #hi hunt is concentrated in 1st quarter
      
      
      qtr2 = qtr1  #could alter later if you want to change the survival rate for infants halfway through
      
      
      if(a==1){
        top.row <- c(Faip, rep(Fap, length=(L-1)))
      }
      if(a==2){
        top.row <- c(0,Fap1, rep(Fap, length=(L-2)))
      }
      if(a>=3){
        top.row <- c(rep(0, length=(a-1)),Fap1,rep(Fap, length=(L-a)))
      }
      
      
      qtr4 = matrix(0,nrow=(L-1), ncol=(L))
      diag(qtr4) = c(pi, rep(p_hi_mort, times= (length(diag(qtr4))-1))) #add in sites
      qtr4[nrow(qtr4), ncol(qtr4)] = p_hi_mort
      
      
      qtr4 = rbind(top.row, qtr4) #survives 4th to next year
      
      
      pop.mat =   qtr1 %*%  qtr2 %*% qtr3 %*% qtr4
      
    } else if (hunt_sim==8){
      Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #make into years again
      Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #first year breeders
      Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
      L = as.numeric(as.character(df$avg_lifespan))
      
      if(dens.dep==TRUE){
        Faip <- Faip*(1-N/K) #only applicable if a==1
        Faip[Faip<0] <- 0 #no negative births allowed
        
        Fap1 <- Fap1*(1-N/K)
        Fap1[Fap1<0] <- 0
        
        #and all other reproduction is acted on less strongly - but we need a way to hold it steady
        Fap <- Fap*(1-(N/(K*2)))
        Fap[Fap<0] <- 0
      }
      
      #build seasonal survival matrices for each quarter - this varies by hunt sim.
      
      qtr1 = matrix(0,nrow=L, ncol=L) 
      qtr3 = qtr1
      diag(qtr1) = c(pi, rep(p_hi_mort, times= (length(diag(qtr1))-1))) #hi hunt is concentrated in 1st quarter
      diag(qtr3) = c(pi, rep(p_low_mort, times= (length(diag(qtr3))-1))) #hi hunt is concentrated in 1st quarter
      
      
      qtr2 = qtr3  #could alter later if you want to change the survival rate for infants halfway through
      
      
      if(a==1){
        top.row <- c(Faip, rep(Fap, length=(L-1)))
      }
      if(a==2){
        top.row <- c(0,Fap1, rep(Fap, length=(L-2)))
      }
      if(a>=3){
        top.row <- c(rep(0, length=(a-1)),Fap1,rep(Fap, length=(L-a)))
      }
      
      
      qtr4 = matrix(0,nrow=(L-1), ncol=(L))
      diag(qtr4) = c(pi, rep(p_hi_mort, times= (length(diag(qtr4))-1))) #add in sites
      qtr4[nrow(qtr4), ncol(qtr4)] = p_hi_mort
      
      
      qtr4 = rbind(top.row, qtr4) #survives 4th to next year
      
      
      pop.mat =   qtr1 %*%  qtr2 %*% qtr3 %*% qtr4
    }
    
  }else if (hunt_sim==9 | hunt_sim==10 | hunt_sim==11 | hunt_sim==12){
    #these tabulations should hold for all 1 quarter hunts
    infant_mort <- (as.numeric(df$IBI_infant_mortality)/as.numeric(df$IBI))/4
    adult_mort_hi <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.9)/3
    adult_mort_low <- ((as.numeric(df$IBI_adult_mortality)/as.numeric(df$IBI))*.1)
    
    p_hi_mort = 1- adult_mort_hi 
    p_low_mort = 1-adult_mort_low 
    pi= 1-infant_mort
    a =  as.numeric(as.character(df$age_1st_rep)) #this time, we do not want IBI units
    
    if(hunt_sim==9){
      #fecundity also affected, depending on the season of hunting (when hunting is in the 4th quarter)
      Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #make into years again
      Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_low_mort #first year breeders
      Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
      L = as.numeric(as.character(df$avg_lifespan))
      
      if(dens.dep==TRUE){
        Faip <- Faip*(1-N/K) #only applicable if a==1
        Faip[Faip<0] <- 0 #no negative births allowed
        
        Fap1 <- Fap1*(1-N/K)
        Fap1[Fap1<0] <- 0
        
        #and all other reproduction is acted on less strongly - but we need a way to hold it steady
        Fap <- Fap*(1-(N/(K*2)))
        Fap[Fap<0] <- 0
      }
      
      #build seasonal survival matrices for each quarter - this varies by hunt sim.
      
      qtr1 = matrix(0,nrow=L, ncol=L) 
      diag(qtr1) = c(pi, rep(p_hi_mort, times= (length(diag(qtr1))-1))) #hi hunt is concentrated in 1st quarter
      
      
      qtr3 = qtr2 = qtr1  #could alter later if you want to change the survival rate for infants halfway through
      
      
      if(a==1){
        top.row <- c(Faip, rep(Fap, length=(L-1)))
      }
      if(a==2){
        top.row <- c(0,Fap1, rep(Fap, length=(L-2)))
      }
      if(a>=3){
        top.row <- c(rep(0, length=(a-1)),Fap1,rep(Fap, length=(L-a)))
      }
      
      
      qtr4 = matrix(0,nrow=(L-1), ncol=(L))
      diag(qtr4) = c(pi, rep(p_low_mort, times= (length(diag(qtr4))-1))) #add in sites
      qtr4[nrow(qtr4), ncol(qtr4)] = p_low_mort
      
      
      qtr4 = rbind(top.row, qtr4) #survives 4th to next year
      
      
      pop.mat =   qtr1 %*%  qtr2 %*% qtr3 %*% qtr4
      
      
    }else if (hunt_sim==10){
      
      #fecundity also affected, depending on the season of hunting (when hunting is in the 4th quarter)
      Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #make into years again
      Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #first year breeders
      Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
      L = as.numeric(as.character(df$avg_lifespan))
      
      if(dens.dep==TRUE){
        Faip <- Faip*(1-N/K) #only applicable if a==1
        Faip[Faip<0] <- 0 #no negative births allowed
        
        Fap1 <- Fap1*(1-N/K)
        Fap1[Fap1<0] <- 0
        
        #and all other reproduction is acted on less strongly - but we need a way to hold it steady
        Fap <- Fap*(1-(N/(K*2)))
        Fap[Fap<0] <- 0
      }
      
      #build seasonal survival matrices for each quarter - this varies by hunt sim.
      
      qtr1 = matrix(0,nrow=L, ncol=L) 
      qtr3 = qtr1
      diag(qtr1) = c(pi, rep(p_low_mort, times= (length(diag(qtr1))-1))) #hi hunt is concentrated in 1st quarter
      diag(qtr3) = c(pi, rep(p_hi_mort, times= (length(diag(qtr3))-1)))
      
      qtr2 = qtr3  #could alter later if you want to change the survival rate for infants halfway through
      
      
      if(a==1){
        top.row <- c(Faip, rep(Fap, length=(L-1)))
      }
      if(a==2){
        top.row <- c(0,Fap1, rep(Fap, length=(L-2)))
      }
      if(a>=3){
        top.row <- c(rep(0, length=(a-1)),Fap1,rep(Fap, length=(L-a)))
      }
      
      
      qtr4 = matrix(0,nrow=(L-1), ncol=(L))
      diag(qtr4) = c(pi, rep(p_hi_mort, times= (length(diag(qtr4))-1))) #add in sites
      qtr4[nrow(qtr4), ncol(qtr4)] = p_hi_mort
      
      
      qtr4 = rbind(top.row, qtr4) #survives 4th to next year
      
      
      pop.mat =   qtr1 %*%  qtr2 %*% qtr3 %*% qtr4
      
    }else if (hunt_sim==11){
      #fecundity also affected, depending on the season of hunting (when hunting is in the 4th quarter)
      Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #make into years again
      Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #first year breeders
      Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
      L = as.numeric(as.character(df$avg_lifespan))
      
      if(dens.dep==TRUE){
        Faip <- Faip*(1-N/K) #only applicable if a==1
        Faip[Faip<0] <- 0 #no negative births allowed
        
        Fap1 <- Fap1*(1-N/K)
        Fap1[Fap1<0] <- 0
        
        #and all other reproduction is acted on less strongly - but we need a way to hold it steady
        Fap <- Fap*(1-(N/(K*2)))
        Fap[Fap<0] <- 0
      }
      
      #build seasonal survival matrices for each quarter - this varies by hunt sim.
      
      qtr1 = matrix(0,nrow=L, ncol=L) 
      qtr2 = qtr1
      diag(qtr1) = c(pi, rep(p_hi_mort, times= (length(diag(qtr1))-1))) #hi hunt is concentrated in 1st quarter
      diag(qtr2) = c(pi, rep(p_low_mort, times= (length(diag(qtr2))-1)))
      
      qtr3 = qtr1  #could alter later if you want to change the survival rate for infants halfway through
      
      
      if(a==1){
        top.row <- c(Faip, rep(Fap, length=(L-1)))
      }
      if(a==2){
        top.row <- c(0,Fap1, rep(Fap, length=(L-2)))
      }
      if(a>=3){
        top.row <- c(rep(0, length=(a-1)),Fap1,rep(Fap, length=(L-a)))
      }
      
      
      qtr4 = matrix(0,nrow=(L-1), ncol=(L))
      diag(qtr4) = c(pi, rep(p_hi_mort, times= (length(diag(qtr4))-1))) #add in sites
      qtr4[nrow(qtr4), ncol(qtr4)] = p_hi_mort
      
      
      qtr4 = rbind(top.row, qtr4) #survives 4th to next year
      
      
      pop.mat =   qtr1 %*%  qtr2 %*% qtr3 %*% qtr4
      
      
    }else if (hunt_sim ==12){
      
      #fecundity also affected, depending on the season of hunting (when hunting is in the 4th quarter)
      Fap= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #make into years again
      Fap1= (as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*p_hi_mort #first year breeders
      Faip=(as.numeric(as.character(df$IBI_fecundity))/as.numeric(as.character(df$IBI)))*pi #first year breeders if only infants (only applicable if a=1)
      L = as.numeric(as.character(df$avg_lifespan))
      
      if(dens.dep==TRUE){
        Faip <- Faip*(1-N/K) #only applicable if a==1
        Faip[Faip<0] <- 0 #no negative births allowed
        
        Fap1 <- Fap1*(1-N/K)
        Fap1[Fap1<0] <- 0
        
        #and all other reproduction is acted on less strongly - but we need a way to hold it steady
        Fap <- Fap*(1-(N/(K*2)))
        Fap[Fap<0] <- 0
      }
      
      #build seasonal survival matrices for each quarter - this varies by hunt sim.
      
      qtr1 = matrix(0,nrow=L, ncol=L) 
      qtr3 = qtr1
      diag(qtr1) = c(pi, rep(p_hi_mort, times= (length(diag(qtr1))-1))) #hi hunt is concentrated in 1st quarter
      diag(qtr3) = c(pi, rep(p_low_mort, times= (length(diag(qtr3))-1)))
      
      qtr2 = qtr1  #could alter later if you want to change the survival rate for infants halfway through
      
      
      if(a==1){
        top.row <- c(Faip, rep(Fap, length=(L-1)))
      }
      if(a==2){
        top.row <- c(0,Fap1, rep(Fap, length=(L-2)))
      }
      if(a>=3){
        top.row <- c(rep(0, length=(a-1)),Fap1,rep(Fap, length=(L-a)))
      }
      
      
      qtr4 = matrix(0,nrow=(L-1), ncol=(L))
      diag(qtr4) = c(pi, rep(p_hi_mort, times= (length(diag(qtr4))-1))) #add in sites
      qtr4[nrow(qtr4), ncol(qtr4)] = p_hi_mort
      
      
      qtr4 = rbind(top.row, qtr4) #survives 4th to next year
      
      
      pop.mat =   qtr1 %*%  qtr2 %*% qtr3 %*% qtr4
      
    }
  }
  
  #and return
  
  return(pop.mat)
}
get_lambda <- function(pop.mat){ 
  lambda <- max(Re(eigen(pop.mat)$value)) #largest real eigenvalue...
  return(lambda)
} 
iterate.one.timestep.lefkovitch <- function(df, Nt_next, dens.dep, use.uci, N, K, nat.mort, apply.hunt,  hunt_sim, add.noise){
  df <- data.frame(df)
  # need to rebuild the population matrix with each timestep
  
  #time series only takes you one timestep further
  times <- seq(1, (2), by = 1) 
  Nt <- matrix(0, nrow=length(Nt_next), length(times)) #create a matrix to hold your outputs. should have same # rows as your age classes and be long enough to hold your time series 
  
  Nt[,1] <- Nt_next
  
  pop.mat <- build.pop.mat(df, nat.mort=nat.mort,use.uci=use.uci, dens.dep=dens.dep, N=N, K=K)   
  
  # if(apply.hunt==TRUE){
  #   #and apply hunt penalty to population matrix
  #   pop.mat <- impose_hunt(pop.mat, hunt.intensity=hunt.intensity, hunt_sim=hunt_sim)
  # }
  
  #apply noise if stochastic simulation desired
  if(add.noise==TRUE){
    pop.mat <- wrap.noise(pop.mat) #shakes pop.mat with each time_step if true
  }
  
  #get lambda and store it.
  lambda.next <- get_lambda(pop.mat)
  
  #then take that matrix and iterate forward in time 
  Nt[,2] <- (pop.mat) %*% Nt[,1] 
  Nt[Nt<0] <- 0 #can't have neg lemurs
  
  Nt_next2 <- Nt[,ncol(Nt)]
  
  return(list(Nt_next2, lambda.next))
}
iterate.one.timestep.annual <- function(df, Nt_next, dens.dep, use.uci, N, K, nat.mort, apply.hunt, hunt_sim, add.noise){
  df <- data.frame(df)
  # need to rebuild the population matrix with each timestep
  
  #time series only takes you one timestep further
  times <- seq(1, (2), by = 1) 
  Nt <- matrix(0, nrow=length(Nt_next), length(times)) #create a matrix to hold your outputs. should have same # rows as your age classes and be long enough to hold your time series 
  
  Nt[,1] <- Nt_next
  
  pop.mat <- build.leslie.annual(df, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, N=N, K=K)   
  

  
  #apply noise if stochastic simulation desired
  if(add.noise==TRUE){
    pop.mat <- wrap.noise(pop.mat) #shakes pop.mat with each time_step if true
  }
  
  #get lambda and store it.
  lambda.next <- get_lambda(pop.mat)
  
  #then take that matrix and iterate forward in time 
  Nt[,2] <- (pop.mat) %*% Nt[,1] 
  Nt[Nt<0] <- 0 #can't have neg lemurs
  
  Nt_next2 <- Nt[,ncol(Nt)]
  
  return(list(Nt_next2, lambda.next))
}
iterate.one.timestep.periodic <- function(df, Nt_next, dens.dep, use.uci, hunt.override, set_hunt, N, K, nat.mort, apply.hunt, hunt_sim, add.noise){
  df <- data.frame(df)
  # need to rebuild the population matrix with each timestep
  
  #time series only takes you one timestep further
  times <- seq(1, (2), by = 1) 
  Nt <- matrix(0, nrow=length(Nt_next), length(times)) #create a matrix to hold your outputs. should have same # rows as your age classes and be long enough to hold your time series 
  
  Nt[,1] <- Nt_next
  
  if(apply.hunt==FALSE){ #no redistributed hunt here
    pop.mat <- build.leslie.periodic(df, nat.mort=nat.mort,use.uci=use.uci, hunt.override = hunt.override, set_hunt = set_hunt, dens.dep=dens.dep, N=N, K=K)   
  }
  
  if(apply.hunt==TRUE){
    #and apply hunt penalty to population matrix
    pop.mat <- build.leslie.periodic.hunt(df, nat.mort=nat.mort,use.uci=use.uci, hunt.override = hunt.override, set_hunt = set_hunt, dens.dep=dens.dep, N=N, K=K,  hunt_sim=hunt_sim)   
  }
  
  
  #if(apply.hunt==TRUE){
    #and apply hunt penalty to population matrix
   # pop.mat <- impose_hunt(pop.mat, hunt.intensity=hunt.intensity, hunt_sim=hunt_sim)
  #}
  
  #apply noise if stochastic simulation desired
  if(add.noise==TRUE){
    pop.mat <- wrap.noise(pop.mat) #shakes pop.mat with each time_step if true
  }
  

  
  #get lambda and store it.
  lambda.next <- get_lambda(pop.mat)
  
  #then take that matrix and iterate forward in time 
  Nt[,2] <- (pop.mat) %*% Nt[,1] 
  Nt[Nt<0] <- 0 #can't have neg lemurs
  
  Nt_next2 <- Nt[,ncol(Nt)]
  
  return(list(Nt_next2, lambda.next))
}
iterate.one.timestep.seasonal <- function(df, Nt_next, dens.dep,use.uci,  N, K, nat.mort, apply.hunt,  hunt_sim, add.noise, seas){
  df <- data.frame(df)
  # need to rebuild the population matrix with each timestep
  
  #time series only takes you one timestep further
  times <- seq(1, (2), by = 1) 
  Nt <- matrix(0, nrow=length(Nt_next), length(times)) #create a matrix to hold your outputs. should have same # rows as your age classes and be long enough to hold your time series 
  
  Nt[,1] <- Nt_next

  
  if(apply.hunt==FALSE){ #no redistributed hunt here
    pop.mat <- build.leslie.seasonal(df, nat.mort=nat.mort,use.uci=use.uci, dens.dep=dens.dep, N=N, K=K, seas=seas)   
  }
  
  if(apply.hunt==TRUE){
    #and apply hunt penalty to population matrix
    pop.mat <- build.leslie.seasonal.hunt(df, nat.mort=nat.mort,use.uci=use.uci, dens.dep=dens.dep, N=N, K=K, seas=seas, hunt_sim=hunt_sim)   
  }
  
  
  
  #apply noise if stochastic simulation desired
  if(add.noise==TRUE){
    pop.mat <- wrap.noise(pop.mat) #shakes pop.mat with each time_step if true
  }
  
  #get lambda and store it.
  lambda.next <- get_lambda(pop.mat)
  
  #then take that matrix and iterate forward in time 
  Nt[,2] <- (pop.mat) %*% Nt[,1] 
  Nt[Nt<0] <- 0 #can't have neg lemurs
  
  Nt_next2 <- Nt[,ncol(Nt)]
  
  return(list(Nt_next2, lambda.next))
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
}#for matrix slicing
slice <- function(input, by){ 
  starts <- seq(1,length(input),by)
  tt <- lapply(starts, function(y) input[y:(y+(by-1))])
  llply(tt, function(x) x[!is.na(x)])
} #for vector slicing
stab.structure <- function(pop.mat){
  lambda <- Re(eigen(pop.mat)$vector[eigen(pop.mat)$values== max(Re(eigen(pop.mat)$values))]) #largest real eigenvalue.
  stab.struct <-  Re(eigen(pop.mat)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  return(stab.struct)
} 
get_Nt_init_annual <- function(df, nat.mort, dens.dep, use.uci, N, K){
  
  pop.mat <- build.leslie.annual(df, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, N=N, K=K) #build pop.mat without density dependence
  Nt_init <- matrix(0, nrow(pop.mat), 1) #create a matrix to hold your outputs. 
  
  eigenvect <- stab.structure(pop.mat) #get the stable age structure from the leslie matrix, and remember this has quarterly classes
  plot(eigenvect,type="b",pch=19,xlab="ages - in years", ylab="proportion", main = unique(df$species), col="red", ylim=range(eigenvect))
  
  #sum the proportion of the pop in each stage and position appropriately
  n.inf <- 1
  n.juv <- as.numeric(df$age_1st_rep)-1
  n.adult <- (as.numeric(df$avg_lifespan)-n.juv - n.inf)
  
  prop.adult <- sum(eigenvect[(n.juv + n.inf + 1):length(eigenvect)])
  adult.pop <-  N
  total.pop <- adult.pop/prop.adult
  prop.inf <- sum(eigenvect[1:n.inf])
  prop.juv <- sum(eigenvect[(n.inf+1):(n.juv + n.inf)])
  inf.pop <- prop.inf*total.pop
  juv.pop <- prop.juv*total.pop
  
  # Now fill distribute those amongst the age classes. We want everyone in the same month at
  # the same time because this is a birth pulse population. The number of rows per age class
  # will vary based on the number of timesteps (explicit months modelled per year)
  # We'll start the simulation at September 15, the with the first new infants in the population,
  # such that the first babies are added a full year later.
  
  
  Nt_init[,1] <- cbind(eigenvect*total.pop)
  
  return(Nt_init)
  
}
get_Nt_init_lefkovitch <- function(df, nat.mort, use.uci, dens.dep, N, K){
  
  pop.mat <- build.pop.mat(df, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, N=N, K=K) #build pop.mat without density dependence
  Nt_init <- matrix(0, nrow(pop.mat), 1) #create a matrix to hold your outputs. 
  
  eigenvect <- stab.structure(pop.mat) #get the stable age structure from the leslie matrix, and remember this has quarterly classes
  plot(eigenvect,type="b",pch=19,xlab="ages - in stages", ylab="proportion", main = unique(df$species), col="red", ylim=range(eigenvect))
  
  
  a = as.numeric(as.character(unique(df$age_1st_rep)))/as.numeric(as.character(unique(df$IBI)))
  
  #sum the proportion of the pop in each stage and position appropriately
  if (a==1){
    prop.adult <- eigenvect[2]
    adult.pop <-  N
    total.pop <- adult.pop/prop.adult
    prop.inf <- eigenvect[1]
    inf.pop <- prop.inf*total.pop
    Nt_init[,1] <- rbind(inf.pop,adult.pop)  
  }
  if (a==2){
    prop.adult <- eigenvect[3]
    adult.pop <-  N
    total.pop <- adult.pop/prop.adult
    prop.inf <- eigenvect[1]
    prop.juv <- eigenvect[2]
    inf.pop <- prop.inf*total.pop
    juv.pop <- prop.juv*total.pop
    Nt_init[,1] <- rbind(inf.pop,juv.pop,adult.pop)  
  }
  if (a>=3){
    prop.adult1 <- eigenvect[3]
    prop.adult2 <- eigenvect[4]
    adult.pop <-  N
    total.pop <- adult.pop/sum(prop.adult1, prop.adult2)
    prop.inf <- eigenvect[1]
    prop.juv <- eigenvect[2]
    inf.pop <- prop.inf*total.pop
    juv.pop <- prop.juv*total.pop
    adult.pop1 <- prop.adult1*total.pop
    adult.pop2 <- prop.adult2*total.pop
    Nt_init[,1] <- rbind(inf.pop,juv.pop,adult.pop1, adult.pop2) 
  }
  
  return(Nt_init)
  
}
get_Nt_init_seasonal <- function(df, nat.mort, use.uci, dens.dep, N, K, seas){
  
  pop.mat <- build.leslie.seasonal(df, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, N=N, K=K, seas=seas) #build pop.mat without density dependence
  Nt_init <- matrix(0, nrow(pop.mat), 1) #create a matrix to hold your outputs. 
  
  eigenvect <- stab.structure(pop.mat) #get the stable age structure from the leslie matrix, and remember this has quarterly classes
  plot(eigenvect,type="b",pch=19,xlab="ages - in quarter years", ylab="proportion", main = unique(df$species), col="red", ylim=range(eigenvect))
  
  #sum the proportion of the pop in each stage and position appropriately
  n.inf <- 1 #where infant starts
  n.juv <- 1 + seas #where juv starts
  n.adult <- (as.numeric(df$age_1st_rep)*seas +1) #where adult starts
  
  prop.adult <- sum(eigenvect[n.adult:length(eigenvect)])
  adult.pop <-  N
  total.pop <- adult.pop/prop.adult
  prop.inf <- sum(eigenvect[1:4])
  prop.juv <- sum(eigenvect[n.juv:(n.adult-1)])
  inf.pop <- prop.inf*total.pop
  juv.pop <- prop.juv*total.pop
  
  # Now fill distribute those amongst the age classes. We want everyone in the same month at
  # the same time because this is a birth pulse population. The number of rows per age class
  # will vary based on the number of timesteps (explicit months modelled per year)
  # We'll start the simulation at September 15, the with the first new infants in the population,
  # such that the first babies are added a full year later.
  
  
  Nt_init[n.inf,1] <- inf.pop 
  Nt_init[n.juv,1] <- juv.pop 
  Nt_init[n.adult,1] <- adult.pop 
  
  return(Nt_init)
  
}
get_Nt_init_periodic <- function(df, nat.mort, use.uci, hunt.override, set_hunt, dens.dep, N, K){
  
  pop.mat <- build.leslie.periodic(df, nat.mort=nat.mort, use.uci=use.uci, hunt.override=hunt.override, set_hunt=set_hunt, dens.dep=dens.dep, N=N, K=K) #build pop.mat without density dependence
  #pop.mat.1 <- build.qtr.periodic(df, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, N=N, K=K)#and try just the quarter pop mat
  Nt_init <- matrix(0, nrow(pop.mat), 1) #create a matrix to hold your outputs. 
  #in this case, we don't have a stable age structure to start our simulations (landscape is periodic) so we want the stable age strucutre
  #of the first quarter (the sub-matrix at which point we start
  
  eigenvect <- stab.structure(pop.mat) #get the stable age structure from the leslie matrix, and remember this has quarterly classes
  plot(eigenvect,type="b",pch=19,xlab="ages", ylab="proportion", main = unique(df$species), col="red", ylim=range(eigenvect))
  
  #sum the proportion of the pop in each stage and position appropriately
  n.inf <- 1
  n.juv <- as.numeric(df$age_1st_rep)-1
  n.adult <- (as.numeric(df$avg_lifespan)-n.juv - n.inf)
  
  prop.adult <- sum(eigenvect[(n.juv + n.inf + 1):length(eigenvect)])
  adult.pop <-  N
  total.pop <- adult.pop/prop.adult
  prop.inf <- sum(eigenvect[1:n.inf])
  prop.juv <- sum(eigenvect[(n.inf+1):(n.juv + n.inf)])
  inf.pop <- prop.inf*total.pop
  juv.pop <- prop.juv*total.pop
  
  
  # Now fill distribute those amongst the age classes. We want everyone in the same month at
  # the same time because this is a birth pulse population. The number of rows per age class
  # will vary based on the number of timesteps (explicit months modelled per year)
  # We'll start the simulation at September 15, the with the first new infants in the population,
  # such that the first babies are added a full year later.
  
  Nt_init[,1] <- cbind(eigenvect*total.pop)
  return(Nt_init)
  
}
sum.by.stage <- function(Nmat){
  n.sites <- dim(Nmat)[3]
  #now assign "a" based on matrix shape
  a <- (nrow(Nmat[,,1]))-1
  
  if(a==1){
    infants <- list()
    adults <- list()
    
    for (i in 1:n.sites){
      infants[[i]] <- Nmat[1,,i] #gives a list for each site that is a string of vectors of infants
      adults[[i]] <- Nmat[2,,i] #gives a list of vectors of adults
    }
    
    infants <- colSums(do.call("rbind", infants))
    adults <- colSums(do.call("rbind", adults))
    #and include 0s for juveniles so that our populations can later merge together
    juveniles = rep(0, length(adults))
    
    
    N.meta.pop <- rbind(infants, juveniles, adults)
  }
  
  
  if(a==2){
    infants <- list()
    juveniles <- list()
    adults <- list()
    
    for (i in 1:n.sites){
      infants[[i]] <- Nmat[1,,i] #gives a list for each site that is a string of vectors of infants
      juveniles[[i]] <- Nmat[2,,i]#gives a list of vectors of juveniles
      adults[[i]] <- Nmat[3,,i] #gives a list of vectors of adults
    }
    
    infants <- colSums(do.call("rbind", infants))
    juveniles <- colSums(do.call("rbind", juveniles))
    adults <- colSums(do.call("rbind", adults))
    
    N.meta.pop <- rbind(infants, juveniles, adults)
  }
  
  if(a>=3){
    infants <- list()
    juveniles <- list()
    adults1 <- list()
    adults2 <- list()
    
    for (i in 1:n.sites){
      infants[[i]] <- Nmat[1,,i] #gives a list for each site that is a string of vectors of infants
      juveniles[[i]] <- Nmat[2,,i]#gives a list of vectors of juveniles
      adults1[[i]] <- Nmat[3,,i] #gives a list of vectors of adults of first breeding year
      adults2[[i]] <- Nmat[4,,i] #gives a list of vectors of adults
    }
    
    infants <- colSums(do.call("rbind", infants))
    juveniles <- colSums(do.call("rbind", juveniles))
    adults1 <- colSums(do.call("rbind", adults1))
    adults2 <- colSums(do.call("rbind", adults2))
    
    #but you want the sums of adults1 and adults2 together as just adults
    adults <- colSums(do.call("rbind", list(adults1, adults2)))
    
    N.meta.pop <- rbind(infants, juveniles, adults)
  }
  
  return(N.meta.pop)
}
make.connect.mat <- function(site.df, ts1, connect.hi){
  
  
  #first make a matrix of connected sites, then expand those to reflect each age class
  n.dim  <- length(unique(site.df$site))
  # gives the number of sites to set the dimensions of our connectivity matrix
  
  site.df$lat <- as.numeric(as.character(site.df$lat))
  site.df$long <- as.numeric(as.character(site.df$long))
  
  
  #and build connectivity matrix from that
  #the rowSums must always add to one. we can insert warnings to check
  euc.dist.matrix<-matrix(NA,nrow=n.dim, ncol=n.dim) #sets up matrix with one space per row per column for each site
  
  for(i in 1:n.dim){
    x.from<-rep(site.df[i,]$long, n.dim) #will take the long value from the row corresponding to i (each site is numbered 1-67)
    y.from<-rep(site.df[i,]$lat, n.dim)
    euc.dist<-sqrt((site.df$long - x.from)^2 + (site.df$lat-y.from)^2)
    #gives a matrix of euclidian distances between sites. now rescale these to 1 
    #unless, we want hi connectivity, in which case, we say dispersal is equal across sites
    if(connect.hi==TRUE){
      small.sites <- length(euc.dist)
      euc.dist =  rep(sum(euc.dist)/small.sites, times=small.sites)
    }
    euc.dist.matrix[i,] =  euc.dist/sum(euc.dist)
  }

  #now make everything a probability of dispersal by taking 1- the matrix
  euc.dist.matrix <- 1-euc.dist.matrix

  #and rescale again as a proportion of one:
  for (i in 1:nrow(euc.dist.matrix)){
    tmp <- sum(euc.dist.matrix[i,])
    past <- euc.dist.matrix[i,] 
    euc.dist.matrix[i,] <- past/tmp
    if(round(sum(euc.dist.matrix[i,]), digits=4) !=1){
      warning(paste("at i=", i, "rowSums probabilities do not add to 1", sep=" "))
    }
  }
  #and once more the other way
  for (i in 1:ncol(euc.dist.matrix)){
    tmp <- sum(euc.dist.matrix[,i])
    past <- euc.dist.matrix[,i] 
    euc.dist.matrix[,i] <- past/tmp
    
    if(round(sum(euc.dist.matrix[,i]), digits=4) !=1){
      warning(paste("at i=", i, "colSums probabilities do not add to 1", sep=" "))
    } 
  }

  #now disallow dispersal between makira and masoala - can add back in if needed later
  
  #and now for the makira-masoala sites
  #if(length(unique(site.df2$region)) >1){
  # no.disperse.to.makira <- c(as.character(site.df2$region))
  #no.disperse.to.makira[no.disperse.to.makira=="makira"] <- 0
  #no.disperse.to.makira[no.disperse.to.makira=="masoala"] <- 1
  #no.disperse.to.makira <- as.numeric(no.disperse.to.makira)
  
  #no.disperse.to.masoala <- c(as.character(site.df2$region))
  #no.disperse.to.masoala[no.disperse.to.masoala=="masoala"] <- 0
  #no.disperse.to.masoala[no.disperse.to.masoala=="makira"] <- 1
  #no.disperse.to.masoala <- as.numeric(no.disperse.to.masoala)
  
  #now multiply each column, element-wise, by the appropriate vector
  #  n.mak <- length(site.df2$region[site.df2$region=="makira"])
  # n.mas <- length(site.df2$region[site.df2$region=="masoala"])
  
  #makira always comes first (factored that way), so first disallow dispersal to masoala
  #for (i in 1:n.mak){
  # tmp <- as.vector(euc.dist.matrix[,i])*no.disperse.to.masoala
  #euc.dist.matrix[,i] <- tmp
  #}
  #then disallow dispersal to makira
  #for (i in (n.mak+1):(n.mak+n.mas)){
  # tmp <- as.vector(euc.dist.matrix[,i])*no.disperse.to.makira
  #euc.dist.matrix[,i] <- tmp
  #}
  #}  
  
  #columns will no longer total to 1 but they will equal each other
  #and a check
  # if(round(sum(colSums(euc.dist.matrix)),4)!=round(sum(rowSums(euc.dist.matrix)),4)){
  #  warning(paste("post-correction #2, row and column sums are unequal"))
  #}
  
  #here's where you rescale:
  #for (i in 1:ncol(euc.dist.matrix)){
  # tmp <- sum(euc.dist.matrix[,i])
  #past <- euc.dist.matrix[,i] 
  #euc.dist.matrix[,i] <- past/tmp
  #}
  
  #now test with warning
  #and a check
  #if(round(sum(euc.dist.matrix[i,]), digits=4) !=1){
  # warning(paste("at i=", i, "rowSums probabilities do not add to 1", sep=" "))
  #}
  
  
  #now replicate each such that each age class interconnects with other (similar)
  #age classes across many sites. first make a blank matrix of 0s and then fill in 
  #where needed with the site-connectivity values
  
  #that blank matrix:
  n.dim2 <- length(ts1)
  
  euc.dist.matrix2 <-matrix(NA,nrow=n.dim2, ncol=n.dim2) 
  
  #take each element of the matrix and make it the diagonal of its own matrix that is 
  #the same length as the number of age classes
  
  #first, find those months
  #n.mon.inf <- unique(site.df$stage_duration[site.df$stage=="infant"])
  #n.mon.juv <- unique(site.df$stage_duration[site.df$stage=="juvenile"])
  #n.mon.adult <- unique(site.df$stage_duration[site.df$stage=="adult"])
  #sum.classes <- sum(n.mon.inf,n.mon.juv,n.mon.adult)
  n.sites <- length(unique(site.df$site))
  sum.classes=n.dim2/n.sites
  
  
  tot.mat.list <- list()
  for (i in 1:n.sites){
    #take site "i" and make a list of dispersal probabilities to the other sites
    tmp <- euc.dist.matrix[i,]
    
    #expand each probability by the number of age classes
    matrix.list <- c()
    for (j in 1:length(tmp)){
      matrix.list[[j]] <- matrix(0, nrow=sum.classes, ncol=sum.classes)
      diag(matrix.list[[j]]) <- tmp[j]
    }
    
    #bind them together end-to-end in a long matrix
    tmp2 <- matrix(unlist(matrix.list), nrow=sum.classes, ncol=sum.classes*n.sites)
    
    #store them in a list and then bind them later by rows
    tot.mat.list[[i]] <- tmp2
  }
  
  euc.dist.matrix2 <- do.call("rbind", tot.mat.list)
  print("euc.dist.matrix 2.  1")
  print(euc.dist.matrix2)
  
  #now write over the columns with infants such that they will not disperse at all
  #first, find the infant columns (every first +3)
  n.inf=1
  inf.columns <- list()
  for (i in 0:(n.sites-1)){
    inf.columns[[i+1]] <- c(seq(1, (n.inf),1) + sum.classes*i)
  }
  inf.columns <- c(unlist(inf.columns))
  
  #now make sure the infants can't disperse
  
  for (i in 1:length(inf.columns)){
    j <- inf.columns[i]
    tmp <- rep(0, ncol(euc.dist.matrix2))
    tmp[j] <- 1
    euc.dist.matrix2[,j] <- tmp #fill in the appropriate columns
  }
  
  #and name its rows and columns the same as our sites - not working, but who cares anyway?
  #name.tmp <- paste(site.df$title, site.df$stage, sep = "_")
  #name.rep <- rep(c(unique(site.df$stage_duration)), times=n.sites)
  #names.tot <- rep(name.tmp, name.rep)
  #rownames(euc.dist.matrix2) <- as.character(names.tot)
  #colnames(euc.dist.matrix2) <- as.character(names.tot)
  
  print("euc.dist.matrix 2.  2")
  print(euc.dist.matrix2)
  
  #and return the named matrix
  return(euc.dist.matrix2)
}
plot.meta.pop <- function(meta.list, df, yrs,seas, do.save, plot.as.metapop, leslie){ #bring in the results of meta.wrap--both your site-specific totals and your summation
  if(leslie=="lefkovitch"){
    timesteps = yrs/as.numeric(as.character(unique(df$IBI)))
    a <- as.numeric(as.character(unique(df$age_1st_rep)))/as.numeric(as.character(unique(df$IBI)))
  }
  if(leslie=="annual" | leslie=="periodic"){
    timesteps = yrs
    a <- as.numeric(as.character(unique(df$age_1st_rep)))
  }
  
  if(leslie=="seasonal"){
    timesteps = yrs*seas
    a <- as.numeric(as.character(unique(df$age_1st_rep)))
  }
  
  
  #main.title <- unique(dat.list$species)
  #first pull apart each part of the list
  
  sum.df <- data.frame(meta.list[1])
  site.df <- data.frame(meta.list[2])
  lambda.list <- as.vector(c(unlist(meta.list[3])))
  
  n.sites=length(unique(df$site))
  
  n.stages <- nrow(sum.df)
  
  rownames(site.df) <- paste(rep(df$site, each=n.stages),rep(rownames(sum.df), times=n.sites), sep="-")
  
  sum.df <- cbind.data.frame(rownames(sum.df), sum.df)
  colnames(sum.df) <- c("stage", paste0("timestep_", seq(1, (timesteps), 1)))
  rownames(sum.df) <- c()
  
  #and for site.df too
  site.df <- cbind.data.frame(rownames(site.df), site.df)
  colnames(site.df) <- c("stage", paste0("timestep_", seq(1, (timesteps), 1)))
  rownames(site.df) <- c()
  
  #now melt them such that we have a LONG table with these column headers: (site, stage, timestep, count)
  #one of the sites can be "metapop"
  site.df <- reshape(site.df, varying=paste0("timestep_",seq(1, (timesteps), 1)), v.names="count", timevar="timestep", times=seq(1,(timesteps),1), direction="long")
  #now split into "stage" and "site"
  site.df$site <- NA
  site.df$site <- sapply(strsplit(as.character(site.df$stage), "-"), "[", 1)
  site.df$stage <- sapply(strsplit(as.character(site.df$stage), "-"), "[", 2)
  
  sum.df <- reshape(sum.df, varying=paste0("timestep_",seq(1, (timesteps), 1)), v.names="count", timevar="timestep", times=seq(1,(timesteps),1), direction="long")
  sum.df$site <- "metapop"
  
  max.y <- max(sum.df$count)
  max.y.site <- max(site.df$count)
  lambda.init <- round(lambda.list[2], 2)
  
  if(plot.as.metapop==TRUE){
    p <- ggplot() +
      geom_line(data=sum.df, aes(x=timestep, y=count, linetype=stage), colour="black", size=1.5) +
      geom_line(data=site.df, aes(x=timestep, y=count, colour=site, linetype=stage), show.legend=FALSE) +
      annotate("text", x=25, y=max.y, label=paste0("initial lambda = ", lambda.init)) +
      theme_bw() + xlab("timesteps (=1000 yrs)") + ylab("# lemurs") + labs(title=paste(unique(df$species), unique(df$scenario), sep="-"))# scale_y_log10() + 
    
    print(p)
    
    if(do.save==TRUE){
      ggsave(file = paste0(paste("/Users/Cara/Documents/R/R_repositories/hunting-lemurs/final_manuscript/figures/connect-plots/TS", leslie, unique(df$species), unique(df$scenario), sep="-"), ".tiff"),
             units="mm",  
             width=60, 
             height=70, 
             scale=3, 
             dpi=300)
      
    }
  }
  
  if (plot.as.metapop==FALSE){ #print just by site
    p <- ggplot() +
      geom_line(data=site.df, aes(x=timestep, y=count, colour=site, linetype=stage)) + scale_colour_hue(guide='none')+
      annotate("text", x=25, y=max.y.site, label=paste0("initial lambda = ", lambda.init)) +
      theme_bw() + xlab("timesteps (=1000 yrs)") + ylab("# lemurs") + 
      labs(title=paste(unique(df$species), unique(df$scenario), sep="-"))
    print(p)
    
    if(do.save==TRUE){
      ggsave(file = paste0(paste("/Users/Cara/Documents/R/R_repositories/hunting-lemurs/final_manuscript/figures/connect-plots/TS-site-plots", leslie, unique(df$species), unique(df$scenario), sep="-"), ".tiff"),
             units="mm",  
             width=60, 
             height=70, 
             scale=3, 
             dpi=300)
    }
  }
}
plot.single.pop <- function(spp.list, yrs, df, do.save, sim_type){ #bring in the results of meta.wrap--both your site-specific totals and your summation
  Nvector <- spp.list[1]
  lambda.list <- c(unlist(spp.list[2]))
  main_title <- paste(df$species, df$scenario, sep="-")
  timesteps = round(yrs/as.numeric(df$IBI),0)
  sum.df <- data.frame(Nvector)
  rownames(sum.df) <- c("infants", "juveniles", "adults")
  sum.df <- cbind.data.frame(rownames(sum.df), sum.df)
  colnames(sum.df) <- c("stage", paste0("timestep_", seq(1, (timesteps), 1)))
  rownames(sum.df) <- c()
  
  sum.df <- reshape(sum.df, varying=paste0("timestep_",seq(1, (timesteps), 1)), v.names="count", timevar="timestep", times=seq(1,(timesteps),1), direction="long")
  max.y <- max(sum.df$count)
  
  #include initial lambda and equilibrium lambda after density-dependence kicks into play
  #init is the first lambda and equilibrium is averaged over the last 100 timesteps
  #only plot init for now
  lambda.init <- round(lambda.list[2], 2)
  #lambda.equil <- mean(lambda.list[(length(lambda.list)-100):length(lambda.list)])
  p <- ggplot() +
    geom_line(data=sum.df, aes(x=timestep, y=count, colour=stage), size=1.5) +
    annotate("text", x=25, y=max.y, label=paste0("initial lambda = ", lambda.init)) +
    theme_bw() + xlab("timesteps (=1000 yrs)") + ylab("# lemurs") + labs(title=main_title)# scale_y_log10() + 
  
  print(p)
  
  filetitle = paste(paste(strsplit(main_title, " ")[[1]][1], strsplit(main_title, " ")[[1]][2], sep="-"), sim_type, sep="-")
  
  if (do.save==TRUE){
    ggsave(file = paste0(paste0("/Users/Cara/Documents/R/R_repositories/hunting-lemurs/final_manuscript/figures/TS", filetitle), ".tiff"),
           units="mm",  
           width=60, 
           height=70, 
           scale=3, 
           dpi=300)
  }
}
collapse.to.stage <- function(Nmat, meta.df){
  n.sites <- length(unique(meta.df$site))
  #n.mon.inf <- unique(meta.df$stage_duration[meta.df$stage=="infant"])
  #n.mon.juv <- unique(meta.df$stage_duration[meta.df$stage=="juvenile"])
  #n.mon.adult <- unique(meta.df$stage_duration[meta.df$stage=="adult"])
  a = as.numeric(unique(meta.df$age_1st_rep))
  L = as.numeric(unique(meta.df$avg_lifespan))
  
  if (a>=2){
    infants <- list()
    juveniles <- list()
    adults <- list()
    Nmeta.pop.list <- list()
    for (i in 1:n.sites){
      infants[[i]] <- Nmat[1,,i] #gives a list of vectors of infants over time
      if(a==2){
        juveniles[[i]] <- Nmat[2,,i]
      }
      if(a>2){
        juveniles[[i]] <- colSums(Nmat[2:a,,i]) #gives a list of vectors of juveniles  
      }
      adults[[i]] <- colSums(Nmat[(a+1):(L),,i]) #gives a list of vectors of adults
      
      Nmeta.pop.list[[i]] <- rbind(infants[[i]],juveniles[[i]], adults[[i]])
    }
    
    #Nmeta.pop is what you want.
    #Now bind it
    N.meta.pop <-  do.call("rbind", Nmeta.pop.list)
  }
  
  if (a==1){
    infants <- list()
    juveniles <- list()
    adults <- list()
    Nmeta.pop.list <- list()
    for (i in 1:n.sites){
      infants[[i]] <- Nmat[1,,i] #gives a list of vectors of infants over time
      adults[[i]] <- colSums(Nmat[2:(L),,i]) #gives a list of vectors of adults
      juveniles[[i]] = 0
      Nmeta.pop.list[[i]] <- rbind(infants[[i]],  juveniles[[i]], adults[[i]])
    }
    
    
    #Nmeta.pop is what you want.
    #Now bind it
    N.meta.pop <-  do.call("rbind", Nmeta.pop.list)
    
  }
  
  
  return(N.meta.pop)
}
collapse.to.stage.seasonal <- function(Nmat, meta.df, seas){
  n.sites <- length(unique(meta.df$site))
  a = as.numeric(unique(meta.df$age_1st_rep))
  L = as.numeric(unique(meta.df$avg_lifespan))
  
  #now we decide based on whether or not there are juveniles
  if(a>=2){
    infants <- list()
    juveniles <- list()
    adults <- list()
    Nmeta.pop.list <- list()
    for (i in 1:n.sites){
      infants[[i]] <- colSums(Nmat[1:(1*seas),,i]) #gives a list of vectors of infants over time
      juveniles[[i]] <- colSums(Nmat[2:(a*seas),,i]) #gives a list of vectors of juveniles  
      adults[[i]] <- colSums(Nmat[(a*seas+1):(L*seas),,i]) #gives a list of vectors of adults
      Nmeta.pop.list[[i]] <- rbind(infants[[i]],juveniles[[i]], adults[[i]])
    }
    
    #Nmeta.pop is what you want.
    #Now bind it
    N.meta.pop <-  do.call("rbind", Nmeta.pop.list)
  }
  if(a==1){
    infants <- list()
    adults <- list()
    Nmeta.pop.list <- list()
    for (i in 1:n.sites){
      infants[[i]] <- colSums(Nmat[1:(1*seas),,i]) #gives a list of vectors of infants over time
      adults[[i]] <- colSums(Nmat[(1*seas+1):(L*seas),,i]) #gives a list of vectors of adults
      Nmeta.pop.list[[i]] <- rbind(infants[[i]], adults[[i]])
    }
    
    #Nmeta.pop is what you want.
    #Now bind it
    N.meta.pop <-  do.call("rbind", Nmeta.pop.list)
  }
  
  
  return(N.meta.pop)
}
connect.leslie.sim.annual <- function(df, N, yrs, nat.mort, use.uci, dens.dep, apply.hunt, hunt.intensity, hunt_sim, add.noise, connect.hi, return.meta, allow.migration){
  
  #first, convert yrs into the timesteps of your model
  timesteps <- yrs
  
  times <- seq(1, timesteps, by = 1) 
  
  #Then make a population vector with columns to cover every timestep and
  #a row for every age class. In the case of the Lefkovitch matrix, this is always just 3
  L = unique(as.numeric(df$avg_lifespan))
  n.sites=length(unique(df$site))
  #set K as double N (will later elevate to metapop)
  K= 2*N
  
  Nmeta <- matrix(NA, nrow=L*n.sites, ncol=timesteps)
  
  a=as.numeric(unique(df$age_1st_rep))
  
  #now split by site.
  df.list <- dlply(df, .(site))
  #now fill the first column of that matrix in with the initial population size
  #we'll start with a stable structure
  
  init_list <- lapply(df.list, get_Nt_init_annual, nat.mort=nat.mort, dens.dep=dens.dep, use.uci=use.uci, N=N, K=K)
  
  Nmeta[,1] <- do.call("rbind", init_list)
  
  #if there are NAs or Infs in this vector, then we need to get rid of them
  Nmeta[,1][is.na(Nmeta[,1])] <- 0
  Nmeta[,1][Nmeta[,1]==Inf] <- 0
  
  #and make an empty vector to store your lambdas
  lambda.chain <- rep(NA, length=timesteps)
  
  #and your proportion of unoccupied patches
  prop.zero.chain <- rep(NA, length=timesteps)
  
  #fill in with your starting prop (should be 0 sites)
  prop.zero.chain[1] <- find.prop.zero(Nmeta[,1], a)
  
  #reset N and K at the metapop level
  # N=N*n.sites
  #K= 2*N
  
  #build your connectivity matrix here:
  #then, to prep for the sims to come, make connectivity matrix between sites
  #this connectivity matrix does not allow for dispersal between makira or 
  #masoala or from infant age classes 
  #(both these qualification can still be shaken slightly with the noise, below)
  if(allow.migration==TRUE){
    connect.mat.origin <- make.connect.mat(site.df=df, ts1=Nmeta[,1], connect.hi=connect.hi) 
  }
  
  #and then we iterate forward
  for (t in 1:(length(times)-1)){
    
    print(paste0("t=", t))
    
    Nt_next <- Nmeta[,t]
    
    #and make sure that "N" changes based on adults in the age class
    index = rep(c(rep(FALSE, times=a-1), rep(TRUE, times = (L-a+1))), times=n.sites)
    # but we want to keep it local, as a list. Density-dependence will sum. all reproducing individuals,
    #including last year juveniles who will reproduce at the end of this year. So in cases where first age of 
    #rep is 1, these will be counted
    #so slice by site
    N.list <- lapply(slice(Nt_next[index], (L-a+1)),sum)
    
    print(paste("metapopulation total into dens.dep survivorship =", sum(c(unlist(N.list))), sep=" "))
    
    #and print "K" so we know it is not changing through the time series
    print(paste0("K =", K))
    
    #we now apply the survival probabilities (and density dependent fecundities) over 
    #one year (and shake as necessary)
    #first split matrix so that you have one matrix per site
    Nt_next.list <- slice(Nt_next, L)
    
    output <- mapply(iterate.one.timestep.annual, df=df.list, Nt_next=Nt_next.list, N=N.list, K=K, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, apply.hunt=apply.hunt, hunt.intensity=hunt.intensity, hunt_sim=hunt_sim, add.noise=add.noise, SIMPLIFY=FALSE)
    #produces a list of by site of vertical vectors by site one timestep later plus a list of lambdas
    Nnext2 <- list()
    for (i in 1:n.sites){
      Nnext2[[i]] <- output[[i]][1]
    }
    Nnext2 <- matrix(unlist(Nnext2), nrow=length(Nt_next), ncol=1, byrow=T)
    lambda_next <- list()
    for(i in 1:n.sites){
      lambda_next[[i]] <- output[[i]][2]
    }
    #take the mean lambda across all sites
    lambda_next <- mean(c(unlist(lambda_next)))
    
    Nnext2[Nnext2<1] <- 0
    #can't have negative lemurs (but partial lemurs OK)
    
    print(paste("population post-survival, pre-dispersal =", sum(Nnext2), sep=" "))
    
    #if the population has reached an infinite sum by the end of the simulation, 
    #we are going to need to constrain it somehow - here we just hold it constant at 
    #whatever the timestep before, and we tell you this. 
    if (sum(Nnext2)==Inf){
      Nnext2 <- Nmeta[,t]
      #then print again to announce correction
      print(paste("population reached infinity. corrected, and reset to previous timestep. population, pre-dispersal =", sum(Nnext2), sep=" ")) 
    }
    
    #now shuffle lemurs based on connectivity matrix--some stay in current
    #location and some move to others
    
    #convert connect.mat to probabilities based on dispersal prob of species (and stages) in question
    #needs to have the same dimensions as entries in the Nnext2 matrix (so one locality per "age" class and a corresponding probability of disperesal)
    if (allow.migration==TRUE){ #& t %% (12/lts) != 0){
      #dispersal only happens in the non-breeding intervals and when it is explicitly turned on. 
      
      # we remake the connectivity matrix with each   
      #add noise to the dispersal probabilities if applicable
      if (add.noise==TRUE){
        connect.mat <-  wrap.noise(connect.mat.origin) #shake it up a little
        #we specifiy that we are shaking the original matrix because we don't want to
        #get stuck on a positivie migration trajectory and just keep building from there
        #this way, we're just shaking the original every time--sometimes this might allow for 
        #introductions and sometimes dispersals in and out of our system
      }
      else{ #no noise just gives the original connectivity matrix
        connect.mat <- connect.mat.origin
      }
      
      Nnext2 <- connect.mat %*% Nnext2 #and we disperse 
    }
    
    
    print(paste("population post-survival, post-dispersal (if applicable) =", sum(Nnext2), sep=" "))
    
    #now track the proportion of 0s in your sub-sites within the metapop
    prop.zero <- find.prop.zero(Nnext2,a)
    prop.zero.chain[t+1] <- prop.zero
    
    #now fill in the next space in the larger matrix and step forward again
    Nmeta[,(t+1)] <-   Nnext2
    #and fill in your lambda chain, as well (will be 1 short of length of time series because it corresponds to each transition)
    #we'll leave slot one as NA to keep the ts matching
    lambda.chain[t+1] <- lambda_next
    
  }
  
  #first slice by site
  Nmeta.site.init <- mat_split(Nmeta, r=L, c=ncol(Nmeta)) 
  
  #then collapse to stage within each site
  Nmeta_pop <- collapse.to.stage(Nmat=Nmeta.site.init, meta.df=df)
  
  #first slice by site
  if(a>=2){
    Nmeta.site <- mat_split(Nmeta_pop, r=3, c=ncol(Nmeta))  
  }
  if(a==1){
    Nmeta.site <- mat_split(Nmeta_pop, r=2, c=ncol(Nmeta))  
  }
  
  
  #then sum across all
  Nmeta2 <- sum.by.stage(Nmat=Nmeta.site)
  
  if(return.meta==TRUE){
    return(list(Nmeta2, Nmeta_pop, lambda.chain, prop.zero.chain)) #returns the sum across all sites by stage (Nmeta2), the joined version by each site as one long matrix (Nmeta) and the chain of averaged lambdas across all sites at the origin
  }
  if(return.meta==FALSE){
    return(list(Nmeta.site, lambda.chain, prop.zero.chain)) #gives the chain by stage by site and the chain of lambdas
  }
  
}
connect.leslie.sim.seasonal <- function(df, N, yrs, nat.mort, dens.dep, use.uci, apply.hunt, hunt.intensity, hunt_sim, add.noise, seas, connect.hi, return.meta, allow.migration){
  
  #first, convert yrs into the timesteps of your model
  timesteps <- yrs*seas
  
  times <- seq(1, timesteps, by = 1) 
  
  #Then make a population vector with columns to cover every timestep and
  #a row for every age class. In the case of the Lefkovitch matrix, this is always just 3
  L = unique(as.numeric(as.character(df$avg_lifespan)))
  n.sites=length(unique(df$site))
  #set K as double N
  K= 2*N
  
  Nmeta <- matrix(NA, nrow=(L*seas*n.sites), ncol=timesteps)
  a=as.numeric(unique(df$age_1st_rep))
  #now fill the first column of that matrix in with the initial population size
  #we'll start with a stable structure
  
  df.list <- dlply(df, .(site))
  
  #split by site
  init_list <- lapply(df.list, get_Nt_init_seasonal, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, N=N, K=K, seas=seas)
  
  
  Nmeta[,1] <- do.call("rbind", init_list)
  
  #if there are NAs or Infs in this vector, then we need to get rid of them
  Nmeta[,1][is.na(Nmeta[,1])] <- 0
  Nmeta[,1][Nmeta[,1]==Inf] <- 0
  
  #and make an empty vector to store your lambdas
  lambda.chain <- rep(NA, length=timesteps)
  
  #and your proportion of unoccupied patches
  prop.zero.chain <- rep(NA, length=timesteps)
  
  #fill in with your starting prop (should be 0 sites)
  prop.zero.chain[1] <- find.prop.zero(Nmeta[,1], a)
  
  
  #reset N and K at the metapop level
  # N=N*n.sites
  #K= 2*N
  
  #build your connectivity matrix here:
  #then, to prep for the sims to come, make connectivity matrix between sites
  #this connectivity matrix does not allow for dispersal between makira or 
  #masoala or from infant age classes 
  #(both these qualification can still be shaken slightly with the noise, below)
  if (allow.migration==TRUE){
    connect.mat.origin <- make.connect.mat(site.df=df, ts1=Nmeta[,1], connect.hi=connect.hi) 
  }
  
  
  #and make an empty vector to store your lambdas
  lambda.chain <- rep(NA, length=timesteps)
  
  #and then we iterate forward
  for (t in 1:(length(times)-1)){
    
    print(paste0("t=", t))
    
    Nt_next <- Nmeta[,t]
    
    
    #and make sure that "N" changes based on adults in the age class
    index = rep(c(rep(FALSE, times=a*seas-1), rep(TRUE, times = ((L-a)*seas+1))), times=n.sites)
    #keep N local
    N.list <- lapply(slice(Nt_next[index], (((L-a)*seas+1))),sum)
    
    
    print(paste("metapopulation total into dens.dep survivorship =", sum(c(unlist(N.list))), sep=" "))
    
    #and print "K" so we know it is not changing through the time series
    print(paste0("K =", K))
    
    #we now apply the survival probabilities (and density dependent fecundities) over 
    #one year (and shake as necessary)
    #first split matrix so that you have one matrix per site
    Nt_next.list <- slice(Nt_next, (L*seas))
    
    output <- mapply(iterate.one.timestep.seasonal, df=df.list, Nt_next=Nt_next.list, N=N.list, K=K, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, apply.hunt=apply.hunt, hunt.intensity=hunt.intensity, hunt_sim=hunt_sim, add.noise=add.noise, seas=seas, SIMPLIFY=FALSE)
    
    Nnext2 <- list()
    for (i in 1:n.sites){
      Nnext2[[i]] <- output[[i]][1]
    }
    Nnext2 <- matrix(unlist(Nnext2), nrow=length(Nt_next), ncol=1, byrow=T)
    lambda_next <- list()
    for(i in 1:n.sites){
      lambda_next[[i]] <- output[[i]][2]
    }
    #take the mean lambda across all sites
    lambda_next <- mean(c(unlist(lambda_next)))
    
    Nnext2[Nnext2<1] <- 0
    #can't have negative lemurs (but partial lemurs OK)
    
    print(paste("population post-survival, pre-dispersal =", sum(Nnext2), sep=" "))
    
    if (sum(Nnext2)==Inf){
      Nnext2 <- Nmeta[,t]
      #then print again to announce correction
      print(paste("population reached infinity. corrected, and reset to previous timestep. population, pre-dispersal =", sum(Nnext2), sep=" ")) 
    }
    
    #now shuffle lemurs based on connectivity matrix--some stay in current
    #location and some move to others
    
    #convert connect.mat to probabilities based on dispersal prob of species (and stages) in question
    #needs to have the same dimensions as entries in the Nnext2 matrix (so one locality per "age" class and a corresponding probability of disperesal)
    if ((allow.migration==TRUE) & (t %% (seas) == 0)){
      #dispersal only happens in the non-breeding intervals and when it is explicitly turned on. 
      
      # we remake the connectivity matrix with each   
      #add noise to the dispersal probabilities if applicable
      if (add.noise==TRUE){
        connect.mat <-  wrap.noise(connect.mat.origin) #shake it up a little
        #we specifiy that we are shaking the original matrix because we don't want to
        #get stuck on a positivie migration trajectory and just keep building from there
        #this way, we're just shaking the original every time--sometimes this might allow for 
        #introductions and sometimes dispersals in and out of our system
      }
      else{ #no noise just gives the original connectivity matrix
        connect.mat <- connect.mat.origin
      }
      
      Nnext2 <- connect.mat %*% Nnext2 #and we disperse 
    }
    
    
    print(paste("population post-survival, post-dispersal (if applicable) =", sum(c(unlist(lapply(slice(Nnext2[index], (((L-a)*seas+1))),sum)))), sep=" "))
    
    #now track the proportion of 0s in your sub-sites within the metapop
    prop.zero <- find.prop.zero(Nnext2,a)
    prop.zero.chain[t+1] <- prop.zero
    
    #now fill in the next space in the larger matrix and step forward again
    Nmeta[,(t+1)] <-   Nnext2
    #and fill in your lambda chain, as well (will be 1 short of length of time series because it corresponds to each transition)
    #we'll leave slot one as NA to keep the ts matching
    lambda.chain[t+1] <- lambda_next
  }
  
  #first slice by site
  Nmeta.site.init <- mat_split(Nmeta, r=(L*seas), c=ncol(Nmeta)) 
  
  #then collapse to stage within each site
  Nmeta_pop <- collapse.to.stage.seasonal(Nmat=Nmeta.site.init, meta.df=df, seas=seas)
  
  #first slice by site
  if(a>=2){
    Nmeta.site <- mat_split(Nmeta_pop, r=3, c=ncol(Nmeta))   
  }
  if(a==1){
    Nmeta.site <- mat_split(Nmeta_pop, r=2, c=ncol(Nmeta))   
  }
  
  
  
  #then sum across all
  Nmeta2 <- sum.by.stage(Nmat=Nmeta.site)
  
  if(return.meta==TRUE){
    return(list(Nmeta2, Nmeta_pop, lambda.chain, prop.zero.chain)) #returns the sum across all sites by stage (Nmeta2), the joined version by each site as one long matrix (Nmeta) and the chain of averaged lambdas across all sites at the origin
  }
  if(return.meta==FALSE){
    return(list(Nmeta.site, lambda.chain, prop.zero,chain)) #gives the chain by stage by site and the chain of lambdas
  }
  
  
}
connect.leslie.sim.periodic <- function(df, N, yrs, nat.mort, hunt.override, set_hunt, dens.dep, use.uci, apply.hunt, hunt_sim, add.noise,  connect.hi, return.meta, allow.migration){
  #first, convert yrs into the timesteps of your model
  timesteps <- yrs
  
  times <- seq(1, timesteps, by = 1) 
  
  #Then make a population vector with columns to cover every timestep and
  #a row for every age class. In the case of the Lefkovitch matrix, this is always just 3
  L = unique(as.numeric(as.character(df$avg_lifespan)))
  n.sites=length(unique(df$site))
  #set K as double N
  K= 2*N
  
  Nmeta <- matrix(NA, nrow=(L*n.sites), ncol=timesteps)
  a=as.numeric(unique(df$age_1st_rep))
  #now fill the first column of that matrix in with the initial population size
  #we'll start with a stable structure
  
  df.list <- dlply(df, .(site))
  
  #split by site
  init_list <- lapply(df.list, get_Nt_init_periodic, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, N=N, K=K, hunt.override=hunt.override, set_hunt=set_hunt)
  
  
  
  Nmeta[,1] <- do.call("rbind", init_list)
  
  #if there are NAs or Infs in this vector, then we need to get rid of them
  Nmeta[,1][is.na(Nmeta[,1])] <- 0
  Nmeta[,1][Nmeta[,1]==Inf] <- 0
  
  #and make an empty vector to store your lambdas
  lambda.chain <- rep(NA, length=timesteps)
  
  #and your proportion of unoccupied patches
  prop.zero.chain <- rep(NA, length=timesteps)
  
  #fill in with your starting prop (should be 0 sites)
  prop.zero.chain[1] <- find.prop.zero(Nmeta[,1], a)
  
  
  #reset N and K at the metapop level
  # N=N*n.sites
  #K= 2*N
  
  #build your connectivity matrix here:
  #then, to prep for the sims to come, make connectivity matrix between sites
  #this connectivity matrix does not allow for dispersal between makira or 
  #masoala or from infant age classes 
  #(both these qualification can still be shaken slightly with the noise, below)
  if (allow.migration==TRUE){
    connect.mat.origin <- make.connect.mat(site.df=df, ts1=Nmeta[,1], connect.hi=connect.hi) 
  }
  
  
  #and make an empty vector to store your lambdas
  lambda.chain <- rep(NA, length=timesteps)
  
  #and then we iterate forward in time -
  #this is running the population model!
  for (t in 1:(length(times)-1)){
    
    print(paste0("t=", t))
    
    Nt_next <- Nmeta[,t]
    
    
    #and make sure that "N" changes based on adults in the age class
    index = rep(c(rep(FALSE, times=a-1), rep(TRUE, times = ((L-a)+1))), times=n.sites)
    #keep N local
    N.list <- lapply(slice(Nt_next[index], (((L-a)+1))),sum)
    
    
    print(paste("metapopulation total into dens.dep survivorship =", sum(c(unlist(N.list))), sep=" "))
    
    #and print "K" so we know it is not changing through the time series
    print(paste0("K =", K))
    
    #we now apply the survival probabilities (and density dependent fecundities) over 
    #one year (and shake as necessary)
    #first split matrix so that you have one matrix per site
    Nt_next.list <- slice(Nt_next, (L))
    
    output <- mapply(iterate.one.timestep.periodic, df=df.list, Nt_next=Nt_next.list, N=N.list, K=K, nat.mort=nat.mort, hunt.override=hunt.override, set_hunt=set_hunt, use.uci=use.uci, dens.dep=dens.dep, apply.hunt=apply.hunt, hunt_sim=hunt_sim, add.noise=add.noise, SIMPLIFY=FALSE)
    
    Nnext2 <- list()
    for (i in 1:n.sites){
      Nnext2[[i]] <- output[[i]][1]
    }
    Nnext2 <- matrix(unlist(Nnext2), nrow=length(Nt_next), ncol=1, byrow=T)
    lambda_next <- list()
    for(i in 1:n.sites){
      lambda_next[[i]] <- output[[i]][2]
    }
    #take the mean lambda across all sites
    lambda_next <- mean(c(unlist(lambda_next)))
    
    Nnext2[Nnext2<1] <- 0
    #can't have negative lemurs (but partial lemurs OK)
    
    print(paste("population post-survival, pre-dispersal =", sum(Nnext2), sep=" "))
    
    if (sum(Nnext2)==Inf){
      Nnext2 <- Nmeta[,t]
      #then print again to announce correction
      print(paste("population reached infinity. corrected, and reset to previous timestep. population, pre-dispersal =", sum(Nnext2), sep=" ")) 
    }
    
    #now shuffle lemurs based on connectivity matrix--some stay in current
    #location and some move to others
    
    #convert connect.mat to probabilities based on dispersal prob of species (and stages) in question
    #needs to have the same dimensions as entries in the Nnext2 matrix (so one locality per "age" class and a corresponding probability of disperesal)
    if (allow.migration==TRUE){
      #dispersal only happens  when it is explicitly turned on. 
      
      # we remake the connectivity matrix with each   
      #add noise to the dispersal probabilities if applicable
      if (add.noise==TRUE){
        connect.mat <-  wrap.noise(connect.mat.origin) #shake it up a little
        #we specifiy that we are shaking the original matrix because we don't want to
        #get stuck on a positivie migration trajectory and just keep building from there
        #this way, we're just shaking the original every time--sometimes this might allow for 
        #introductions and sometimes dispersals in and out of our system
      }
      else{ #no noise just gives the original connectivity matrix
        connect.mat <- connect.mat.origin
      }
      
      Nnext2 <- connect.mat %*% Nnext2 #and we disperse 
    }
    
    
    print(paste("population post-survival, post-dispersal (if applicable) =", sum(c(unlist(lapply(slice(Nnext2[index], (((L-a)+1))),sum)))), sep=" "))
    
    #now track the proportion of 0s in your sub-sites within the metapop
    prop.zero <- find.prop.zero(Nnext2,a)
    prop.zero.chain[t+1] <- prop.zero
    
    #now fill in the next space in the larger matrix and step forward again
    Nmeta[,(t+1)] <-   Nnext2
    #and fill in your lambda chain, as well (will be 1 short of length of time series because it corresponds to each transition)
    #we'll leave slot one as NA to keep the ts matching
    lambda.chain[t] <- lambda_next
   }
  #first slice by site
  Nmeta.site.init <- mat_split(Nmeta, r=(L), c=ncol(Nmeta)) 
  
  #then collapse to stage within each site
  Nmeta_pop <- collapse.to.stage(Nmat=Nmeta.site.init, meta.df=df)
  
  #now, slice again by site (each site has all three stages and is chained together in a list)
  #if(a>=2){
    Nmeta.site <- mat_split(Nmeta_pop, r=3, c=ncol(Nmeta))   
  #}
  # if(a==1){
  #   Nmeta.site <- mat_split(Nmeta_pop, r=2, c=ncol(Nmeta))   
  # }
  # 
  
  
  #then sum across all - this is the entire metapopulation added up
  Nmeta2 <- sum.by.stage(Nmat=Nmeta.site)
  
  if(return.meta==TRUE){
    return(list(Nmeta2, Nmeta_pop, lambda.chain, prop.zero.chain)) #returns the sum across all sites by stage (Nmeta2), the joined version by each site as one long matrix (Nmeta) and the chain of averaged lambdas across all sites at the origin
  }
  if(return.meta==FALSE){
    return(list(Nmeta.site, lambda.chain, prop.zero.chain)) #gives the chain by stage by site and the chain of lambdas
  }
  
  
}
connect.leslie.sim.lefkovitch <- function(df, N, yrs, use.uci, nat.mort, dens.dep, apply.hunt, hunt.intensity, hunt_sim, add.noise, connect.hi, return.meta, allow.migration) {
  
  #first, convert yrs into the timesteps of your model
  timesteps <- round(yrs/unique(as.numeric(df$IBI)),0)
  
  times <- seq(1, timesteps, by = 1) 
  
  #Then make a population vector with columns to cover every timestep and
  #a row for every age class. In the case of the Lefkovitch matrix, this is always just 3
  n.sites=length(unique(df$site))
  
  a = as.numeric(as.character(unique(df$age_1st_rep)))/as.numeric(as.character(unique(df$IBI)))
  
  if (a==1){
    Nmeta <- matrix(NA, nrow=2*n.sites, ncol=timesteps)  
  }
  if (a==2){
    Nmeta <- matrix(NA, nrow=3*n.sites, ncol=timesteps)  
  }
  if (a>=3){
    Nmeta <- matrix(NA, nrow=4*n.sites, ncol=timesteps)  
  }
  
  #set K as double N for now (will reset at the metapop below)
  K= 2*N
  
  #split into sites
  df.list <- dlply(df, .(site))
  
  init_list <- lapply(df.list, get_Nt_init_lefkovitch, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, N=N, K=K)
  
  #now bind these together and fill space 
  #now fill the first column of that matrix in with the initial population size
  #N refers just to adults. we'll fill according to the stable age structure
  
  Nmeta[,1] <- do.call("rbind", init_list)
  
  #if there are NAs or Infs in this vector, then we need to get rid of them
  Nmeta[,1][is.na(Nmeta[,1])] <- 0
  Nmeta[,1][Nmeta[,1]==Inf] <- 0
  
  #and make an empty vector to store your lambdas
  lambda.chain <- rep(NA, length=timesteps)
  
  #and your proportion of unoccupied patches
  prop.zero.chain <- rep(NA, length=timesteps)
  
  #fill in with your starting prop (should be 0 sites)
  prop.zero.chain[1] <- find.prop.zero(Nmeta[,1], a)
  
  #reset N and K at the metapop level
  #N=N*n.sites
  #K= 2*N
  
  #build your connectivity matrix here:
  #then, to prep for the sims to come, make connectivity matrix between sites
  #this connectivity matrix does not allow for dispersal between makira or 
  #masoala or from infant age classes 
  #(both these qualification can still be shaken slightly with the noise, below)
  if (allow.migration==TRUE){
    connect.mat.origin <- make.connect.mat(site.df=df, ts1=Nmeta[,1], connect.hi=connect.hi)  
  }
  
  
  
  #and then we iterate forward
  for (t in 1:(length(times)-1)){
    
    print(paste0("t=", t))
    
    Nt_next <- Nmeta[,t]
    
    #print(paste0("Ntnext=",Nt_next))
    
    #and make sure that "N" changes based on the # adults in the age class
    #that's every 3rd individual.
    #we'll keep N local
    if (a==1){
      index <- c(FALSE,TRUE)  
      N.list <- slice(Nt_next, 2)
      N.list <- sapply(N.list, '[', index)
      N.list <- as.list(N.list)
    }
    if (a==2){
      index <- c(FALSE,FALSE,TRUE)
      N.list <- slice(Nt_next, 3)
      N.list <- sapply(N.list, '[', index)
      N.list <- as.list(N.list)
    }
    if (a>=3){
      index <- c(FALSE,FALSE,TRUE,TRUE) #rep(, times=n.sites)  
      #this one is a little trickier because we want the sum of every two true entries
      N.list <- slice(Nt_next, 4) #first slice by site
      N.list <- sapply(N.list, '[', index)
      #take column sums
      N.list <- as.list(colSums(N.list))
    }
    
    #and let's track the proportion of Ns (adults only) which are 0 across all sites
    
    print(paste("total metapopulation into dens.dep survivorship =", sum(c(unlist(N.list))), sep=" "))
    
    #and print "K" so we know it is not changing through the time series
    print(paste0("K =", K))
    
    #we now apply the survival probabilities (and density dependent fecundities) over 
    #one year (and shake as necessary)
    
    #first split matrix by site
    if (a==1){
      Nt_next.list <- slice(Nt_next, 2)
    }
    if (a==2){
      Nt_next.list <- slice(Nt_next, 3)
    }
    if (a>=3){
      Nt_next.list <- slice(Nt_next, 4)
    }
    
    output <- mapply(iterate.one.timestep.lefkovitch, df=df.list, Nt_next=Nt_next.list, N=N.list, K=K, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, apply.hunt=apply.hunt, hunt.intensity=hunt.intensity, hunt_sim=hunt_sim, add.noise=add.noise, SIMPLIFY=FALSE)
    #produces a list of by site of vertical vectors by site one timestep later plus a list of lambdas
    Nnext2 <- list()
    for (i in 1:n.sites){
      Nnext2[[i]] <- output[[i]][1]
    }
    if (a==1){
      Nnext2 <- matrix(unlist(Nnext2), nrow=2*n.sites, ncol=1, byrow=T)
    }
    if (a==2){
      Nnext2 <- matrix(unlist(Nnext2), nrow=3*n.sites, ncol=1, byrow=T)
    }
    if (a>=3){
      Nnext2 <- matrix(unlist(Nnext2), nrow=4*n.sites, ncol=1, byrow=T)
    }
    lambda_next <- list()
    for(i in 1:n.sites){
      lambda_next[[i]] <- output[[i]][2]
    }
    #take the mean lambda across all sites
    lambda_next <- mean(c(unlist(lambda_next)))
    
    Nnext2[Nnext2<1] <- 0
    #can't have negative lemurs (but partial lemurs OK)
    
    print(paste("population post-survival, pre-dispersal =", sum(Nnext2), sep=" "))
    
    #if the population has reached an infinite sum by the end of the simulation, 
    #we are going to need to constrain it somehow - here we just hold it constant at 
    #whatever the timestep before, and we tell you this. 
    if (sum(Nnext2)==Inf){
      Nnext2 <- Nmeta[,t]
      #then print again to announce correction
      print(paste("population reached infinity. corrected, and reset to previous timestep. population, pre-dispersal =", sum(Nnext2), sep=" ")) 
    }
    
    #now shuffle lemurs based on connectivity matrix--some stay in current
    #location and some move to others
    
    #convert connect.mat to probabilities based on dispersal prob of species (and stages) in question
    #needs to have the same dimensions as entries in the Nnext2 matrix (so one locality per "age" class and a corresponding probability of disperesal)
    if (allow.migration==TRUE){ #& t %% (12/lts) != 0){
      #dispersal only happens in the non-breeding intervals and when it is explicitly turned on. 
      
      # we remake the connectivity matrix with each   
      #add noise to the dispersal probabilities if applicable
      if (add.noise==TRUE){
        connect.mat <-  wrap.noise(connect.mat.origin) #shake it up a little
        #we specifiy that we are shaking the original matrix because we don't want to
        #get stuck on a positivie migration trajectory and just keep building from there
        #this way, we're just shaking the original every time--sometimes this might allow for 
        #introductions and sometimes dispersals in and out of our system
      }
      else{ #no noise just gives the original connectivity matrix
        connect.mat <- connect.mat.origin
      }
      
      Nnext2 <- connect.mat %*% Nnext2 #and we disperse 
    }
    
    print(paste("population post-survival, post-dispersal (if applicable) =", sum(Nnext2), sep=" "))
    
    #fill in vector of the proportion of unoccupied sites (for all age class)
    
    #now track the proportion of 0s in your sub-sites within the metapop
    prop.zero <- find.prop.zero(Nnext2,a)
    
    prop.zero.chain[t+1] <- prop.zero
    
    #now fill in the next space in the larger matrix and step forward again
    Nmeta[,(t+1)] <-   Nnext2
    #and fill in your lambda chain, as well (will be 1 short of length of time series because it corresponds to each transition)
    #we'll leave slot one as NA to keep the ts matching
    lambda.chain[t+1] <- lambda_next
  }
  #need to sum across sites and add in connectivity
  #Nmeta gives a vertical vector with 3 stages per site
  #make a summation vector too
  
  #first slice by site
  if(a==1){
    Nmeta.site <- mat_split(Nmeta, r=2, c=ncol(Nmeta))   
  }
  if(a==2){
    Nmeta.site <- mat_split(Nmeta, r=3, c=ncol(Nmeta)) 
  }
  if(a>=3){
    Nmeta.site <- mat_split(Nmeta, r=4, c=ncol(Nmeta)) 
  }
  #then sum across all
  Nmeta2 <- sum.by.stage(Nmat=Nmeta.site)
  
  #once you've run through the time series, return the matrix and the chain of lambdas that corresponds
  
  #then, because for plotting by site you would like just 3 age classes, as well, we need to adjust Nmeta for some to reduce it
  
  if(a>=3){
    Nmeta.site.tmp <- mat_split(Nmeta, r=4, c=ncol(Nmeta)) 
    Nmeta.site.new <- collapse.adults(Nmeta.site.tmp, n.sites=n.sites)
    Nmeta <- do.call("rbind", Nmeta.site.new)
    }
  
  #you return depends on whether this is a metapop sim or just a fake replication of a single pop
  if(return.meta==TRUE){
    return(list(Nmeta2, Nmeta, lambda.chain, prop.zero.chain)) #returns the sum across all sites by stage (Nmeta2), the joined version by each site as one long matrix (Nmeta) and the chain of averaged lambdas across all sites at the origin
  }
  if(return.meta==FALSE){
    return(list(Nmeta.site, lambda.chain, prop.zero.chain)) #gives the chain by stage by site and the chain of lambdas
  }
  
}
connect.leslie.sim.wrap <- function(df, simtype, N, yrs, leslie, use.uci, nat.mort, dens.dep, apply.hunt, hunt.intensity, hunt_sim, add.noise, seas, connect.hi, return.meta, allow.migration, do.plot,plot.as.metapop, do.save){
  
  #first, pick whether this is a "best" or worst sim
   
  df1 <- subset(df, scenario==simtype)
  
  #make a list of species -- each species will travel within all sites
  
  #build a separate population matrix for every site and scenario
  df.list <- dlply(df1, .(species))
  
  #now apply the sim over each species (based on leslie type)
  
  #now apply leslie sim to each combo of species and scenario
  #produces your population simulation
  if(leslie=="lefkovitch"){
    leslie.sim.list <- lapply(df.list, connect.leslie.sim.lefkovitch, N=N, yrs=yrs, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, apply.hunt=apply.hunt, hunt.intensity=hunt.intensity, hunt_sim=hunt_sim, add.noise=add.noise, connect.hi=connect.hi, return.meta=return.meta, allow.migration=allow.migration)  
  }
  if(leslie=="annual"){
    leslie.sim.list <- lapply(df.list, connect.leslie.sim.annual, N=N, yrs=yrs, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, apply.hunt=apply.hunt, hunt.intensity=hunt.intensity, hunt_sim=hunt_sim, add.noise=add.noise, connect.hi=connect.hi, return.meta=return.meta, allow.migration=allow.migration)  
  }
  if(leslie=="seasonal"){
    leslie.sim.list <- lapply(df.list, connect.leslie.sim.seasonal, N=N, yrs=yrs, nat.mort=nat.mort,use.uci=use.uci, dens.dep=dens.dep, apply.hunt=apply.hunt, hunt.intensity=hunt.intensity, hunt_sim=hunt_sim, add.noise=add.noise, seas=seas, connect.hi=connect.hi, return.meta=return.meta, allow.migration=allow.migration)  
  }
  
  if(leslie=="periodic"){
    leslie.sim.list <- lapply(df.list, connect.leslie.sim.periodic, N=N, yrs=yrs, nat.mort=nat.mort,use.uci=use.uci, dens.dep=dens.dep, apply.hunt=apply.hunt,  hunt_sim=hunt_sim, add.noise=add.noise, connect.hi=connect.hi, return.meta=return.meta, allow.migration=allow.migration)  
  }
  
  
  
  if(do.plot==TRUE){
    
    mapply(plot.meta.pop, meta.list=leslie.sim.list, df=df.list, yrs=yrs, seas=seas, do.save=do.save, plot.as.metapop=plot.as.metapop, leslie=leslie, SIMPLIFY=FALSE)
  }
  
  #and store and return
  return(leslie.sim.list)
}
choose.par<- function(df){
  #take your min/max values for the two pars that vary and select
  df %>% mutate_if(is.factor, as.character) -> df
  df$litter_size = (df$IBI_fecundity*2)
  age.1st.rep = as.numeric(c(max(df$age_1st_rep), min(df$age_1st_rep)))
  litter.size = as.numeric(c(max(df$litter_size), min(df$litter_size)))
  
  if ((length(unique(age.1st.rep))>1) & (length(unique(litter.size))>1)){
     mat = matrix(rbind(age.1st.rep, litter.size), ncol=2)  
     par = mvrnorm(1, c(mean(age.1st.rep), mean(litter.size)), Sigma = cor(mat))
     par = round(par, 0)
     par[par<1] = 1
     #and we can't allow littersize or age 1st rep under 1
     
     
     #then add to data and return
     df$age_1st_rep = par[1]
     df$IBI_fecundity = par[2]/as.numeric(df$IBI)
     
  } else if ((length(unique(age.1st.rep))>1) & (length(unique(litter.size))==1)){
    age.rand = rnorm(n=1, mean = mean(age.1st.rep), sd= sd(age.1st.rep))
    
    age.rand = round(age.rand,0)
    age.rand[age.rand<1] = 1
    
    #then add to data and return
    df$age_1st_rep = age.rand
    df$litter_size = unique(df$litter_size)
  
  
  }else if((length(unique(age.1st.rep))==1) & (length(unique(litter.size))>1)){
    litter.rand = rnorm(n=1, mean = mean(litter.size), sd= sd(litter.size))
    litter.rand = round(litter.rand,0)
    litter.rand[litter.rand<1] = 1
    
    
    #then add to data and return
    df$age_1st_rep = unique( df$age_1st_rep)
    df$IBI_fecundity =  litter.rand/as.numeric(df$IBI)
   
    }else{
    
    df$age_1st_rep = unique(df$age_1st_rep)
    df$litter_size = unique(df$litter_size)
    
  }
  df = dplyr::select(df, -(litter_size))
  
  #and sample in correlation.
  
  df$scenario = "avg"
  df = df[!duplicated(df),]
  df$IBI_fecundity = as.numeric(df$IBI_fecundity)
  df$age_1st_rep= as.numeric(df$age_1st_rep)
  df$avg_lifespan = as.numeric(df$avg_lifespan)
  df$IBI = as.numeric(df$IBI)
  df$IBI_infant_mortality = as.numeric(df$IBI_infant_mortality)
  df$IBI_adult_mortality = as.numeric(df$IBI_adult_mortality)
  df$lat = as.numeric(df$lat)
  df$long = as.numeric(df$long)
  df$site = as.numeric(df$site)
  
  print(unique(df$species))
  print(paste0("IBI_fecundity = ", unique(df$IBI_fecundity)))
  print(paste0("age 1st rep = ", unique(df$age_1st_rep)))
  return(df)
}
connect.leslie.sim.wrap.draw <- function(df, N, yrs, hunt.override, set_hunt, leslie, use.uci, nat.mort, dens.dep, apply.hunt, hunt.intensity, hunt_sim, add.noise, seas, connect.hi, return.meta, allow.migration, do.plot,plot.as.metapop, do.save){
  print("connect.leslie.sim.wrap.draw")
  print(df)
  #make a list of species, then choose your parameters for each
  df.list <- dlply(df, .(species))
  print(df.list)
  #randomly draw your parameters from the range of parameter values...and run enough sims to encapsulate that range distribution
  #parameters that have ranges are: litter size and age at 1st reproduction. everything else is fixed
  
  
  df.list1 <- lapply(df.list, choose.par) #each time, will randomly draw from normal distribution using the range that we fed it in the data table
  print(df.list1)

  
  #now apply the sim over each species (based on the chosen structure of the transition matrix--this is coded as "leslie")
  
  #now apply leslie sim to each combo of species and scenario
  #produces your population simulation
  if(leslie=="lefkovitch"){
    leslie.sim.list <- lapply(df.list1, connect.leslie.sim.lefkovitch, N=N, yrs=yrs, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, apply.hunt=apply.hunt, hunt.intensity=hunt.intensity, hunt_sim=hunt_sim, add.noise=add.noise, connect.hi=connect.hi, return.meta=return.meta, allow.migration=allow.migration)  
  }
  if(leslie=="annual"){
    leslie.sim.list <- lapply(df.list1, connect.leslie.sim.annual, N=N, yrs=yrs, nat.mort=nat.mort, use.uci=use.uci, dens.dep=dens.dep, apply.hunt=apply.hunt, hunt.intensity=hunt.intensity, hunt_sim=hunt_sim, add.noise=add.noise, connect.hi=connect.hi, return.meta=return.meta, allow.migration=allow.migration)  
  }
  if(leslie=="seasonal"){
    leslie.sim.list <- lapply(df.list1, connect.leslie.sim.seasonal, N=N, yrs=yrs, nat.mort=nat.mort,use.uci=use.uci, dens.dep=dens.dep, apply.hunt=apply.hunt, hunt.intensity=hunt.intensity, hunt_sim=hunt_sim, add.noise=add.noise, seas=seas, connect.hi=connect.hi, return.meta=return.meta, allow.migration=allow.migration)  
  }
  
  if(leslie=="periodic"){
    leslie.sim.list <- lapply(df.list1, connect.leslie.sim.periodic, N=N, yrs=yrs, hunt.override=hunt.override, set_hunt=set_hunt, nat.mort=nat.mort,use.uci=use.uci, dens.dep=dens.dep, apply.hunt=apply.hunt,  hunt_sim=hunt_sim, add.noise=add.noise, connect.hi=connect.hi, return.meta=return.meta, allow.migration=allow.migration)  
  }
  
  
  
  if(do.plot==TRUE){
    
    mapply(plot.meta.pop, meta.list=leslie.sim.list, df=df.list, yrs=yrs, seas=seas, do.save=do.save, plot.as.metapop=plot.as.metapop, leslie=leslie, SIMPLIFY=FALSE)
  }
  
  #and store and return
  return(leslie.sim.list)
}
merge.sites <- function(meta.list){
  time.series <- meta.list[[1]] #now a strange vector thing with a time series for every site  
  
  n.sites=dim(meta.list[[1]])[3]
  
  #average all stages across all sites
  n.stage=nrow(time.series[,,1])
  
  if(n.stage==3){
    infants <- list()
    juveniles <- list()
    adults <- list()
    
    for (i in 1:n.sites){
      infants[[i]] <- time.series[1,,i] #gives a list for each site that is a string of vectors of infants
      juveniles[[i]] <- time.series[2,,i]#gives a list of vectors of juveniles
      adults[[i]] <- time.series[3,,i] #gives a list of vectors of adults
    }
    
    infants <- do.call("rbind", infants) #gives TS of count per class for all sites. later stick together and take mean of them all
    juveniles <- do.call("rbind", juveniles)
    adults <- do.call("rbind", adults)
    
    return(list(infants, juveniles, adults)) 
  }
  
  if(n.stage==2){
    infants <- list()
    adults <- list()
    
    for (i in 1:n.sites){
      infants[[i]] <- time.series[1,,i] #gives a list for each site that is a string of vectors of infants
      adults[[i]] <- time.series[2,,i] #gives a list of vectors of adults
    }
    
    infants <- do.call("rbind", infants) #gives TS of count per class for all sites. later stick together and take mean of them al
    adults <- do.call("rbind", adults)
    
    return(list(infants, adults)) 
  }
  
}
merge.sites.metapop <- function(meta.list){
  sum.ts <- meta.list[[1]] #time series across the metapopulation 
  site.ts <- meta.list[[2]] 
  n.stage <- nrow(sum.ts)
  n.sites <- nrow(site.ts)/n.stage
  
return(list(sum.ts, site.ts))
} 
group.spp <- function(spp.list, reps, n.species){
  #every n.sites the species repeat
  #make an index for each species
  #if cleverer, would here use paste to make 7 empty index vectors
  index1 <- rep(c(TRUE, rep(FALSE, times=(n.species-1))), reps)
  index2 <- rep(c(FALSE, TRUE, rep(FALSE, times=(n.species-2))), reps)
  index3 <- rep(c(FALSE, FALSE, TRUE, rep(FALSE, times=(n.species-3))), reps)
  index4 <- rep(c(FALSE, FALSE, FALSE, TRUE, rep(FALSE, times=(n.species-4))), reps)
  index5 <- rep(c(rep(FALSE,4), TRUE, FALSE, FALSE), reps)
  index6 <- rep(c(rep(FALSE,5), TRUE,  FALSE), reps)
  index7 <- rep(c(rep(FALSE,6), TRUE), reps)
  
  index.list <- list(index1,index2,index3,index4,index5,index6,index7)
  spp.list.new <- list()
  
  for (i in 1:n.species){
    spp.list.new[[i]] <- spp.list[index.list[[i]]]  
  }
  
  return(spp.list.new)
}
merge.stage <- function(single.spp.list){
  #get all of every class
  n.stage <- length(single.spp.list[[1]])
  if(n.stage==3){
    infants <- do.call("rbind", sapply(single.spp.list, '[', 1)) #gets all infants
    juveniles <- do.call("rbind", sapply(single.spp.list, '[', 2)) #gets all juveniles
    adults <- do.call("rbind", sapply(single.spp.list, '[', 3)) #gets all adults
    names.vector <- c("infants", "juveniles", "adults")
    
    #and then report the average of each timestep and the max and min
    mean.ts.inf <- colMeans(infants)
    mean.ts.juv <- colMeans(juveniles)
    mean.ts.adults <- colMeans(adults)
    mean.ts <- rbind(mean.ts.inf,mean.ts.juv, mean.ts.adults)
    rownames(mean.ts) <- paste("mean", names.vector, sep="-")
    
    max.ts.inf <- colMaxs(infants)
    max.ts.juv <- colMaxs(juveniles)
    max.ts.adults <- colMaxs(adults)
    max.ts <- rbind(max.ts.inf,max.ts.juv, max.ts.adults)
    rownames(max.ts) <- paste("max", names.vector, sep="-")
    
    min.ts.inf <- colMins(infants)
    min.ts.juv <- colMins(juveniles)
    min.ts.adults <- colMins(adults)
    min.ts <- rbind(min.ts.inf,min.ts.juv, min.ts.adults)
    rownames(min.ts) <- paste("min", names.vector, sep="-")
  }
  
  if(n.stage==2){
    infants <- do.call("rbind", sapply(single.spp.list, '[', 1)) #gets all infants
    adults <- do.call("rbind", sapply(single.spp.list, '[', 2)) #gets all adults
    names.vector <- c("infants", "adults")
    
    mean.ts.inf <- colMeans(infants)
    mean.ts.adults <- colMeans(adults)
    mean.ts <- rbind(mean.ts.inf, mean.ts.adults)
    rownames(mean.ts) <- paste("mean", names.vector, sep="-")
    
    max.ts.inf <- colMaxs(infants)
    max.ts.adults <- colMaxs(adults)
    max.ts <- rbind(max.ts.inf, max.ts.adults)
    rownames(max.ts) <- paste("max", names.vector, sep="-")
    
    min.ts.inf <- colMins(infants)
    min.ts.adults <- colMins(adults)
    min.ts <- rbind(min.ts.inf, min.ts.adults)
    rownames(min.ts) <- paste("min", names.vector, sep="-")   
  }
  
  #bind and return
  sum.ts <- rbind(mean.ts,max.ts,min.ts)
  return(sum.ts)
}
merge.stage.metapop <- function(single.spp.list){
  
  metapops <- sapply(single.spp.list, '[', 1) #list of all the different metapop trials
  site.df <- sapply(single.spp.list, '[', 2)
  
  n.stage <- nrow(metapops[[1]])
  n.sites <- (nrow(site.df[[1]]))/n.stage 
  
  #get the mean, min, max
  metapops.mean <- Reduce("+", metapops) / length(metapops)
  metapops.max <- do.call("pmax", metapops)
  metapops.min <- do.call("pmin", metapops)
  rownames(metapops.mean) <- paste(rownames(metapops.mean), "mean", sep="-")
  rownames(metapops.max) <- paste(rownames(metapops.max), "max", sep="-")
  rownames(metapops.min) <- paste(rownames(metapops.min), "min", sep="-")
  
  sum.ts <- rbind(metapops.mean,metapops.max,metapops.min) #this is the mean, max, and min for the whole metapop
  
  #now we also want the same for all of the little sites
  #take the mean element-wise across your reps
  site.df.mean <- Reduce("+", site.df) / length(site.df)
  site.df.max <- do.call("pmax", site.df)
  site.df.min <- do.call("pmin", site.df)
  
  if(n.stage==3){
  rownames(site.df.mean) <- paste(rep(c("infants", "juveniles", "adults"), times=n.sites), "mean", "site", rep(seq(1,n.sites,1), each=3), sep="-")
  rownames(site.df.max) <- paste(rep(c("infants", "juveniles", "adults"), times=n.sites), "max", "site", rep(seq(1,n.sites,1), each=3), sep="-")
  rownames(site.df.min) <- paste(rep(c("infants", "juveniles", "adults"), times=n.sites), "min", "site", rep(seq(1,n.sites,1), each=3), sep="-")
  }
  
  if(n.stage==2){
  rownames(site.df.mean) <- paste(rep(c("infants", "adults"), times=n.sites), "mean", "site", rep(seq(1,n.sites,1), each=2), sep="-")
  rownames(site.df.max) <- paste(rep(c("infants",  "adults"), times=n.sites), "max", "site", rep(seq(1,n.sites,1), each=2), sep="-")
  rownames(site.df.min) <- paste(rep(c("infants",  "adults"), times=n.sites), "min", "site", rep(seq(1,n.sites,1), each=2), sep="-")
  }
  site.ts <- rbind(site.df.mean,site.df.max, site.df.min)
  
  return(list(sum.ts, site.ts))
  }
average.reps <- function(replicate.list, meta.df, reps, yrs, return.avg){
  n.species=length(unique(meta.df$species))
  
  #first we average all the sites within each run (specific only to a non-metapop context)
  sim.group <- lapply(replicate.list, merge.sites)
  #now each entry is a separate species and within that list, there is a list of all reps
  #then lists within that for each stage class. 
  
  #we reorganize by species - note that this code requires a 7 species entry. would need to 
  #edit for flexibility
  list.by.spp <- group.spp(spp.list=sim.group, reps=reps, n.species=n.species)
  
  if(return.avg==TRUE){
    #now we summarize to produce the mean, max and min per timestep per stage per species
    sum.lists <- lapply(list.by.spp, merge.stage)
    
    #and name by species
    names(sum.lists) <- unique(meta.df$species)
    
    #and return - to plot later
    return(sum.lists)
  }
  if(return.avg==FALSE){
    #you want a list of dfs for all the species with every run a separate trial
    trial.lists <- mapply(trial.group, one.spp.list=list.by.spp,name=as.list(unique(meta.df$species)), yrs=yrs, SIMPLIFY=FALSE)
    
    names(trial.lists) <- unique(meta.df$species)
    #then you return the list by species of dataframes of all the trials and plot them together
    return(trial.lists)
  }   
} 
average.reps.metapop <- function(replicate.list, meta.df, reps, yrs, return.avg){
  n.species=length(unique(meta.df$species))
  sim.group <- lapply(replicate.list, merge.sites.metapop)
  #we reorganize by species - note that this code requires a 7 species entry. would need to 
  #edit for flexibility
  list.by.spp <- group.spp(spp.list=sim.group, reps=reps, n.species=n.species)
  #now have a 7 item list, and within each item, we have a list of time series
  
  if(return.avg==TRUE){
    #now we summarize to produce the mean, max and min per timestep per stage per species
    
    #first group them together to produce an average of each metapop estimate per timestep
    #and the same for each subsite
    sum.lists <- lapply(list.by.spp, merge.stage.metapop)
    #gives a list the length of species, each with two time series: the metapop, and the site list
    
    #and name by species
    names(sum.lists) <- unique(meta.df$species)
    
    #and return - to plot later
    return(sum.lists)
  }
  if(return.avg==FALSE){
    
    names(list.by.spp) <- unique(meta.df$species)
    
    #and return - to plot later
    return(list.by.spp)
  }   
} 
make.df <- function(df, name, yrs){
  timesteps <- ncol(df)
  yr.per.timestep <- yrs/timesteps
  sum.df <- data.frame(df)
  sum.df <- cbind.data.frame(rownames(sum.df), sum.df)
  colnames(sum.df) <- c("stage", paste0("timestep_", seq(1, (timesteps), 1)))
  rownames(sum.df) <- c()
  
  sum.df2 <- reshape(sum.df, varying=paste0("timestep_",seq(1, (timesteps), 1)), v.names="count", timevar="timestep", times=seq(1,(timesteps),1), direction="long")
  #add type column
  sum.df2$type <- sapply(strsplit(as.character(sum.df2$stage), "-"), '[', 1)
  sum.df2$stage <- sapply(strsplit(as.character(sum.df2$stage), "-"), '[', 2)
  
  sum.df3 <- subset(sum.df2, type=="mean")
  sum.df4 <- subset(sum.df2, type=="min")
  sum.df5 <- subset(sum.df2, type=="max")
  names(sum.df3)[3] <- "mean_count"
  names(sum.df4)[3] <- "min_count"
  names(sum.df5)[3] <- "max_count"
  #drop the "type" and "id" columns
  sum.df3[,ncol(sum.df3)] <- c()
  sum.df3[,ncol(sum.df3)] <- c()
  
  #and remerge
  sum.df3$min_count <- sum.df4[,3]
  sum.df3$max_count <- sum.df5[,3]
  
  #and attach species name
  sum.df3$species <- name
  
  #and attach the year for each timestep
  sum.df3$timestep <- as.numeric(as.character(sum.df3$timestep)) - 1
  sum.df3$year <- as.numeric(as.character(sum.df3$timestep))*yr.per.timestep
  
  
  return(sum.df3)
}
make.df.meta <- function(df, name, yrs){
  timesteps <- ncol(df)
  yr.per.timestep <- yrs/timesteps
  
  sum.df <- data.frame(df)
  sum.df <- cbind.data.frame(rownames(sum.df), sum.df)
  colnames(sum.df) <- c("title", paste0("timestep_", seq(1, (timesteps), 1)))
  rownames(sum.df) <- c()
  
  sum.df2 <- reshape(sum.df, varying=paste0("timestep_",seq(1, (timesteps), 1)), v.names="count", timevar="timestep", times=seq(1,(timesteps),1), direction="long")
  #add type column
  sum.df2$type <- sapply(strsplit(as.character(sum.df2$title), "-"), '[', 2)
  sum.df2$stage <- sapply(strsplit(as.character(sum.df2$title), "-"), '[', 1)
  
  sum.df3 <- subset(sum.df2, type=="mean")
  sum.df4 <- subset(sum.df2, type=="min")
  sum.df5 <- subset(sum.df2, type=="max")
  names(sum.df3)[3] <- "mean_count"
  names(sum.df4)[3] <- "min_count"
  names(sum.df5)[3] <- "max_count"
  #drop the "type" and "id" columns
  sum.df3[,(ncol(sum.df3)-1)] <- c()
  sum.df3[,(ncol(sum.df3)-1)] <- c()
  
  #and remerge
  sum.df3$min_count <- sum.df4[,3]
  sum.df3$max_count <- sum.df5[,3]
  
  #and attach species name
  sum.df3$species <- name
  
  #and attach the year for each timestep
  sum.df3$timestep <- as.numeric(as.character(sum.df3$timestep)) - 1
  sum.df3$year <- as.numeric(as.character(sum.df3$timestep))*yr.per.timestep
  sum.df3$title <- paste(sum.df3$species, sum.df3$title, sep="-")
  
  
  
  return(sum.df3)
}
make.df.meta.site <- function(df, name, yrs){
  timesteps <- ncol(df)
  yr.per.timestep <- yrs/timesteps
  sum.df <- data.frame(df)
  
  #does not always have rownames, so add them here
  sum.df <- cbind.data.frame(rownames(sum.df), sum.df)
  colnames(sum.df) <- c("title", paste0("timestep_", seq(1, (timesteps), 1)))
  rownames(sum.df) <- c()
  
  sum.df2 <- reshape(sum.df, varying=paste0("timestep_",seq(1, (timesteps), 1)), v.names="count", timevar="timestep", times=seq(1,(timesteps),1), direction="long")
  #add type and site column
  sum.df2$type <- sapply(strsplit(as.character(sum.df2$title), "-"), '[', 2)
  sum.df2$site <- sapply(strsplit(as.character(sum.df2$title), "-"), '[', 4)
  sum.df2$stage <- sapply(strsplit(as.character(sum.df2$title), "-"), '[', 1)
  
  
  sum.df3 <- subset(sum.df2, type=="mean")
  sum.df4 <- subset(sum.df2, type=="min")
  sum.df5 <- subset(sum.df2, type=="max")
  names(sum.df3)[3] <- "mean_count"
  names(sum.df4)[3] <- "min_count"
  names(sum.df5)[3] <- "max_count"
  #drop the "type" and "id" columns
  sum.df3[,5] <- c()
  sum.df3[,4] <- c()
  
  #and remerge
  sum.df3$min_count <- sum.df4[,3]
  sum.df3$max_count <- sum.df5[,3]
  
  #and attach species name
  sum.df3$species <- unlist(name)
  
  #and attach the year for each timestep
  sum.df3$timestep <- as.numeric(as.character(sum.df3$timestep)) - 1
  sum.df3$year <- as.numeric(as.character(sum.df3$timestep))*yr.per.timestep
  sum.df3$title <- paste(sum.df3$species, sum.df3$title, sep="-")
  
  
  return(sum.df3)
}
make.df.meta.trial <- function(df, name, yrs){
  meta.list <- sapply(df, '[', 1)
  site.list <- sapply(df, '[', 2)
  
  meta.df <- do.call("rbind", meta.list)
  site.df <- do.call("rbind", site.list)
  
  n.stage=nrow(meta.list[[1]])
  n.sites=nrow(site.list[[1]])/n.stage
  reps=length(meta.list)
  
  timesteps <- ncol(meta.df)
  yr.per.timestep <- yrs/timesteps
  
  rownames(meta.df) <- paste(rownames(meta.df), "trial", rep(seq(1,reps,1), each=n.stage), sep="-")
  meta.df <- cbind.data.frame(rownames(meta.df), meta.df)
  colnames(meta.df) <- c("title", paste0("timestep_", seq(1, (timesteps), 1)))
  rownames(meta.df) <- c()
  
  meta.df2 <- reshape(meta.df, varying=paste0("timestep_",seq(1, (timesteps), 1)), v.names="count", timevar="timestep", times=seq(1,(timesteps),1), direction="long")
  #add trial column
  meta.df2$trial <- sapply(strsplit(as.character(meta.df2$title), "-"), '[', 3)
  meta.df2$stage <- sapply(strsplit(as.character(meta.df2$title), "-"), '[', 1)
  meta.df2$species <- name
  meta.df2$timestep <- as.numeric(as.character(meta.df2$timestep)) -1
  meta.df2$year <- as.numeric(as.character(meta.df2$timestep))*yr.per.timestep
  meta.df2$title <- paste(meta.df2$species, meta.df2$title, sep="-")
  
  #and the same for the site.df
  rownames(site.df) <- paste(rep(c(unique(meta.df2$stage)), times=n.sites),"site", rep(seq(1,n.sites,1), times=n.stage), "trial", rep(seq(1,reps,1), each=(n.stage*n.sites)), sep="-")
  site.df <- cbind.data.frame(rownames(site.df), site.df)
  colnames(site.df) <- c("title", paste0("timestep_", seq(1, (timesteps), 1)))
  rownames(site.df) <- c()
  
  site.df2 <- reshape(site.df, varying=paste0("timestep_",seq(1, (timesteps), 1)), v.names="count", timevar="timestep", times=seq(1,(timesteps),1), direction="long")
  #add trial column
  site.df2$trial <- sapply(strsplit(as.character(site.df2$title), "-"), '[', 5)
  site.df2$site <- sapply(strsplit(as.character(site.df2$title), "-"), '[', 3)
  site.df2$stage <- sapply(strsplit(as.character(site.df2$title), "-"), '[', 1)
  site.df2$species <- name
  site.df2$timestep <- as.numeric(as.character(site.df2$timestep)) - 1
  site.df2$year <- as.numeric(as.character(site.df2$timestep))*yr.per.timestep 
  site.df2$title <- paste(site.df2$species, site.df2$title, sep="-")
  
  meta.df2 <- arrange(meta.df2, species, year)
  site.df2 <- arrange(site.df2, species, year)
  
  return(c(list(meta.df2, site.df2)))
}
trial.group <- function(one.spp.list, name, yrs){
  n.stage <- length(one.spp.list[[1]])
  if(n.stage==3){
    infants <- data.frame(do.call("rbind", sapply(one.spp.list, "[", 1)))
    juveniles <- data.frame(do.call("rbind", sapply(one.spp.list, "[", 2)))
    adults <- data.frame(do.call("rbind", sapply(one.spp.list, "[", 3)))
    
    timesteps <- ncol(infants)
    trial.num <- as.numeric(as.character(nrow(infants)))
    
    
    rownames(infants) <- paste(rep("infants", trial.num), "trial", seq(1,trial.num, by=1), sep="-")
    rownames(juveniles) <- paste(rep("juveniles", trial.num), "trial", seq(1,trial.num, by=1), sep="-")
    rownames(adults) <- paste(rep("adults", trial.num), "trial", seq(1,trial.num, by=1), sep="-")
    
    tot.df <- rbind(infants, juveniles, adults)
  }
  
  if(n.stage==2){
    infants <- data.frame(do.call("rbind", sapply(one.spp.list, "[", 1)))
    adults <- data.frame(do.call("rbind", sapply(one.spp.list, "[", 2)))
    
    timesteps <- ncol(infants)
    trial.num <- as.numeric(as.character(nrow(infants)))
    
    rownames(infants) <- paste(rep("infants", trial.num), "trial", seq(1,trial.num, by=1), sep="-")
    rownames(adults) <- paste(rep("adults", trial.num), "trial", seq(1,trial.num, by=1), sep="-")
    
    tot.df <- rbind(infants, adults)
  }
  
  tot.df <- cbind.data.frame(rownames(tot.df), tot.df)
  colnames(tot.df) <- c("stage", paste0("timestep_", seq(1, (timesteps), 1)))
  rownames(tot.df) <- c()
  
  tot.df2 <- reshape(tot.df, varying=paste0("timestep_",seq(1, (timesteps), 1)), v.names="count", timevar="timestep", times=seq(1,(timesteps),1), direction="long")
  #remove ID
  tot.df2[, ncol(tot.df2)] <- c()
  
  #add columns
  tot.df2$trial <- as.character(sapply(strsplit(as.character(tot.df2$stage), "-"), '[', 3))
  tot.df2$stage <- sapply(strsplit(as.character(tot.df2$stage), "-"), '[', 1)
  
  #and attach species name
  tot.df2$species <- name
  
  #and year
  yr.per.timestep <- yrs/timesteps
  tot.df2$timestep <- as.numeric(as.character(tot.df2$timestep)) - 1
  tot.df2$year <- as.numeric(as.character(tot.df2$timestep))*(yr.per.timestep)
  
  rownames(tot.df2) <- c()
  
  #and return
  return(tot.df2)  
}
plot.avg <- function(sum.list, yrs, do.save, filename){
  #bind list
  #some have uneven lengths because we fed it years and it converted to timesteps
  #first group them all as a data.frame
  
  df.list<-  mapply(make.df, df=sum.list, name=as.list(names(sum.list)), yrs=yrs, SIMPLIFY=FALSE)
  #and bind
  sum.df <- do.call("rbind", df.list)
  rownames(sum.df) <- c()
  
  #sort by species, then timestep
  sum.df <- arrange(sum.df, species, timestep)
  
  #and make a total
  
  #and assign colors
  colourCount <- 11 #length minus the titles
  getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  color.vect <- getPalette(colourCount)
  #and how do these 8 fit into the other color vector?
  spp.index <- c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
  color.vect <- color.vect[spp.index]
  
  p <- ggplot() +
    geom_line(data=sum.df, aes(x=year, y=mean_count, linetype=stage, color=species), size=1.5) + scale_colour_manual(values=color.vect) + 
    geom_ribbon(data=sum.df, aes(x=year, ymin=min_count,ymax=max_count, linetype=stage, fill=species), size=.5, alpha=.3) + scale_fill_manual(values=color.vect) + 
    theme_bw() + xlab("years") + ylab("# lemurs") #
  
  print(p)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=60, 
           height=70, 
           scale=3, 
           dpi=300)
    
  }
  
  
}  
plot.avg.metapop <- function(sum.list, yrs, plot.meta, spp.facet, do.save, filename){
  #bind list
  #some have uneven lengths because we fed it years and it converted to timesteps
  #first group them all as a data.frame
  meta.list <- sapply(sum.list, '[', 1)
  site.list <- sapply(sum.list, '[', 2)
  
  df.list <-  mapply(make.df.meta, df=meta.list, name=as.list(names(sum.list)), yrs=yrs, SIMPLIFY=FALSE)
  #and bind
  sum.df <- do.call("rbind", df.list)
  rownames(sum.df) <- c()
  #and add "site"
  sum.df$site <- "metapop"
  
  #and add in the sub-sites
  site.df.list <-  mapply(make.df.meta.site, df=site.list, name=as.list(names(sum.list)), yrs=yrs, SIMPLIFY=FALSE)
  #and bind
  site.df <- do.call("rbind", site.df.list)
  rownames(site.df) <- c()
  
  #sort by species, then year
  sum.df <- arrange(sum.df, species, year)
  site.df <- arrange(site.df, species, year)
  #and make a total
  
  #and assign colors
  colourCount <- 11 #length minus the titles
  getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  color.vect <- getPalette(colourCount)
  #and how do these 8 fit into the other color vector?
  spp.index <- c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
  color.vect <- color.vect[spp.index]
  
  max.y <- max(sum.df$max_count)
  max.y.site <- max(site.df$max_count)
  
  new_names <- c(
    'Avahi laniger'="Avahi\nlaniger",
    'Cheirogaleus sp.'="Cheirogaleus\nsp.",
    'Daubentonia madagascariensis'="Daubentonia\nmadagascariensis",
    'Eulemur albifrons'="Eulemur\nalbifrons",
    'Eulemur rubriventer'="Eulemur\nrubriventer",
    'Hapalemur griseus'="Hapalemur\ngriseus",
    'Indri indri'="Indri\nindri",
    'Lepilemur sp.'="Lepilemur\nsp.",
    'Microcebus sp.'="Microcebus\nsp.",
    'Varecia rubra'="Varecia\nrubra",
    'Varecia variegata'="Varecia\nvariegata"
  )
  
  
  if(plot.meta==TRUE){
    p <- ggplot() +
      geom_line(data=sum.df, aes(x=year, y=mean_count, linetype=stage, color=species, group=title), size=1.2) + scale_colour_manual(values=color.vect) + 
      geom_ribbon(data=sum.df, aes(x=year, ymin=min_count,ymax=max_count, linetype=stage, fill=species, group=title), size=.5, alpha=.3) + scale_fill_manual(values=color.vect) + 
      geom_line(data=site.df, aes(x=year, y=mean_count, linetype=stage, color=species, group=title), size=.5, show.legend=FALSE) + scale_alpha_discrete(guide="none") + 
      geom_ribbon(data=site.df, aes(x=year, ymin=min_count,ymax=max_count, linetype=stage, fill=species, group=title), alpha=.2, show.legend=FALSE) +
      theme_classic() + xlab("years") + ylab("# lemurs") #
    
    print(p) 
    
    if (spp.facet==TRUE){
      p <- ggplot() +
        geom_line(data=sum.df, aes(x=year, y=mean_count, linetype=stage, group=title), color="black", size=1.2) +# scale_colour_manual(values=color.vect) + 
        geom_ribbon(data=sum.df, aes(x=year, ymin=min_count,ymax=max_count, linetype=stage,  group=title), fill="black", alpha=.3) + #scale_fill_manual(values=color.vect) + 
        geom_line(data=site.df, aes(x=year, y=mean_count, linetype=stage, color=site, group=title), size=.5) + scale_alpha_discrete(guide="none") + 
        geom_ribbon(data=site.df, aes(x=year, ymin=min_count,ymax=max_count, linetype=stage, fill=site, group=title), alpha=.2, show.legend=FALSE) +
        xlab("years") + ylab("# lemurs") + facet_grid(species~., labeller = labeller(species=new_names)) + theme_classic() + theme(strip.text = element_text(face = "italic")) +
        geom_segment(data=sum.df, aes(x=0,xend=0,y=0, yend=max.y), colour="black", size=.1) +
        geom_segment(data=sum.df, aes(x=0,xend=yrs,y=0, yend=0), colour="black", size=.1) 
      
      print(p) 
    }
  }
  if(plot.meta==FALSE){
    site.df$stage <- factor(site.df$stage, levels=c("infants", "juveniles", "adults"))
    p <- ggplot() +
      geom_line(data=site.df, aes(x=year, y=mean_count, linetype=stage, color=species, group=title), size=.6) + scale_colour_manual(values=color.vect) +
      geom_ribbon(data=site.df, aes(x=year, ymin=min_count,ymax=max_count, fill=species, group=title), alpha=.3, show.legend=FALSE) +
      theme_classic() + xlab("years") + ylab("# lemurs") #
    
    print(p) 
    
    if (spp.facet==TRUE){
      p <- ggplot() +
        geom_line(data=site.df, aes(x=year, y=mean_count, linetype=stage, color=site, group=title), size=.6) + #scale_colour_manual(values=color.vect) +
        geom_ribbon(data=site.df, aes(x=year, ymin=min_count,ymax=max_count, fill=site, group=title), alpha=.3, show.legend=FALSE) +
        xlab("years") + ylab("# lemurs") + facet_grid(species~., labeller = labeller(species=new_names)) + theme_classic() + theme(strip.text = element_text(face = "italic")) +
        geom_segment(data=sum.df, aes(x=0,xend=0,y=0, yend=max.y.site), colour="black", size=.1) +
        geom_segment(data=sum.df, aes(x=0,xend=yrs,y=0, yend=0), colour="black", size=.1) 
      
      print(p)
    }
    
    
  }
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=60, 
           height=70, 
           scale=3, 
           dpi=300)
    
  }
  
  
}  
plot.trial <- function(trial.list, do.save, filename){
  trial.df <- do.call("rbind",trial.list)
  rownames(trial.df) <- c()
  
  colourCount <- 11 #length minus the titles
  getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  color.vect <- getPalette(colourCount)
  #and how do these 8 fit into the other color vector?
  spp.index <- c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
  color.vect <- color.vect[spp.index]
  
  p <- ggplot() +
    geom_line(data=trial.df, aes(x=year, y=count, linetype=stage, color=species,  alpha=trial), size=.8) + scale_colour_manual(values=color.vect) + scale_alpha_discrete(guide="none")+
    theme_bw() + xlab("years") + ylab("# lemurs") #
  
  print(p)
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=60, 
           height=70, 
           scale=3, 
           dpi=300)
  }
}
get.avg <- function(df, name, yrs){
  mean <- Reduce('+', df)/length(df)
  min <- do.call("pmin", df)
  max <- do.call("pmax", df)
  
  tmp <- data.frame(rbind(mean, min, max))
  rownames(tmp) <- c("mean", "min", "max")
  colnames(tmp) <- paste0("timestep_",seq(1,length(mean),1))
  tmp <- cbind.data.frame(rownames(tmp), tmp)
  names(tmp)[1] <- "type"
  
  #then reshape
  timesteps <- length(mean)
  tmp2 <- reshape(tmp, varying=paste0("timestep_",seq(1, (timesteps), 1)), v.names="prop", timevar="timestep", times=seq(1,(timesteps),1), direction="long")
  
  #drop id and rownames
  tmp2[,ncol(tmp2)] <- c()
  rownames(tmp2) <- c()
  tmp2$species=name
  tmp2$timestep <- as.numeric(as.character(tmp2$timestep)) - 1
  tmp2$year <- tmp2$timestep*(yrs/timesteps)
  
  tmp3 <- subset(tmp2, type=="mean")
  tmp4 <- subset(tmp2, type=="min")
  tmp5 <- subset(tmp2, type=="max")
  tmp3[,1] <- c()
  names(tmp3)[2] <- "mean_prop"
  tmp3$min_prop <- tmp4$prop
  tmp3$max_prop <- tmp5$prop
  
  return(tmp3)
}
plot.trial.metapop <- function(trial.list, yrs, do.save, filename, plot.meta, spp.facet){
  
  df.list <- mapply(make.df.meta.trial,df=trial.list,name=as.list(names(trial.list)),yrs=yrs, SIMPLIFY=FALSE)
  
  sum.list <- sapply(df.list, '[',1)
  site.list <- sapply(df.list, '[',2)
  sum.df <- do.call("rbind", sum.list)
  site.df <- do.call("rbind", site.list)
  rownames(sum.df) <- rownames(site.df) <- c()
  sum.df$site <- "metapop"
  
  colourCount <- 11 #length minus the titles
  getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  color.vect <- getPalette(colourCount)
  #and how do these 8 fit into the other color vector?
  spp.index <- c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
  color.vect <- color.vect[spp.index]
  
  #first make a vector to rename things
  new_names <- c(
    'Avahi laniger'="Avahi\nlaniger",
    'Cheirogaleus sp.'="Cheirogaleus\nsp.",
    'Daubentonia madagascariensis'="Daubentonia\nmadagascariensis",
    'Eulemur albifrons'="Eulemur\nalbifrons",
    'Eulemur rubriventer'="Eulemur\nrubriventer",
    'Hapalemur griseus'="Hapalemur\ngriseus",
    'Indri indri'="Indri\nindri",
    'Lepilemur sp.'="Lepilemur\nsp.",
    'Microcebus sp.'="Microcebus\nsp.",
    'Varecia rubra'="Varecia\nrubra",
    'Varecia variegata'="Varecia\nvariegata"
  )
  
  max.y <- max(sum.df$count)
  max.y.site <- max(site.df$count)
  
  
  if(plot.meta==TRUE){
  p <- ggplot() +
    geom_line(data=sum.df, aes(x=year, y=count, linetype=stage, color=species,  alpha=trial, group=title), size=.8) + scale_colour_manual(values=color.vect) + scale_alpha_discrete(guide="none")+
    geom_line(data=site.df, aes(x=year, y=count, linetype=stage, color=species, alpha=trial,  group=title), size=.4) +
    theme_classic() + xlab("years") + ylab("# lemurs") #
  print(p)
  
  if (spp.facet==TRUE){
    
    p <- ggplot() +
      geom_line(data=sum.df, aes(x=year, y=count, linetype=stage, alpha=trial, group=title), color="black", size=.8) + scale_alpha_discrete(guide="none")+
      geom_line(data=site.df, aes(x=year, y=count, linetype=stage, color=site, alpha=trial,  group=title), size=.4) +
      xlab("years") + ylab("# lemurs") + facet_grid(species~., labeller = labeller(species=new_names)) + theme_classic() + theme(strip.text = element_text(face = "italic")) +
      geom_segment(data=sum.df, aes(x=0,xend=0,y=0, yend=max.y), colour="black", size=.1) + geom_segment(data=sum.df, aes(x=0,xend=yrs,y=0, yend=0), colour="black", size=.1) 
      print(p)
  }
  }
  
  if(plot.meta==FALSE){
    p <- ggplot() +
      geom_line(data=site.df, aes(x=year, y=count, linetype=stage, color=species,  alpha=trial,  group=title), size=.4) +
      scale_colour_manual(values=color.vect) + scale_alpha_discrete(guide="none") +
      theme_classic() + xlab("years") + ylab("# lemurs") #
    print(p)
    if (spp.facet==TRUE){
      p <- ggplot() +
        geom_line(data=site.df, aes(x=year, y=count, linetype=stage, color=site,  alpha=trial,  group=title), size=.4) +
        scale_alpha_discrete(guide="none") +  xlab("years") + ylab("# lemurs") + #scale_colour_manual(values=color.vect) +
        facet_grid(species~., labeller = labeller(species=new_names)) + theme_classic() + theme(strip.text = element_text(face = "italic")) +
        geom_segment(data=sum.df, aes(x=0,xend=0,y=0, yend=max.y.site), colour="black", size=.1) +
        geom_segment(data=sum.df, aes(x=0,xend=yrs,y=0, yend=0), colour="black", size=.1) 
      
      print(p)
    
    }
  }
  
  if(do.save==TRUE){
    ggsave(file = filename,
           units="mm",  
           width=60, 
           height=70, 
           scale=3, 
           dpi=300)
  }
}
pick.zero.vector <- function(replicate.list){
  zero.vect <- replicate.list[[4]]
  return(zero.vect)
}
plot.open.patches <- function(replicate.list, reps, yrs, meta.df, do.plot, do.save, filename, return.data){
  n.species=length(unique(meta.df$species))
  list.patches <- lapply(replicate.list, pick.zero.vector) #now have a list of all the vectors
  #now sort by species
  list.by.spp <- group.spp(spp.list=list.patches, reps=reps, n.species=n.species)
  #now have a list by species, each of which contains a list of the proportion of unoccupied patches
  names(list.by.spp) <- unique(meta.df$species)
  
  #now get avg, max, min
  avg.list <- mapply(get.avg, df=list.by.spp, name=as.list(names(list.by.spp)), yrs=yrs, SIMPLIFY=FALSE)
  
  avg.df <- do.call("rbind", avg.list)
  rownames(avg.df) <- c()
  
  colourCount <- 11 #length minus the titles
  getPalette <- colorRampPalette(brewer.pal(11, "Spectral"))
  color.vect <- getPalette(colourCount)
  #and how do these 8 fit into the other color vector?
  spp.index <- c(TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE)
  color.vect <- color.vect[spp.index]
  
  #first make a vector to rename things
  new_names <- c(
    'Avahi laniger'="Avahi\nlaniger",
    'Cheirogaleus sp.'="Cheirogaleus\nsp.",
    'Daubentonia madagascariensis'="Daubentonia\nmadagascariensis",
    'Eulemur albifrons'="Eulemur\nalbifrons",
    'Eulemur rubriventer'="Eulemur\nrubriventer",
    'Hapalemur griseus'="Hapalemur\ngriseus",
    'Indri indri'="Indri\nindri",
    'Lepilemur sp.'="Lepilemur\nsp.",
    'Microcebus sp.'="Microcebus\nsp.",
    'Varecia rubra'="Varecia\nrubra",
    'Varecia variegata'="Varecia\nvariegata"
  )
  
  #for each list, take
  if (do.plot==TRUE){
    
    p <- ggplot() +
      geom_line(data=avg.df, aes(x=year,y=mean_prop, color=species), size=1.2) + scale_colour_manual(values=color.vect) +
      geom_ribbon(data=avg.df, aes(x=year,ymin=min_prop, ymax=max_prop, fill=species), alpha=.3) +  scale_fill_manual(values=color.vect) +
      facet_grid(species~., labeller = labeller(species=new_names)) + theme_classic() + theme(strip.text = element_text(face = "italic")) + 
      xlab("years") + ylab("proportion of sites occuppied within metapopulation")
    
    print(p)
    
    if(do.save==TRUE){
      ggsave(file = filename,
             units="mm",  
             width=60, 
             height=70, 
             scale=3, 
             dpi=300)
    }
  }
  
  if(return.data==TRUE){
    return(avg.df)  
  }
}
wrap_save_plot<- function(reps, yrs, simtype, df, seas, leslie, allow.migration,  avg.plot.name, trial.plot.name, prop.plot.name, data.save.name){

  data.name <- replicate(n=reps, connect.leslie.sim.wrap(df=df, simtype=simtype, use.uci=FALSE, N=200, yrs=yrs, 
                                                       leslie=leslie, nat.mort=0, dens.dep=TRUE, seas=seas,
                                                       apply.hunt=FALSE, hunt.intensity=NA, hunt_sim=NA, 
                                                       add.noise=TRUE,  connect.hi=FALSE, return.meta=TRUE,
                                                       allow.migration=allow.migration, do.plot=FALSE, plot.as.metapop=FALSE, do.save=FALSE))

save(data.name, file=data.save.name)
metapop.avg <- average.reps.metapop(replicate.list=data.name, meta.df=df, reps=reps,yrs=yrs, return.avg=TRUE)
plot.avg.metapop(sum.list=metapop.avg, yrs=yrs, do.save=TRUE, plot.meta=TRUE, spp.facet=TRUE, filename=avg.plot.name)
metapop.trial <- average.reps.metapop(replicate.list=data.name, meta.df=df, reps=reps, yrs=yrs, return.avg=FALSE)
plot.trial.metapop(trial.list=metapop.trial, yrs=yrs, do.save=TRUE, plot.meta=TRUE, spp.facet=TRUE, filename=trial.plot.name)
plot.open.patches(replicate.list=data.name, reps=reps, yrs=yrs, meta.df=df, do.plot=TRUE, do.save=TRUE, filename=prop.plot.name, return.data=FALSE)
}
wrap_save_no_plot<- function(reps, yrs, simtype, df, leslie, allow.migration, seas, apply.hunt, hunt_sim, filename){
  data.name <- replicate(n=reps, connect.leslie.sim.wrap.draw(df=df, simtype=simtype, use.uci=FALSE, N=200, yrs=yrs, 
                                                         leslie=leslie, nat.mort=0, dens.dep=TRUE, 
                                                         apply.hunt=apply.hunt, hunt_sim=hunt_sim, #only applicable to leslie =seasonal
                                                         add.noise=TRUE, seas=seas, connect.hi=FALSE, return.meta=TRUE,
                                                         allow.migration=allow.migration, do.plot=FALSE, plot.as.metapop=FALSE, do.save=FALSE))
  save(data.name, file=filename)
}
wrap_save_no_plot_draw<- function(reps, yrs, df, hunt.override, set_hunt, leslie, nat.mort, allow.migration, seas, apply.hunt, hunt_sim, filename){
  data.name <- replicate(n=reps, connect.leslie.sim.wrap.draw(df=df, use.uci=FALSE, N=200, yrs=yrs, hunt.override=hunt.override, set_hunt=set_hunt,
                                                              leslie=leslie, nat.mort=nat.mort, dens.dep=TRUE, 
                                                              apply.hunt=apply.hunt, hunt_sim=hunt_sim, #only applicable to leslie =seasonal
                                                              add.noise=TRUE, seas=seas, connect.hi=FALSE, return.meta=TRUE,
                                                              allow.migration=allow.migration, do.plot=FALSE, plot.as.metapop=FALSE, do.save=FALSE))
  save(data.name, file=filename)
}


#now run the model for 100 years, with 100 random iterations for each run
wrap_save_no_plot_draw(reps=100, yrs=100, hunt.override=FALSE, set_hunt=NA, nat.mort=0, df=lambda.df, leslie="periodic", apply.hunt =FALSE, hunt_sim= NA, seas= NA,
                  allow.migration=TRUE, filename=paste0(homewd, "/prior-scripts/ConsBioLemur/model-output/periodic_no_hunt_yes_mig.Rdata"))

#and repeat with natural background mortality rate of 10% on top of the observed survival
#wrap_save_no_plot_draw(reps=100, yrs=100, hunt.override=FALSE, set_hunt=NA, nat.mort=.1, df=lambda.df, leslie="periodic", apply.hunt =FALSE, hunt_sim= NA, seas= NA,
               #        allow.migration=TRUE, filename=paste0(homewd, "/prior-scripts/ConsBioLemur/model-output/periodic_no_hunt_yes_mig_background_mortality_10.Rdata"))

#you can change these input parameters for each iteration (set hunt, nat.mort, etc.) to get different simulations
#that correspond to those explored in the paper. 

#you should start by working backwards from these wrapper functions to understand the population model at its core.

#then you can simulate with fewer iterations to test the model in its composite form (reps=10 or so instead of 100)
#eventually you can run many iterations of simulations like this on the computing cluster, then collect the output for plotting

#these runs both produce data in the form of an enormous matrix, with 6 rows, each corresponding to a different species
#and 100 columns, each correspinding to a different run of the model

#within each of those cells, you will find a list with 4 different components:
#part 1 of each list has the total infants, juveniles, and adults in the population over time, from 1-100
#part 2 of the list has the tht total infants, juveniles, and adults by site, stacked on top of each other in a matrix
#over time, from 1-100. since there are 7 sites, this means there are 21 rows
#part 3 has the chain of lambda values for the population from timestep 1-100
#part 4 returns the proportion of the 7 sites with 0 lemurs (so the proportion of sites going locally extinct)
#from timestep 1-100

#all of these outputs across a range of different simulation qualities were used to create Figure 4 and 5
#of the Conservation Biology paper


