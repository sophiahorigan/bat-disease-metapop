rm(list=ls())
library(plyr)
library(dplyr)
library(ggplot2)


#first, load in the seropositive/negative data and combine.
#then, join with age and plot age-seroprevalence.
#then, run the model

homewd = "/Users/carabrook/Developer/bat-sero-dynamics"
madabatwd = "/Users/carabrook/Developer/Madagascar-Bat-Data"
luminexdat = "/Users/carabrook/Developer/luminex-serology"

#and load the datasets
lumdat <- read.csv(paste0(luminexdat, "/data/seropos_data_70perc.csv"),header = T, stringsAsFactors = F )
head(lumdat)
histdat  <- read.csv(paste0(madabatwd, "/Serological-Data/Luminex_Henipa_Filo_Brook_2019.csv"),header = T, stringsAsFactors = F )
head(histdat)
catchdat  <- read.csv(paste0(madabatwd, "/Catching-Data/catching_data.csv"),header = T, stringsAsFactors = F )
head(catchdat)
agedat  <- read.csv(paste0(madabatwd, "/Age-Data/aged_bats.csv"),header = T, stringsAsFactors = F )
head(agedat)

#now subset out just NiV for historical and GhV for modern and E. dupreanum
lumdat = subset(lumdat, bat_species=="Eidolon dupreanum" & antigen=="GhV" | bat_species=="Eidolon dupreanum" & antigen=="NiV")
head(lumdat)
histdat <- dplyr::select(histdat, sampleid, species, NIVG_pos)
histdat = subset(histdat, species=="Eidolon dupreanum")
head(histdat)
histdat$antigen <- "NiV"
names(histdat) <- c("sampleid", "bat_species", "seropos", "antigen")
lumdat <- dplyr::select(lumdat, -(MFI), -(log10MFI))
lumdat <- dplyr::select(lumdat, names(histdat))
lumdat <- rbind(lumdat, histdat)


#and attach to age
head(lumdat)#770
head(agedat)
#slim to what you need
agedat <- dplyr::select(agedat, sampleid, age)
#and merge
merge.dat <- merge(lumdat, agedat, by="sampleid")
#merge.dat <- merge(histdat, agedat, by="sampleid")
head(merge.dat)#243 bats that have age and serological info

#and subset only by those specific to GhV
#merge.dat= subset(merge.dat, antigen=="GhV")
merge.dat= subset(merge.dat, antigen=="NiV")

#and bin by age class to plot - bin by .5 and 1, then by year
age_bins <- c(.1,.8,2,4,6,9,11,14,16)
#age_bins <- c(.5,1,2,3,5,7,11,14) #this for the historical data
merge.dat$age_bin <-NA
for(i in 1:length(merge.dat$age)){
merge.dat$age_bin[i] <- min(age_bins[merge.dat$age[i]<=age_bins])
}
sero.sum <- ddply(merge.dat, .(age_bin), summarise, N= length(age), N_seropos = sum(seropos))
sero.sum$seroprevalence <- sero.sum$N_seropos/sero.sum$N
sero.sum$N_seroneg = sero.sum$N - sero.sum$N_seropos
get.binom.CI <- function(df){
  out = prop.test(x=df$N_seropos, n= df$N, alternative="two.sided", conf.level = .95)
  df$seroprev_lci = out$conf.int[1]
  df$seroprev_uci = out$conf.int[2]
  return(df)
}
sero.split <- dlply(sero.sum, .(age_bin))
sero.sum <- data.table::rbindlist(lapply(sero.split, get.binom.CI))


#and fill in the holes
#tmp.df <- cbind.data.frame(age_bin = c(.5, seq(1,16, by=1)))
sero.sum$N[is.na(sero.sum$N)] <- 0
sero.sum$N_seropos[is.na(sero.sum$N_seropos)] <- 0
sero.sum$seroprevalence[is.na(sero.sum$seroprevalence)] <- 0

p1 <- ggplot(data = sero.sum) + theme_bw() + theme(panel.grid = element_blank()) +
      geom_point(aes(x=age_bin, y=seroprevalence, size=N), shape=1) +
      geom_linerange(aes(x=age_bin, ymin=seroprev_lci, ymax=seroprev_uci), linetype=3) +
      geom_line(aes(x=age_bin, y=seroprevalence))
p1


#add date and biweek
catch.merge <- dplyr::select(catchdat, sampleid, collection_date)


merge.dat <- merge(merge.dat, catch.merge, by = "sampleid")
head(merge.dat)
merge.dat$collection_date <- as.Date(merge.dat$collection_date, format="%m/%d/%y")

library(lubridate)
merge.dat$week <- week(merge.dat$collection_date)

#make dates
biwk.dat <-cbind.data.frame(week=1:52, biwk = rep(1:26, each=2)) 

merge.dat <- merge(merge.dat, biwk.dat, by="week")
head(merge.dat)
#now save this dataset 
write.csv(merge.dat, file = paste0(homewd, "/data/EidNiVfit.csv"), row.names = F)
#and run the model on these data
