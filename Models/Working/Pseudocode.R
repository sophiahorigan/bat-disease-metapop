########################
## PSEUDOCODE ##
#######################

# Created by Sophia Horigan (shorigan@uchicago.edu)
# Started: 01-17-24
# Last updated: 01-17-24

# Goal: Draft pseudocode for bat disease metapop model

# Project Overview: This project seeks to use existing demographic data and the current best-fit transmission model
# to explore persistence threshholds for hypothetical pathogens under data-based metapopulation structure of Pteropus
# rufus in Madagascar. 

## PART 1 : SIMULATE MSIRN MODEL
# set model parameters
# function : run model
# plotting : 
# statistics : AP (annual persistence: with prob > 50% infection will persist in population for 1 year)
#           : LP (long-term persisence: with prob > 50% infection will persist in population for 100 years)

## PART 2 : MAKE IT STOCHASTIC 
# instead of fixed variables, draw value from a distribution

## PART 3 : MAKE IT A METAPOPULATION MODEL
# make SIR model into an array with variable 'patch' numbers
# make K patch specific
# add in dispersal, that can seasonally change -- use buildFmatrix as guide - that changes over the year

## PART 4 : ADD IN METAPOP POP ESTIMATES AND CONNECTIVITY ESTIMATES
# use telemetry data to estimate
# use roost count pop estimates? use stochastic range of estimates? 

## PART 5 : TEST UNDER DIFFERENT PATHOGENS
# we know there are henipaviruses in bats - can we guess their parameters?
# use estimates from serology of Brook 2019 paper
# or just try out a range of parameters that represent different pathogens 


