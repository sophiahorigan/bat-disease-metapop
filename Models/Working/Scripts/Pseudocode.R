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

## PART 4 : ADD IN CONNECTIVITY STRUCTURE
# potentially two different things here
# 1. literal dispersal
# 2. probability of interacting i.e. transmitting

## PART 5 : DETERMINE AND CODE IN FINAL OUTPUT
# proportion of simulations that go extinct? 
# AP vs LP?
# print param sets - capture stochasticity

## PART 6 : DETERMINE DIFFERENT LOOPS
# different pathogen parameter values, different init pop sizes
# different introduction times

## PART 7 : ADD IN TELEMETETRY DATA

## PART 8 : SIMULATE!! 


