
setwd("/Users/shorigan/Documents/GitHub/Bat-disease-metapop/bat-disease-metapop/Models/WIP/")


LOCAL = 1
CLUSTER = 0
do.save.data = FALSE



## set values
yrs = 1
numsims = 1 
prop_occupied_patches = 1
inf_low = 0.1 # used for fixed also
inf_high = 0.45
burnin = 0
stoch_demo = TRUE
stoch_disease = TRUE
model_vect <- as.vector(c(1))
grid_vect = as.vector(c(100, 1000))
patch_vect = as.vector(c(2))
pop_size = as.vector(c(10000))
intermingle_vect <- as.vector(c("TRUE"))
dispersal_vect <- as.vector(c("TRUE"))
simtype <- as.vector(c("Int_Disp"))
#intermingle_vect <- as.vector(c("TRUE", "FALSE", "FALSE", "TRUE"))
#dispersal_vect <- as.vector(c("TRUE", "FALSE", "TRUE", "FALSE"))
#simtype <- as.vector(c("Int_Disp", "No int_No disp", "No int_Disp", "Int_No disp"))
param.array <- read.csv("Input/fit_parameters.csv")

COUNT = 1
## generate IC list
for (i in 1:length(model_vect)){ # models
  for (j in 1:length(grid_vect)){ # grid sizes
    for (k in 1:length(patch_vect)){ # patch dim
      for (m in 1:length(pop_size)){ # init population size
        for (n in 1:length(dispersal_vect)){
          
          model.list = c(rep(i, numsims))
          num_patches.list = c(rep(patch_vect[k]^2, numsims))
          patch_dim.list = c(rep(patch_vect[k], numsims))
          grid_size.list = c(rep(grid_vect[j], numsims))
          prop_occupied_patches.list = c(rep(prop_occupied_patches, numsims))
          yrs.list = c(rep(yrs, numsims))
          pop_low.list = rep(list(rep(pop_size[m]/patch_vect[k]^2, patch_vect[k]^2)), numsims)
          pop_high.list = rep(list(rep(pop_size[m]/patch_vect[k]^2, patch_vect[k]^2)), numsims)
          sim_pop.df = GeneratePops(num_patches = patch_vect[k]^2, prop_occupied_patches = prop_occupied_patches, grid_mode = "stochastic", pop_mode = "static", pop_low = pop_low.list, pop_high = pop_high.list)
          intermingling.list = c(rep(intermingle_vect[n], numsims))
          dispersal.list = c(rep(dispersal_vect[n], numsims))
          do.save.data.list = c(rep(do.save.data, numsims))
          sim.list = c(rep(1:numsims))
          burnin.list = c(rep(burnin, numsims))
          ntyr.list = c(rep(26, numsims))
          s.list = c(rep(20, numsims))
          age.brk.list = c(rep(20, numsims))
          init_fracI.df = GenerateFracI(num_patches = patch_vect[k]^2, num_sims = numsims, inf_mode = "stochastic", inf_low = inf_low, inf_high = inf_high)
          add.inf.mort.list = c(rep(FALSE, numsims))
          printpop.list = c(rep(FALSE, numsims))
          stoch_demo.list = c(rep(stoch_demo, numsims))
          stoch_disease.list = c(rep(stoch_disease, numsims))
          simtype.list = c(rep(simtype[n], numsims))
          param.array.list = c(rep(param.array), 1)
          iter = 1
          
          IC_list <- c(list(sim.list, model.list, num_patches.list, patch_dim.list, grid_size.list, prop_occupied_patches.list, yrs.list, intermingling.list, dispersal.list, do.save.data.list, burnin.list, ntyr.list, s.list, age.brk.list, add.inf.mort.list, printpop.list, stoch_demo.list, stoch_disease.list, simtype.list, iter, COUNT))
          IC_pops <- c(list(sim_pop.df, init_fracI.df))
          param.array.list <- c(list(param.array.list))

          if (CLUSTER == 1){
            write.table(IC_list, file = paste0("Input/IC_list_", COUNT), col.names = FALSE)
            write.table(IC_pops, file = paste0("Input/IC_pops_", COUNT), col.names = FALSE)
            write.table(param.array.list, file = paste("Input/param.array.list_", COUNT), col.names = FALSE)
          }
          
          if (LOCAL == 1){
            RunModel(IC_list = IC_list, IC_pops = IC_pops, param.array.list = param.array.list)
          }
          
          COUNT = COUNT + 1
        }
      }
    }
  }
}


          
          
          
          