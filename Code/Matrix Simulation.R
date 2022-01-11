# ---------------------------------------------------------------------------- #
# - FILE NAME:   MatrixSimulation.R         
# - DATE:        10/10/2021
# - DESCRIPTION: Code to simulate different matrix population models and calculate 
#                resilience and life history traits
# - AUTHORS:     Pol Capdevila Lanzaco
# ---------------------------------------------------------------------------- #

rm(list = ls(all=T))

# Librearies

library(doSNOW)
library(parallel)
library(tidyverse)
library(popbio)
library(popdemo)
library(Matrix)
library(Rage)
library(Rcompadre)
library(plyr)
library(mpmtools)

#Working directories

path <- gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath <- paste0(path,"/Data")
ResultPath <- paste0(path, "/Results") 
CodePath <- paste0(path, "/Code") 

# Load functions to simulate the data 

source(paste0(CodePath, "/SimulationFunctions.R"))
source(paste0(CodePath, "/Functions.R"))

# Load the species that we used 

load(paste0(DataPath,"/ResData.RData"))

# Store the range of matrix sizes

dimensions <- range(c(smallandata$Dimension, smallplandata$Dimension))

# Delete the data we are not going to use

rm(smallandata, smallantree, smallplandata, smalltree)

# Simulate random matrices #####################################################

# Random matrices --------------------------------------------------------------
# Here we simulate matrix population models with random non-negative elements, 
# keeping the Umat column sums to <= 1 (=survival cannot be higher than 1). 

# We will use a for each loop 

clus  <- makeCluster(detectCores() - 1)

# register the cluster with doParallel

registerDoSNOW(clus)

# Build random matrices with stasis, progression and retrogression

random <- foreach(i = c(3:max(dimensions)),
                  .combine = "rbind") %dopar% {
                    matrices <- replicate(n=100, 
                                          random_matrices(dimension=i),
                                          simplify = F)
                    return(matrices)
                  }
# Estimate the life history traits 

random_m <- foreach(i = c(1:length(random)),
                    .combine = "rbind",
                    .packages = c("tidyverse", "popdemo", "popbio", "Rage")) %dopar% {
                      matrix <- random[[i]]
                      if(isErgodic(matrix$matrix_A)){
                        tryCatch(lam <- lambda(matrix$matrix_A), 
                                 error=function(e) NA)
                        tryCatch(GenT <- generation.time(matrix$matrix_A), 
                                 error=function(e) NA)
                        tryCatch(Fec <- vitalRates(matrix$matrix_U,
                                                   matrix$matrix_F)$fec,
                                 error=function(e) NA)
                        tryCatch(rupr <- reac(A = matrix$matrix_A,bound="upper"),
                                 error=function(e) NA)
                        tryCatch(rlwr<- reac(A = matrix$matrix_A,
                                             bound="lower"),
                                 error=function(e) NA)
                        tryCatch(xt <- return.time(matrix$matrix_A), 
                                 error=function(e) NA)
                        tryCatch(Dimension <- dim(matrix$matrix_A)[1], 
                                 error=function(e) NA)
                        
                      }else{
                        lam <- NA
                        GenT <- NA
                        Fec <- NA
                        rupr<- NA
                        rlwr<- NA
                        xt <- NA
                        Dimension <- NA
                      }
                      return(tibble(mat=list(matrix), lam,
                                    GenT, Fec, rupr, rlwr, xt, Dimension))
                    }

# Build random matrices only with stasis and progression

random_p <- foreach(i = c(3:max(dimensions)),
                  .combine = "rbind") %dopar% {
                    matrices <- replicate(n=100, 
                                          random_matrices(dimension=i, upd=TRUE),
                                          simplify = F)
                    return(matrices)
                  }

# Estimate the life history traits 

random_mp <- foreach(i = c(1:length(random)),
                    .combine = "rbind",
                    .packages = c("tidyverse", "popdemo", "popbio", "Rage")) %dopar% {
                      matrix <- random_p[[i]]
                      if(isErgodic(matrix$matrix_A)){
                        tryCatch(lam <- lambda(matrix$matrix_A), 
                                 error=function(e) NA)
                        tryCatch(GenT <- generation.time(matrix$matrix_A), 
                                 error=function(e) NA)
                        tryCatch(Fec <- vitalRates(matrix$matrix_U,
                                                   matrix$matrix_F)$fec,
                                 error=function(e) NA)
                        tryCatch(rupr <- reac(A = matrix$matrix_A,bound="upper"),
                                 error=function(e) NA)
                        tryCatch(rlwr<- reac(A = matrix$matrix_A,
                                             bound="lower"),
                                 error=function(e) NA)
                        tryCatch(xt <- return.time(matrix$matrix_A), 
                                 error=function(e) NA)
                        tryCatch(Dimension <- dim(matrix$matrix_A)[1], 
                                 error=function(e) NA)
                        
                      }else{
                        lam <- NA
                        GenT <- NA
                        Fec <- NA
                        rupr<- NA
                        rlwr<- NA
                        xt <- NA
                        Dimension <- NA
                      }
                      return(tibble(mat=list(matrix), lam,
                                    GenT, Fec, rupr, rlwr, xt, Dimension))
                    }

#Stop the cluster 

stopCluster(clus)

# Data cleaning and preparation ################################################
# First, we constrain the matrices to certain lambda values

leslie_matrices <- leslie_matrices %>% filter(lam>-1.5&lam<1.5)
lefkovitch_matrices <- lefkovitch_matrices %>% filter(lam>-1.5&lam<1.5)


# We now clean the data as we did in the original analyses 

random_m <- random_m %>% 
  filter(!is.infinite(GenT)) %>% 
  mutate(xt=log(xt+1),
         rupr=log(rupr+1),
         rlwr=abs(log((1 - rlwr)+1)),
         Fec=log(Fec+1),
         GenT=log(GenT+1))

random_mp <- random_mp %>% 
  filter(!is.infinite(GenT)) %>% 
  mutate(xt=log(xt+1),
         rupr=log(rupr+1),
         rlwr=abs(log((1 - rlwr)+1)),
         Fec=log(Fec+1),
         GenT=log(GenT+1))

# Save them ####################################################################

setwd(DataPath)

save(random_m, 
     random_mp, file="SimData.RData")

