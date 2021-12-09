# ---------------------------------------------------------------------------- #
# - FILE NAME:   MatrixSimulation.R         
# - DATE:        10/10/2021
# - DESCRIPTION: Code to randomise matrices and calculate resilience and 
#                life history traits
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

# Simulation ----

# Store a list of values for the different survival and reproductive curves 

params <- list(values= seq(3,max(dimensions),1), #Matrix sizes 
               types =  c(1.2,1.4, 1.6, 1.8, 1.9, 2.1, 2.3, 2.5, 2.7), #Survival curve types
               exponents = c(5, 10, 15), # Different exponents rerpoductive curve
               cutoffs = c(0.2, 0.4, 0.6, 0.8), #Reproductive curve cutoffs
               slopes = c(0.1,1, 2), # Slopes of the curves
               starting_values = c(0.1,1, 5))

# Make a data.frame of every possible combination of the parameters

all_var <- expand.grid(params)

# Convert that into a list, where each object is 1 row of the data frame

all_var_list <- as.list(as.data.frame(t(all_var)))

# Calculate the survival and reproductive curves

surv_curvs <- lapply(all_var_list, function(x){
  generate_survival_curve(type = x[2], 
                          values = x[1], 
                          exp=x[3])})
rep_curvs <- lapply(all_var_list, function(x){
  generate_reproduction_curve(values = x[1], 
                              cutoff = x[4], 
                              slope = x[5],
                              starting_value = x[6])})

# Transform into a data.frame 

surv_data <- ldply(surv_curvs, data.frame)
rep_data <- ldply(rep_curvs, data.frame)

# Join them 

life_tables <- surv_data %>% cbind(rep_data[-c(1,2)])

# We remove one matrix that gives error with the isPrimitive function
# 
# life_tables <- life_tables %>% filter(.id!=c("V51187"))

# Now we convert the curves into matrix population models ###################### 

# make a cluster, defining the number of cores to use

clus  <- makeCluster(detectCores() - 1)

# register the cluster with doParallel

registerDoSNOW(clus)

# Build the Leslie matrices

leslie_matrices <- foreach(i = unique(life_tables$.id),
                        .combine = "rbind",
                        .packages = c("tidyverse", "popdemo", "popbio", "Rage")) %dopar% {
                          x <- life_tables %>% filter(.id==i)
                          matrix <- build_leslie(survival = x$px, 
                                                 reproduction = x$reproduction)
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
    return(tibble(mat=list(matrix),lam, 
                      GenT, Fec, rupr, rlwr, xt, Dimension))
  }


# Now for lefkovitch matrices
lefkovitch_matrices <- foreach(i = unique(life_tables$.id),
                           .combine = "rbind",
                           .packages = c("tidyverse", "popdemo", "popbio", "Rage")) %dopar% {
                             x <- life_tables %>% filter(.id==i)
                             tryCatch(matrix <-build_lefkovitch(x=x$x,
                                                       survival = x$px, 
                                                       reproduction = x$reproduction),
                                      error=function(e) NA)
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

# We set resitance measures to inverse

#Remove infinites

leslie_matrices <- leslie_matrices[!is.infinite(leslie_matrices$GenT),]
lefkovitch_matrices <- lefkovitch_matrices[!is.infinite(lefkovitch_matrices$GenT),]

#We logtransform the transient data 

leslie_matrices <- leslie_matrices %>% 
  mutate(xt=log(xt+1),
         rupr=log(rupr+1),
         rlwr=abs(log((1 - rlwr)+1)),
         Fec=log(Fec+1),
         GenT=log(GenT+1))

lefkovitch_matrices <- leslie_matrices %>% 
  mutate(xt=log(xt+1),
         rupr=log(rupr+1),
         rlwr=abs(log((1 - rlwr)+1)),
         Fec=log(Fec+1),
         GenT=log(GenT+1))

# Save them ####################################################################

setwd(DataPath)

save(lefkovitch_matrices, leslie_matrices, file="SimData.RData")

