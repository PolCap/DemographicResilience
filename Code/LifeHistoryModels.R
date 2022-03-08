# --------------------------------------------------------------------------------------- #
# - FILE NAME:   GenTModdels.R         
# - DATE:        04/05/2020
# - DESCRIPTION: Code to perform brms analyses of generation time vs demographic
#                resilience. 
# - AUTHORS:     Pol Capdevila Lanzaco
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

library(plyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ape)
library(caper)
library(geiger)
library(phytools)
library(brms)
library(bayestestR)
library(tidybayes)
library(shinystan)
library(reshape)
library(boot)
library(dplyr)
library(rstan)

#Working directories

path = gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath = paste0(path,"/Data")
ResultPath = paste0(path, "/Results") 

#Load  and prepare data----

setwd(DataPath)

# Load animals and plants

load("ResData.RData")

#These options help Stan run faster:

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

# Plants generation time ###################################

smallplandata <- smallplandata %>% 
  drop_na(GenT, Fec)

sp <- unique(smallplandata$species)
names(sp) <- unique(smallplandata$species)
(chk<- name.check(smalltree, sp))
smallplantree <- drop.tip(smalltree,chk$tree_not_data)
smallplandata <- smallplandata[smallplandata$species%in%smalltree$tip.label,]

# Define the correlation structure

A <- ape::vcv.phylo(smalltree)

#  Analyses ####

#Set modelling parameters

iter = 8000
thin = 0.0005 * iter
warmup = 0.1 * iter

# Define priors 

priors <- c(prior(normal(0, 10), class = Intercept),
            prior(normal(0, 1), class = b))


# Define the correlation structure

A <- ape::vcv.phylo(smalltree)

# The three of them #################

mGenTp <- brm(mvbind(scale(xt), scale(rupr), scale(rlwr)) ~ scale(GenT) + scale(Fec) + 
                scale(GenT):scale(Fec) + scale(Dimension)+ (1|gr(animal, cov = A)) +(1|popID)+(1|typeM),
              iter = iter, thin = thin, warmup = warmup,
              prior= priors, 
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = smallplandata,
              family = gaussian(), 
              data2 = list(A = A),
              cores = 20)


# Animals ###################################

smallandata <- smallandata %>% drop_na(GenT)

# Define the correlation structure
sp <- unique(smallandata$species)
names(sp) <- unique(smallandata$species)
(chk<- name.check(smallantree, sp))
smallantree <- drop.tip(smallantree,chk$tree_not_data)
smallandata <- smallandata[smallandata$species%in%smallantree$tip.label,]

A <- ape::vcv.phylo(smallantree)

mGenTa <- brm(mvbind(scale(xt), scale(rupr), scale(rlwr)) ~ scale(GenT) + scale(Fec) + scale(GenT):scale(Fec)+ scale(Dimension)+ (1|gr(animal, cov = A)) +(1|popID)+(1|typeM),
              iter = iter, thin = thin, warmup = warmup,
              prior= priors, 
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = smallandata,
              family = gaussian(), 
              data2 = list(A = A),
              cores = 20)

# Save the models

setwd(ResultPath)
save(mGenTa, mGenTp,
     file = "ModelsGenTime2.RData")