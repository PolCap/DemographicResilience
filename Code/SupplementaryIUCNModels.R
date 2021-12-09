# --------------------------------------------------------------------------------------- #
# - FILE NAME:   7.IUCNModels.R         
# - DATE:        21/05/2020
# - DESCRIPTION: Code to descrive the resilience patterns between species and groups of
#                organisms. 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

library(boot)
library(plyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ape)
library(caper)
library(geiger)
library(phytools)
library(reshape)
library(brms)
library(bayestestR)
library(tidybayes)
library(shinystan)

#Working directories

path = gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath = paste0(path,"/Data")
ResultPath = paste0(path, "/Results") 

options(mc.cores = parallel::detectCores()-1)

n_cores <- parallel::detectCores()
n_cores
cluster <- parallel::makeCluster(n_cores-1, type="PSOCK")
# Register this cluster with do snow
doSNOW::registerDoSNOW(cluster)

#Load  and prepare data----

setwd(DataPath)

# Load animals and plants

load("ResData.RData") 

# Transfor IUCN category

smallplandata$category <- ordered(smallplandata$category,
                             c("LC", "NT", "VU", "EN", "CR"))
smallandata$category <- ordered(smallandata$category,
                           c("LC", "NT", "VU", "EN", "CR"))

# Animals ##############################################################################

# Remove the NAs

smallandata <- subset(smallandata, !is.na(category))

# Analyses --------------------------------------------------------------------------------------------------------------
#Set modelling parameters

iter = 8000
thin = 0.0005 * iter
warmup = 0.1 * iter

# Define priors 

priors <- c(#prior(normal(0, 10), class = Intercept),
            prior(normal(0, 1), class = b))


# Define the correlation structure


sp <- unique(smallandata$species)
names(sp) <- unique(smallandata$species)
(chk<- name.check(smallantree, sp))
smallantree <- drop.tip(smallantree,chk$tree_not_data)
smallandata <- smallandata[smallandata$species%in%smallantree$tip.label,]

A <- ape::vcv.phylo(smallantree)

# The three of them #################

mcon_a <- brm(mvbind(scale(xt), scale(rupr), scale(rlwr)) ~ category-1 + scale(Dimension)+ (1|gr(animal, cov = A)) + (1|popID)+(1|typeM),
              #iter = iter, thin = thin, warmup = warmup,
              prior= priors, 
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = smallandata,
              family = gaussian(), 
              data2 = list(A = A),
              cores = 8)

# Plants ##############################################################################

# Remove the NAs

smallplandata <- subset(smallplandata, !is.na(category))
sp <- unique(smallplandata$species)
names(sp) <- unique(smallplandata$species)
(chk<- name.check(smalltree, sp))
smallplantree <- drop.tip(smalltree,chk$tree_not_data)
smallplandata <- smallplandata[smallplandata$species%in%smalltree$tip.label,]

A <- ape::vcv.phylo(smalltree)

mcon_p <- brm(mvbind(scale(xt), scale(rupr), scale(rlwr)) ~ category-1 + scale(Dimension)+ (1|gr(animal, cov = A)) + (1|popID)+(1|typeM),
              iter = iter, thin = thin, warmup = warmup,
              prior= priors, 
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = smallplandata,
              family = gaussian(), 
              data2 = list(A = A),
              cores = 8)

setwd(ResultPath)
save(mcon_a, mcon_p, file="ConservationModels.RData")

