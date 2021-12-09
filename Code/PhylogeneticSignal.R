# --------------------------------------------------------------------------------------- #
# - FILE NAME:   PhyloSig.R         
# - DATE: 07/12/2020
# - DESCRIPTION: 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything
options(mc.cores = parallel::detectCores())

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
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(tidytree)
library(ggnewscale)
library(ggpmisc)
library(ggplot2)
library(ggimage)
library(patchwork)
library(tidyverse)
library(viridis)

#Working directories

path <-  gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath <- paste0(path,"/Data")
ResultPath <- paste0(path, "/Results") 

#Load and prepare data----

setwd(DataPath)

# Load animals and plants

load("ResData.RData")
#load("Explore.RData")

# Phylogenetic signal ##########################################################
# Set modelling parameters

iter <- 8000
thin <- 0.0005 * iter
warmup <- 0.1 * iter

# Set priors

priors <- c(prior(normal(0, 10), class = Intercept),
            #prior(normal(0, 1), class = b),
            prior(exponential(1), class = sd))

# Animals ----------------------------------------------------------------------

# Define the correlation structure

A <- ape::vcv.phylo(smallantree)

# Model resistance -------------------------------------------------------------

res_a <- brm(rlwr ~ 1 + (1|gr(animal, cov = A)) + (1|species)+(1|typeM),
             iter = iter, thin = thin, warmup = warmup,
             prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
             data = smallandata,
             family = gaussian(), 
             data2 = list(A = A),
             cores = 8)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_species__Intercept^2 + sd_typeM__Intercept^2 + sigma^2) = 0")

(psign_res <- hypothesis(res_a, hyp, class = NULL))

# Model recovery ---------------------------------------------------------------

rec_a <- brm(xt ~ 1 + (1|gr(animal, cov = A)) + (1|species)+(1|typeM),
             iter = iter, thin = thin, warmup = warmup,
             prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
             data = smallandata,
             family = gaussian(), 
             data2 = list(A = A),
             cores = 8)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_species__Intercept^2 + sd_typeM__Intercept^2 + sigma^2) = 0")

(psign_rec <- hypothesis(rec_a, hyp, class = NULL))

# Model compensation -----------------------------------------------------------

com_a <- brm(rupr ~ 1 + (1|gr(animal, cov = A)) + (1|species)+(1|typeM),
             iter = iter, thin = thin, warmup = warmup,
             prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
             data = smallandata,
             family = gaussian(), 
             data2 = list(A = A),
             cores = 8)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_species__Intercept^2 + sd_typeM__Intercept^2 + sigma^2) = 0")

(psign_com <- hypothesis(com_a, hyp, class = NULL))

# Plants ----------------------------------------------------------------------

# Define the correlation structure

A2 <- ape::vcv.phylo(smalltree)

# Model resistance -------------------------------------------------------------

res_p <- brm(rlwr ~ 1 + (1|gr(animal, cov = A2)) + (1|species)+(1|typeM),
             iter = iter, thin = thin, warmup = warmup,
             prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
             data = smallplandata,
             family = gaussian(), 
             data2 = list(A2 = A2),
             cores = 8)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_species__Intercept^2 + sd_typeM__Intercept^2 + sigma^2) = 0")

(psign_res_p <- hypothesis(res_p, hyp, class = NULL))

# Model recovery ---------------------------------------------------------------

rec_p <- brm(xt ~ 1 + (1|gr(animal, cov = A2)) + (1|species)+(1|typeM),
             iter = iter, thin = thin, warmup = warmup,
             prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
             data = smallplandata,
             family = gaussian(), 
             data2 = list(A2 = A2),
             cores = 8)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_species__Intercept^2 + sd_typeM__Intercept^2 + sigma^2) = 0")

(psign_rec_p <- hypothesis(rec_p, hyp, class = NULL))

# Model compensation -----------------------------------------------------------

com_p <- brm(rupr ~ 1 + (1|gr(animal, cov = A2)) + (1|species)+(1|typeM),
             iter = iter, thin = thin, warmup = warmup,
             prior= priors, control = list(adapt_delta = .975, max_treedepth = 20),
             data = smallplandata,
             family = gaussian(), 
             data2 = list(A2 = A2),
             cores = 8)

# Now we measure the phylogenetic signal

hyp <-  paste("sd_animal__Intercept^2 /", 
              "(sd_animal__Intercept^2 + sd_species__Intercept^2 + sd_typeM__Intercept^2 + sigma^2) = 0")

(psign_com_p <- hypothesis(com_p, hyp, class = NULL))


# Save the results

setwd(ResultPath)
save(res_a, rec_a,com_a, psign_res, psign_rec, psign_com, 
     res_p, rec_p,com_p, psign_res_p, psign_rec_p, psign_com_p, 
     file="PhyloSingal.RData")