# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Supplementary Analyses.R         
# - DATE:        27/05/2020
# - DESCRIPTION: Code to explore the influence of body mass, Raunkier's growth
#                form and conservation status on the three components of 
#                demographic resilience. 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything
options(mc.cores = parallel::detectCores())

library(plyr)
library(dplyr)
library(tidyr)
library(ape)
library(caper)
library(geiger)
library(phytools)
library(brms)
library(bayestestR)

#Working directories

path <-  gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath <- paste0(path,"/Data")
ResultPath <- paste0(path, "/Results") 

#Load  and prepare data----

setwd(DataPath)

# Load animals and plants

load("ResData.RData")

# Raunkier analyses ############################################################

raunkier <- read.csv("RaunkiaerGrowthForms.csv")

# Join with the data

smallplandata <- smallplandata %>% 
  left_join(raunkier, by=c("species"="Species"))

# Correct some species manually 

smallplandata$Growth.form.Raunkiaer[which(smallplandata$SpeciesAccepted=="Callitris columellaris")]="Megaphanerophyte"
smallplandata$Growth.form.Raunkiaer[which(smallplandata$SpeciesAccepted=="Sonchus pustulatus")]="Hemicryptophyte"
smallplandata$Growth.form.Raunkiaer[which(smallplandata$SpeciesAccepted=="Thymus vulgaris")]="Chamaephyte"
smallplandata$Growth.form.Raunkiaer[which(smallplandata$SpeciesAccepted=="Trollius europaeus")]="Hemicryptophyte"
smallplandata$Growth.form.Raunkiaer[which(smallplandata$SpeciesAccepted=="Vitaliana primuliflora")]="Hemicryptophyte"
smallplandata$Growth.form.Raunkiaer[which(smallplandata$SpeciesAccepted=="Chaerophyllum aureum")]="Hemicryptophyte"
smallplandata$Growth.form.Raunkiaer[which(smallplandata$SpeciesAccepted=="Colchicum autumnale")]="Geophyte"
smallplandata$Growth.form.Raunkiaer[which(smallplandata$SpeciesAccepted=="Cornus florida")]="Mesophanerophyte"
smallplandata$Growth.form.Raunkiaer[which(smallplandata$SpeciesAccepted=="Dracocephalum austriacum")]="Chamaephyte"
smallplandata$Growth.form.Raunkiaer[which(smallplandata$SpeciesAccepted=="Linum flavum")]="Chamaephyte"
smallplandata$Growth.form.Raunkiaer[which(smallplandata$SpeciesAccepted=="Linum tenuifolium")]="Chamaephyte"
smallplandata$Growth.form.Raunkiaer[which(smallplandata$SpeciesAccepted=="Saussurea medusa")]="Hemicryptophyte"
smallplandata$Growth.form.Raunkiaer[which(smallplandata$SpeciesAccepted=="Trollius laxus")]="Hemicryptophyte"

# Let's have a look at it

table(smallplandata$Growth.form.Raunkiaer)

# Analyses ---------------------------------------------------------------------

#Set modelling parameters

iter <- 8000
thin <- 0.0005 * iter
warmup <- 0.1 * iter

priors <-  c( prior(normal(0, 10), "b"),
              prior(normal(0, 1e6), "Intercept"))

# Subset the dataset

plants <- smallplandata %>% 
  filter(!is.na(Growth.form.Raunkiaer))

# Correct the tree

sp <- unique(plants$species)
names(sp) <- unique(plants$species)
(chk<- name.check(smalltree, sp))
plantree <- drop.tip(smalltree,chk$tree_not_data)
plants <- plants[plants$species%in%plantree$tip.label,]

# Define the correlation structure

A <- ape::vcv.phylo(plantree)

# The three of them #################

mraunkier <- brm(mvbind(scale(xt), scale(rupr), scale(rlwr)) ~ Growth.form.Raunkiaer + scale(Dimension)+ 
                   (1|gr(animal, cov = A)) +(1|popID)+(1|typeM),
              iter = iter, thin = thin, warmup = warmup,
              prior= priors, 
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = plants,
              family = gaussian(), 
              data2 = list(A = A),
              cores = 8)

# Plants body size #############################################################

# Remove nas from the dataset ### 

plants <- smallplandata %>% 
  filter(!is.na(max_height))

# Match the tree

sp <- unique(plants$species)
names(sp) <- unique(plants$species)
(chk<- name.check(smalltree, sp))

# Prune the tree

plantree <- drop.tip(smalltree,chk$tree_not_data)

# Sort dataframe

plants <- plants[plants$species%in%plantree$tip.label,]

# Analyses ####

#Set modelling priors

priors <-  c( prior(normal(0, 10), "b"),
              prior(normal(0, 1), "Intercept"))

# Define the correlation structure

A <- ape::vcv.phylo(smalltree)

# Run the model 

smallplandata <- smallplandata %>% drop_na(GenT,max_height)

# Define the correlation structure

sp <- unique(smallplandata$species)
names(sp) <- unique(smallplandata$species)
(chk<- name.check(smalltree, sp))
smallplantree <- drop.tip(smalltree,chk$tree_not_data)
smallplandata <- smallplandata[smallplandata$species%in%smalltree$tip.label,]

A <- ape::vcv.phylo(smalltree)

mbodyplants <- brm(mvbind(scale(xt), scale(rupr), scale(rlwr)) ~  + scale(log10(max_height+1)) + scale(Dimension)+ (1|gr(animal, cov = A)) +(1|popID)+(1|typeM),
              iter = iter, thin = thin, warmup = warmup,
              prior= priors, 
              control = list(adapt_delta = .975, max_treedepth = 20),
              data = smallplandata,
              family = gaussian(), 
              data2 = list(A = A),
              cores = 8)

# Animals

andata <- smallandata %>% drop_na(GenT,adult_body_mass_g)

# Define the correlation structure
sp <- unique(andata$species)
names(sp) <- unique(andata$species)
(chk<- name.check(smallantree, sp))
smallantree2 <- drop.tip(smallantree,chk$tree_not_data)
andata2 <- andata[andata$species%in%smallantree2$tip.label,]

A <- ape::vcv.phylo(smallantree2)

mbodyanimals <- brm(mvbind(scale(xt), scale(rupr), scale(rlwr)) ~ scale(log10(adult_body_mass_g+1)) + scale(Dimension)+ (1|gr(animal, cov = A)) +(1|popID)+(1|typeM),
                   iter = iter, thin = thin, warmup = warmup,
                   prior= priors, 
                   control = list(adapt_delta = .975, max_treedepth = 20),
                   data = andata2,
                   family = gaussian(), 
                   data2 = list(A = A),
                   cores = 8)

# Save models

setwd(ResultPath)
save(plants, mbodyanimals, mbodyplants, mraunkier, file="Supplementary.RData")
