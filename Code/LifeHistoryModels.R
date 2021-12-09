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

#plot(conditional_effects(mGenTp))

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

# Plots 

dat <- mGenTa %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(.variable = gsub("b_", "", .variable),
         .variable = gsub("scale", "", .variable),
         .variable = gsub("rlwr", "Resistance", .variable),
         .variable = gsub("rupr", "Compensation", .variable),
         .variable = gsub("xt", "Recovery time", .variable)) %>% 
  filter(.variable!="Intercept",
         .variable!="Dimension") 

(g1a <- mGenTa %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    mutate(.variable = gsub("b_", "", .variable),
           .variable = gsub("scale", "", .variable),
           .variable = gsub("rlwr", "Resistance", .variable),
           .variable = gsub("rupr", "Compensation", .variable),
           .variable = gsub("xt", "Recovery time", .variable)) %>% 
    filter(.variable!="Intercept",
           .variable!="Dimension") %>% 
    median_qi(.width=c(.5,.8,.95)) %>% 
    ggplot(aes(x = .variable, y = .value)) +
    geom_pointinterval(aes(ymin = .lower, 
                           ymax = .upper),
                       interval_size_range = c(0.5, 2), 
                       colour="grey35") +
    geom_hline(yintercept = 0, linetype="dashed", colour="grey40")+
    stat_slab(data=dat, alpha=0.2, aes(fill=.variable))+
    labs(x="", y = "Posterior estimate \n slope with Generation time") +
    theme(legend.position = "none"))

dat <- mGenTp %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(.variable = gsub("b_", "", .variable),
         .variable = gsub("scale", "", .variable),
         .variable = gsub("log10max_heightP1", "Body size", .variable),
         .variable = gsub("rlwr", "Resistance", .variable),
         .variable = gsub("rupr", "Compensation", .variable),
         .variable = gsub("xt", "Recovery time", .variable)) %>% 
  filter(.variable!="Intercept",
         .variable!="Dimension") 

(g1b <- mGenTp %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    mutate(.variable = gsub("b_", "", .variable),
           .variable = gsub("scale", "", .variable),
           .variable = gsub("log10max_heightP1", "Body size", .variable),
           .variable = gsub("rlwr", "Resistance", .variable),
           .variable = gsub("rupr", "Compensation", .variable),
           .variable = gsub("xt", "Recovery time", .variable)) %>% 
    filter(.variable!="Intercept",
           .variable!="Dimension") %>% 
    median_qi(.width=c(.5,.8,.95)) %>% 
    ggplot(aes(x = .variable, y = .value)) +
    geom_pointinterval(aes(ymin = .lower, 
                           ymax = .upper),
                       interval_size_range = c(0.5, 2), 
                       colour="grey35") +
    geom_hline(yintercept = 0, linetype="dashed", colour="grey40")+
    stat_slab(data=dat, alpha=0.2, aes(fill=.variable))+
    labs(x="", y = "Posterior estimate \n slope with Generation time") +
    theme(legend.position = "none"))

# All of them 

# mgenTa <- brm(scale(GenT) ~ scale(xt) + (1|gr(animal, cov = A)) + (1|popID),
#           #iter = iter, thin = thin, warmup = warmup,
#           data = smallandata,
#           family = gaussian(), 
#           data2 = list(A = A),
#           cores = 8)


(g1 <- smallandata %>% 
ggplot(aes(xt, GenT, fill=Class, colour=Class))+
  geom_point(shape=21, alpha=.4, size=5) + 
  stat_smooth(method = "lm", se=F) +
  xlab("Recovery time")+
  ylab("Generation time") +
  ggtitle("Animals"))

library(ggrepel)

(g2 <- smallandata %>% 
  filter(Class=="Mammalia") %>% 
  ggplot(aes(xt, GenT, fill=Order, colour=Order))+
  geom_point(shape=21, alpha=.4, size=5) + 
  stat_smooth(method = "lm", se=F) +
  geom_text_repel(aes(label=species), colour="black", size=1.2)+
  xlab("Recovery time")+
  ylab("Generation time")+
  ggtitle("Mammals"))

final <- gridExtra::grid.arrange(g1, g2, ncol = 2)

ggsave(final,
       filename = "GenTRecov.pdf",
       width = 12, height = 6,
       path = ResultPath)

# Old 

# # Generation time vs time of recovery ----
# # Run the model 
# 
# m1.1 <- brm(xt ~ GenT+Dimension + (1|gr(animal, cov = A)) + (1|popID),
#             iter = iter, thin = thin, warmup = warmup,
#             data = smallplandata,
#             family = gaussian(), 
#             data2 = list(A = A),
#             cores = 8)
# 
# # Model diagnostics
# 
# #launch_shinystan(m1.1)
# 
# # Phylogenetic signal
# 
# hyp <-  paste("sd_animal__Intercept^2 /", 
#               "(sd_animal__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0")
# 
# (hyp <- hypothesis(m1.1, hyp, class = NULL))
# 
# # Look to the overall model
# 
# summary(m1.1)
# 
# # Resistance ----
# # Run the model 
# 
# m1.2 <- brm(rlwr ~ GenT+Dimension + (1|gr(animal, cov = A)) + (1|species),
#             iter = iter, thin = thin,warmup = warmup,
#             data = smallplandata,
#             family = gaussian(), 
#             data2 = list(A = A),
#             cores = 8)
# 
# # Model diagnostics
# 
# #launch_shinystan(m1.1)
# 
# # Phylogenetic signal
# 
# hyp <-  paste("sd_animal__Intercept^2 /", 
#               "(sd_animal__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0")
# 
# (hyp <- hypothesis(m1.2, hyp, class = NULL))
# 
# # Look to the overall model
# 
# summary(m1.2)
# 
# # Compensation ----
# 
# m1.3 <- brm(rupr ~ GenT+Dimension + (1|gr(animal, cov = A)) + (1|species),
#             iter = iter, thin = thin, warmup = warmup,
#             data = smallplandata,
#             family = gaussian(), 
#             data2 = list(A = A),
#             cores = 8)
# 
# # Model diagnostics
# 
# #launch_shinystan(m1.1)
# 
# # Phylogenetic signal
# 
# # Look to the overall model
# 
# summary(m1.3)
# Time to recover ----

# m2.1.1 <- brm(xt ~ GenT+log(adult_body_mass_g+1)+ Dimension+ (1|gr(animal, cov = A)) + (1|species),
#               #iter = iter, thin = thin, warmup = warmup,
#               data = smallandata2,
#               family = gaussian(), 
#               data2 = list(A = A),
#               cores = 8)
# 
# m2.1.2 <- brm(xt ~ GenT*log(adult_body_mass_g+1)+ Dimension+ (1|gr(animal, cov = A)) + (1|species),
#               #iter = iter, thin = thin, warmup = warmup,
#               data = smallandata2,
#               family = gaussian(), 
#               data2 = list(A = A),
#               cores = 8)
# 
# m2.1.3 <- brm(xt ~ log(adult_body_mass_g+1)+ Dimension+ (1|gr(animal, cov = A)) + (1|species),
#               #iter = iter, thin = thin, warmup = warmup,
#               data = smallandata2,
#               family = gaussian(), 
#               data2 = list(A = A),
#               cores = 8)
# 
# m2.1.4 <- brm(xt ~ GenT+ Dimension+ (1|gr(animal, cov = A)) + (1|species),
#               #iter = iter, thin = thin, warmup = warmup,
#               data = smallandata2,
#               family = gaussian(), 
#               data2 = list(A = A),
#               cores = 8)
# 
# 
# loo(m2.1.1, m2.1.2, m2.1.3, m2.1.4)
# 
# # Phylogenetic signal
# 
# hyp <-  paste("sd_animal__Intercept^2 /", 
#               "(sd_animal__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0")
# 
# (hyp <- hypothesis(m2.1, hyp, class = NULL))
# 
# # Look to the overall model
# 
# summary(m2.1)
# 
# # Resistance ----
# # Run the model 
# 
# m2.2 <- brm(rlwr ~ GenT+Dimension + (1|gr(animal, cov = A)) + (1|species),
#             iter = iter, thin = thin, warmup = warmup,
#             data = smallandata,
#             family = gaussian(), 
#             data2 = list(A = A),
#             cores = 8)
# 
# # Model diagnostics
# 
# #launch_shinystan(m1.1)
# 
# # Phylogenetic signal
# 
# hyp <-  paste("sd_animal__Intercept^2 /", 
#               "(sd_animal__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0")
# 
# (hyp <- hypothesis(m2.2, hyp, class = NULL))
# 
# # Look to the overall model
# 
# summary(m2.2)
# 
# # Compensation ----
# 
# m2.3 <- brm(rupr ~ GenT+Dimension + (1|gr(animal, cov = A)) + (1|species),
#             iter = iter, thin = thin, warmup = warmup,
#             data = smallandata,
#             family = gaussian(), 
#             data2 = list(A = A),
#             cores = 8)
# 
# # Model diagnostics
# 
# #launch_shinystan(m1.1)
# 
# # Phylogenetic signal
# 
# hyp <-  paste("sd_animal__Intercept^2 /", 
#               "(sd_animal__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0")
# 
# (hyp <- hypothesis(m2.3, hyp, class = NULL))
# 
# # Look to the overall model
# 
# summary(m2.3)
