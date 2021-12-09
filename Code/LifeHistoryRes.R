# --------------------------------------------------------------------------------------- #
# - FILE NAME:   LifeHistoryRes.R         
# - DATE:        04/05/2020
# - DESCRIPTION: Code to perform MCMCglmm analyses of generation time vs demographic
#                resilience. 
# - AUTHORS:     Pol Capdevila Lanzaco
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

library(plyr)
library(tidyr)
library(tidyverse)
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

#Working directories

path = gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath = paste0(path,"/Data")
ResultPath = paste0(path, "/Results") 

#Load  and prepare data----

setwd(DataPath)

# Load animals and plants

load("ResData.RData")

# Load PCA

load(paste0(ResultPath, "/PCA.RData"))

# Join datasets

smallandata$PC1 <- pcaA$x[,"PC1"]
smallandata$PC2 <- pcaA$x[,"PC2"]
smallplandata$PC1 <- pcaP$x[,"PC1"]
smallplandata$PC2 <- pcaP$x[,"PC2"]

# Plants ###################################

#Set modelling parameters

iter = 10000
thin = 0.0005 * iter
warmup = 0.1 * iter

# Set the variance covariance for the phylogeny

A <- ape::vcv.phylo(smalltree)

# Generation time vs time of recovery ----
# Run the model 

mpgen1 <- brm(PC1~ GenT + Dimension + (1|species),
            iter = iter, thin = thin, warmup = warmup,
            data = subset(smallplandata, !is.na(GenT)),
            family = gaussian(), 
            #data2 = list(A = A),
            cores = 8)

mpgen2 <- brm(PC2~ GenT + Dimension + (1|species),
              iter = iter, thin = thin, warmup = warmup,
              data = subset(smallplandata, !is.na(GenT)),
              family = gaussian(), 
              #data2 = list(A = A),
              cores = 8)

mpfec1 <- brm(PC1~Fec+Dimension + (1|species),
              iter = iter, thin = thin, warmup = warmup,
              data = subset(smallplandata, !is.na(Fec)),
              family = gaussian(), 
              #data2 = list(A = A),
              cores = 8)
mpfec2 <- brm(PC2~Fec+Dimension + (1|species),
              iter = iter, thin = thin, warmup = warmup,
              data = subset(smallplandata, !is.na(Fec)),
              family = gaussian(), 
              #data2 = list(A = A),
              cores = 8)
mpbs1 <- brm(PC1~max_height+Dimension + (1|species),
            iter = iter,  thin = thin, warmup = warmup,
            data = subset(smallplandata, !is.na(max_height)),
            family = gaussian(), 
            #data2 = list(A = A),
            cores = 8)

mpbs2 <- brm(PC2~max_height+Dimension + (1|species),
             iter = iter,  thin = thin, warmup = warmup,
             data = subset(smallplandata, !is.na(max_height)),
             family = gaussian(), 
             #data2 = list(A = A),
             cores = 8)

# Animals

magen1 <- brm(PC1~GenT+Dimension + (1|species),
              iter = iter, thin = thin, warmup = warmup,
              data = subset(smallandata, !is.na(GenT)),
              family = gaussian(), 
              #data2 = list(A = A),
              cores = 8)
magen2 <- brm(PC2~GenT+Dimension + (1|species),
             iter = iter, thin = thin, warmup = warmup,
             data = subset(smallandata, !is.na(GenT)),
             family = gaussian(), 
             #data2 = list(A = A),
             cores = 8)

mafec1 <- brm(PC1~Fec+Dimension + (1|species),
              iter = iter, thin = thin, warmup = warmup,
              data = subset(smallandata, !is.na(Fec)),
              family = gaussian(), 
              #data2 = list(A = A),
              cores = 8)
mafec2 <- brm(PC2~Fec+Dimension + (1|species),
              iter = iter, thin = thin, warmup = warmup,
              data = subset(smallandata, !is.na(Fec)),
              family = gaussian(), 
              #data2 = list(A = A),
              cores = 8)

mabs1 <- brm(PC1~adult_body_mass_g+Dimension+ (1|species),
             iter = iter, thin = thin, warmup = warmup,
             data = subset(smallandata, !is.na(adult_body_mass_g)),
             family = gaussian(), 
             #data2 = list(A = A),
             cores = 8)
mabs2 <- brm(PC2~adult_body_mass_g+Dimension+ (1|species),
             iter = iter, thin = thin, warmup = warmup,
             data = subset(smallandata, !is.na(adult_body_mass_g)),
             family = gaussian(), 
             #data2 = list(A = A),
             cores = 8)

# Save the models
setwd(ResultPath)
save(mpgen1, mpfec1, mpbs1,mpgen2, mpfec2, mpbs2,
     magen1, mafec1, mabs1,magen2, mafec2, mabs2,
     file = "LifeHistoryModels.RData")

# Join both

datat <- plyr::rbind.fill(smallandata, smallplandata)

# Load the models for the trade-offs

colors <- c("#187F99", "#E58E3C")

# Panel 1: 

#Predict plants
fe_only <- tibble(GenT = seq(min(smallplandata$GenT, na.rm = T),
                           max(smallplandata$GenT, na.rm = T),
                           length.out=100),
                  Dimension = mean(smallplandata$Dimension) ) %>%
  add_fitted_draws(mpgen1,
                   re_formula = NA,
                   scale = "response", n = 250)

fe_only_mean <- fe_only %>% 
  group_by(GenT) %>%
  summarize(.value = mean(.value))


#Predict animals

fe_only2 <- tibble(GenT = seq(min(smallandata$GenT, na.rm=T),
                            max(smallandata$GenT, na.rm=T),
                            length.out=100),
                   Dimension = mean(smallandata$Dimension) ) %>%
  add_fitted_draws(magen1,
                   re_formula = NA,
                   scale = "response", n = 250)

fe_only_mean2 <- fe_only2 %>% 
  group_by(GenT) %>%
  summarize(.value = mean(.value))

# this is the predicted line of multiple linear regression

(g1 <- ggplot(datat, aes(x = GenT, y = PC1)) +
    geom_point(size=2, alpha=.25, shape=21, 
               colour="black", aes(fill=Kingdom)) +
    geom_line(data=fe_only, aes(x = GenT, y = .value,
                                group = .draw),
              alpha = 0.05, colour=colors[2]) +
    geom_line(data = fe_only_mean,
              aes(x = GenT, y = .value,
                  group = .draw),
              color = colors[2], lwd = 2, group = 1) +
    geom_line(data=fe_only2, aes(x = GenT, y = .value,
                                 group = .draw),
              alpha =0.05, colour=colors[1]) +
    geom_line(data = fe_only_mean2,
              aes(x = GenT, y = .value,
                  group = .draw),
              color = colors[1], lwd = 2, group = 1) +
    scale_fill_manual(values = colors) +
    labs(x="Generation time", y="PC1"))


#Predict plants
fe_only <- tibble(GenT = seq(min(smallplandata$GenT, na.rm = T),
                             max(smallplandata$GenT, na.rm = T),
                             length.out=100),
                  Dimension = mean(smallplandata$Dimension) ) %>%
  add_fitted_draws(mpgen2,
                   re_formula = NA,
                   scale = "response", n = 250)

fe_only_mean <- fe_only %>% 
  group_by(GenT) %>%
  summarize(.value = mean(.value))


#Predict animals

fe_only2 <- tibble(GenT = seq(min(smallandata$GenT, na.rm=T),
                              max(smallandata$GenT, na.rm=T),
                              length.out=100),
                   Dimension = mean(smallandata$Dimension) ) %>%
  add_fitted_draws(magen2,
                   re_formula = NA,
                   scale = "response", n = 250)

fe_only_mean2 <- fe_only2 %>% 
  group_by(GenT) %>%
  summarize(.value = mean(.value))

# this is the predicted line of multiple linear regression

(g2 <- ggplot(datat, aes(x = GenT, y = PC2)) +
    geom_point(size=2, alpha=.25, shape=21, 
               colour="black", aes(fill=Kingdom)) +
    geom_line(data=fe_only, aes(x = GenT, y = .value,
                                group = .draw),
              alpha = 0.05, colour=colors[2]) +
    geom_line(data = fe_only_mean,
              aes(x = GenT, y = .value,
                  group = .draw),
              color = colors[2], lwd = 2, group = 1) +
    geom_line(data=fe_only2, aes(x = GenT, y = .value,
                                 group = .draw),
              alpha =0.05, colour=colors[1]) +
    geom_line(data = fe_only_mean2,
              aes(x = GenT, y = .value,
                  group = .draw),
              color = colors[1], lwd = 2, group = 1) +
    scale_fill_manual(values = colors) +
    labs(x="Generation time", y="PC2"))


#Predict plants
fe_only <- tibble(Fec = seq(min(smallplandata$Fec, na.rm = T),
                             max(smallplandata$Fec, na.rm = T),
                             length.out=100),
                  Dimension = mean(smallplandata$Dimension) ) %>%
  add_fitted_draws(mpfec1,
                   re_formula = NA,
                   scale = "response", n = 250)

fe_only_mean <- fe_only %>% 
  group_by(Fec) %>%
  summarize(.value = mean(.value))


#Predict animals

fe_only2 <- tibble(Fec = seq(min(smallandata$Fec, na.rm=T),
                              max(smallandata$Fec, na.rm=T),
                              length.out=100),
                   Dimension = mean(smallandata$Dimension) ) %>%
  add_fitted_draws(mafec1,
                   re_formula = NA,
                   scale = "response", n = 250)

fe_only_mean2 <- fe_only2 %>% 
  group_by(Fec) %>%
  summarize(.value = mean(.value))

# this is the predicted line of multiple linear regression

(g1 <- ggplot(datat, aes(x = Fec, y = PC1)) +
    geom_point(size=2, alpha=.25, shape=21, 
               colour="black", aes(fill=Kingdom)) +
    geom_line(data=fe_only, aes(x = Fec, y = .value,
                                group = .draw),
              alpha = 0.05, colour=colors[2]) +
    geom_line(data = fe_only_mean,
              aes(x = Fec, y = .value,
                  group = .draw),
              color = colors[2], lwd = 2, group = 1) +
    geom_line(data=fe_only2, aes(x = Fec, y = .value,
                                 group = .draw),
              alpha =0.05, colour=colors[1]) +
    geom_line(data = fe_only_mean2,
              aes(x = Fec, y = .value,
                  group = .draw),
              color = colors[1], lwd = 2, group = 1) +
    scale_fill_manual(values = colors) +
    labs(x="Reproductive output", y="PC1"))

#Predict plants
fe_only <- tibble(Fec = seq(min(smallplandata$Fec, na.rm = T),
                            max(smallplandata$Fec, na.rm = T),
                            length.out=100),
                  Dimension = mean(smallplandata$Dimension) ) %>%
  add_fitted_draws(mpfec2,
                   re_formula = NA,
                   scale = "response", n = 250)

fe_only_mean <- fe_only %>% 
  group_by(Fec) %>%
  summarize(.value = mean(.value))


#Predict animals

fe_only2 <- tibble(Fec = seq(min(smallandata$Fec, na.rm=T),
                             max(smallandata$Fec, na.rm=T),
                             length.out=100),
                   Dimension = mean(smallandata$Dimension) ) %>%
  add_fitted_draws(mafec2,
                   re_formula = NA,
                   scale = "response", n = 250)

fe_only_mean2 <- fe_only2 %>% 
  group_by(Fec) %>%
  summarize(.value = mean(.value))

# this is the predicted line of multiple linear regression

(g2 <- ggplot(datat, aes(x = Fec, y = PC2)) +
    geom_point(size=2, alpha=.25, shape=21, 
               colour="black", aes(fill=Kingdom)) +
    geom_line(data=fe_only, aes(x = Fec, y = .value,
                                group = .draw),
              alpha = 0.05, colour=colors[2]) +
    geom_line(data = fe_only_mean,
              aes(x = Fec, y = .value,
                  group = .draw),
              color = colors[2], lwd = 2, group = 1) +
    geom_line(data=fe_only2, aes(x = Fec, y = .value,
                                 group = .draw),
              alpha =0.05, colour=colors[1]) +
    geom_line(data = fe_only_mean2,
              aes(x = Fec, y = .value,
                  group = .draw),
              color = colors[1], lwd = 2, group = 1) +
    scale_fill_manual(values = colors) +
    labs(x="Reproductive output", y="PC2"))
