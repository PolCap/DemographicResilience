# --------------------------------------------------------------------------------------- #
# - FILE NAME:   SupplementaryFigures.R         
# - DATE:        27/05/2020
# - DESCRIPTION: Supplementary figures
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything
options(mc.cores = parallel::detectCores())

library(ggplot2) #Now I load the latest version of ggplot2
library(tidyverse)
library(geiger)
library(phytools)
library(caper)
library(tidybayes)
library(rstan)
library(rstanarm)
library(brms)
library(bayesplot)
library(bayestestR)
library(ggdist)
library(cowplot)
library(dplyr)
library(wesanderson)
library(ggrepel)
library(patchwork)

#Personal theme for ggplot

theme_p <-  theme(panel.background=element_blank(), 
                  strip.background=element_blank(), 
                  axis.line=element_line("black"), 
                  axis.ticks=element_line("black"), 
                  axis.text=element_text(colour="black", size=12), 
                  axis.title=element_text(size=16),
                  panel.grid=element_blank(),
                  strip.text=element_text(size=14),
                  legend.key=element_blank(), 
                  axis.line.x=element_line("black"),
                  axis.line.y=element_line("black"),
                  axis.text.x = element_text(size = 14),
                  axis.text.y = element_text(size = 14),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.major.y = element_blank(),
                  plot.margin = unit(c(0, 0, 0, 0), "cm"),
                  plot.title = element_text(size = 15, vjust = 1, 
                                            hjust = -0.10, 
                                            face = "bold"),
                  legend.text = element_text(size = 12, face = "italic"),
                  legend.title = element_blank(),
                  legend.background = element_rect(color = NA))

theme_set(theme_p) 


#Working directories

path <-  gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath <- paste0(path,"/Data")
ResultPath <- paste0(path, "/Results") 


# Body size effects ############################################################

# Load body size and raunkier analyses

load(paste0(ResultPath, "/Supplementary.RData"))

# Create colours for animals and plants that are going to be used in the
# graphs

colors <- c("#187F99", "#E58E3C")

# Panel a: Animal body size, compensation --------------------------------------

post_com  <- posterior_samples(mbodyanimals, 
                               pars = "scalerupr", subset=1:250) %>% 
  mutate(intercept=b_scalerupr_Intercept,
         bodysize=b_scalerupr_scalelog10adult_body_mass_gP1)

meanpost_com  <- post_com %>%
  summarise(intercept=mean(intercept),
            bodysize=mean(bodysize))

(gg1a <-ggplot(smallandata, 
              aes(x = scale(log10(adult_body_mass_g)), 
                  y = scale(rupr))) +
    geom_point(size=2, alpha=.25, colour=colors[1]) + 
    geom_abline(intercept = post_com$intercept, 
                slope     = post_com$bodysize,
                alpha = 0.05, colour=colors[1]) +
    geom_abline(intercept = meanpost_com$intercept, 
                slope     = meanpost_com$bodysize,
                colour=colors[1], size=2) +
    labs(x="Adult body weight (g)", y="Compensation"))

# Panel b: Animal body size, resistance --------------------------------------

post_com  <- posterior_samples(mbodyanimals, 
                               pars = "scalerlwr", subset=1:250) %>% 
  mutate(intercept=b_scalerlwr_Intercept,
         bodysize=b_scalerlwr_scalelog10adult_body_mass_gP1)

meanpost_com  <- post_com %>%
  summarise(intercept=mean(intercept),
            bodysize=mean(bodysize))

(gg1b <-ggplot(smallandata, 
               aes(x = scale(log10(adult_body_mass_g)), 
                   y = scale(rlwr))) +
    geom_point(size=2, alpha=.25, colour=colors[1]) + 
    geom_abline(intercept = post_com$intercept, 
                slope     = post_com$bodysize,
                alpha = 0.05, colour=colors[1]) +
      geom_abline(intercept = meanpost_com$intercept, 
                  slope     = meanpost_com$bodysize,
                 colour=colors[1], size=2) +
    labs(x="Adult body weight (g)", y="Resistance"))

# Panel c: Animal body size, recovery time --------------------------------------

post_com  <- posterior_samples(mbodyanimals, 
                               pars = "scalext", subset=1:250) %>% 
  mutate(intercept=b_scalext_Intercept,
         bodysize=b_scalext_scalelog10adult_body_mass_gP1)

meanpost_com  <- post_com %>%
  summarise(intercept=mean(intercept),
            bodysize=mean(bodysize))

(gg1c <-ggplot(smallandata, 
               aes(x = scale(log10(adult_body_mass_g)), 
                   y = scale(xt))) +
    geom_point(size=2, alpha=.25, colour=colors[1]) + 
    geom_abline(intercept = post_com$intercept, 
                slope     = post_com$bodysize,
                alpha = 0.05, colour=colors[1]) +
    geom_abline(intercept = meanpost_com$intercept, 
                slope     = meanpost_com$bodysize,
                colour=colors[1], size=2) +
    labs(x="Adult body weight (g)", y="Recovery time"))

# Panel a: Plants body size, compensation --------------------------------------

post_com  <- posterior_samples(mbodyplants, 
                               pars = "scalerupr", subset=1:250) %>% 
  mutate(intercept=b_scalerupr_Intercept,
         bodysize=b_scalerupr_scalelog10max_heightP1)

meanpost_com  <- post_com %>%
  summarise(intercept=mean(intercept),
            bodysize=mean(bodysize))

(gg2a <-ggplot(smallplandata, 
               aes(x = scale(log10(max_height)), 
                   y = scale(rupr))) +
    geom_point(size=2, alpha=.25, colour=colors[2]) + 
    geom_abline(intercept = post_com$intercept, 
                slope     = post_com$bodysize,
                alpha = 0.05, colour=colors[2]) +
    geom_abline(intercept = meanpost_com$intercept, 
                slope     = meanpost_com$bodysize,
                colour=colors[2], size=2) +
    labs(x="Maximum body height (m)", y="Compensation"))

# Panel b: Plant body size, resistance --------------------------------------

post_com  <- posterior_samples(mbodyplants, 
                               pars = "scalerlwr", subset=1:250) %>% 
  mutate(intercept=b_scalerlwr_Intercept,
         bodysize=b_scalerlwr_scalelog10max_heightP1)

meanpost_com  <- post_com %>%
  summarise(intercept=mean(intercept),
            bodysize=mean(bodysize))

(gg2b <-ggplot(smallplandata, 
               aes(x = scale(log10(max_height)), 
                   y = scale(rlwr))) +
    geom_point(size=2, alpha=.25, colour=colors[2]) + 
    geom_abline(intercept = post_com$intercept, 
                slope     = post_com$bodysize,
                alpha = 0.05, colour=colors[2]) +
    geom_abline(intercept = meanpost_com$intercept, 
                slope     = meanpost_com$bodysize,
                colour=colors[2], size=2) +
    labs(x="Maximum body height (m)", y="Resistance"))


# Panel c: Plant body size, recovery time --------------------------------------

post_com  <- posterior_samples(mbodyplants, 
                               pars = "scalext", subset=1:250) %>% 
  mutate(intercept=b_scalext_Intercept,
         bodysize=b_scalext_scalelog10max_heightP1)

meanpost_com  <- post_com %>%
  summarise(intercept=mean(intercept),
            bodysize=mean(bodysize))

(gg2c <-ggplot(smallplandata, 
               aes(x = scale(log10(max_height)), 
                   y = scale(xt))) +
    geom_point(size=2, alpha=.25, colour=colors[2]) + 
    geom_abline(intercept = post_com$intercept, 
                slope     = post_com$bodysize,
                alpha = 0.05, colour=colors[2]) +
    geom_abline(intercept = meanpost_com$intercept, 
                slope     = meanpost_com$bodysize,
                colour=colors[2], size=2) +
    labs(x="Maximum body height (m)", y="Recovery time"))


# General plot   ---------------------------------------------------------------

(figureS2 <- (gg1a+xlab(""))+(gg2a+labs(x="",y=""))+
  (gg1b+xlab(""))+(gg2b+labs(x="",y=""))+
  gg1c+(gg2c+ylab(""))+ 
   plot_layout(ncol = 2)+ 
   plot_annotation(tag_levels = "a"))


ggsave(figureS2,
       filename = "Figure S2.pdf",
       width = 8, height = 10,
       path = ResultPath)


# Raukiaer classification ########################################

# Calculate sample size 

sample_size <-plants %>% 
  group_by(Growth.form.Raunkiaer) %>% 
  dplyr::summarize(num=n()) %>% 
  filter(Growth.form.Raunkiaer!="") %>% 
  mutate(raun = ordered(Growth.form.Raunkiaer,
                        levels=c("Hydrophyte",
                                 "Helophyte", 
                                 "Geophyte",
                                 "Hemicryptophyte",
                                 "Chamaephyte",
                                 "Nanophanerophyte",
                                 "Mesophanerophyte",
                                 "Macrophanerophyte")),
         resil="Recovery time",
         resil = ordered(resil,
                         levels=c("Compensation", 
                                  "Resistance", "Recovery time"))) %>% 
  drop_na(raun)


# Now we generate the values for raunkier

(g3 <-  mraunkier %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    mutate(raun = gsub("b_", "", .variable),
           resil = gsub("_.*", "", raun),
           resil = gsub("scalext", "Recovery time", resil),
           resil = gsub("scalerlwr", "Resistance", resil),
           resil = gsub("scalerupr", "Compensation", resil),
           raun = gsub("scalext_", "", raun),
           raun = gsub("scalerlwr_", "", raun),
           raun = gsub("scalerupr_", "", raun),
           raun = gsub("scale", "", raun),
           raun = gsub("Growth.form.Raunkiaer", "", raun)) %>% 
    filter(raun!="Dimension") %>%
    filter(raun!="Intercept") %>%
    filter(raun!="Hydrophyte") %>% 
    mutate(resil = ordered(resil,
                           levels=c("Compensation", 
                                    "Resistance", "Recovery time")),
           raun = ordered(raun,
                           levels=c("Hydrophyte",
                                    "Helophyte", 
                                    "Geophyte",
                                    "Hemicryptophyte",
                                    "Chamaephyte",
                                    "Nanophanerophyte",
                                    "Mesophanerophyte",
                                    "Macrophanerophyte"))) %>%
    filter(!is.na(raun)) %>% 
    ggplot(aes(y = .value, 
               x = raun)) +
    geom_hline(yintercept = 0, linetype="dashed", colour="grey50")+
    geom_violin(aes(colour=raun, 
                    fill=raun), alpha=0.5, colour="white")+
    geom_boxplot(aes(colour=raun), fill="white", width=0.1)+
    facet_wrap(~resil, ncol=1) +
    labs(x="Raunkiaer life form", y = "Posterior estimate") +
    scale_colour_manual("", values=wes_palette("Cavalcanti1",
                                               8,
                                               type = "continuous"))+
    scale_fill_manual("", values = wes_palette("Cavalcanti1",
                                               8,
                                               type = "continuous"))+
    geom_text(data = subset(sample_size, raun!="Hydrophyte"),
                aes(y=-4.5, x=raun,
                    label=paste0("n=", num))) +
    theme(legend.position = "none",
          strip.text = element_text(hjust = 0),
          axis.text.x = element_text(angle = 25, 
                                     hjust=1,
                                     vjust=1)))
# Save the figure

ggsave(g3,
       filename = "Figure S3.pdf",
       width = 8, height = 10,
       path = ResultPath)
# Conservation status ##########################################################

# Load conservation models

load(paste0(ResultPath, "/ConservationModels.RData"))

# Load animals and plants data

load(paste0(DataPath, "/ResData.RData"))

# Join animals and plants

datat <- plyr::rbind.fill(smallandata, smallplandata)

# Transform IUCN category into an ordered factor

datat$category <- ordered(datat$IUCN,
                          c("LC", "NT", "VU", "EN", "CR"))

# Remove species we don't have information about conservation status 

subdatat <- subset(datat, !is.na(category))

# Find the sample size per category and per Kingdom

sample_size <-subdatat %>% 
  mutate(resil="Recovery time",
         threat = ordered(category,
                          levels=c("LC", "NT", "VU", "EN", "CR")),
         resil = ordered(resil,
                         levels=c("Compensation", 
                                  "Resistance", "Recovery time"))) %>%  
  group_by(Kingdom, threat) %>% 
  dplyr::summarize(num=n()) 

# Define colours

iucncolors <- c("#52C148","#C5E51E","#FAEA12","#FD6A32","#CF000A")

# Animals


(g1a <-  mcon_a %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    mutate(threat = gsub("b_", "", .variable),
           resil = gsub("_.*", "", threat),
           resil = gsub("scalext", "Recovery time", resil),
           resil = gsub("scalerlwr", "Resistance", resil),
           resil = gsub("scalerupr", "Compensation", resil),
           threat = gsub("scalext_", "", threat),
           threat = gsub("scalerlwr_", "", threat),
           threat = gsub("scalerupr_", "", threat),
           threat = gsub("category", "", threat)) %>% 
    filter(threat!="b_Dimension") %>%
    mutate(threat = ordered(threat,
                            levels=c("LC", "NT", "VU", "EN", "CR")),
           resil = ordered(resil,
                           levels=c("Compensation", 
                                    "Resistance", "Recovery time"))) %>%
    filter(threat!="<NA>") %>% 
    ggplot(aes(y = .value, 
               x = threat)) +
    geom_hline(yintercept = 0, linetype="dashed", colour="grey50")+
    geom_violin(aes(colour=threat, 
                    fill=threat), alpha=0.5, colour="white")+
    geom_boxplot(aes(colour=threat), 
                 fill="white", width=0.1)+
    facet_wrap(~resil, ncol=1) +
    labs(x="Threat status", y = "Posterior estimate") +
    scale_colour_manual("", values=iucncolors)+
    scale_fill_manual("", values=iucncolors)+
    geom_text(data = subset(sample_size, Kingdom=="Animalia"),
              aes(y=-2, x=threat,
                  label=paste0("n=", num))) +
    theme(legend.position = "none",
          strip.text = element_text(hjust = 0)))

(g1b <-  mcon_p %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    mutate(threat = gsub("b_", "", .variable),
           resil = gsub("_.*", "", threat),
           resil = gsub("scalext", "Recovery time", resil),
           resil = gsub("scalerlwr", "Resistance", resil),
           resil = gsub("scalerupr", "Compensation", resil),
           threat = gsub("scalext_", "", threat),
           threat = gsub("scalerlwr_", "", threat),
           threat = gsub("scalerupr_", "", threat),
           threat = gsub("category", "", threat)) %>% 
    filter(threat!="b_Dimension") %>%
    mutate(threat = ordered(threat,
                            levels=c("LC", "NT", "VU", "EN", "CR")),
           resil = ordered(resil,
                           levels=c("Compensation", 
                                    "Resistance", "Recovery time"))) %>%
    filter(threat!="<NA>") %>% 
    ggplot(aes(y = .value, 
               x = threat)) +
    geom_hline(yintercept = 0, linetype="dashed", colour="grey50")+
    geom_violin(aes(colour=threat, 
                    fill=threat),alpha=0.5, colour="white")+
    geom_boxplot(aes(colour=threat), fill="white", width=0.1)+
    facet_wrap(~resil, ncol=1) +
    labs(x="Threat status", y = "Posterior estimate") +
    scale_colour_manual("", values=iucncolors)+
    scale_fill_manual("", values=iucncolors)+
    geom_text(data = subset(sample_size, Kingdom=="Plantae"),
              aes(y=-2, x=threat,
                  label=paste0("n=", num))) +
    theme(legend.position = "none",
          strip.text = element_text(hjust = 0)))

# Combine figure


(FigS4 <- (g1a|g1b) + plot_annotation(tag_levels = 'a'))


ggsave(FigS4,
       filename = "Figure S4.pdf",
       width = 10, height = 12,
       path = ResultPath)

# Free up some space 


