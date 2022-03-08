# ---------------------------------------------------------------------------- #
# - FILE NAME:   SimulationResults.R         
# - DATE:        10/10/2021
# - DESCRIPTION: Code to re-analyse the data with 
# - AUTHORS:     Pol Capdevila Lanzaco
# ---------------------------------------------------------------------------- #

rm(list = ls(all=T))

library(tidyr)
library(brms)
library(bayestestR)
library(tidybayes)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(Rcompadre)

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
                  plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
                  plot.title = element_text(size = 15, vjust = 1, 
                                            hjust = -0.10, 
                                            face = "bold"),
                  legend.text = element_text(size = 12, face = "italic"),
                  legend.title = element_blank(),
                  legend.background = element_rect(color = NA))

theme_set(theme_p) 


#Working directories

path = gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath = paste0(path,"/Data")
ResultPath = paste0(path, "/Results") 

#Load  and prepare data----

setwd(ResultPath)

# Load models

load("SimulationModels.RData")

# Figur S1: Life histories and demographic resilience ########################## 

# Random matrices --------------------------------------------------------------

data_random <- model_random %>%  
  gather_draws(`b_.*`, regex = TRUE) %>% 
  median_qi(.width = c(.95, .9, .8)) %>%
  mutate(.variable=gsub("b_", "", .variable),
         .variable=gsub("scale", "", .variable),
         .variable=gsub("_", "", .variable),
         resil=gsub("xt", "Recovery time", .variable),
         resil=gsub("rlwr", "Resistance", resil),
         resil=gsub("rupr", "Compensation", resil),
         resil = gsub("(Recovery time).*", "\\1", resil),
         resil = gsub("(Compensation).*", "\\1", resil),
         resil = gsub("(Resistance).*", "\\1", resil),
         resil = gsub("Intercept", "", resil),
         resil = gsub("GenT", "", resil),
         resil = gsub("Fec", "", resil),
         resil = gsub("Dimension", "", resil),
         resil = gsub(":", "", resil),
         resil = gsub("GenT:Fec", "", resil),
         .variable=gsub("xt", "", .variable),
         .variable=gsub("rlwr", "", .variable),
         .variable=gsub("rupr", "", .variable),
         .variable=gsub("Fec", "Mean reproductive output", .variable),
         .variable=gsub("GenT", "Generation time", .variable),
         .variable=gsub("Generation time:Mean reproductive output", 
                        "Interaction", .variable),
         .variable=factor(.variable,levels = c("Interaction",
                                               "Generation time", 
                                               "Mean reproductive output")), 
         dataset = "Random")%>% 
  filter(.variable!="Dimension", .variable!="Intercept")

# Random matrices without shrinkage --------------------------------------------

data_randomp <- model_random_p %>%  
  gather_draws(`b_.*`, regex = TRUE) %>% 
  median_qi(.width = c(.95, .9, .8)) %>%
  mutate(.variable=gsub("b_", "", .variable),
         .variable=gsub("scale", "", .variable),
         .variable=gsub("_", "", .variable),
         resil=gsub("xt", "Recovery time", .variable),
         resil=gsub("rlwr", "Resistance", resil),
         resil=gsub("rupr", "Compensation", resil),
         resil = gsub("(Recovery time).*", "\\1", resil),
         resil = gsub("(Compensation).*", "\\1", resil),
         resil = gsub("(Resistance).*", "\\1", resil),
         resil = gsub("Intercept", "", resil),
         resil = gsub("GenT", "", resil),
         resil = gsub("Fec", "", resil),
         resil = gsub("Dimension", "", resil),
         resil = gsub(":", "", resil),
         resil = gsub("GenT:Fec", "", resil),
         .variable=gsub("xt", "", .variable),
         .variable=gsub("rlwr", "", .variable),
         .variable=gsub("rupr", "", .variable),
         .variable=gsub("Fec", "Mean reproductive output", .variable),
         .variable=gsub("GenT", "Generation time", .variable),
         .variable=gsub("Generation time:Mean reproductive output", 
                        "Interaction", .variable),
         .variable=factor(.variable,levels = c("Interaction",
                                               "Generation time", 
                                               "Mean reproductive output")), 
         dataset = "Random without retrogression")%>% 
  filter(.variable!="Dimension", .variable!="Intercept")

# Empirical data ############################################################### 

load(paste0(ResultPath, "/ModelsGenTime2.RData"))

# Animals ----------------------------------------------------------------------

data_an <- mGenTa %>%  
  gather_draws(`b_.*`, regex = TRUE) %>% 
  median_qi(.width = c(.95, .9, .8)) %>%
  mutate(.variable=gsub("b_", "", .variable),
         .variable=gsub("scale", "", .variable),
         .variable=gsub("_", "", .variable),
         resil=gsub("xt", "Recovery time", .variable),
         resil=gsub("rlwr", "Resistance", resil),
         resil=gsub("rupr", "Compensation", resil),
         resil = gsub("(Recovery time).*", "\\1", resil),
         resil = gsub("(Compensation).*", "\\1", resil),
         resil = gsub("(Resistance).*", "\\1", resil),
         resil = gsub("Intercept", "", resil),
         resil = gsub("GenT", "", resil),
         resil = gsub("Fec", "", resil),
         resil = gsub("Dimension", "", resil),
         resil = gsub(":", "", resil),
         resil = gsub("GenT:Fec", "", resil),
         .variable=gsub("xt", "", .variable),
         .variable=gsub("rlwr", "", .variable),
         .variable=gsub("rupr", "", .variable),
         .variable=gsub("Fec", "Mean reproductive output", .variable),
         .variable=gsub("GenT", "Generation time", .variable),
         .variable=gsub("Generation time:Mean reproductive output", 
                        "Interaction", .variable),
         .variable=factor(.variable,levels = c("Interaction",
                                               "Generation time", 
                                               "Mean reproductive output")),
         dataset="Animals")%>% 
  filter(.variable!="Dimension", .variable!="Intercept") 

# Plants -----------------------------------------------------------------------

data_pl <- mGenTp %>%  
  gather_draws(`b_.*`, regex = TRUE) %>% 
  median_qi(.width = c(.95, .9, .8)) %>%
  mutate(.variable=gsub("b_", "", .variable),
         .variable=gsub("scale", "", .variable),
         .variable=gsub("_", "", .variable),
         resil=gsub("xt", "Recovery time", .variable),
         resil=gsub("rlwr", "Resistance", resil),
         resil=gsub("rupr", "Compensation", resil),
         resil = gsub("(Recovery time).*", "\\1", resil),
         resil = gsub("(Compensation).*", "\\1", resil),
         resil = gsub("(Resistance).*", "\\1", resil),
         resil = gsub("Intercept", "", resil),
         resil = gsub("GenT", "", resil),
         resil = gsub("Fec", "", resil),
         resil = gsub("Dimension", "", resil),
         resil = gsub(":", "", resil),
         resil = gsub("GenT:Fec", "", resil),
         .variable=gsub("xt", "", .variable),
         .variable=gsub("rlwr", "", .variable),
         .variable=gsub("rupr", "", .variable),
         .variable=gsub("Fec", "Mean reproductive output", .variable),
         .variable=gsub("GenT", "Generation time", .variable),
         .variable=gsub("Generation time:Mean reproductive output", 
                        "Interaction", .variable),
         .variable=factor(.variable,levels = c("Interaction",
                                               "Generation time", 
                                               "Mean reproductive output")),
         dataset="Plants")%>% 
  filter(.variable!="Dimension", 
         .variable!="Intercept") 

# Combine the datasets 

data_total <- rbind(data_random, data_randomp,
                    data_an, data_pl)

# Change the order 

data_total <- data_total %>% 
  mutate(dataset=ordered(dataset, levels=c("Animals", 
                                           "Plants", 
                                           "Random",
                                           "Random without retrogression")))

# Define colours 

colors <- c("#187F99", "#E58E3C",   
            "grey45", "grey25")

# Now we plot them 

(figureS1 <- data_total %>% 
    mutate(resil=factor(resil, levels=c("Compensation", 
                                        "Resistance", 
                                        "Recovery time"))) %>% 
    ggplot(aes(y = .variable, x = .value, 
               group=dataset, colour=dataset)) +
    facet_grid(~resil)+
    geom_vline(xintercept = 0, 
               linetype = "dashed", 
               colour="grey35") +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       position =  position_dodge(0.5)) +
    scale_color_manual("", values = colors) +
    labs(x="Effect sizes", 
         y = "") +
    theme(legend.position = "bottom", 
          strip.text = element_text(size=18)))

# Save 

ggsave(figureS1,
       filename = "Figure S1.pdf",
       width = 11, height = 6,
       path = ResultPath)
