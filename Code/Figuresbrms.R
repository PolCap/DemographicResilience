# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Figures.R         
# - DATE:        17/06/2020
# - DESCRIPTION:  
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# -------------------------------------------------------------------------------------- #
rm(list=ls(all=TRUE)) #remove everything

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
library(magrittr)
library(forcats)
library(modelr)
library(emmeans)
library(cowplot)
library(rphylopic)
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
                  plot.margin = margin(0.1, 0, 0, 0, "cm"),
                  plot.title = element_text(size = 15, vjust = 1, 
                                            hjust = -0.10, 
                                            face = "bold"),
                  legend.text = element_text(size = 12, face = "italic"),
                  legend.title = element_blank(),
                  legend.background = element_rect(color = NA))

theme_set(theme_p) 

#Working directories

path <- gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))
DataPath  <-  paste0(path,"/Data")
ResultPath <-  paste0(path, "/Results") 

#Load and prepare resilience data ----------------------------------------------

load(paste0(DataPath,"/ResData.RData"))

# Join animals and plants

datat <- plyr::rbind.fill(smallandata, smallplandata)

# Create colours for animals and plants that are going to be used in the
# graphs

colors <- c("#187F99", "#E58E3C")

# We create a summary of the data that we store as Table S1

TableS1 <- datat %>% 
  group_by(Kingdom, Class,Order, SpeciesAccepted) %>% 
  dplyr::summarise(Npop=n())

# TableS1 <- datat %>% 
#   group_by(Kingdom, Class, Order) %>%
#   distinct(species) %>% 
#   dplyr::summarise(Nsp=n()) %>% 
#   left_join(TableS1)

setwd(ResultPath)
write.csv(TableS1, file="TableS1.csv")

# Figure 3: Trade-offs ------------------------------------------------------------

# Load the models 

load(paste0(ResultPath, "/Trade-offModels2.RData"))

# We will store the residual correlations

posts <-
  tibble(model = str_c("mcor_", c("a","p"))) %>% 
  mutate(mcor  = purrr::map(model, get)) %>% 
  mutate(post = purrr::map(mcor, posterior_samples)) %>% 
  unnest(post) %>% 
  dplyr::select(model, 
                rescor__scalext__scalerupr,
                rescor__scalext__scalerlwr, 
                rescor__scalerupr__scalerlwr) %>% 
  mutate(Kingdom=ifelse(model=="mcor_a", "Animals", "Plants"))

# Panel a: Resistance and recovery time ---------------------------------------- 

# Raw values

(g1a <- smallandata %>%  
   ggplot(aes(x = xt, y = rlwr)) +
   geom_point(size=3, alpha=.35, colour=colors[1]) +
   xlim(min(datat$xt),6.8)+
    labs(x="", y="Resistance"))


# Raw values

(g1b <- smallplandata %>%  
    ggplot(aes(x = xt, y = rlwr)) +
    geom_point(size=3, alpha=.35, colour=colors[2]) +
    xlim(min(datat$xt),6.8)+
    labs(x="Recovery time", y="Resistance"))

# Distribution of the correlation estimate

#Calculate the mean value of the correlation 
label <- posts %>% 
  group_by(Kingdom) %>% 
  summarise(value=round(mean(rescor__scalext__scalerlwr), 2))

# First we plot the distribution

(g1a_ins <- posts %>% filter(Kingdom=="Animals") %>%  
    ggplot(aes(x = rescor__scalext__scalerlwr)) +
    stat_density(alpha=.25, position = "identity",
                 fill=colors[1]) +
    geom_vline(xintercept = 0, linetype="dashed", colour="grey50")+
    labs(x=expression(paste("Residual correlation (",rho,")")), 
         y="") + 
    xlim(-0.55, 0.95)+
    theme(legend.position = "none",
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=10),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA)))

# First we plot the distribution

(g1b_ins <- posts %>% filter(Kingdom=="Plants") %>%  
    ggplot(aes(x = rescor__scalext__scalerlwr)) +
    stat_density(alpha=.25, position = "identity",
                 fill=colors[2]) +
    geom_vline(xintercept = 0, linetype="dashed", colour="grey50")+
    labs(x=expression(paste("Residual correlation (",rho,")")), 
         y="") + 
    xlim(-0.55, 0.95)+
    theme(legend.position = "none",
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=10),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA)))

#Then, we build the data displayed on the plo

p.data_a <- ggplot_build(g1a_ins)$data[[1]]
p.data_b <- ggplot_build(g1b_ins)$data[[1]]

# We extract the max density row for each group 
# Change the factor and we join with label

p.text_a <- p.data_a %>% 
  summarise(max=max(ymax)) %>% 
  mutate(Kingdom= "Animals") %>% 
  left_join(label)

p.text_b <- p.data_b %>% 
  summarise(max=max(ymax)) %>% 
  mutate(Kingdom= "Plants") %>% 
  left_join(label)

# now add the text layer to the plot
inset_a <- g1a_ins + geom_text_repel(data = p.text_a,
                      aes(label = paste("rho==", value), x = value, y = max), 
                      colour=colors[1], size=5, parse = TRUE) 

inset_b <- g1b_ins + geom_text_repel(data = p.text_b,
                               aes(label = paste("rho==", value), x = value, y = max), 
                               colour=colors[2], size=5, parse = TRUE) 

# Combine both

(ga <- ggdraw(g1a+theme(legend.position = "none")) +
    draw_plot(inset_a, .58, .15, .4, .4))

(gd <- ggdraw(g1b+theme(legend.position = "none")) +
    draw_plot(inset_b, .58, .15, .4, .4))

# Panel b: Resistance and compensation  ---------------------------------------- 

# Raw values

(g2a <- smallandata %>%  
   ggplot(aes(x = rupr, y = rlwr)) +
   geom_point(size=3, alpha=.35, colour=colors[1]) +
   xlim(min(datat$rupr),6.8)+
   labs(x="", y="Resistance"))

# Raw values

(g2b <- smallplandata %>%  
    ggplot(aes(x = rupr, y = rlwr)) +
    geom_point(size=3, alpha=.35, colour=colors[2]) +
    xlim(min(datat$rupr),6.8)+
    labs(x="Compensation", y="Resistance"))

# Distribution of the correlation estimate

#Calculate the mean value of the correlation 

label <- posts %>% 
  group_by(Kingdom) %>% 
  summarise(value=round(mean(rescor__scalerupr__scalerlwr), 2))

# First we plot the distribution

(g2a_ins <- posts %>% filter(Kingdom=="Animals") %>%  
    ggplot(aes(x = rescor__scalerupr__scalerlwr)) +
    stat_density(alpha=.25, position = "identity",
                 fill=colors[1]) +
    geom_vline(xintercept = 0, linetype="dashed", colour="grey50")+
    labs(x=expression(paste("Residual correlation (",rho,")")), 
         y="") + 
    xlim(-0.55, 0.95)+
    theme(legend.position = "none",
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=10),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA)))

# First we plot the distribution

(g2b_ins <- posts %>% filter(Kingdom=="Plants") %>%  
    ggplot(aes(x = rescor__scalerupr__scalerlwr)) +
    stat_density(alpha=.25, position = "identity",
                 fill=colors[2]) +
    geom_vline(xintercept = 0, linetype="dashed", colour="grey50")+
    labs(x=expression(paste("Residual correlation (",rho,")")), 
         y="") + 
    xlim(-0.55, 0.95)+
    theme(legend.position = "none",
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=10),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA)))

#Then, we build the data displayed on the plo

p.data_a <- ggplot_build(g2a_ins)$data[[1]]
p.data_b <- ggplot_build(g2b_ins)$data[[1]]

# We extract the max density row for each group 
# Change the factor and we join with label

p.text_a <- p.data_a %>% 
  summarise(max=max(ymax)) %>% 
  mutate(Kingdom= "Animals") %>% 
  left_join(label)

p.text_b <- p.data_b %>% 
  summarise(max=max(ymax)) %>% 
  mutate(Kingdom= "Plants") %>% 
  left_join(label)

# now add the text layer to the plot
inset_a <- g2a_ins + geom_text_repel(data = p.text_a,
                               aes(label = paste("rho==", value), x = value, y = max), 
                               colour=colors[1], size=5, parse = TRUE) 

inset_b <- g2b_ins + geom_text_repel(data = p.text_b,
                               aes(label = paste("rho==", value), x = value, y = max), 
                               colour=colors[2], size=5, parse = TRUE) 

# Combine both

(gb <- ggdraw(g2a+theme(legend.position = "none")) +
    draw_plot(inset_a, .58, .15, .4, .4))

(ge <- ggdraw(g2b+theme(legend.position = "none")) +
    draw_plot(inset_b, .58, .15, .4, .4))

# Panel c: Compensation and recovery time--------------------------------------- 
# Raw values

(g3a <- smallandata %>%  
   ggplot(aes(x = xt, y = rupr)) +
   geom_point(size=3, alpha=.35, colour=colors[1]) +
   xlim(min(datat$xt),6.8)+
   labs(x="", y="Compensation"))


# Raw values

(g3b <- smallplandata %>%  
    ggplot(aes(x = xt, y = rupr)) +
    geom_point(size=3, alpha=.35, colour=colors[2]) +
    xlim(min(datat$xt),6.8)+
    labs(x="Recovery time", y="Compensation"))

# Distribution of the correlation estimate

#Calculate the mean value of the correlation 
label <- posts %>% 
  group_by(Kingdom) %>% 
  summarise(value=round(mean(rescor__scalext__scalerupr), 2))

# First we plot the distribution

(g3a_ins <- posts %>% filter(Kingdom=="Animals") %>%  
    ggplot(aes(x = rescor__scalext__scalerupr)) +
    stat_density(alpha=.25, position = "identity",
                 fill=colors[1]) +
    geom_vline(xintercept = 0, linetype="dashed", colour="grey50")+
    labs(x=expression(paste("Residual correlation (",rho,")")), 
         y="") + 
    xlim(-0.55, 0.95)+
    theme(legend.position = "none",
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=10),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA)))

# First we plot the distribution

(g3b_ins <- posts %>% filter(Kingdom=="Plants") %>%  
    ggplot(aes(x = rescor__scalext__scalerupr)) +
    stat_density(alpha=.25, position = "identity",
                 fill=colors[2]) +
    geom_vline(xintercept = 0, linetype="dashed", colour="grey50")+
    labs(x=expression(paste("Residual correlation (",rho,")")), 
         y="") + 
    xlim(-0.55, 0.95)+
    theme(legend.position = "none",
          axis.text.x = element_text(size=10),
          axis.title.x = element_text(size=10),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.line.y = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA)))

#Then, we build the data displayed on the plo

p.data_a <- ggplot_build(g3a_ins)$data[[1]]
p.data_b <- ggplot_build(g3b_ins)$data[[1]]

# We extract the max density row for each group 
# Change the factor and we join with label

p.text_a <- p.data_a %>% 
  summarise(max=max(ymax)) %>% 
  mutate(Kingdom= "Animals") %>% 
  left_join(label)

p.text_b <- p.data_b %>% 
  summarise(max=max(ymax)) %>% 
  mutate(Kingdom= "Plants") %>% 
  left_join(label)

# now add the text layer to the plot

inset_a <- g3a_ins + geom_text_repel(data = p.text_a,
                               aes(label = paste("rho==", value), x = value, y = max), 
                               colour=colors[1], size=5, parse = TRUE) 

inset_b <- g3b_ins + geom_text_repel(data = p.text_b,
                               aes(label = paste("rho==", value), x = value, y = max), 
                               colour=colors[2], size=5, parse = TRUE) 

# Combine both

(gc <- ggdraw(g3a+theme(legend.position = "none")) +
    draw_plot(inset_a, .58, .15, .4, .4))

(gf <- ggdraw(g3b+theme(legend.position = "none")) +
    draw_plot(inset_b, .58, .15, .4, .4))

# General plot ---------------------------------------
# Combine the plots
(figure3 <- (ga|gb|gc)/(gd|ge|gf)+ plot_annotation(tag_levels = 'a'))

# Save

ggsave(figure3,
       filename = "Figure 3.pdf",
       width = 12, height = 8,
       path = ResultPath)


# Check the confidence intervals 

posts %>% 
  group_by(Kingdom) %>% 
  median_qi(rescor__scalext__scalerupr, .width = c(.5, .95)) 

# Figure 4: Life histories vs resilience -----------------------------------------------------------

# Load the models 

load(paste0(ResultPath, "/ModelsGenTime2.RData"))

# Compensation 

#Store the lines 

post_com  <- posterior_samples(mGenTa, pars = "scalerupr", subset=1:250) %>% 
  mutate(intercept=b_scalerupr_Intercept,
         GenT=b_scalerupr_scaleGenT,
         Fec=b_scalerupr_scaleFec)

meanpost_com  <- post_com %>%
  summarise(intercept=mean(intercept),
            GenT=mean(GenT),
            Fec=mean(Fec))

post_comp  <- posterior_samples(mGenTp, pars = "scalerupr", subset=1:250) %>% 
  mutate(intercept=b_scalerupr_Intercept,
         GenT=b_scalerupr_scaleGenT,
         Fec=b_scalerupr_scaleFec)

meanpost_comp  <- post_comp %>%
  summarise(intercept=mean(intercept),
            GenT=mean(GenT),
            Fec=mean(Fec))

# Plots 

(gg1 <-ggplot(datat, aes(x = scale(GenT), y = scale(rupr))) +
    geom_point(size=3, alpha=.35, aes(colour=Kingdom)) + 
    geom_abline(intercept = post_com$intercept, 
                slope     = post_com$GenT,
                alpha = 0.05, colour=colors[1]) +
     geom_abline(intercept = meanpost_com$intercept, 
                 slope     = meanpost_com$GenT,
                 colour=colors[1], size=2) +
     geom_abline(intercept = post_comp$intercept, 
                 slope     = post_comp$GenT,
                 alpha = 0.05, colour=colors[2]) +
     geom_abline(intercept = meanpost_comp$intercept, 
                 slope     = meanpost_comp$GenT,
                 colour=colors[2], size=2) +
     scale_colour_manual(values = colors) +
    labs(x="scaled(Generation time)", y="scaled(Compensation)"))

(gg4 <-ggplot(datat, aes(x = scale(Fec), y = scale(rupr))) +
    geom_point(size=3, alpha=.35, aes(colour=Kingdom)) + 
    geom_abline(intercept = post_com$intercept, 
                slope     = post_com$Fec,
                alpha = 0.05, colour=colors[1]) +
    geom_abline(intercept = meanpost_com$intercept, 
                slope     = meanpost_com$Fec,
                colour=colors[1], size=2) +
    geom_abline(intercept = post_comp$intercept, 
                slope     = post_comp$Fec,
                alpha = 0.05, colour=colors[2]) +
    geom_abline(intercept = meanpost_comp$intercept, 
                slope     = meanpost_comp$Fec,
                colour=colors[2], size=2) +
    scale_colour_manual(values = colors) +
    labs(x="scaled(Mean reproductive /noutput)", y="scaled(Compensation)"))

# Resistance -----

post_res  <- posterior_samples(mGenTa, pars = "scalerlwr", subset=1:250) %>% 
  mutate(intercept=b_scalerlwr_Intercept,
         GenT=b_scalerlwr_scaleGenT,
         Fec=b_scalerlwr_scaleFec)

meanpost_res  <- post_res %>%
  summarise(intercept=mean(intercept),
            GenT=mean(GenT),
            Fec=mean(Fec))

post_resp  <- posterior_samples(mGenTp, pars = "scalerlwr", subset=1:250) %>% 
  mutate(intercept=b_scalerlwr_Intercept,
         GenT=b_scalerlwr_scaleGenT,
         Fec=b_scalerlwr_scaleFec)

meanpost_resp  <- post_resp %>%
  summarise(intercept=mean(intercept),
            GenT=mean(GenT),
            Fec=mean(Fec))


# this is the predicted line of multiple linear regression

(gg2 <-ggplot(datat, aes(x = scale(GenT), y = scale(rlwr))) +
    geom_point(size=3, alpha=.35, aes(colour=Kingdom)) + 
    geom_abline(intercept = post_res$intercept, 
                slope     = post_res$GenT,
                alpha = 0.05, colour=colors[1]) +
    geom_abline(intercept = meanpost_res$intercept, 
                slope     = meanpost_res$GenT,
                colour=colors[1], size=2) +
    geom_abline(intercept = post_resp$intercept, 
                slope     = post_resp$GenT,
                alpha = 0.05, colour=colors[2]) +
    geom_abline(intercept = meanpost_resp$intercept, 
                slope     = meanpost_resp$GenT,
                colour=colors[2], size=2) +
    scale_colour_manual(values = colors) +
    labs(x="scaled(Generation time)", y="scaled(Resistance)"))

(gg5 <-ggplot(datat, aes(x = scale(Fec), y = scale(rlwr))) +
    geom_point(size=3, alpha=.35, aes(colour=Kingdom)) + 
    geom_abline(intercept = post_res$intercept, 
                slope     = post_res$Fec,
                alpha = 0.05, colour=colors[1]) +
    geom_abline(intercept = meanpost_res$intercept, 
                slope     = meanpost_res$Fec,
                colour=colors[1], size=2) +
    geom_abline(intercept = post_resp$intercept, 
                slope     = post_resp$Fec,
                alpha = 0.05, colour=colors[2]) +
    geom_abline(intercept = meanpost_resp$intercept, 
                slope     = meanpost_resp$Fec,
                colour=colors[2], size=2) +
    scale_colour_manual(values = colors) +
    labs(x="scaled(Mean reproductive \noutput)", y="scaled(Resistance)"))

# Recovery time -----

post_rec  <- posterior_samples(mGenTa, pars = "scalext", subset=1:250) %>% 
  mutate(intercept=b_scalext_Intercept,
         GenT=b_scalext_scaleGenT,
         Fec=b_scalext_scaleFec)

meanpost_rec  <- post_rec %>%
  summarise(intercept=mean(intercept),
            GenT=mean(GenT),
            Fec=mean(Fec))

post_recp  <- posterior_samples(mGenTp, pars = "scalext", subset=1:250) %>% 
  mutate(intercept=b_scalext_Intercept,
         GenT=b_scalext_scaleGenT,
         Fec=b_scalext_scaleFec)

meanpost_recp  <- post_recp %>%
  summarise(intercept=mean(intercept),
            GenT=mean(GenT),
            Fec=mean(Fec))

# Plots 

(gg3 <-ggplot(datat, aes(x = scale(GenT), y = scale(xt))) +
    geom_point(size=3, alpha=.35, aes(colour=Kingdom)) + 
    geom_abline(intercept = post_rec$intercept, 
                slope     = post_rec$GenT,
                alpha = 0.05, colour=colors[1]) +
    geom_abline(intercept = meanpost_rec$intercept, 
                slope     = meanpost_rec$GenT,
                colour=colors[1], size=2) +
    geom_abline(intercept = post_recp$intercept, 
                slope     = post_recp$GenT,
                alpha = 0.05, colour=colors[2]) +
    geom_abline(intercept = meanpost_recp$intercept, 
                slope     = meanpost_recp$GenT,
                colour=colors[2], size=2) +
    scale_colour_manual(values = colors) +
    labs(x="scaled(Generation time)", y="scaled(Recovery time)"))

(gg6 <-ggplot(datat, aes(x = scale(Fec), y = scale(xt))) +
    geom_point(size=3, alpha=.35, aes(colour=Kingdom)) + 
    geom_abline(intercept = post_rec$intercept, 
                slope     = post_rec$Fec,
                alpha = 0.05, colour=colors[1]) +
    geom_abline(intercept = meanpost_rec$intercept, 
                slope     = meanpost_rec$Fec,
                colour=colors[1], size=2) +
    geom_abline(intercept = post_recp$intercept, 
                slope     = post_recp$Fec,
                alpha = 0.05, colour=colors[2]) +
    geom_abline(intercept = meanpost_recp$intercept, 
                slope     = meanpost_recp$Fec,
                colour=colors[2], size=2) +
    scale_colour_manual(values = colors) +
    labs(x="scaled(Mean reproductive output)", y="scaled(Recovery time)"))

# General plot ---------------------------------------

#Create first rows

(first_rows <- plot_grid(gg1+ theme(legend.position = "none")+ 
                           ggtitle("a") + xlab(""),
                         gg4+ theme(legend.position = "none")+
                           ggtitle("d")+ labs(x="", y=""),
                         gg2+theme(legend.position = "none")+
                           ggtitle("b")+xlab(""), 
                         gg5+theme(legend.position = "none") + 
                           ggtitle("e") + labs(x="", y=""),
                         gg3+theme(legend.position = "none")+ 
                           ggtitle("c"),
                         gg6+theme(legend.position = "none")+ 
                           ggtitle("f") + ylab(""),
                         ncol = 2, axis = "b"))


# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  gg1 + geom_line() + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom",
          legend.text=element_text(size=16))
)

# add the legend 

(final_plot <- plot_grid(first_rows,
                         legend, rel_heights = c(3,.1), ncol = 1))

ggsave(final_plot,
       filename = "Figure 4.pdf",
       width = 8, height = 10,
       path = ResultPath)

# Suppplementary tables -------------------------------------------------------- 

tb1 <- describe_posterior(mGenTa,ci = 0.95)
tb1$Kingdom <- "Animals"
tb2 <- describe_posterior(mGenTp, ci=0.95)
tb2$Kingdom <- "Plants"

# We join them 

TableS2 <- rbind(tb1, tb2)

# Now we add some variables 

TableS2 <- TableS2 %>%
  mutate(Parameter=gsub("b_", "", Parameter),
         Response=gsub("(scalext_).*", "\\1", Parameter),
         Response=gsub("scalext_", "Recovery time", Response),
         Response=gsub("(scalerlwr_).*", "\\1", Response),
         Response=gsub("scalerlwr_", "Resistance", Response),
         Response=gsub("(scalerupr_).*", "\\1", Response),
         Response=gsub("scalerupr_", "Compensation", Response),
         Parameter=gsub("scale", "", Parameter),
         Parameter=gsub("Fec", "Mean reproductive output", Parameter),
         Parameter=gsub("GenT", "Generation time", Parameter),
         Parameter=gsub("xt_", "", Parameter),
         Parameter=gsub("rlwr_", "", Parameter),
         Parameter=gsub("rupr_", "", Parameter)) %>% 
  dplyr::select(Kingdom, Response, Parameter, Median, CI_low, CI_high, pd, Rhat)

# Save it

setwd(ResultPath)
write.csv(TableS2, file = "TableS2.csv",row.names = F)
