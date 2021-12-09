# --------------------------------------------------------------------------------------- #
# - FILE NAME:   PhylogenyPlot.R         
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

#Working directories

path <-  gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath <- paste0(path,"/Data")
ResultPath <- paste0(path, "/Results") 

#Load and prepare data----

setwd(DataPath)

# Load animals and plants

load("ResData.RData")
#load("Explore.RData")

# Load the phylogenetic signal results

setwd(ResultPath)
load("PhyloSingal.RData")

# Animal Phylogeny ----

antree <- smallantree

# Change the animal data

data <- smallandata %>% 
  group_by(species,Class) %>% 
  summarise(Resistance=mean(rlwr, na.rm = T),
            Recovery=mean(xt, na.rm = T),
            Compensation=mean(rupr, na.rm = T)) 

data$species <- as.factor(data$species)

# Same for plants

plandata <- smallplandata %>%
  group_by(species) %>% 
  summarise(Resistance=mean(rlwr, na.rm = T),
            Recovery=mean(xt, na.rm = T),
            Compensation=mean(rupr, na.rm = T)) 

plandata$species <- as.factor(plandata$species)

# Format resilience data so that it lines up with species names in tree

data$node <- NA
tipnode <- seq_along(antree$tip.label)
names(tipnode) <- antree$tip.label
data$node <- tipnode[data$species] ## convert the tip label to tip node number
nodes <- data$Class
names(nodes) <- data$species
antree$node.label2 <- nodes[antree$tip.label] ## convert the tip label to tip node number

antree$node.label2 <- data$Class[data$species%in%antree$tip.label] #Add node.labels to the tree

antree$tip.label[antree$node.label2=="Mammalia"]

# Start the tree

(tree <- ggtree(antree, layout = "circular")) 

# Re-arrenge the data for the heatmap

matrix <- as.matrix(data[, "Compensation"])
matrix.names <- data$species
rownames(matrix) <- matrix.names

# Plot it

(heat.tree <- gheatmap(tree, matrix, offset = 0, 
                       width = 0.1,
                        colnames = F) + 
    scale_fill_gradient("Compensation", low = "grey95",
                        high = "#1F84A3",
                        limits = range(plandata$Compensation, 
                        plandata$Compensation),
                        breaks=c(seq(0, max(plandata$Compensation), by=1.5)))+
    theme(legend.position = "bottom"))

# Same for time to recover 

matrix <- as.matrix(data[, "Resistance"])
matrix.names <- data$species
rownames(matrix) <- matrix.names

# New plot  

heat.tree <- heat.tree + new_scale_fill()
(heat.tree <- gheatmap(heat.tree, matrix, offset = 0.1, 
                       width = 0.1,
                        colnames = F) + 
    scale_fill_gradient("Resistance", low = "grey95",
                        high = "#B33815", 
                        limits = c(0, 1), 
                        breaks=c(0, 0.5, 1))+  
    theme(legend.position = "bottom"))

# Compensation

matrix <- as.matrix(data[, "Recovery"])
matrix.names <- data$species
rownames(matrix) <- matrix.names

# New plot  

heat.tree <- heat.tree + new_scale_fill()

(heat.tree <- gheatmap(heat.tree, matrix, offset = 0.2, width = 0.1,
                       colnames = F) + 
    scale_fill_gradient("Recovery time", low = "grey95",
                        high = "#808A7E", 
                        limits=range(plandata$Recovery),
                        breaks=seq(0.5, max(plandata$Recovery), by=1.0))+ 
    theme(legend.position = "bottom"))

# Get the images for the plot from phylopic
# 
# sp <- data.frame()
#  
# for(i in 1:length(antree$tip.label)){
#   tryCatch(x <- phylopic_uid(antree$tip.label[i]),
#            error = function(c) "error")
#   tryCatch(sp <- rbind(sp, x),
#            error = function(c) "error")
# }

load(paste0(DataPath,"/SpPhyloID.RData"))

# Remove the species we are not interested in

d2 <- sp %>% left_join(smallandata[, c("species", "Order", "Class")], 
                      by=c("name"="species")) %>% 
  distinct(name,.keep_all=T) %>%
  filter(!name%in%c("Cercopithecus mitis",
                    "Crocodylus johnsoni",
                    "Papio cynocephalus",
                    "Brachyteles hypoxanthus",
                    "Propithecus verreauxi",
                    "Hippocamelus bisulcus",
                    "Odocoileus virginianus",
                    "Rangifer tarandus",
                    "Ovis canadensis",
                    "Panthera pardus",
                    "Mirounga leonina",
                    "Ursus maritimus",
                    "Ursus americanus",
                    "Ursus arctos",
                    "Lycalopex culpaeus",
                    "Eumetopias jubatus",
                    "Eidolon helvum",
                    "Lontra canadensis",
                    "Alces alces",
                    "Urocyon littoralis",
                    "Elephas maximus",
                    "Vulpes vulpes",
                    "Macaca mulatta",
                    "Phacochoerus aethiopicus",
                    "Halichoerus grypus",
                    "Tamiasciurus hudsonicus",
                    "Emydura macquarii",
                    "Chelodina expansa", 
                    "Malaclemys terrapin",
                    "Chrysemys picta",
                    "Chelydra serpentina",
                    "Anser anser",
                    "Centrocercus minimus",
                    "Certhia americana",
                    "Picoides borealis",
                    "Fulmarus glacialis",
                    "Phoebastria immutabilis",
                    "Anthropoides paradiseus",
                    "Haliaeetus albicilla",
                    "Gyps coprotheres",
                    "Crocodylus johnsoni",
                    "Onychogalea fraenata",
                    "Lagopus muta",
                    "Puffinus auricularis",
                    "Crocodylus johnsoni",
                    "Milvus migrans",
                    "Chrysemys picta",
                    "Gyps coprotheres",
                    "Gorilla beringei",
                    "Cebus capucinus",
                    "Caretta caretta",
                    "Mustela erminea",
                    "Isurus oxyrinchus",
                    "Canis lupus",
                    "Phrynosoma cornutum")) %>% 
  drop_na(uid)

#Modify one of the entrances ###

sp <- rbind(sp,data.frame(name="Epinephelus morio", uid="86c40d81-2613-4bb4-ad57-fb460be56ae5") )

d3 <- sp %>% left_join(smallandata[, c("species", "Order", "Class")], 
                      by=c("name"="species")) %>% 
  distinct(name,.keep_all=T) %>%
  filter(name%in%c("Ursus maritimus")) %>% 
  drop_na(uid)

d4 <- sp %>% left_join(smallandata[, c("species", "Order", "Class")], 
                      by=c("name"="species")) %>% 
  distinct(name,.keep_all=T) %>%
  filter(name%in%c("Epinephelus morio")) %>% 
  drop_na(uid)

# Add pictures

(figA <- heat.tree +
    geom_fruit(data=d2,
               geom=geom_phylopic,
    mapping=aes(y=name, image=uid),
    size=0.06,
    offset=0.52,
    colour="grey15",
    position=position_identityx())+
    geom_fruit(data=d3,
               geom=geom_phylopic,
               mapping=aes(y=name, image=uid),
               size=0.08,
               offset=0.03,
               colour="grey15",
               position=position_identityx())+
     geom_fruit(data=d4,
                geom=geom_phylopic,
                mapping=aes(y=name, image=uid),
                size=0.12,
                offset=0.11,
                colour="grey15",
                position=position_identityx())+
    theme(legend.position = "none")) 

# Subtract the legend  

legend <- get_legend(heat.tree+ theme(legend.title = element_text(size=14),
                                      legend.text = element_text(size=10)))

# Plant phylogeny ----

plantree <- smalltree

# Correct some names (we look for synonims given the scarcity of silouettes)

plantree$tip.label[which(plantree$tip.label=="Callitris columellaris")] <- "Pseudotsuga"
plantree$tip.label[which(plantree$tip.label=="Chamaecrista lineata")] <- "Spartium junceum"
plantree$tip.label[which(plantree$tip.label=="Lathyrus Vernus")] <-"Pisum sativum"
plantree$tip.label[which(plantree$tip.label=="Alnus incana")] <-"Fagales"
plantree$tip.label[which(plantree$tip.label=="Magnolia macrophylla")] <-"Magnolia"
plantree$tip.label[which(plantree$tip.label=="Eryngium alpinum")] <-"Umbelliferae"

smallplandata$species[which(smallplandata$species=="Callitris columellaris")] <- "Pseudotsuga"
smallplandata$species[which(smallplandata$species=="Chamaecrista lineata")] <- "Spartium junceum"
smallplandata$species[which(smallplandata$species=="Lathyrus Vernus")] <-"Pisum sativum"
smallplandata$species[which(smallplandata$species=="Alnus incana")] <-"Fagales"
smallplandata$species[which(smallplandata$species=="Magnolia macrophylla")] <-"Magnolia"
smallplandata$species[which(smallplandata$species=="Eryngium alpinum")] <-"Umbelliferae"


# Format resilience data so that it lines up with species names in tree

# plandata$node <- NA
# tipnode <- seq_along(plantree$tip.label)
# names(tipnode) <- plantree$tip.label
# plandata$node <- tipnode[plandata$species] ## convert the tip label to tip node number
# nodes <- plandata$AngioGymno
# names(nodes) <- plandata$species
# plantree$node.label2 <- nodes[plantree$tip.label] ## convert the tip label to tip node number

# Start the tree

(tree2 <- ggtree(plantree, layout = "circular")) 

# Re-arrenge the data for the heatmap

matrix <- as.matrix(plandata[, "Compensation"])
matrix.names <- plandata$species
rownames(matrix) <- matrix.names

# Plot it

(heat.tree <- gheatmap(tree2, matrix, offset = 0, width = 0.1,
                       colnames = F) + 
    scale_fill_gradient("Compensation", low = "grey95",
                        high = "#1F84A3",
                        limits = range(plandata$Compensation, 
                                       plandata$Compensation),
                        breaks=c(seq(0, max(plandata$Compensation), by=1.5)))+
    theme(legend.position = "bottom"))

# Same for time to recover 

matrix <- as.matrix(plandata[, "Resistance"])
matrix.names <- plandata$species
rownames(matrix) <- matrix.names

# New plot  

heat.tree <- heat.tree + new_scale_fill()

(heat.tree <- gheatmap(heat.tree, matrix, offset = 0.1, width = 0.1,
                       colnames = F) + 
    scale_fill_gradient("Resistance", low = "grey95",
                        high = "#B33815", 
                        limits = c(0, 1), 
                        breaks=c(0, 0.5, 1))+ 
    theme(legend.position = "bottom"))

# Recovery time

matrix <- as.matrix(plandata[, "Recovery"])
matrix.names <- plandata$species
rownames(matrix) <- matrix.names

# New plot  
heat.tree <- heat.tree + new_scale_fill()

(heat.tree <- gheatmap(heat.tree, matrix, offset = 0.2, width = 0.1,
                       colnames = F) + 
    scale_fill_gradient("Recovery time", low = "grey95",
                        high = "#808A7E", 
                        limits=range(plandata$Recovery),
                        breaks=seq(0.5, max(plandata$Recovery), by=1.0))+ 
    
    theme(legend.position = "none"))

# Get the images for the plot from phylopic

# sp2 <- data.frame()
# rm(x)

# Run the loop

# for(i in 1:length(plantree$tip.label)){
#   tryCatch(x <- phylopic_uid(plantree$tip.label[i]),
#            error = function(c) "error") 
#   tryCatch(sp2 <- rbind(sp2, x),
#            error = function(c) "error") 
# }

d2 <- sp2 %>% left_join(smallplandata[, c("species", "Order", "Class")], 
                      by=c("name"="species")) %>% 
  distinct(name,.keep_all=T) %>%
  drop_na(uid) %>% 
  filter(!name%in%c("Ipomopsis aggregata", "Cirsium acaule", 
                    "Magnolia", "Narcissus pseudonarcissus", "Pisum sativum","Spartium junceum")) 

# Add pictures

(figB <- heat.tree +
    #geom_tiplab()+
    geom_fruit(data=d2,
               geom=geom_phylopic,
               mapping=aes(y=name, image=uid),
               size=0.05,
               offset=0.5,
               colour="grey15",
               position=position_identityx())) 

# Final plot #########
# Combine figures

(fig2 <- (figA + theme(plot.margin = margin(-10, -0.28, -10, -0.28, "cm"))|
            figB + theme(plot.margin = margin(-10, -0.28, -10, -0.28, "cm"))) + 
            plot_layout(guides = "collect") & 
            plot_annotation(tag_levels = 'a') & 
    theme(legend.position = "bottom",
          plot.tag = element_text(face = 'bold')))

# Save it 

ggsave(filename = "Figure 2.pdf",
       plot = fig2,
       width = 8, height = 6,
       path = ResultPath)
  
# Save name of the species 
setwd(DataPath)
save(sp, sp2, file="SpPhyloID.RData")
