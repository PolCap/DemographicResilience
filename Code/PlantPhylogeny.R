# --------------------------------------------------------------------------------------- #
# - FILE NAME:   PlantPhylogeny.R         
# - DATE:        15/03/2020
# - DESCRIPTION: Built plant phylogeny and sort trait data 
# - AUTHORS:     Pol Capdevila Lanzaco
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

#Libraries

library(V.PhyloMaker)
library(ape)
library(caper)
library(phytools)
library(rotl)
library(geiger)
library(dplyr)
library(tidyr)
require(data.table)
library(ggtree)

#Working directories

ResultPath <- "C:/Users/zool2260/Dropbox/2018 Projects/Global patterns of resilience/Results"
DataPath <- "C:/Users/zool2260/Dropbox/2018 Projects/Global patterns of resilience/Data"

#Load  data----

setwd(DataPath)

#Load transiets
load("TransData.RData") 

# Plants ######

#First we transform the Database to a dataframe

plants <- cdb_flatten(compadreUse)

#We subset species, genus and family for the phylogeny

plants <- subset(plants,
                 select = c("SpeciesAccepted", "Genus", "Family"))

#we put the right names in the columns

colnames(plants) <- c("species", "genus", "family")

#Correct species names and remove duplicated species
plants <- plants %>%
  mutate(species= gsub("_"," ", species))%>%
  mutate(species= gsub("[0-9]+","", species))%>%
  mutate(species= gsub("\\svar.+","", species))%>%
  mutate(species= gsub("\\(.*"," ", species))%>%
  mutate(species= sub("^(\\S*\\s+\\S+).*","\\1", species)) %>%  
  distinct(species, .keep_all = TRUE)

#Put the tree in a data.frame format

plants <- as.data.frame(plants)

# Generate the phylogeny using phylo.maker

tree <- phylo.maker(sp.list=plants, 
                      scenarios="S3")

#Explore the tree

(tree <- tree$scenario.3)

#Check the tree
#is it binary?

is.binary(tree)

#Resolve Di/Trichotomies.

tree <- multi2di(tree, random=FALSE)

# random=FALSE collapse the zero length branches into 'true 
# multichotomies' 

is.binary(tree)

#Root the tree resolving the basal the basal Di/Trichotomies.

is.rooted(tree)

#Node labels must be unique

tree <- makeNodeLabel(tree)

#Another possible issue could be that the disntance between taxons
#Could be 0, and this is not enabled for some analyses 
#Thus we can solve this by adding an 0.001 to all tips

tree$edge.length <-  tree$edge.length + 0.001

#read in phylogenetic tree 
tree$tip.label <- gsub("_", " ", tree$tip.label)
tree$tip.label <- gsub("[0-9]+", "", tree$tip.label)
tree$tip.label <- gsub(" ott", "", tree$tip.label)

#Step 4: Prepare trait data ##################################################

compadreUse$species <- gsub("_"," ", compadreUse$SpeciesAccepted)
compadreUse$species <- gsub("\\svar.+", "", compadreUse$species)#remove varieties because of no matches in several points
compadreUse$species <- gsub("[0-9]+", "", compadreUse$species)
compadreUse$species<-gsub("\\(.*","",compadreUse$species)
compadreUse$species <- sub("^(\\S*\\s+\\S+).*", "\\1", compadreUse$species)

#Transform database into data.frame

compadreUse <- cdb_flatten(compadreUse)

# Visualize the tree

ggtree(tree, layout = "fan")

# save the traits and the tree

setwd(DataPath)
write.tree(plantree,"plantree.tre")
