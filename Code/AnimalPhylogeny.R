# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Phylogeny.R         
# - DATE:        15/03/2020
# - DESCRIPTION: Code to explore the demographic resilience and its association 
#                with life history data 
# - AUTHORS:     Pol Capdevila Lanzaco
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything


#Libraries

library(data.table)
library(datelife)
library(dplyr)
library(ape)
library(caper)
library(phytools)
library(rotl)
library(geiger)
library(taxize)

#Working directories

path = gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath = paste0(path,"/Data")
ResultPath = paste0(path, "/Results") 

#Load  data----

setwd(DataPath)

#Load transiets
load("transDataPost.RData") 

# Animals ######

taxa <- unique(comadreUse$SpeciesAccepted) #subset the unique names
taxa <- na.omit(taxa)
taxa <- gsub("_"," ", taxa)
taxa <- gsub("[0-9]+", "", taxa)
taxa <- gsub("\\svar.+", "", taxa)#remove varieties because of no matches in several points
taxa<- gsub("\\(.*","",taxa)
taxa <- unique(taxa) #subset the unique names

#Search for a match in the internal names of OTL

resolved_names <- tnrs_match_names(taxa)

#Step 2: Get the tree corresponding to our taxa----

resolved_names <- subset(resolved_names, ott_id!=1065138)

#drop one more species

rn<- na.omit(resolved_names$ott_id)

tree <- tol_induced_subtree(ott_ids = rn)
tree$tip.label <- gsub("_ott.*$", "", tree$tip.label)
tree$tip.label <- gsub("_", " ", tree$tip.label)

plot(tree)

#Root the tree resolving the basal the basal Di/Trichotomies.

tree <- multi2di(tree, random=FALSE)
is.rooted(tree)

# random=FALSE collapse the zero length branches into 'true 
# multichotomies' 

#Node labels must be unique

any(duplicated(tree$node.label))

#Another possible issue could be that the disntance between taxons
#Could be 0, and this is not enabled for some analyses 
#Thus we can solve this by adding an 0.001 to all tips

tree$edge.length <-  tree$edge.length + 0.001


#We need that the tree must be binary

is.binary(tree)

#Transform the tree into an ultrametric one

is.ultrametric(tree)
tree <- chronos(tree) 
is.ultrametric(tree)
class(tree) <- "phylo"


#read in phylogenetic tree 
tree$tip.label <- gsub("_", " ", tree$tip.label)
tree$tip.label <- gsub("[0-9]+", "", tree$tip.label)
tree$tip.label <- gsub(" ott", "", tree$tip.label)

# Step 3: Time calibration ################################################## 

phylo <- datelife_search(input = tree, summary_format = "phylo_median") # Fails with sdm

#Change name 

antree <- phylo

# Save the phlogeny

setwd(DataPath)
save(antree, file="animalsData2.RData")

write.tree(antree, "antree.tre")
