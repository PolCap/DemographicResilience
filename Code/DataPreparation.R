# --------------------------------------------------------------------------------------- #
# - FILE NAME:   DataPreparation.R         
# - DATE:        21/05/2020
# - DESCRIPTION: Code to descrive the resilience patterns between species and groups of
#                organisms. 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Libraries

library(Rcompadre)
library(dlookr)
library(tidyr)
library(plyr)
library(phytools)
library(dplyr)
library(geiger)

#Working directories

path = gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath = paste0(path,"/Data")
ResultPath = paste0(path, "/Results") 

# Load animals and plants

#load("IUCNData.RData") 
load(paste0(DataPath,"/IUCNPopDataPost.RData"))

# Data cleaning and preparation ###########################################################

#Remove infinites

andata <- andata[!is.infinite(andata$GenT),]
plandata <- plandata[!is.infinite(plandata$GenT),]

#Lets have a look at the normality of the data 

plot_normality(andata,rlwr, xt, rupr)
plot_normality(plandata, rlwr, xt, rupr)

#We logtransform the transient data 

plandata <- plandata %>% 
  mutate(xt=log(xt+1),
         rupr=log(rupr+1),
         rlwr=abs(log((1 - rlwr)+1)),
         Fec=log(Fec+1),
         GenT=log(GenT+1))
andata <- andata %>% 
  mutate(xt=log(xt+1),
         rupr=log(rupr+1),
         rlwr=abs(log((1 - rlwr)+1)),
         Fec=log(Fec+1),
         GenT=log(GenT+1))

# Remove outliers

#we create a function to go quicker
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# Now we apply the function to the set of transients

plandata$GenT <- remove_outliers(plandata$GenT)
plandata$Fec <- remove_outliers(plandata$Fec)
plandata$xt <- remove_outliers(plandata$xt)
plandata$rlwr <- remove_outliers(plandata$rlwr)
plandata$rupr <- remove_outliers(plandata$rupr)
plandata <- plandata[!is.na(plandata$xt),]
plandata <- plandata[!is.na(plandata$rlwr),]
plandata <- plandata[!is.na(plandata$rupr),]

andata$GenT <- remove_outliers(andata$GenT)
andata$Fec <- remove_outliers(andata$Fec)
andata$xt <- remove_outliers(andata$xt)
andata$rlwr <- remove_outliers(andata$rlwr)
andata$rupr <- remove_outliers(andata$rupr)
andata <- andata[!is.na(andata$xt),]
andata <- andata[!is.na(andata$rlwr),]
andata <- andata[!is.na(andata$rupr),]

#Prepare the trees ###################################################################

#Plants phylogeny

setwd(DataPath)
plantree <- read.tree("plantree.tre")

#Check the tree
#is it binary?

is.binary(plantree)

# is it rooted?

is.rooted(plantree)

#Node labels must be unique

plantree <- makeNodeLabel(plantree)

#Is ultrametric

is.ultrametric(plantree)
plantree <- chronos(plantree,model = "correlated") 
is.ultrametric(plantree)
class(plantree) <- "phylo"

# Correct the names of the species

plantree$tip.label <- gsub("_", " ", plantree$tip.label)
plantree$tip.label <- gsub("[0-9]+", "", plantree$tip.label)
plantree$tip.label <- gsub(" ott", "", plantree$tip.label)
plantree$tip.label <- sub("^(\\S*\\s+\\S+).*", "\\1", plantree$tip.label)

# Animals

antree <- read.tree("antree.tre")

#Check the tree
#is it binary?

is.binary(antree)

# is it rooted?

is.rooted(antree)

#Node labels must be unique

antree <- makeNodeLabel(antree)

#Is ultrametric

is.ultrametric(antree)
#plantree <- chronos(plantree,model = "correlated") 
#is.ultrametric(plantree)
#class(tree) <- "phylo"


# Correct the names of the species

antree$tip.label <- gsub("_", " ", antree$tip.label)
antree$tip.label <- gsub("[0-9]+", "", antree$tip.label)
antree$tip.label <- gsub(" ott", "", antree$tip.label)
antree$tip.label <- sub("^(\\S*\\s+\\S+).*", "\\1", antree$tip.label)


# Match the tree

sp <- unique(plandata$species)
names(sp) <- unique(plandata$species)
(chk<- name.check(plantree, sp))

# Prune the tree

smalltree <- drop.tip(plantree,chk$tree_not_data)

# Sort dataframe

smallplandata <- plandata[plandata$species%in%smalltree$tip.label,]


#Change to a data frame

smallplandata <- as.data.frame(smallplandata)
smallplandata$animal <- smallplandata$species

# Match the tree

sp <- unique(andata$species)
names(sp) <- unique(andata$species)
(chk<- name.check(antree, sp))

# Prune the tree

smallantree <- drop.tip(antree,chk$tree_not_data)

# Sort dataframe

smallandata <- andata[andata$species%in%smallantree$tip.label,]


#Change to a data frame

smallandata <- as.data.frame(smallandata)
smallandata$animal <- smallandata$species

# Create a factor distinguishing size/stage-based matrices #####################

smallandata <- smallandata %>% 
  mutate(typeM=ifelse(MatrixCriteriaSize=="Yes", "Size-based",
                      ifelse(MatrixCriteriaOntogeny=="Yes"|
                               MatrixCriteriaAge=="Yes", "Stage-based", NA)))

smallplandata <- smallplandata %>% 
  mutate(typeM=ifelse(MatrixCriteriaSize=="Yes", "Size-based",
                      ifelse(MatrixCriteriaOntogeny=="Yes"|
                               MatrixCriteriaAge=="Yes", "Stage-based", NA)))

# Save the data ----

setwd(DataPath)
save(smallplandata, smallandata, smalltree, smallantree, 
     file = "ResData.RData")
