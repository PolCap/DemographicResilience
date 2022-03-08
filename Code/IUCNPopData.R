# --------------------------------------------------------------------------------------- #
# - FILE NAME:   IUCNPopData.R         
# - DATE:        22/04/2020
# - DESCRIPTION: Code to create IUCN data for resilience analyses  
# - AUTHORS:     Pol Capdevila Lanzaco
# -------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

library(rredlist)
library(tidyr)
library(dplyr)
library(Rcompadre)
library(tidyverse)
library(magrittr)

#Working directories

path = gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath = paste0(path,"/Data")
ResultPath = paste0(path, "/Results") 
CodePath = paste0(path, "/Code") 

#Load  data----

setwd(DataPath)

#Load transiens

load("transDataPost.RData") 

# Clean the data

comadreUse$species <- gsub("_"," ", comadreUse$SpeciesAccepted)
comadreUse$species <- gsub("\\svar.+", "", comadreUse$species)#remove varieties because of no matches in several points
comadreUse$species <- gsub("[0-9]+", "", comadreUse$species)
comadreUse$species<-gsub("\\(.*","",comadreUse$species)

#Now we remove the subspecies 
comadreUse$species <- sub("^(\\S*\\s+\\S+).*", "\\1", comadreUse$species)

#Transform database into data.frame

andata <- cdb_flatten(comadreUse)

#Now plants

compadreUse$species <- gsub("_"," ", compadreUse$SpeciesAccepted)
compadreUse$species <- gsub("\\svar.+", "", compadreUse$species)#remove varieties because of no matches in several points
compadreUse$species <- gsub("[0-9]+", "", compadreUse$species)
compadreUse$species<-gsub("\\(.*","",compadreUse$species)
compadreUse$species <- sub("^(\\S*\\s+\\S+).*", "\\1", compadreUse$species)

#Transform database into data.frame

plandata <- cdb_flatten(compadreUse)

#Match the IUCN criteria ------------------------------------------------------

IUCN_REDLIST_KEY #Provide your own key"
#get Red List version
rl_version(key=IUCN_REDLIST_KEY)
#February, 2020

#Remove duplicates temporarly 

plan <- plandata %>% 
  distinct(species)

#Loop for IUCN

plan <- plan %>% 
  mutate(IUCN=NA,
         trend= NA,
         criteria=NA,
         threat=NA,
         threats=NA)

for(i in 1:length(plan$species)){
  x <- tryCatch(rl_search(plan$species[i],
                          key = IUCN_REDLIST_KEY,
                          parse = T),
                error=function(e) NULL)
  tryCatch(plan$IUCN[i] <- x$result$category,
           error=function(e) NULL)
  print(plan$species[i])
}

# Join with plandata per population

plandata <- merge(plandata,plan, by="species",all.x = T,all.y = F)

# Animals

#Remove duplicates temporarly 

an <- andata %>% 
  distinct(species)

# Loop

an <- an %>% 
  mutate(IUCN=NA,
         trend= NA,
         criteria=NA,
         threat=NA,
         threats=NA)

for(i in 1:length(an$species)){
  x <- tryCatch(rl_search(an$species[i], 
                          key = IUCN_REDLIST_KEY, 
                          parse = T),
                error=function(e) NULL)
  tryCatch(an$IUCN[i] <- x$result$category,
           error=function(e) NULL)
  print(an$species[i])
}

# Join with andata per population

andata <- merge(andata,an, by="species", all.x=T)

andata$category <- as.character(andata$IUCN)
andata$category[andata$category=="LR/cd"] <- "LC"
andata$category[andata$category=="LR/nt"] <- "LC"
andata$category[andata$category=="LR/lc"] <- "LC"

plandata$category <- as.character(plandata$IUCN)
plandata$category[plandata$category=="LR/cd"] <- "LC"
plandata$category[plandata$category=="LR/nt"] <- "LC"
plandata$category[plandata$category=="LR/lc"] <- "LC"

# save the data

setwd(DataPath)
save(andata, plandata, file = "IUCNPopDataPost.RData")
