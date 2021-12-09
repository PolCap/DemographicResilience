# ---------------------------------------------------------------------------- #
# - FILE NAME:   Transients.R         
# - DATE:        12/08/2019
# - DESCRIPTION:  Code to estimate the transient dynamics of populations in 
#                 compadre and comadre. 
# - AUTHORS:     Iain Stott & Pol Capdevila Lanzaco
# ---------------------------------------------------------------------------- #

rm(list=ls(all=TRUE))

library(tidyverse)
library(fields)
library(scales)
library(MASS)
library(popbio)
library(popdemo)
library(Matrix)
library(Rage)
library(Rcompadre)
library(readr)
library(data.table)

# Data, results and code paths

path = gsub("/Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

DataPath = paste0(path,"/Data")
ResultPath = paste0(path, "/Results") 
CodePath = paste0(path, "/Code") 

# Load functions 

source(paste0(CodePath, "/Functions.R"))

### COMPADRE ###################################################################

# Estimate life history traits --------------------------------------------

# Load compadre

compadreUse <- cdb_fetch("compadre")

# Subset the database according to the criteria expressed in the methods

compadreUse <- compadreUse %>% 
  cdb_flag() %>% 
  filter(check_NA_A == FALSE, 
         check_NA_U == FALSE, 
         check_NA_F == FALSE, 
         check_component_sum ==TRUE, 
         check_ergodic == TRUE, 
         check_primitive == TRUE, 
         MatrixDimension > 2,
         ProjectionInterval == 1,
         MatrixTreatment %in% "Unmanipulated",
         MatrixCaptivity %in% "W", 
         MatrixFec=="Yes", 
         SurvivalIssue<=1, 
         OrganismType!="Epiphyte", 
         OrganismType!="Palm", 
         OrganismType!="Fern")

# Create a unique identifier for different populations

compadreUse$popID <- cdb_id(compadreUse, columns = c("SpeciesAuthor", "MatrixPopulation"))

# Subset mean matrices

compadreUse1 <- compadreUse %>% 
  filter(MatrixComposite=="Mean")

# Look for the individual matrices where there is no mean matrices

compadreUse2 <- compadreUse %>% 
  filter(MatrixComposite=="Individual",
         popID%in%setdiff(compadreUse$popID, 
                          compadreUse1$popID))  
# We now merge them

compadreUse <- cdb_rbind(compadreUse1, compadreUse2)

# Now we download the data from Jelbert et al. Nat Comms to check the pre-reproductive 
# matrices. 

jelbert <- read.csv(paste0(DataPath, "/41467_2019_13556_MOESM6_ESM.csv"))

# We create unique identifiers for the census type distinguishing by DOI

jelbert <- jelbert %>% distinct(publication_DOI_ISBN, .keep_all=T)

# Now we join both data frames to check which are post/pre reproductives 

compadreUse <- compadreUse %>% left_join(jelbert[c("publication_DOI_ISBN", 
                                                   "fixed_census_timing")], 
                                         by=c("DOI_ISBN"="publication_DOI_ISBN")) 

# Now we introduce manually the values for the remaining studies 

compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.1111/j.1600-0706.2011.19946.x"] <- "pre"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.1111/j.1365-2745.2009.01611.x"] <- "pre"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.1023/A:1015506019670"] <- "pre"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.1023/A:1010017510471"] <- "pre"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.1007/978-3-662-09389-4_12"] <- "pre"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.1086/297246"] <- "pre"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.2307/1941656"] <- "pre"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.1111/j.1756-1051.2004.tb01647.x"] <- "pre"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.1111/j.0030-1299.2007.15705.x"] <- "pre"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.2307/2261568"] <- "post"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.1111/1365-2664.12057"] <- "post"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.1111/j.1600-0587.2012.07425.x"] <- "pre"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.1890/07-1908.1"] <- "pre"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.1046/j.1365-2745.1998.00280.x"] <- "pre"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.2307/2446413"] <- "pre"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.1016/j.foreco.2010.11.007"] <- "post"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.2307/177080"] <- "post"
compadreUse$fixed_census_timing[compadreUse$DOI_ISBN=="10.1016/j.jnc.2011.05.005"] <- "pre"

# Correct the post/pre reproductive census (I could not do this with dplyr)

compadreUse$A <- NA

for(i in 1:length(compadreUse$MatrixID)){
  compadreUse$A[i] <- convert2pre(matA(compadreUse[i]), 
                       matF(compadreUse[i]), 
                       matC(compadreUse[i]), 
                       matU(compadreUse[i]),
                       compadreUse$fixed_census_timing[i])
  
  }


# Subset when we have matrices

compadreUse <- compadreUse %>% 
  subset(!is.na(A)) 

# I feel forced to use a for loop because tidyverse is not working for a reason
# beyond my comprehension

for(i in 1:length(compadreUse$MatrixID)){
  if(isErgodic(compadreUse$A[[i]])&
     isPrimitive(compadreUse$A[[i]])&
     isIrreducible(compadreUse$A[[i]])){
    stage <- tryCatch(mpm_first_active(compadreUse$mat[[i]]), 
                      error=function(e) NA)
    if(dim(matU(compadreUse$mat[[i]]))[1]==dim(compadreUse$A[[i]])[1]){
      matR <- tryCatch(compadreUse$A[[i]]-matU(compadreUse$mat[[i]]),
                       error=function(e) NA)
      matU <- compadreUse$A[[i]]-matR
      }else{
      pos <- which(rowSums(matF(compadreUse[i])[[1]])!=0)
      matR <- compadreUse$A[[i]]
      matR[-pos,]=0 
      matU <- compadreUse$A[[i]]-matR
    }
    
    compadreUse$GenT[i]=  tryCatch(unlist(lapply(compadreUse$A[i], 
                                                generation.time)), 
                                  error=function(e) NA)
    compadreUse$Fec[i] = tryCatch(vitalRates(matU,
                                   matR)$fec, 
                                  error=function(e) NA)
    compadreUse$rupr[i] = tryCatch(mapply(FUN = reac,
                                         A = compadreUse$A[i],
                                         bound="upper"),
                                   error=function(e) NA)
    
    compadreUse$rlwr[i] = tryCatch(mapply(FUN = reac,
                                         A = compadreUse$A[i],
                                         bound="lower"),
                                  error=function(e) NA)
    compadreUse$xt[i] = tryCatch(unlist(lapply(compadreUse$A[i], 
                                              return.time)), 
                                error=function(e) NA)
    compadreUse$Dimension[i]=dim(compadreUse$A[[i]])[1]
    
  }else{
    compadreUse$GenT[i]=NA
    compadreUse$rupr[i] = NA
    compadreUse$rlwr[i] = NA
    compadreUse$xt[i] =NA
    compadreUse$Dimension[i]=NA
    }
  
  print(i)
}

# Now we will estimate the life history traits and resilience components

# compadreUse <- compadreUse %>% 
#   mutate(GenT=  tryCatch(unlist(lapply(A, generation.time)), 
#                          error=function(e) NA),
#          rupr = tryCatch(mapply(FUN = reac,
#                 A = A,
#                 bound="upper"),error=function(e) NA),
#          rlwr = tryCatch(mapply(FUN = reac,
#                        A = A,
#                        bound="lower"),error=function(e) NA),
#          dr = tryCatch(unlist(lapply(A, dr)),error=function(e) NA),
#          xt = tryCatch(unlist(lapply(A, return.time)), error=function(e) NA),
#          Dimension=dim(A)[1])

### COMADRE ####################################################################

# Load database

comadreUse <- cdb_fetch("comadre")

# Subset the database according to the criteria expressed in the methods

comadreUse <- comadreUse %>% 
  cdb_flag() %>% 
  filter(check_NA_A == FALSE, 
         check_NA_U == FALSE, 
         check_NA_F == FALSE, 
         check_component_sum ==TRUE, 
         check_ergodic == TRUE, 
         check_primitive == TRUE, 
         MatrixDimension > 2,
         ProjectionInterval == 1,
         MatrixTreatment %in% "Unmanipulated",
         MatrixCaptivity %in% "W", 
         MatrixFec=="Yes", 
         SurvivalIssue<=1,
         Phylum=="Chordata")

# Create a unique identifier for different populations

comadreUse$popID <- cdb_id(comadreUse, columns = c("SpeciesAuthor", "MatrixPopulation"))

# Subset mean matrices

comadreUse1 <- comadreUse %>% 
  filter(MatrixComposite=="Mean")

# Look for the individual matrices where there is no mean matrices

comadreUse2 <- comadreUse %>% 
  filter(MatrixComposite=="Individual",
         popID%in%setdiff(comadreUse$popID, 
                          comadreUse1$popID))  
# We now merge them

comadreUse <- cdb_rbind(comadreUse1, comadreUse2)

# Enter manually the census type 

comadreUse$fixed_census_timing <- NA

comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/j.ecolmodel.2010.06.026"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1890/07-0892.1"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/j.1600-0633.2005.00084.x"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/j.ecolmodel.2012.09.022"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1139/cjfas-2012-0520"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1046/j.1523-1739.2003.01535.x"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1371/jourNAl.pone.0085464"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1080/20018091094835"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1007/s00227-012-1933-6"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1890/07-0018.1"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1163/157075610X523260"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1890/0012-9615(2001)071[0377:RAODRT]2.0.CO;2"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/j.1365-2435.2009.01563.x"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/1365-2664.12476"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1007/s10336-011-0758-2"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2307/3803072"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/j.0908-8857.2008.04189.x"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/j.1365-2664.2012.02163.x"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1002/ece3.1290"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/j.1523-1739.2005.00276.x"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1890/03-5340"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/ibi.12125"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1890/06-1090.1"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1007/s10144-012-0306-9"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/acv.12054"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/j.biocon.2009.12.010"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2981/13-021"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/j.baae.2010.11.004"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.5253/078.100.0208"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1002/jwmg.628"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/j.1469-1795.2009.00311.x"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2307/1369293"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1007/s10336-011-0745-7"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1642/0004-8038(2004)121[1056:TVITVR]2.0.CO;2"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1093/icesjms/fsu056"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/evo.12952"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1086/657443"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1371/jourNAl.pone.0034379"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/S0304-3800(01)00493-8"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1007/BF01237655"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/j.1600-0706.2011.19436.x"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2193/2005-608"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1017/S0950268812000167"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/j.biocon.2014.05.026"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2307/2402601"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1002/ajp.22177"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1007/s10764-010-9461-z"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1002/ajp.22323"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1086/597225"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.3354/meps10547"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/j.1365-2656.2007.01274.x"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1002/jwmg.835"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2307/2641054"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2307/3803144"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/j.1365-2664.2010.01846.x"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/j.biocon.2009.06.020"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/S0006-3207(02)00421-4"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/j.biocon.2007.10.006"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2193/2006-349"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1890/0012-9658(2001)082[1921:EODROU]2.0.CO;2"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2980/1195-6860(2007)14[362:LHOFRS]2.0.CO;2"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1007/s00442-010-1761-7"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/S0006-3207(02)00341-5"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2192/1537-6176-20.2.77"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/j.ecolmodel.2014.08.021"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1890/09-1641"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/j.1600-0706.2012.20706.x"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2307/1941948"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/j.ecolmodel.2014.11.025"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/j.1523-1739.2005.00487.x"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/S0304-3800(01)00433-1"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="978-0-19-803726-2"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1093/icb/34.3.397"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2307/1941498"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/j.biocon.2008.04.001"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="978-0-949324-89-4"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1071/wr9850541"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1046/j.1440-1703.2002.00463.x"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1007/s10144-011-0292-3"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2307/1941572"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1111/1365-2664.12194"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1007/s10144-014-0450-5"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2744/CCB-0778.1"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1655/HERPETOLOGICA-D-12-00038R2"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1643/0045-8511(2007)7[324:AOTPDO]2.0.CO;2"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1163/156853808784124992"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="https://www.jstor.org/stable/3830713"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1577/M08-034.1"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1890/08-0305.1"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/j.ecolmodel.2015.12.002"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/j.ecoenv.2012.01.019"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1016/j.jas.2012.04.037"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.1007/s00227-015-2695-8"] <- "post"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2307/1381947"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.2307/1935300"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$DOI_ISBN=="10.3398/1527-0904(2007)67[492:GADOOP]2.0.CO;2"] <- "post"
comadreUse$fixed_census_timing[comadreUse$MatrixID==240502] <- "post"
comadreUse$fixed_census_timing[comadreUse$MatrixID==248046] <- "post"
comadreUse$fixed_census_timing[comadreUse$SpeciesAuthor=="Notropis_anogenus_2"] <- "post"
comadreUse$fixed_census_timing[comadreUse$SpeciesAuthor=="Percina_copelandi_2"] <- "post"
comadreUse$fixed_census_timing[comadreUse$SpeciesAuthor=="Notropis_photogenis_2"] <- "post"
comadreUse$fixed_census_timing[comadreUse$SpeciesAuthor=="Alces_alces_6"] <- "post"
comadreUse$fixed_census_timing[comadreUse$SpeciesAuthor=="Kinosternon_integrum"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$SpeciesAuthor=="Eumetopias_jubatus_3"] <- "pre"
comadreUse$fixed_census_timing[comadreUse$Authors=="Peoples"] <- "pre"

# Correct the post/pre reproductive census (I could not do this with dplyr)

comadreUse$A <- NA

for(i in 1:length(comadreUse$MatrixID)){
  comadreUse$A[i] <- convert2pre(matA(comadreUse[i]), 
                                  matF(comadreUse[i]), 
                                  matC(comadreUse[i]), 
                                  matU(comadreUse[i]), 
                                 comadreUse$fixed_census_timing[i])
  
}

# Subset when we have matrices

comadreUse <- comadreUse %>% 
  subset(!is.na(A))

# I feel forced to use a for loop because tidyverse is not working for a reason
# beyond my comprehension


for(i in 1:length(comadreUse$MatrixID)){
  if(isErgodic(comadreUse$A[[i]])&
     isPrimitive(comadreUse$A[[i]])&
     isIrreducible(comadreUse$A[[i]])){
    stage <- tryCatch(mpm_first_active(comadreUse$mat[[i]]), 
                      error=function(e) NA)
    if(dim(matU(comadreUse$mat[[i]]))[1]==dim(comadreUse$A[[i]])[1]){
      matR <- tryCatch(comadreUse$A[[i]]-matU(comadreUse$mat[[i]]),
                       error=function(e) NA)
      matU <- comadreUse$A[[i]]-matR
    }else{
      pos <- which(rowSums(matF(comadreUse[i])[[1]])!=0)
      matR <- comadreUse$A[[i]]
      matR[-pos,]=0 
      matU <- comadreUse$A[[i]]-matR
      }
      
    comadreUse$GenT[i]=  tryCatch(unlist(lapply(comadreUse$A[i], 
                                                 generation.time)), 
                                   error=function(e) NA)
    comadreUse$Fec[i] = vitalRates(matU,
                                    matR)$fec
    comadreUse$rupr[i] = tryCatch(mapply(FUN = reac,
                                          A = comadreUse$A[i],
                                          bound="upper"),error=function(e) NA)
    
    comadreUse$rlwr[i] = tryCatch(mapply(FUN = reac,
                                          A = comadreUse$A[i],
                                          bound="lower"),
                                   error=function(e) NA)
    comadreUse$xt[i] = tryCatch(unlist(lapply(comadreUse$A[i], 
                                               return.time)), 
                                 error=function(e) NA)
    comadreUse$Dimension[i]=dim(comadreUse$A[[i]])[1]
  
    }else{
    comadreUse$GenT[i]=NA
    comadreUse$rupr[i] = NA
    comadreUse$rlwr[i] = NA
    comadreUse$xt[i] =NA
    comadreUse$Dimension[i]=NA
  }
  
  print(i)
}

# Add body size information #### 

# Plants
# Load data from TRY

setwd(DataPath)

tryData <- fread("8609.txt",header = T, sep = "\t", dec = ".", 
                 quote = "", data.table = T, showProgress = FALSE)

# Estimate the mean 

tryData <- tryData %>%
  filter(AccSpeciesName%in%compadreUse$SpeciesAccepted,
         TraitName=="Plant height vegetative") %>% 
  group_by(AccSpeciesName) %>%
  summarise(mean_height=mean(StdValue, na.rm=TRUE),
            max_height=max(StdValue, na.rm = T)) %>% 
  left_join(tryData, ., by = c('AccSpeciesName')) %>%
  group_by(AccSpeciesName) %>%
  filter(StdValue==max_height) %>%
  dplyr::select(AccSpeciesName, mean_height, max_height, Reference)%>%
  distinct(AccSpeciesName, .keep_all = T)

# Join with comadre

compadreUse <- left_join(compadreUse, tryData, 
                   by =c( "SpeciesAccepted" ="AccSpeciesName"))

# Animals

bm <- read.csv("Amniote_Database_Aug_2015.csv")

# Create a species name

bm <- bm %>% 
  mutate(SpeciesAccepted = paste(genus, species)) %>% 
  dplyr::select(SpeciesAccepted, adult_body_mass_g)

# Join with comadre

comadreUse <- left_join(comadreUse, bm)

#Save both data

setwd(DataPath)
save(comadreUse, compadreUse, file="transDataPost.RData") #transData2 
