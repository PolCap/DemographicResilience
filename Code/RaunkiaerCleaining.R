

setwd("/Users/rob/Desktop/")

d1 <- read.csv("PlantData.csv")
d2 <- read.csv("RaunkiaerGrowthForms.csv")

names(d1)
head(d1)
dim(d1)

d1=d1[which(!duplicated(d1$SpeciesAccepted)),]
dim(d1)

names(d2)
head(d2)
d2$Species_=gsub(" ","_",d2$Species)

length(which(d1$SpeciesAccepted%in%d2$Accepted.name))
length(which(d1$SpeciesAuthor%in%d2$Species_))
length(which(d1$SpeciesAccepted%in%d2$Species))
d1$SpeciesAccepted[which(!d1$SpeciesAccepted%in%d2$Species)]
#"Alnus incana subsp. rugosa" =--> "Alnus indicata subspp. rugosa"
d2$Species[which(d2$Species=="Alnus indicata subspp. rugosa")]="Alnus incana subsp. rugosa"
#"Callitris columellaris"   <- not in COMPADRE Raunkiaer file, but yes Callitris intratropica 
#"Sonchus pustulatus" <- not in compadre data
#"Thymus vulgaris" --> "Thymus loscosii" OR "Thymus webbianus"            
#"Trollius europaeus" --> "Trollius laxus"
#"Vitaliana primuliflora" --> Nope

d2$Growth.form.author[which(!d1$SpeciesAccepted%in%d2$Species)]

d2.1=d2[,c("Species","Growth.form.Raunkiaer","Growth.form.author")]

d3=merge(d1,d2.1,by.x="SpeciesAccepted",by.y="Species", all.x=TRUE)
dim(d3)
head(d3)

d3[which(is.na(d3$Growth.form.Raunkiaer)),c("SpeciesAccepted","max_height")]

d3$Growth.form.Raunkiaer[which(d3$SpeciesAccepted=="Callitris columellaris")]="Megaphanerophyte"
d3$Growth.form.Raunkiaer[which(d3$SpeciesAccepted=="Sonchus pustulatus")]="Hemicryptophyte"
d3$Growth.form.Raunkiaer[which(d3$SpeciesAccepted=="Thymus vulgaris")]="Chamaephyte"
d3$Growth.form.Raunkiaer[which(d3$SpeciesAccepted=="Trollius europaeus")]="Hemicryptophyte"
d3$Growth.form.Raunkiaer[which(d3$SpeciesAccepted=="Vitaliana primuliflora")]="Hemicryptophyte"

d3[which(d3$Growth.form.Raunkiaer==""),c("SpeciesAccepted","max_height")]

d3$Growth.form.Raunkiaer[which(d3$SpeciesAccepted=="Chaerophyllum aureum")]="Hemicryptophyte"
d3$Growth.form.Raunkiaer[which(d3$SpeciesAccepted=="Colchicum autumnale")]="Geophyte"
d3$Growth.form.Raunkiaer[which(d3$SpeciesAccepted=="Cornus florida")]="Mesophanerophyte"
d3$Growth.form.Raunkiaer[which(d3$SpeciesAccepted=="Dracocephalum austriacum")]="Chamaephyte"
d3$Growth.form.Raunkiaer[which(d3$SpeciesAccepted=="Linum flavum")]="Chamaephyte"
d3$Growth.form.Raunkiaer[which(d3$SpeciesAccepted=="Linum tenuifolium")]="Chamaephyte"
d3$Growth.form.Raunkiaer[which(d3$SpeciesAccepted=="Saussurea medusa")]="Hemicryptophyte"
d3$Growth.form.Raunkiaer[which(d3$SpeciesAccepted=="Trollius laxus")]="Hemicryptophyte"

table(d3$Growth.form.Raunkiaer)

#Tidying things up
d3$Growth.form.Raunkiaer[which(d3$Growth.form.Raunkiaer=="Macrophanerophyte")]="Megaphanerophyte"

table(d3$Growth.form.Raunkiaer)


write.csv(d3,"finalRaunk.csv")




  
  