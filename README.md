# Life history mediates the trade-offs among different components of demographic resilience

Pol Capdevila<sup>1,2</sup>*, Iain Stott<sup>3</sup>, James Cant<sup>4</sup>, Maria Beger<sup>4,5</sup>, Gwilym Rowlands<sup>1</sup>, Molly Grace<sup>1</sup>, Roberto Salguero-Gómez<sup>1,5,6</sup>
 
<sup>1</sup>Zoology Department, Oxford University, Zoology Research and Administration Building, 11a Mansfield Rd, Oxford OX1 3SZ, UK

 <sup>2</sup>School of Biological Sciences, University of Bristol, 24 Tyndall Ave, BS8 1TQ, Bristol, UK 
 
 <sup>3</sup>School of Life and Environmental Sciences, University of Lincoln, Brayford Pool, Lincoln LN6 7TS, UK
 
 <sup>4</sup>School of Biology, Faculty of Biological Sciences, University of Leeds, UK, LS2 9JT
 
 <sup>5</sup>Centre for Biodiversity and Conservation Science, School of Biological Sciences, University of Queensland, Brisbane, 4072, Australia
 
 <sup>6</sup>Max Planck Institute for Demographic Research, Konrad Zuse Straße 1, Rostock 18057, Germany

#### Contact: pcapdevila.pc[at]gmail.com

---

## Abstract

_Accelerating rates of biodiversity loss underscore the need to understand how species achieve resilience – the ability to resist and recover from a/biotic disturbances. Yet, the factors determining the resilience of species remain poorly understood, due to disagreements on its definition and the lack of large-scale analyses. Here, we investigate how the life history of 910 natural populations of animals and plants predict their intrinsic ability to be resilient. We show that demographic resilience can be achieved through different combinations of compensation, resistance, and recovery after a disturbance. We demonstrate that these resilience components are highly correlated with life history traits related to the species’ pace of life and reproductive strategy. Species with longer generation times require longer recovery times post-disturbance, while those with greater reproductive capacity have greater resistance and compensation. Our findings highlight the key role of life history traits to understand species resilience, improving our ability to predict how natural populations cope with disturbance regimes._

---

## Data

- __`ResData`__: final data set containing the components of resilience and life history traits data. 
- __`transDataPost.RData`__: initial calculations of the components of resilience and life history traits without transformation and normalisation. 
- __`SimData.RData`__: simulated matrices to create Figure S1.  
- __`RaunkiaerGrowthForms.csv`__: list of species with different Raunkiaer classification.
- __`finalRaunk.csv`__: final Raunkiaer classification for Figure S3.
- __`IUCNPopDataPost.RData`__: data set containing the components of resilience and life history traits data and their conservation status used to produce Figure S4.  

---

# Code

To run the statistical analyses we used different R scripts: 

- __`AnimalPhylogeny.R`__: code to produce the phylogeny of animals.
- __`DataPreparation.R`__: code to prepare the data for the analyses.
- __`DemographicCalculations.R`__: code to calculate the demographic resilience components and life history traits. 
- __`Figuresbrms.R`__: code to produce the figures. 
- __`Functions.R`__: generic functions used in the code.
- __`IUCNPopData.R`__: code to get the conservation status of the species.
- __`LifeHistoryModels.R`__: code to model the correlations among the demographic resilience components and life history traits. 
- __`Matrix Simulation.R`__: code to simulate the matrices. 
- __`PhylogeneticSignal.R`__: code to calculate the phylogenetic signal.
- __`PhylogenyPlot.R`__: code to produce Figure 2.
- __`PlantPhylogeny.R`__: code to create the plant phylogeny. 
- __`RaunkiaerCleaining.R`__: code to obtain the Raunkiaer classification. 
- __`SimulationFunctions.R`__: functions used in the matrix simulation to produce Figure S1.
- __`SimulationModelling.R`__: code to replicate the same models than in the main manuscript but using the simulated matrices.
- __`SimulationResults.R`__: code to produce Figure S1. 
- __`SupplementaryAnalyses.R`__: code to produce the supplementary analyses Figures S2-S3. 
- __`SupplementaryFigures.R`__: code to produce Figures S2-S4. 
- __`SupplementaryIUCNModels.R`__: code to produce the supplementary analyses showed in Figure S4. 

 
---

# Software

_R version 4.0.2 or greater_

To download `R`, go to https://www.r-project.org and for `RStudio`, go to https://www.rstudio.com/products/rstudio/download/ .


