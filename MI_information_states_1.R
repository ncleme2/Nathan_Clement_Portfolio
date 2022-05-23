#Computes the Mutual Information between each combination of residues in a 
#protein.
#This script is designed to run in parallel in the Palmetto Cluster.

#Allows simple tasks to be performed for each of some index value.
library(foreach)

#Enables a parallel backend for foreach. This makes computations over multiple cores possible.
library(doParallel)

#This package is used for its entropy function. Entropy is necessary to compute for calculation mutual information.
library(infotheo)

#This package is used to force diagonal matrices to be symmetric.
library(Matrix)

library(cluster)

library(bio3d)

################################################################################
#paste this into every script
#turn it into an argument for the script?
#include production number if necessary
#for example: 6m0j_1, 6m0j_2...
################################################################################

Protein_of_Interest <- "1v4s"
setwd(paste("/scratch1/ncleme2/Info_States/", Protein_of_Interest, "/charmm-gui/namd/MI_Analysis_Data/10k", sep = ""))


#this will be used to add residue names to the data later
#use this as an index with a foreach loop to create a name vector. 
pdb <- bio3d::read.pdb(file = paste("/scratch1/ncleme2/Info_States/", Protein_of_Interest, "/charmm-gui/namd/step3_input.pdb", sep = ""))
ca_inds <- atom.select(pdb, elety = "CA")


production_clusters_input <- as.matrix(read.csv(file = paste(Protein_of_Interest, "_production_clusters.csv", sep = "")))
#replace 1:ncol() with the created name vector. run whole script and check how the names work with the calculations. might have to fix later?

#colnames(production_clusters_input) <- paste(pdb$atom$resid[ca_inds$atom], ca_inds$atom, sep = "")


################################################################################
#Initialize parralelization backend
#The cores availible for computation are evaluated and registered to the backend.
#One core is left out to ensure that there is a core to run RStudio.
################################################################################

co <- detectCores()-1
cl <- makeCluster(co)
registerDoParallel(cl)


################################################################################
#Alpha carbon individual entropy calculation and list formation
#Mutual information calculations use individual and joint entropies.
#The individual entropies are calculated here to that they can be utilized later.
################################################################################

ind.entropy<-foreach (k = 1:(ncol(production_clusters_input)), .packages='infotheo') %dopar% entropy(production_clusters_input[,k], method = "emp")

################################################################################
#Positional Mutual Information calculation
#First an empty matrix is generated so that the foreach loop can write to it during calculations.
#We are writing to a matrix during this calculation so that a diagonal matrix can be made from the calculation.
#Mutual inforamtion heatmaps are symmetric and this cuts out half of the calculations.
################################################################################

mut<-matrix(data=NA, nrow=ncol(production_clusters_input), ncol=ncol(production_clusters_input))

#The mutual information function calculations the mutual information for residue i and j.
#The goal of the mutual information calculations is to calculate the relationship of each combination of residues.

mutual_information <- function(i,j)
{
  ind.entropy[[i]] + ind.entropy[[j]] - infotheo::entropy(cbind(production_clusters_input[,((1)+(i-1))],production_clusters_input[,((1)+(j-1))]), method = "emp")
}

#This foreach loop indexes over each combination of residues i and j and writes the result to the corresponding location in the diagonal matrix.

foreach(i=1:(ncol(production_clusters_input)), .packages='infotheo') %:% foreach(j=1:i) %do% (mut[i,j]<-mutual_information(i,j))

#The matrix is then forced to be symmetric for analysis.

mut<-forceSymmetric(mut, uplo = "L")
#colnames(mut) <- paste(pdb$atom$resid[ca_inds$atom], ca_inds$atom, sep = "")
#rownames(mut) <- paste(pdb$atom$resid[ca_inds$atom], ca_inds$atom, sep = "")

################################################################################
#Normalize the Mutual Information
#Mutual Information is calculated in bits and it is easier to analyze on a scale from 0 to 1.
#The normalize function performs this normalization based upon the mutual information value and the individual entropies of the residues.
################################################################################

normalize <- function(calculated_mutual_information, Entropy_X, Entropy_Y) 
{
  2*calculated_mutual_information/(Entropy_X+Entropy_Y)
}

#These next lines of code function exactly like those in the Mutual Information section, but calculate the normalized information instead.

nmut<-matrix(data=NA, nrow=ncol(mut), ncol=ncol(mut))

foreach(ii=1:ncol(mut)) %:% foreach(jj=1:ii) %do% (nmut[ii,jj]<-normalize(mut[ii,jj], ind.entropy[[ii]], ind.entropy[[jj]]))
   
nmut<-forceSymmetric(nmut, uplo = "L")
#colnames(nmut) <- paste(pdb$atom$resid[ca_inds$atom], ca_inds$atom, sep = "")
#rownames(nmut) <- paste(pdb$atom$resid[ca_inds$atom], ca_inds$atom, sep = "")

mut<-as.matrix(mut)
nmut<-as.matrix(nmut)

################################################################################
#Stop the cluster
################################################################################

stopCluster(cl)

################################################################################
#Output .csv of Mutual Information
#Shiny apps do not work well on Palmetto.
# A .csv is required so that the data can be sent to a local computer and graphed.
#nmut must be evaluated as a matrix because the Matrix package turns it into a dsyMatrix which is not compatible to export with many commands.
################################################################################

write.table(mut, file = paste(Protein_of_Interest, "_mut.csv", sep = ""), row.names = F, col.names = F)
write.table(nmut, file = paste(Protein_of_Interest, "_nmut.csv", sep = ""), row.names = F, col.names = F)