#MI Prerun Script
#Load trajectory file, remove indices of interest, output trajectory file of data to be analyzed.
#Clean the data (remove unequilibrated frames using RMSD
#Perform several calculations to determine the optimal number of clusters for each column of the data set.

library(foreach)
library(doParallel)
library(fpc)
library(cluster)

################################################################################
#paste this into every script
################################################################################

Protein_of_Interest <- "1v4s"

#Change Test_Data to Information_States for real production!
setwd(paste("/scratch1/ncleme2/Info_States/", Protein_of_Interest, "/charmm-gui/namd/MI_Analysis_Data", sep = ""))



#alpha_carbons_input <- read.csv(file = paste(Protein_of_Interest, "_alpha_trj.csv", sep = ""))
#load master energy file
master <- read.table(file = paste(Protein_of_Interest, "_net_master.nrg", sep = ""))



################################################################################
#Clustering and Validation
################################################################################

co <- detectCores()-1
cl <- makeCluster(co)
registerDoParallel(cl)

################################################################################
#cluster each residue with an increasing number of test clusters.
#the removal of the variable names may need to be moved into the for loop for larger systems
################################################################################
residue_index <- 1:ncol(master)
k_test_clusters <- 5

for (r in residue_index) {
  
  cluster_set <- paste("cluster_set_", r, sep = "")
  assign(cluster_set, foreach(i = 2:k_test_clusters, .combine = cbind, .packages = 'cluster') %dopar% pam(master[,r], i, pamonce=6, metric="euclidean",   cluster.only = T))
  
}

cluster.lookup <- foreach(r = residue_index) %do% get(paste("cluster_set_", r, sep = ""))
rm(list = ls(pattern = "cluster_set_", envir = .GlobalEnv), envir = .GlobalEnv)

################################################################################
#search the results of each cluster value for each residue to compute the optimal number of energetic microstates
#the removal of the variable names may need to be moved into the for loop for larger systems
################################################################################

for (r in residue_index) {
  
  CH_Index <- paste("CH_Index_", r, sep = "")
  assign(CH_Index, foreach(k = 1:(k_test_clusters-1), .combine = c, .packages = 'fpc') %dopar% round(calinhara(master[,r], cluster.lookup[[r]][,k]),digits=2))
  
}  

CH.lookup <- foreach(r = residue_index) %do% get(paste("CH_Index_", r, sep = ""))
rm(list = ls(pattern = "CH_Index_", envir = .GlobalEnv), envir = .GlobalEnv)
production_index <- foreach (r = residue_index, .combine = c) %do% which(CH.lookup[[r]] == max(CH.lookup[[r]]))+1  

###################################################################
#the production index is a vector containing the optimal number of microstates for each residue
#if a vector has the maximum microstates test number rerun this script with a higher number of test clusters
###################################################################

write.csv(production_index, file = paste(Protein_of_Interest, "_production_index.csv", sep = ""), col.names = F, row.names = F)

###################################################################
#production_clusters_output so the data does not need to be clustered again in next script
###################################################################

production_clusters <- foreach(r = residue_index) %do% cluster.lookup[[r]][,(production_index[r]-1)]
write.csv(production_clusters, file = paste(Protein_of_Interest, "_production_clusters.csv", sep = ""), col.names = F, row.names = F)

###################################################################
#compute the silhouette for each residue's cluster value to interpret cluster quality
###################################################################

silhouette_run<-foreach(r = residue_index, .combine = cbind, .packages = 'cluster') %do% pam(master[, r], production_index[r], pamonce=6, metric="euclidean")
avg_sil<-foreach(i = residue_index, .combine = c, .packages = 'cluster') %do% silhouette_run[,i]$silinfo$avg.width

stopCluster(cl)


