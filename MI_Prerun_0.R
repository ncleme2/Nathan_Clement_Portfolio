#MI Prerun Script
#Load trajectory file, remove indices of interest, output trajectory file of data to be analyzed.
#Clean the data (remove unequilibrated frames using RMSD
#Perform several calculations to determine the optimal number of clusters for each column of the data set.

library(bio3d)
#library(parallel)

################################################################################
#Protein selection 
#Production names
#Allows for automatic loading of nrg from its appropriate directory
################################################################################

Protein_of_Interest <- "1v4s"
#Vector of all production names
Production_Name <- c("step5_production", "step6_restart1", "step7_restart2", "step8_restart3", "step9_restart4", "step10_restart5")

#Change Test_Data to Information_States for real production!
setwd(paste("/scratch1/ncleme2/Info_States/", Protein_of_Interest, "/charmm-gui/namd", sep = ""))

################################################################################
#Load trajectory and protein data bank files.
################################################################################

pdb <- bio3d::read.pdb(file = "step3_input.pdb")
pdf(paste(Protein_of_Interest,"_rmsd.pdf", sep = "")) 

################################################################################
#Trajectory frame superposition on alpha carbons atoms
################################################################################

for (job in Production_Name) {
  
  full_trj <- bio3d::read.dcd(trjfile = paste(job, "_x.dcd", sep = ""))
  
  ca_inds <- atom.select(pdb, elety = "CA")
  xyz <- fit.xyz(fixed = pdb$xyz, mobile = full_trj, fixed.inds = ca_inds$xyz, mobile.inds = ca_inds$xyz)
  cleaned_alphas <- cbind(xyz[,ca_inds$xyz])
  
  current_rmsd <- paste(job, "_rmsd", sep = "")
  assign(current_rmsd, rmsd(xyz[1,ca_inds$xyz], xyz[,ca_inds$xyz]))
  
  plot(get(paste(job, "rmsd", sep = "")), typ="l", ylab="RMSD", xlab="Frame No.", main = paste(job, "_rmsd", sep = ""))
  points(lowess(current_rmsd), typ="l", col="red", lty=2, lwd=2)
  
  rm(full_trj)
  rm(xyz)
  
}

################################################################################
#End Loop and Write Output
################################################################################

#cleaned_alphas <- cleaned_alphas[500:nrow(cleaned_alphas),]

dev.off()

#write.csv(cleaned_alphas, paste(Protein_of_Interest, "_alpha_trj.csv", sep = ""))
#write.csv(ca_inds$atom, paste(Protein_of_Interest, "_ca_inds.csv", sep = ""))