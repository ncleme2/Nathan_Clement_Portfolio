library(bio3d)

################################################################################
#Protein selection 
#Trajectory selection
#Allows for automatic loading of pdb and nrg from their appropriate directories
################################################################################

Protein_of_Interest <- "1v4s"
Production_Name <- "step5_production"

setwd(paste("/scratch1/ncleme2/Info_States/", Protein_of_Interest, "/charmm-gui/namd", sep = ""))

pdb <- bio3d::read.pdb(file = paste("/scratch1/ncleme2/Info_States/", Protein_of_Interest, "/charmm-gui/namd/step3_input.pdb", sep = ""))

rd <- rmsd(xyz[1,ca_inds$xyz], xyz[,ca_inds$xyz])

pdf(paste(Protein_of_Interest,"_rmsd.pdf", sep = "")) 

plot(rd, typ="l", ylab="RMSD", xlab="Frame No.")
points(lowess(rd), typ="l", col="red", lty=2, lwd=2)

dev.off()
