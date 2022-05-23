################################################################################
#Protein selection 
#Production names
#Allows for automatic loading of nrg from its appropriate directory
################################################################################

Protein_of_Interest <- "1v4s"
#Vector of all production names
Production_Name <- c("Production", "Restart1", "Restart2")

#Change Test_Data to Information_States for real production!
setwd(paste("/scratch1/ncleme2/Info_States/", Protein_of_Interest, "/charmm-gui/namd/MI_Analysis_Data", sep = ""))

for (job in Production_Name) {
  
  current_job <- paste(Protein_of_Interest, "_", job, ".nrg", sep = "")
  assign(current_job, read.table(file = paste(Protein_of_Interest, "_", job, "_master.nrg", sep = "")))
  
}



job_count <- 1
net_master <- matrix(data = NA, ncol = ncol(get(paste(Protein_of_Interest, "_", Production_Name[1], ".nrg", sep = ""))))
colnames(net_master) <- colnames(get(paste(Protein_of_Interest, "_", Production_Name[1], ".nrg", sep = "")))

for (job in Production_Name) {
  
net_master <- rbind(net_master, get(paste(Protein_of_Interest, "_", Production_Name[job_count], ".nrg", sep = "")))
job_count <- job_count+1  

}

net_master <- net_master[2:nrow(net_master), ]

write.table(net_master, file = paste(Protein_of_Interest, "_net_master.nrg", sep = ""))
