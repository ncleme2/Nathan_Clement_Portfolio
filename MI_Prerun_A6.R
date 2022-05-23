################################################################################
#MI Prerun Script
#Load the .nrg files for each residue in the protein and collate them into a master file
################################################################################

library(bio3d)

library(foreach)

library(doParallel)

library(stringr)

library(tidyverse)

################################################################################
#Protein selection 
#Trajectory selection
#Allows for automatic loading of pdb and nrg from their appropriate directories
################################################################################

Protein_of_Interest <- "1v4s"
Production_Name <- "Restart2"

setwd(paste("/scratch1/ncleme2/Info_States/", Protein_of_Interest, "/charmm-gui/namd/NRGs/", Production_Name, sep = ""))

#choose num with
#test <- read.table(file = "100.nrg")
number_of_frames <- 10135

#check this number every time because shennanigans can happen
start_res <- 14

################################################################################
#Load pdb and pull the protein's number of residues
################################################################################

pdb <- bio3d::read.pdb(file = paste("/scratch1/ncleme2/Info_States/", Protein_of_Interest, "/charmm-gui/namd/step3_input.pdb", sep = ""))

ca_inds <- atom.select(pdb, elety = "CA")
resno <- length(ca_inds$atom)
#resno <- 4

################################################################################
#Open residue nrg files and collate the total energy into a master enegy file
################################################################################


residue_vector <- start_res:(resno-1)
#residue_vector <- 22:25
master <- data.frame(matrix(data=NA, nrow = number_of_frames, ncol=length(residue_vector)))
name_index <- c()

assignment_counter <- 1

for (r in residue_vector) {
  
  current_residue <- paste(r, ".nrg", sep = "")
  assign(current_residue, read.table(file = paste(residue_vector[assignment_counter], ".nrg", sep = ""), stringsAsFactors = F))
  master[, assignment_counter] <- get(paste(r, ".nrg", sep = ""))[11]
  name_index[assignment_counter] <- paste(r, ".nrg", sep = "")
  assignment_counter <- assignment_counter + 1
  
  #cleanup the excessive number of variables in the environment to save memory
  rm(list = ls(pattern = ".nrg", envir = .GlobalEnv), envir = .GlobalEnv)
  
}

master <- as.tibble(master[2:nrow(master),])
master <- master %>% mutate_all(funs(str_replace(., "\\+", "")))
master <- as.data.frame(sapply(master, as.numeric)) 

################################################################################
#route and save master nrg to MI_Analysis_Data in the protein's directory
################################################################################

setwd(paste("/scratch1/ncleme2/Info_States/", Protein_of_Interest, "/charmm-gui/namd/MI_Analysis_Data", sep = ""))

write.table(master, paste(Protein_of_Interest, "_", Production_Name, "_master.nrg", sep = ""))

