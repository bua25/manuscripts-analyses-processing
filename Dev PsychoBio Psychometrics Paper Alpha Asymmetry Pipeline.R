# ASM reliability pipeline
# In this script, we take alpha power exports and compute asymmetry scores second-by-second, in order to compute reliability metrics for each
# participant, and average these metrics across the sample. We do a similar analysis for alpha power. While we include only one example of each 
# script for simplicity, our actual analyses repeated these scripts on alpha power for eyes-open files, eyes-closed files, and files that 
# included both conditions. We also repeated these scripts on files for each time visit (8-, 12-, and 18-months).

list.of.packages <- c("tidyverse", "psych", "ggplot2", "stats", "reshape2", "dplyr", "plyr", "multicon")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# Loading packages
lapply(list.of.packages, require, character.only = TRUE)


#Setting the Working directory for Alpha power files
path = "Y:/LANTS EEG/8MO/ASM reliability/"
setwd(path)

#Reading files in and merging them
ASMfiles <- dir(path, pattern=".txt")

# looping over data to clean files, log transform alpha power, create a time columns since alpha is currently in second-by-second format,
# computing second-by-second asymmetry scores between electrode pairs as well as cluster, in order to compute reliability metrics per participant 
for(i in 1:length(ASMfiles)){
  
  #creating print statement to track loop progress
  dataname <- ASMfiles[i]
  print(paste0("Now in ", dataname))
  
  # read in file 
  df <- read.table(ASMfiles[i], sep = "\t", header = TRUE)
  
  # running log function
  log(df[,-1]) -> df
  
  
  # create an id variable
  df$id  <- substr(ASMfiles[i], 17, 21)
  
  # adding in time variable
  df$time <- seq_along(df$id)
  
  
  #computing  specific ASM scores
  if("X32" %in% colnames(df) && "X1" %in% colnames(df)){
    df$Fp1Fp2ASM <- df$X32 - df$X1
  } else {
    df$Fp1Fp2ASM <- NA
  }
  
  if("X30" %in% colnames(df) && "X3" %in% colnames(df)){
    df$F3F4ASM <- df$X30 - df$X3
  } else {
    df$F3F4ASM <- NA
  }
  
  if("X31" %in% colnames(df) && "X4" %in% colnames(df)){
    df$F7F8ASM <- df$X31 - df$X4
  } else {
    df$F7F8ASM <- NA
  }
  
  #computing region cluster averages
  df$left <- rowMeans(df[, c("X1","X3","X4")], na.rm =T)
  
  df$right <- rowMeans(df[, c("X30","X31","X32")], na.rm = T)
  
  #Computing cluster based ASM scores
  df$ASMcl <- df$right - df$left
  
  
  # saving cleaned version of the files 
  write.csv(df, file.path("Y:/LANTS EEG/8MO/ASM reliability/CSV", 
                          paste(ASMfiles[i], ".csv")), row.names = FALSE)
}


## Reliability Analysis for ASM scores

# create list of processed files
path = "Y:/LANTS EEG/8MO/ASM reliability/CSV"
setwd(path)

asmfiles <- dir(path, pattern=".csv")

# looping over data compute split-half reliability and cronbach's alpha on each participant's second by second ASM score
for(i in 1:length(asmfiles)){
  
  #creating print statement to track loop progress
  dataname <- asmfiles[i]
  print(paste0("Now in ", dataname))
  
  # read in file 
  df <- read.table(asmfiles[i], sep = ",", header = TRUE)
  
  #sub-setting the files to only have ASM scores
  df1<- df[, c("id","time","Fp1Fp2ASM","F3F4ASM","F7F8ASM","ASMcl")]
  
  # computing reliability and exporting to another object
  ASMrelF<- splithalf.r(df1[, c(3:6)], sims = 100, graph = F, seed = 45)
  
  
  write.csv(ASMrelF, file.path("Y:/LANTS EEG/8MO/ASM reliability/CSV/RelResults", 
                               paste(asmfiles[i], ".csv")), row.names = FALSE)
}

# Binding all participants to get a full dataset for alpha power and asm across full time series
path = "Y:/LANTS EEG/8MO/ASM reliability/CSV/"
setwd(path)
ASMpowerFULL <- dir(path, pattern=".csv")
ASMpower8<- lapply(ASMpowerFULL, read.csv, header = TRUE)
ASMpower8df = do.call(rbind.fill, ASMpower8)
write.table(ASMpower8df, file = "Y:/LANTS EEG/Processed Datasets/8MO/asm8_all_segments.csv", sep = ",", row.names = T)


# Binding reliability results from each participant to get average reliability metrics across the sample.
# create list of processed files
path = "Y:/LANTS EEG/8MO/ASM reliability/CSV/RelResults"
setwd(path)

ASMRresults <- dir(path, pattern=".csv")
RelASM8<- lapply(ASMRresults, read.csv, header = TRUE)

# binding all the reliability result files together 
RelASM8All = do.call(rbind.fill, RelASM8)
write.table(RelASM8All, file = "Y:/LANTS EEG/Processed Datasets/8MO/ReliabilityASM8F.csv", sep = ",", row.names = T)

# Getting averages
RelResultsTable8<-psych::describe(RelASM8All)

# Exporting results
write.table(RelResultsTable8, file = "Y:/LANTS EEG/Processed Datasets/8MO/RelResultsTable8F.csv", sep = ",", row.names = T)


## Summary
# This processing was repeated on files from the eyes-open and eyes-closed conditions separately.


#########################################################################################

## Reliability analysis for alpha power in the Eyes-open condition

# create list of processed files of alpha power
path = "Y:/LANTS EEG/12MO/EOP/CSV"
setwd(path)

asmfiles <- dir(path, pattern=".csv")

# looping over data to clean file and create id
for(i in 1:length(asmfiles)){
  
  #creating print statement to track loop progress
  dataname <- asmfiles[i]
  print(paste0("Now in ", dataname))
  
  # read in file 
  df <- read.table(asmfiles[i], sep = ",", header = TRUE)
  
  #subsetting the files to only have electrodes and averages 
  if("X1" %in% colnames(df) && "X3" %in% colnames(df)&& "X4" %in% colnames(df)){
    dfleft<- df[, c("id","time","X1","X3","X4","left")]
  } else {
    dfleft$X1 <- NA
    dfleft$X3 <- NA
    dfleft$X4 <- NA
  }
  
  if("X30" %in% colnames(df) && "X31" %in% colnames(df)&& "X32" %in% colnames(df)){
    dfright<- df[, c("id","time","X30","X31","X32","right")]
  } else {
    dfleft$X1 <- NA
    dfleft$X3 <- NA
    dfleft$X4 <- NA
  }
  
  # computing reliability for alpha power at left and right hemispheres and exporting to another object
  ASMrelLeft<- splithalf.r(dfleft[,c(3:6)], sims = 100, graph = F, seed = 45)
  
  ASMrelRight<- splithalf.r(dfright[,c(3:6)], sims = 100, graph = F, seed = 45)
  
  write.csv(ASMrelLeft, file.path("Y:/LANTS EEG/12MO/EOP/CSV/RelPowerL", 
                                  paste(asmfiles[i], ".csv")), row.names = FALSE)
  write.csv(ASMrelRight, file.path("Y:/LANTS EEG/12MO/EOP/CSV/RelPowerR", 
                                   paste(asmfiles[i], ".csv")), row.names = FALSE)
  
}

# Bindings all the reliability results to compute mean reliability metrics across the sample
#Right 
path = "Y:/LANTS EEG/12MO/EOP/CSV/RelPowerR"
setwd(path)

PowerR <- dir(path, pattern=".csv")
RelPowerR12<- lapply(PowerR, read.csv, header = TRUE)

# binding all the reliability result files together 
RelPR12All = do.call(rbind.fill, RelPowerR12)

# averaging metrics across all participants
RelPRTable12<-psych::describe(RelPR12All)
write.table(RelPRTable12, file = "Y:/LANTS EEG/Processed Datasets/12MO/RelpowerTable12EOPr.csv", sep = ",", row.names = T)


#Left 
path = "Y:/LANTS EEG/12MO/EOP/CSV/RelPowerL"
setwd(path)

PowerL <- dir(path, pattern=".csv")
RelPowerL12<- lapply(PowerL, read.csv, header = TRUE)

# binding all the reliability result files together 
RelPL12All = do.call(rbind.fill, RelPowerL12)

# averaging metrics across all participants
RelPLTable12<-psych::describe(RelPL12All)
write.table(RelPLTable12, file = "Y:/LANTS EEG/Processed Datasets/12MO/RelpowerTable12EOPl.csv", sep = ",", row.names = T)

# Summary
# This script for processing reliability of alpha power covers the eyes-open condition. For our analyses, we applied the same script to files
# for the eyes-closed condition. 
