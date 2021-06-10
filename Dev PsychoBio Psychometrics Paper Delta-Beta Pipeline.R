#Deltab-beta reduction and correlation pipeline
# This pipeline takes second-by-second exports (.txt format) of delta and beta power for each participant, and concatenates them to be 
# analyzed at the second-by-second level. It also creates a dataset with individual delta-beta coupling scores for each participant, 
# to be analyzed in that way as well. Finally, at the end of the script we also include reliability analyses of delta-beta coupling scores and 
# reliability of average delta and beta power. 


# Installing packages 
list.of.packages <- c("tidyverse", "stats", "tidyr", "reshape2", "psych")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# Loading packages
lapply(list.of.packages, require, character.only = TRUE)


#Setting the Working directory for Delta files (this is to the folder where all the delta power exports for participants are)
path = "Y:/LANTS EEG/8MO/Delta/"
setwd(path)

#Reading files in and merging them
deltafiles <- dir(path, pattern=".txt")
datlist <- lapply(deltafiles, read.table, sep="\t", header = TRUE)
long = do.call(rbind.fill, datlist)

#Creating an id column  that will go inside the merged file (this will include both BLN and OCM)
long$id  <- substr(deltafiles[], 17, 21) #change these numbers to match where your id number starts and ends

#Log transformation on the data set
LOG <- function(x){log(x)}
long[1:16] <- lapply(long[1:16], LOG) #applying it to only the channel columns  and not the id column


#Creating averaged power variables for delta from the frontal, central, and parietal electrodes of interest
long$DeltaFrontal<- rowMeans(long[c('X3','X30')], na.rm = TRUE)
long$DeltaCentral<- rowMeans(long[c('X8','X24', 'X25')], na.rm = TRUE)
long$DeltaParietal<- rowMeans(long[c('X14','X13', 'X19')], na.rm = TRUE)

#Sub-setting data to only export ID column and Average Power variables
DeltaAve<- subset(long,
                  select = c("id", "DeltaFrontal","DeltaCentral","DeltaParietal"))

###########################################################################################
# Repeating the process above with beta exports

#Setting the Working directory for Beta files
path = "C:/Users/bua25/Documents/Free Play/SRCD 2019/EEG/8MO/Beta/"
setwd(path)

#Reading files in and merging them
betafiles <- dir(path, pattern=".txt")
betalist <- lapply(betafiles, read.table, sep="\t", header = TRUE)
longB = do.call(rbind.fill, betalist)


#Creating an id column  that will go inside the merged file 
longB$id  <- substr(betafiles[],17, 21)

#Log transformation on the data set
LOG <- function(x){log(x)}
longB[1:16] <- lapply(longB[1:16], LOG)


#Creating the averaged power for beta including fronto-central electrodes of interest
longB$BetaFrontal<- rowMeans(longB[c('X3','X30')], na.rm = TRUE)
longB$BetaCentral<- rowMeans(longB[c('X8','X24', 'X25')], na.rm = TRUE)
longB$BetaParietal<- rowMeans(longB[c('X14','X13', 'X19')], na.rm = TRUE)


#Subsetting data to only export ID column and BetaPower averaged
BetaAve<- subset(longB,
                 select = c("id", "BetaFrontal", "BetaCentral","BetaParietal"))


# Binding Delta and Beta files together and getting rid of the extra id column (column 3)
DBcoupling<- cbind(DeltaAve, BetaAve)
DBcoupling[3] <- NULL

#Exporting long dataset (this is the second-by-second data set)
write.table(DBcoupling, file="C:/Users/bua25/Documents/Free Play/SRCD 2019/EEG/8DBcoupling.csv", sep=",",row.names=F)

# Computing delta-beta correlation for each participant
xx <- data.frame(group = rep(1:4, 100), a = rnorm(400) , b = rnorm(400) )
head(xx)

library(plyr)
func <- function(DBcoupling)
{
  return(data.frame(FrontalSynchrony = cor(DBcoupling$DeltaFrontal, DBcoupling$BetaFrontal),
                    CentralSynchrony = cor(DBcoupling$DeltaCentral, DBcoupling$BetaCentral),
                    ParietalSynchrony = cor(DBcoupling$DeltaParietal, DBcoupling$BetaParietal)))
}

#running function to produce correlations
Synchrony<-ddply(DBcoupling, .(id), func)

#exporting wide data set with delta-beta coupling scores (correlations) for each participant across frontal, central, and parietal regions
write.csv(Synchrony, file="C:/Users/bua25/Documents/Free Play/SRCD 2019/EEG/8DBSynchronyScores.csv",sep=",",row.names=F)


###################################################################################

#Delta-Beta coupling Reliability

#reading second by second dataset
epoch8<-read.csv("8moDBsecbysec.csv", header = T)

#splitting the second-by-second dataset by even and odd rows
epoch8%>%
  group_by(id)%>%
  arrange(time)%>%
  dplyr::filter(row_number() %% 2 == 0) ->df8e

epoch8%>%
  group_by(id)%>%
  arrange(time)%>%
  dplyr::filter(row_number() %% 2 == 1) ->df8o

# computing delt-beta coupling scores from these two datasets, to reflect coupling from even and odd rows
func <- function(x)
{
  return(data.frame(FrontalSync = cor(x$DeltaFrontal, x$BetaFrontal),
                    CentralSync = cor(x$DeltaCentral, x$BetaCentral),
                    ParietalSync = cor(x$DeltaParietal, x$BetaParietal)))
}

SyncE<-ddply(df8e, .(id), func)
SyncO<-ddply(df8o, .(id), func)


# Merging the even and odd synchrony scores computed above
EOr8<-merge(SyncE, SyncO, by="id")

# Computing correlations between odd and even coupling scores
EO8cor <- cor(EOr8[-1],use="complete.obs", method = "pearson")

# Applying the Spearman-Brown correction
2*EO8cor[4]/ 1+ EO8cor[4]

2*EO8cor[11]/ 1+EO8cor[11]

2*EO8cor[18]/ 1+EO8cor[18]



## Computing random half, split half reliability on the average power data

epoch8M<-epoch8%>%
  group_by(id)%>%
  dplyr::summarise_each(funs(mean))

#Delta power
splithalf.r(epoch8M[, c(3:5)], sims = 100, graph = F, seed = 2)

#Beta power
splithalf.r(epoch8M[, c(6:8)], sims = 100, graph = TRUE, seed = 2)


## These analyses were repeated on the 12- and 18-month data.
