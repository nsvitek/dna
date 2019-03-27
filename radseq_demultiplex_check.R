#check results of demultiplex step, see if any specimen has too few retained reads.

# Dependencies -----
library(dplyr)

# set file paths ----- 
path_in<-"C://Users/N.S/Dropbox/Documents/Dissertation/modern/onychomys_leucogaster/radseq_analysis/02-demultiplex"
# path_in<-"D://Dropbox/Documents/Dissertation/modern/onychomys_leucogaster"
# path_scripts<-"C://cygwin/home/N.S/scripts/cataloging"

name_samples<-"retained_reads.txt"

# read, clean ---------
specimen_df<-paste(path_in,name_samples,sep="/") %>% 
  read.table(., header = TRUE, sep = "\t")

#look at % of reads retained for each specimen
specimen_df<-specimen_df %>% mutate(.,percent = Retained/Total) 


specimen_df$Retained %>% hist()
specimen_df$percent %>% hist()

specimen_df[which(specimen_df$percent<0.4),] #still at least 400k reads.

specimen_df<-specimen_df[order(specimen_df$Retained),]

specimen_df$Filename[78:82]
