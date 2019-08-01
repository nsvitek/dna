#check results of demultiplex step, see if any specimen has too few retained reads.

# Dependencies -----
library(dplyr)

# set file paths ----- 
path_in<-"C://Users/N.S/Dropbox/Documents/Dissertation/modern/peromyscus_maniculatus/radseq"
# path_in<-"D://rs"

name_samples<-"retained_reads.txt"

# read, clean ---------
specimen_df<-paste(path_in,name_samples,sep="/") %>% 
  read.table(., header = TRUE, sep = "\t")

#look at % of reads retained for each specimen
specimen_df<-specimen_df %>% mutate(.,percent = Retained/Total) 


specimen_df$Retained %>% hist()
specimen_df$percent %>% hist()

#look at which specimens have a very low % retained reads
specimen_df[which(specimen_df$percent<0.4),] #still at least 400k reads.

#order the specimens by # reads retained
specimen_df<-specimen_df[order(specimen_df$Retained),]

#look at specimens with fewest retained reads
specimen_df$Retained[1:20] %>% hist()
specimen_df$Filename[1:5]


#look at lowest number of reads
specimen_df[which(specimen_df$Filename=="UWBM80320"),]
specimen_df[which(specimen_df$Filename=="MSB150769"),]

specimen_df[c(95:97),]

