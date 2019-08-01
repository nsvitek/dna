data<-read.csv("C:/Users/N.S/Documents/Dissertation/isd_window_reps.csv")
(loci<-apply(data,2,mean))
(readnum<-loci*100)
(locnum<-readnum*7)
(newguess<-reads/locnum)
# coverage<-20


samples_per_locality<-6
reads<-110000000 #110,000,000; NextSeq500mid

# localities<-round(reads/(coverage*samples_per_locality*loci))
# x-coverage * #-loci * #-samples < or =  #-reads 

str(data)
hist(data$X250.400bp,xlab="number of loci")
mean(data$X250.400bp)

#notes from Peterson 2012. How did they do it? 
# they aimed for: 1.5-2k loci/individual; window 275-325
# for 10x of that, they thought they needed 20-50k reads/individual
# when they sampled 2-3k regions, they got  1886 SNPs in 1638 loci, with about 
# 10x coverage and ~250,000 reads/individual

#I'd like about the same numbers of SNPs, 2k seems pretty good to me.
#Let's say that means 3k loci sampled.
#But here's the rub, they got saturation after about 200k reads per critter in sim
#with narrow size selection (30 bp); so if you want ~2k good SNPs, you need about
#100x more reads. 

15*6
# to try again: http://www.ncbi.nlm.nih.gov/guide/howto/dwn-genome/
