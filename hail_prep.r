#Hail Preparer: 
#Reformat fastGWA files to .tsv readable 
#by Hail (for QQ plots etc).

#Developer: Felix O'Farrell
#May 2020


library(data.table)

#command line args to plug into pipe
args <- commandArgs(trailingOnly=TRUE)

#read in master fastGWA
DT <- fread(args[1])

#add 'rs' string to SNP collumn
DT[,SNP:=paste0("rs",SNP)]

#paste CHR and SNP to new locus collumn
#allows to be read by Hail
DT[,locus:=paste0(CHR, sep = ":", SNP)]

#write to h_master tsv file
fwrite(DT, args[2], sep = '\t')




