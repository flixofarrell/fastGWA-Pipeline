#Concatenater: 
#Concatenate dataframe to master dataframe
#Needs to be executed with one thread to make sure
#no rows are duplicated/misssed out.

#Developer: Felix O'Farrell
#May 2020


library('data.table')
require(data.table)

#command line args to plug into pipe
args <- commandArgs(trailingOnly=TRUE)

#read in master
m_dt <- fread(args[1])

#read in Chrx+1 .fastgwa file
n_dt <- fread(args[2])

#concat Chrx+1 to master 
m_dt <- merge.data.table(m_dt, n_dt, all = TRUE)

#write to master file
fwrite(m_dt, args[1])

