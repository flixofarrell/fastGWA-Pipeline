library('data.table')
require(data.table)

args <- commandArgs(trailingOnly=TRUE)


system.time(m_dt <- fread(args[1]))

system.time(n_dt <- fread(args[2]))


m_dt <- merge(m_dt, n_dt, all = TRUE)



fwrite(m_dt, args[1])
