# Split the master data file
hsiStocks <- read.csv("C:/Users/yaoz/Documents/FACDapp/blocktrade/20140224data/20140224data.csv",header=F)
setwd("C:/Users/yaoz/Documents/FACDapp/blocktrade/20140224data/stocksdata")

for (i in 1: 50) {
  temp <- na.exclude(hsiStocks[,(5*(i-1)+1):(5*(i-1)+3)])
  temp <- cbind(temp, rep(i, nrow(temp)))
  write.table(temp, file = paste("stock", i, ".csv", sep=""), col.names=F, row.names=F)
}

rm(temp)
rm(hsiStocks)

# Import data of all 50 stocks and combine them into one data frame
for (i in 1:50) {
  assign(paste("stock",i,sep=""), read.table(paste("C:/Users/yaoz/Documents/FACDapp/blocktrade/20140224data/stocksdata/stock", i, ".csv", sep="")))
}

a <- NULL
for (i in lapply(ls(pattern="stock"), get)) {
  a <- rbind(a, i)
}

rm(i, list=ls(pattern="stock"))
save(a, file = "a.RData")