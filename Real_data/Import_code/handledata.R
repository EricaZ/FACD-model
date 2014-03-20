handledata <- function(a) {
  
  threshold <- 500000
  
  amount <- a$V2 * a$V3
  a1 <- a[which(amount > threshold),]
  
  a1$V1 <- mdy_hms(a1$V1)  
  a2 <- a1[order(a1$V1),]
  
  # compute time of the day (in seconds)
  TofD <- as.numeric(difftime(a2$V1, floor_date(a2$V1, "day"), unit="secs")) 
  
  # Retain trades that occurred during 9:30 - 12:00 & 13:00 - 16:00
  a3 <- cbind(a2, TofD)
  a3 <- subset(a3, (a3$TofD >= 34200) & (a3$TofD < 57600) & ((a3$TofD < 43200)|(a3$TofD >= 46800)))
  
  # Obtain indicator of the first trade of a day (nonzero if true, zero if false)
  ## a vector of the current date
  getdate1 <- floor_date(a3$V1[-nrow(a3)], "day") 
  ## a vector of the next date
  getdate2 <- floor_date(a3$V1[-1], "day") 
  ## a vector of the difference between the previous date and the current date
  ## (nonzero iff it is the first trade of a day)
  datediff <- c(1, as.numeric(difftime(getdate2, getdate1, unit="days"))) 
  
  # Obtain indicator of morning(1) or afternoon trade(0)
  ismorning <- rep(0, nrow(a3))
  ismorning[a3$TofD < 43200] <- 1
  morningdiff <- c(0, ismorning[-1]-ismorning[-nrow(a3)])
  ## indicator of the first afternoon trade(1) of a day, or else(0)
  isfirst <- rep(0, nrow(a3))
  isfirst[(datediff==0) & (morningdiff== -1)] <- 1
  
  # compute duration and add duration column into data set
  dur <- c(0,int_length(int_diff(a3$V1)))
  a3 <- cbind(a3, dur)
  
  # remove the first trade of each day and the first trade of each afternoon
  a4 <- a3[-which((datediff != 0)|(isfirst == 1)), ]
  
  # remove obs with zero duration
  a4 <- a4[-which(a4$dur == 0), ]
  
  # Retain trades that occurred during 10:00 - 12:00 & 13:00 - 15:30
  a5 <- subset(a4, (a4$TofD >= 36000) & (a4$TofD < 55800))
  
  return(a5)
}

