library(lubridate)
library(evd)
setwd("~/FACDapp/Newtons_code")
source("handledata.R")

a5 <- handledata(a)

plot(a5$dur,type="l")
summary(a5$TofD)

a1 <- subset(a5, a5$TofD <= 43200)
plot(a1$dur,type="l")
lines(a1$TofD/1000)
summary(a1$TofD)

a2 <- subset(a5, a5$TofD > 43200)
plot(a2$dur,type="l")
lines(a2$TofD/1000)
summary(a2$TofD)

# n <- 5
# plot(a1$TofD, a1$dur, main = paste("spline[fun](.) through", n, "points"))
# lines(spline(a1$TofD, a1$dur, n = n), col = 2)
timeofday1 <- a1$TofD
m1 <- smooth.spline(timeofday1, a1$dur, nknots=5)
plot(m1, type="l", xlab="sec", ylab=expression(phi))
timeofday2 <- a2$TofD
m2 <- smooth.spline(timeofday2, a2$dur, nknots=6)
plot(m2, type="l", xlab="sec", ylab=expression(phi))

# m2 <- smooth.spline(timeofday2, a2$dur, nknots=6)
# plot(m2, type="l", xlab="sec", ylab=expression(phi))

# library(lubridate)
# hms("10:00:00") + minutes(round(m1$fit$knot*60*2,0))
# m2$fit$knot*60*2.5

dur.smooth <- rep(0, nrow(a5))
for (i in 1: nrow(a5)) {
  if (a5$TofD[i] <= 43200) {
    dur.smooth[i] <- predict(m1, a5$TofD[i])$y
  } else {
    dur.smooth[i] <- predict(m2, a5$TofD[i])$y
  } 
}

dur.adj <- a5$dur / dur.smooth
plot(dur.adj,type="l")

save(dur.adj, file="dur.RData")

