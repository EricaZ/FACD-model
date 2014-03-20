library(evd)
r <- 1.2
xdata <- seq(from=0, to=7, by=0.01)
theo.den <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)

plot(theo.den~xdata, type="l", col="black", lwd=2, ylim=c(0, max(theo.den)), xlim=c(0,3.5))

r <- 1.6
theo.den <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
lines(theo.den~xdata, col="gray30", lwd=2)

r <- 2
theo.den <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
lines(theo.den~xdata, col="gray60", lwd=2)

r <- 2.4
theo.den2 <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
lines(theo.den2~xdata, col="gray90", lwd=2)

r <- 3
theo.den <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
lines(theo.den~xdata, col="darkblue", lwd=2)

r <- 4
theo.den <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
lines(theo.den~xdata, col="blue3", lwd=2)

r <- 5
theo.den2 <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
lines(theo.den2~xdata, col="dodgerblue", lwd=2)

r <- 10
theo.den2 <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
lines(theo.den2~xdata, col="skyblue", lwd=2)




r <- 1.2
xdata <- seq(from=0, to=7, by=0.01)
theo.den <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)

plot(theo.den~xdata, type="l", col="black", lwd=2, ylim=c(0, 0.03), xlim=c(3,7))

r <- 1.6
theo.den <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
lines(theo.den~xdata, col="gray30", lwd=2)

r <- 2
theo.den <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
lines(theo.den~xdata, col="gray60", lwd=2)

r <- 2.4
theo.den2 <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
lines(theo.den2~xdata, col="gray90", lwd=2)

# r <- 3
# theo.den <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
# lines(theo.den~xdata, col="darkblue", lwd=2)
# r <- 4
# theo.den <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
# lines(theo.den~xdata, col="blue3", lwd=2)
# r <- 5
# theo.den2 <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
# lines(theo.den2~xdata, col="dodgerblue", lwd=2)
# r <- 10
# theo.den2 <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
# lines(theo.den2~xdata, col="skyblue", lwd=2)

r <- 1.2
xdata <- seq(from=0, to=9, by=0.01)
theo.den <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)

plot(theo.den~xdata, type="l", col="black", lwd=2, ylim=c(0, 0.01), xlim=c(5.5,9))

r <- 1.6
theo.den <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
lines(theo.den~xdata, col="gray30", lwd=2)

r <- 2
theo.den <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
lines(theo.den~xdata, col="gray60", lwd=2)

r <- 2.4
theo.den2 <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
lines(theo.den2~xdata, col="gray90", lwd=2)


