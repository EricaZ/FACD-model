setwd("~/FACDapp/Newtons_code")
source("ComputeIS_repar.R")
source("ComputeGrad_repar.R")
source("ComputeCondDur_repar.R")
source("ComputeC_repar.R")
source("FitFACD_repar.R")
source("SimulateACDgeneral.R")

id.str <- "test"
xd <- SimulateACDgeneral(param = list(w = 0.1, a.vec = 0.2, b.vec = 0.6, r = 2), offset = 200, num.n = 1000, num.rep = 1, distrib = "frechet") 
simfit <- FitFACD_repar(x = xd[,1], id.str)

r <- simfit$param[1]
s <- simfit$param[2]
a <- simfit$param[3]
b <- simfit$param[4]
w <- s * gamma(1 - 1 / r)
alpha <- a * w 

xdata <- seq(from=0, to=7, by=0.01)
par.str <- paste(round(r,digits=3),round(w,digits=3), round(alpha,digits=3), round(b,digits=3), sep=", ")

res.den1 <- density(simfit$res, adj=1)
theo.den1 <- dfrechet(xdata, loc = 0, scale = s, shape = r, log = FALSE)
emp.q1 <- quantile(simfit$res, probs = seq(0.01, 0.99, 0.01))
theo.q1 <- qfrechet(p=seq(0.01, 0.99, 0.01), loc= 0, scale= s, shape = r, lower.tail = TRUE)

pdf(file=paste(id.str, '_plot.pdf', sep=""), width = 8, height = 4)

par(mfrow=c(1,2))
plot(theo.den1~xdata, type="l", col="blue", lty="dashed", main = paste("FACD(1,1):", par.str))
lines(res.den1)
plot(theo.q1, emp.q1, ylim=c(0,1), xlim=c(0,1))
lines(0:20, 0:20)

#################################################
setwd("~/FACDapp/blocktrade/20140106data")

id.str <- "test_repar"
testfit <- FitFACD_repar(x = dur.adj, id.str)

# Output plots
r <- testfit$param[1]
s <- testfit$param[2]
a <- testfit$param[3]
b <- testfit$param[4]
w <- s * gamma(1 - 1 / r)
k <- mean(testfit$res)
m0.alpha <- a * w 
m2.alpha <- a * k

mean(testfit$res)

xdata <- seq(from=0, to=7, by=0.01)
m0par.str <- paste(round(r,digits=3),round(w,digits=3), round(m0.alpha,digits=3), round(b,digits=3), sep=", ")
m2par.str <- paste(round(r,digits=3),round(k,digits=3), round(m2.alpha,digits=3), round(b,digits=3), sep=", ")

res.den1 <- density(testfit$res, adj=1)
theo.den1 <- dfrechet(xdata, loc = 0, scale = s, shape = r, log = FALSE)
emp.q1 <- quantile(testfit$res, probs = seq(0.01, 0.99, 0.01))
theo.q1 <- qfrechet(p=seq(0.01, 0.99, 0.01), loc=0, scale= s, shape = r, lower.tail = TRUE)

scaled.res <- testfit$res / k
res.den2 <- density(scaled.res, adj=1)
theo.den2 <- dfrechet(xdata, loc = 0, scale = s / k, shape = r, log = FALSE)
emp.q2 <- quantile(scaled.res, probs = seq(0.01, 0.99, 0.01))
theo.q2 <- qfrechet(p=seq(0.01, 0.99, 0.01), loc=0, scale= s / k, shape = r, lower.tail = TRUE)


pdf(file=paste(id.str, '_plot.pdf', sep=""), width = 8, height = 4)
par(mfrow=c(1,2))
plot(theo.den1~xdata, type="l", col="blue", xlim=c(0,1), ylim=c(0, max(theo.den1, res.den1$y)), lty="dashed", main = paste("FACD(1,1):", m0par.str))
lines(res.den1)
plot(theo.q1, emp.q1)
lines(0:20, 0:20)

plot(theo.den2~xdata, type="l", col="blue", ylim = c(0, max(theo.den2, res.den2$y)), lty="dashed", main = paste("FACD(1,1):", m2par.str))
lines(res.den2)
plot(theo.q2, emp.q2)
lines(0:20, 0:20)

# plot(theo.q1, emp.q1/mean(testfit$res), ylim=c(0,1), xlim=c(0,0.2))
# lines(0:20, 0:20)
dev.off()
