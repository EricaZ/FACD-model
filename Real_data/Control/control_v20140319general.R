library(lubridate)
library(evd)
setwd("~/FACDapp/Newtons_code")
source("ComputeIS.R")
source("ComputeGrad.R")
source("ComputeCondDur.R")
source("ComputeC.R")
source("FitFACD.R")
# source("FitWACD.R")
source("FitWACD2.R")
source("FitEACD.R")
source("ComputeOmega.R")
source("ComputeOmega2.R")
source("ComputeOmega3.R")
source("handledata.R")
source("ComputeIS_dual.R")
source("ComputeC_dual.R")
source("ComputeLoglik_dual.R")
source("FitACD_2-stage_correct.R")
# source("FitFACD_fix.R")
source("FitFACD_fix2.R")
source("FitACDgeneral.R")
source("ComputeISgeneral.R")
source("ComputeCondDurgeneral.R")
source("ComputeGradgeneral.R")

setwd("C:/Users/yaoz/Documents/FACDapp/blocktrade/20140106data")

# init.param = list(w = 0.1, a.vec = c(0.1,0.1), b.vec = 0.5, r = 1)

a.vec <- 0.1
b.vec <- c(0.3,0.3)

# a.vec <- c(0.1,0.1)
# b.vec <- 0.6

# a.vec <- c(0.1,0.1)
# b.vec <- c(0.3,0.3)

init.param = list(w = 0.1, a.vec = a.vec, b.vec = b.vec, r = 1)
p <- length(a.vec)
q <- length(b.vec)

id.str <- paste("wb_", p, q, sep="")

wbfit <- FitACDgeneral(dur.adj, portmanteau = TRUE, distrib = "weibull", p = p, q = q, maxit = 100, id.str = id.str, init.param = init.param) 

# Output plots
xdata <- seq(from=0, to=9, by=0.01)

res.den2 <- density(wbfit$res, adj=1)
r <- wbfit$param[1]
theo.den2 <- dweibull(xdata, shape = r, scale = 1 / gamma(1 + 1 / r))
emp.q2 <- quantile(wbfit$res, probs = seq(0.01, 0.99, 0.01))
theo.q2 <- qweibull(p=seq(0.01, 0.99, 0.01), shape = r, scale = 1 / gamma(1 + 1 / r), lower.tail = TRUE)

# pdf(file=paste(id.str, '_plot.pdf', sep=""), width = 16, height = 8)
# par(mfrow=c(2,4), mar=c(5,5,5,5))
plot(theo.den2~xdata, type="l", col="blue", ylim=c(0, 1.5), lty="dashed", main = paste("WACD(", p, ",", q, ")"))
lines(res.den2)
plot(theo.den2~xdata, type="l", col="blue", xlim=c(5, 8), ylim=c(0, 0.015), lty="dashed", main = paste("mean residual = ", round(mean(wbfit$res), digits=3)))
lines(res.den2)
plot(theo.q2, emp.q2, ylim=c(0,6), xlim=c(0,5))
lines(0:20, 0:20)
acf(wbfit$res, ylim=c(-0.03,0.03), main=paste("WACD(", p, ",", q, ")"))
# dev.off()

wbfit0106p1q2 <- wbfit
save(wbfit0106p1q2, file="wbfit0106p1q2.RData")
