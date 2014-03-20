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
source("FitACDgeneral_fix.R")
source("ComputeISgeneral_dual.R")

setwd("C:/Users/yaoz/Documents/FACDapp/blocktrade/20140210data")

# init.param = list(w = 0.1, a.vec = c(0.1,0.1), b.vec = 0.5, r = 1)

a.vec <- 0.1
b.vec <- c(0.3,0.3)

a.vec <- c(0.1,0.1)
b.vec <- 0.6

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

# =====================================================
a.vec <- 0.1
b.vec <- c(0.3,0.3)

init.param = list(w = 0.1, a.vec = a.vec, b.vec = b.vec)
p <- length(a.vec)
q <- length(b.vec)

# mfr2fit <- FitACDgeneral_fix(dur.adj, portmanteau = TRUE, distrib = "modifrechet", p = p, q = q, r = 2, maxit = 100, id.str = id.str, init.param = init.param) 

id.str <- paste("mfr2_", p, q, sep="")
rseq <- seq(from = 1.5, to = 2.5, by = 0.1)
for (r in rseq) {
  try({
    assign(paste("mfr2fit_", r, sep=""), FitACDgeneral_fix(dur.adj, portmanteau = TRUE, distrib = "modifrechet", p = p, q = q, r = r, maxit = 100, id.str = id.str, init.param = init.param))
  }, TRUE)
}

pdf(file=paste(id.str, '_plot.pdf', sep=""), width = 16, height = 5)
xdata <- seq(from=0, to=9, by=0.01)
rseq <- seq(from = 1.5, to = 2.5, by = 0.1)
for (r in rseq) {
  # Output plots
  try({       
    mfr2fit <- get(paste("mfr2fit_", r, sep=""))
    res.den5 <- density(mfr2fit$res, adj=1)
    p <- 1 - 1 / r
    g <- gamma(p)
    q <- 9.21^(-1 / r)  # 9.21 = -log(delta) if delta = 0.0001
    bga <- 1 / (g - q)
    aga <- -q * bga
    theo.den5 <- dfrechet(xdata, loc = aga, scale = bga, shape = r, log = FALSE)
    emp.q5 <- quantile(mfr2fit$re, probs = seq(0.01, 0.99, 0.01))
    theo.q5 <- qfrechet(p=seq(0.01, 0.99, 0.01), loc = aga, scale = bga, shape = r, lower.tail = TRUE)
    
    par(mfrow=c(1,4), mar=c(5,5,5,5))
    plot(theo.den5~xdata, type="l", col="blue", ylim=c(0, 3.5), lty="dashed", main = paste("trunFACD(1,2) fix r = ", r))
    lines(res.den5)
    
    plot(theo.den5~xdata, type="l", col="blue", xlim=c(3, 7), ylim=c(0, 0.1), lty="dashed", main = paste("mean residual = ", round(mean(mfr2fit$res), digits=3)))
    lines(res.den5)
    
    plot(theo.q5, emp.q5/mean(mfr2fit$res), ylim=c(0,12), xlim=c(0,9))
    lines(0:20, 0:20)
    
    acf(mfr2fit$res, ylim=c(-0.03,0.03))
  }, TRUE)
}
dev.off()

save(mfr2fit_1.9, mfr2fit_2, file="mfr2fit0210p1q2.RData")


#===============================================

wbfit <- wbfit0210p1q2
xdata <- seq(from=0, to=9, by=0.01)
# wbpar.str <- paste(round(wbfit$param[1],digits=3),round(wbfit$param[2],digits=3), round(wbfit$param[3],digits=3), round(wbfit$param[4],digits=3), sep=", ")
res.den2 <- density(wbfit$res, adj=1)
r <- wbfit$param[1]
theo.den2 <- dweibull(xdata, shape = r, scale = 1 / gamma(1 + 1 / r))
emp.q2 <- quantile(wbfit$res, probs = seq(0.01, 0.99, 0.01))
theo.q2 <- qweibull(p=seq(0.01, 0.99, 0.01), shape = r, scale = 1 / gamma(1 + 1 / r), lower.tail = TRUE)

r <- 1.9
mfr2fit <- get(paste("mfr2fit_", r, sep=""))
res.den5 <- density(mfr2fit$res, adj=1)
p <- 1 - 1 / r
g <- gamma(p)
q <- 9.21^(-1 / r)  # 9.21 = -log(delta) if delta = 0.0001
bga <- 1 / (g - q)
aga <- -q * bga
theo.den5 <- dfrechet(xdata, loc = aga, scale = bga, shape = r, log = FALSE)
emp.q5 <- quantile(mfr2fit$re, probs = seq(0.01, 0.99, 0.01))
theo.q5 <- qfrechet(p=seq(0.01, 0.99, 0.01), loc = aga, scale = bga, shape = r, lower.tail = TRUE)
# mfr2par.str <- paste(round(mfr2fit$param[2],digits=3), round(mfr2fit$param[3],digits=3), round(mfr2fit$param[4],digits=3), sep=", ")


pdf(file=paste(id.str, '_compare.pdf', sep=""), width = 12, height = 8)

par(mfrow=c(2,3), mar=c(4,4,4,0.5))
plot(theo.den5~xdata, type="l", col="blue", ylim=c(0, 1.5), lty="dashed", xlab="", ylab="", main="")
mtext("Density function", side = 3, line = 2, cex = 0.9, font=2)
mtext("FACD", side = 2, line = 2.5, cex = 0.9, font=2)
lines(res.den5)
plot(theo.den5~xdata, type="l", col="blue", xlim=c(5, 8), ylim=c(0, 0.015), xlab="", ylab="", lty="dashed", main="")
mtext("Density function (right tail)", side = 3, line = 2, cex = 0.9, font=2)
lines(res.den5)
# plot(theo.q5, emp.q5/mean(mfr2fit$res), ylim=c(0,12), xlim=c(0,9))
# lines(0:20, 0:20)
acf(mfr2fit$res, ylim=c(-0.03,0.03), lag.max=19, main="", ann=FALSE)
mtext("Residual ACF", side = 3, line = 2, cex = 0.9, font=2)
mtext("Lag", side = 1, line = 2.5, cex = 0.7)

plot(theo.den2~xdata, type="l", col="blue", ylim=c(0, 1.5), lty="dashed", xlab="", ylab="", main="")
mtext("WACD", side = 2, line = 2.5, cex = 0.9, font=2)
lines(res.den2)
plot(theo.den2~xdata, type="l", col="blue", xlim=c(5, 8), ylim=c(0, 0.015), xlab="", ylab="", lty="dashed", main="")
lines(res.den2)
# plot(theo.q2, emp.q2, ylim=c(0,12), xlim=c(0,9))
# lines(0:20, 0:20)
acf(wbfit$res, ylim=c(-0.03,0.03),lag.max=19, main="", ann=FALSE)
mtext("Lag", side = 1, line = 2.5, cex = 0.7)

dev.off()



