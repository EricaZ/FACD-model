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
setwd("C:/Users/yaoz/Documents/FACDapp/blocktrade/20131104data")

init.param <- c(0.1, 0.2, 0.6)
id.str <- "fix_gamma_mFACD_new"
rseq <- seq(from = 1.1, to = 3, by = 0.05)
for (r in rseq) {
  try({
#     assign(paste("fr2fit_", r, sep=""), FitFACD_fix(dur.adj, id.str = id.str, r = r, init.param = init.param, distrib="frechet"))
    assign(paste("mfr2fit_", r, sep=""), FitFACD_fix(dur.adj, id.str = id.str, r = r, portmanteau = TRUE, init.param = init.param, distrib="modifrechet"))
  }, TRUE)
}

# 4-param estimtaion of the modified FACD --> unstable, not converge 
# init.frparam <- c(2, 0.1, 0.2, 0.6)
# id.str <- "mFACD"
# rseq <- seq(from = 1.05, to = 3, by = 0.05)
# for (r in rseq) {
#   try({
#     assign(paste("mfrfit_", r, sep=""), FitFACD(dur.adj, id.str = id.str, portmanteau = TRUE, init.param = c(2, 0.1, 0.2, 0.6), distrib="modifrechet"))
# #     assign(paste("mfr2fit_", r, sep=""), FitFACD_fix(dur.adj, id.str = id.str, r = r, portmanteau = TRUE, init.param = init.param, distrib="modifrechet"))
#   }, TRUE)
# }

pdf(file=paste(id.str, '_plot.pdf', sep=""), width = 16, height = 5)
xdata <- seq(from=0, to=7, by=0.01)
for (r in rseq) {
  # Output plots
  try({
#        fr2fit <- get(paste("fr2fit_", r, sep=""))
#        res.den4 <- density(fr2fit$res, adj=1)
#        theo.den4 <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
#        emp.q4 <- quantile(fr2fit$res, probs = seq(0.01, 0.99, 0.01))
#        theo.q4 <- qfrechet(p=seq(0.01, 0.99, 0.01), loc=0, scale=1 / gamma(1 - 1 / r), shape = r, lower.tail = TRUE)
#        
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
       
#        fr2par.str <- paste(round(fr2fit$param[2],digits=3), round(fr2fit$param[3],digits=3), round(fr2fit$param[4],digits=3), sep=", ")
       mfr2par.str <- paste(round(mfr2fit$param[2],digits=3), round(mfr2fit$param[3],digits=3), round(mfr2fit$param[4],digits=3), sep=", ")
       
       par(mfrow=c(1,4), mar=c(5,5,5,5))
#        plot(theo.den4~xdata, type="l", col="blue", ylim=c(0, 3.5), lty="dashed", main = paste("FACD(1,1)","fix gamma = ", r, "\n", fr2par.str))
#        lines(res.den4)
#        
#        plot(theo.den4~xdata, type="l", col="blue", xlim=c(3, 7), ylim=c(0, 0.1), lty="dashed", main = paste("mean residual = ", round(mean(fr2fit$res), digits=3)))
#        lines(res.den4)
#        
#        plot(theo.q4, emp.q4/mean(fr2fit$res), ylim=c(0,12), xlim=c(0,9))
#        lines(0:20, 0:20)
#        
#        acf(fr2fit$res, ylim=c(-0.03,0.03))
       
       plot(theo.den5~xdata, type="l", col="blue", ylim=c(0, 3.5), lty="dashed", main = paste("trunFACD(1,1)","fix gamma = ", r, "\n", mfr2par.str))
       lines(res.den5)
       
       plot(theo.den5~xdata, type="l", col="blue", xlim=c(3, 7), ylim=c(0, 0.1), lty="dashed", main = paste("mean residual = ", round(mean(mfr2fit$res), digits=3)))
       lines(res.den5)
       
       plot(theo.q5, emp.q5/mean(mfr2fit$res), ylim=c(0,12), xlim=c(0,9))
       lines(0:20, 0:20)
       
       acf(mfr2fit$res, ylim=c(-0.03,0.03))
  }, TRUE)
}
dev.off()

#===================================================

init.wbparam <- c(1, 0.1, 0.2, 0.6)
id.str <- "wb_new"
wbfit <- FitWACD(dur.adj, portmanteau = TRUE, id.str = id.str, init.param = init.wbparam) 

# Output plots
xdata <- seq(from=0, to=7, by=0.01)
wbpar.str <- paste(round(wbfit$param[1],digits=3),round(wbfit$param[2],digits=3), round(wbfit$param[3],digits=3), round(wbfit$param[4],digits=3), sep=", ")

res.den2 <- density(wbfit$res, adj=1)
r <- wbfit$param[1]
theo.den2 <- dweibull(xdata, shape = r, scale = 1 / gamma(1 + 1 / r))
emp.q2 <- quantile(wbfit$res, probs = seq(0.01, 0.99, 0.01))
theo.q2 <- qweibull(p=seq(0.01, 0.99, 0.01), shape = r, scale = 1 / gamma(1 + 1 / r), lower.tail = TRUE)

pdf(file=paste(id.str, '_plot.pdf', sep=""), width = 16, height = 8)
par(mfrow=c(2,4), mar=c(5,5,5,5))
plot(theo.den2~xdata, type="l", col="blue", ylim=c(0, 3.5), lty="dashed", main = paste("WACD(1,1):", wbpar.str))
lines(res.den2)
plot(theo.den2~xdata, type="l", col="blue", xlim=c(3, 7), ylim=c(0, 0.1), lty="dashed", main = paste("mean residual = ", round(mean(wbfit$res), digits=3)))
lines(res.den2)
plot(theo.q2, emp.q2, ylim=c(0,12), xlim=c(0,9))
lines(0:20, 0:20)
acf(wbfit$res, ylim=c(-0.03,0.03))

r <- 1.7
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

mfr2par.str <- paste(round(mfr2fit$param[2],digits=3), round(mfr2fit$param[3],digits=3), round(mfr2fit$param[4],digits=3), sep=", ")

plot(theo.den5~xdata, type="l", col="blue", ylim=c(0, 3.5), lty="dashed", main = paste("trunFACD(1,1)","fix gamma = ", r, "\n", mfr2par.str))
lines(res.den5)

plot(theo.den5~xdata, type="l", col="blue", xlim=c(3, 7), ylim=c(0, 0.1), lty="dashed", main = paste("mean residual = ", round(mean(mfr2fit$res), digits=3)))
lines(res.den5)

plot(theo.q5, emp.q5/mean(mfr2fit$res), ylim=c(0,12), xlim=c(0,9))
lines(0:20, 0:20)

acf(mfr2fit$res, ylim=c(-0.03,0.03))

dev.off()


save(list=ls(pattern="fr2fit_"), wbfit, file="fit_results.RData")

#===================================================

id.str <- "fix_gamma_FACD"
rseq <- seq(from = 1.1, to = 3, by = 0.05)
pdf(file=paste(id.str, '_compare_res.pdf', sep=""), width = 8, height = 8)
xdata <- seq(from=-0.5, to=7, by=0.01)
# par(mfrow=c(1,2), mar=c(5,5,5,5))
for (r in rseq) {
  # Output plots
  #   try({
  fr2fit <- get(paste("fr2fit_", r, sep=""))
  theo.den4 <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
  res.den4 <- density(fr2fit$res, adj=1)
  
  mfr2fit <- get(paste("mfr2fit_", r, sep=""))
  res.den5 <- density(mfr2fit$res, adj=1)
  p <- 1 - 1 / r
  g <- gamma(p)
  q <- 9.21^(-1 / r)  # 9.21 = -log(delta) if delta = 0.0001
  bga <- 1 / (g - q)
  aga <- -q * bga
  theo.den5 <- dfrechet(xdata, loc = aga, scale = bga, shape = r, log = FALSE)
  
  fr2par.str <- paste(round(fr2fit$param[2],digits=3), round(fr2fit$param[3],digits=3), round(fr2fit$param[4],digits=3), sep=", ")
  mfr2par.str <- paste(round(mfr2fit$param[2],digits=3), round(mfr2fit$param[3],digits=3), round(mfr2fit$param[4],digits=3), sep=", ")
  
  plot(res.den5, col="blue", ylim=c(0,5.5), xlim=c(0,5), main = paste("FACD(1,1)","fix gamma = ", r, "\n", "frechet_param: ", fr2par.str, ", avg(res) = ", round(mean(fr2fit$res),digits=2), 
                                                                      "\n", "trunfre_param: ", mfr2par.str, ", avg(res) = ", round(mean(mfr2fit$res),digits=2)))
  lines(res.den4)
  lines(theo.den5~xdata, col="blue", lty = "dashed")
  lines(theo.den4~xdata, lty = "dashed")
  lines(rep(mean(fr2fit$res),2), c(5,0))
  lines(rep(mean(mfr2fit$res),2), c(5,0), col="blue")
  
  points(3:5, 5*fr2fit$param[2:4])
  points(3:5, 5*mfr2fit$param[2:4], pch=19, col="blue")
  #   }, TRUE)
}
dev.off()



