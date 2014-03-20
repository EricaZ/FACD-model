num.rep <- 1
true.coeff <- c(0.1, 0.2, 0.6)

true.frparam <- c(2, true.coeff)
true.wbparam <- c(1.3, true.coeff)
true.expparam <- true.coeff

num.n <- 1000
xfr.mat <- SimulateACD(param = true.frparam, num.n = num.n, num.rep = num.rep, distrib = "frechet")
write.csv(xfr.mat, file = "dgp1_1k_data.csv")

xwb.mat <- SimulateACD(param = true.wbparam, num.n = num.n, num.rep = num.rep, distrib = "weibull")
write.csv(xwb.mat, file = "wbdata.csv")

xexp.mat <- SimulateACD(param = true.expparam, num.n = num.n, num.rep = num.rep, distrib = "exp")
write.csv(xexp.mat, file = "dgp3_1k_data.csv")

# Estimation

init.frparam <- c(2, 0.1, 0.2, 0.6)
init.param <- c(0.1, 0.2, 0.6)
id.str <- "dgpwb_fix"

# frfit <- FitFACD(dur.adj, id.str = id.str, init.param = init.frparam, portmanteau = FALSE)

rseq <- seq(from = 1.2, to = 2.7, by = 0.03)

for (r in rseq) {
  assign(paste("fr2fit_", r, sep=""), FitFACD_fix(wbdata[,2], id.str = id.str, r = r, init.param = init.param))
}

# save(list=ls(pattern="fr2fit_"), file="fr2fit.RData")

pdf(file=paste(id.str, '_plot.pdf', sep=""), width = 9, height = 5)
xdata <- seq(from=0, to=7, by=0.001)
for (r in rseq) {
  # Output plots
  fr2fit <- get(paste("fr2fit_", r, sep=""))
  res.den4 <- density(fr2fit$res, adj=1)
  theo.den4 <- dfrechet(xdata, loc = 0, scale = 1 / gamma(1 - 1 / r), shape = r, log = FALSE)
  emp.q4 <- quantile(fr2fit$res, probs = seq(0.01, 0.99, 0.01))
  theo.q4 <- qfrechet(p=seq(0.01, 0.99, 0.01), loc=0, scale=1 / gamma(1 - 1 / r), shape = r, lower.tail = TRUE)
  
  par(mfrow=c(1,3), mar=c(5,5,2.5,2.5))
  #   plot(theo.den4~xdata, type="l", col="blue", ylim=c(0, max(c(theo.den4, res.den4$y))), lty="dashed", main = paste("fix gamma = ", r))
  plot(theo.den4~xdata, type="l", col="blue", ylim=c(0, 3.5), lty="dashed", ylab=NA, main = paste("fix gamma = ", r))
  lines(res.den4)
  
  plot(theo.den4~xdata, type="l", col="blue", xlim=c(3, 7), ylim=c(0, 0.1), lty="dashed", ylab=NA)
  lines(res.den4)
  
  plot(theo.q4, emp.q4, ylim=c(0,12), xlim=c(0,9))
  lines(0:20, 0:20)
}
dev.off()


pdf(file=paste(id.str, '_acf.pdf', sep=""), width = 9, height = 5)
par(mfrow=c(1,1), mar=c(5,5,5,2.5))
for (r in rseq) {
  fr2fit <- get(paste("fr2fit_", r, sep=""))
  acf(fr2fit$res, main=paste("fix gamma = ", r), ylim=c(0,0.1))
}
dev.off()