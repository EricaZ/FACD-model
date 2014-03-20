library(evd)
# Output plots
xdata <- seq(from=0, to=9, by=0.01)

wbfit <- wbfit1021
wbpar.str1 <- paste(round(wbfit$param[1],digits=3),round(wbfit$param[2],digits=3), round(wbfit$param[3],digits=3), round(wbfit$param[4],digits=3), sep=", ")
res.den1 <- density(wbfit$res, adj=1)
r <- wbfit$param[1]
theo.den1 <- dweibull(xdata, shape = r, scale = 1 / gamma(1 + 1 / r))
emp.q1 <- quantile(wbfit$res, probs = seq(0.01, 0.99, 0.01))
theo.q1 <- qweibull(p=seq(0.01, 0.99, 0.01), shape = r, scale = 1 / gamma(1 + 1 / r), lower.tail = TRUE)
wbfit1 <- wbfit

wbfit <- wbfit1104
wbpar.str2 <- paste(round(wbfit$param[1],digits=3),round(wbfit$param[2],digits=3), round(wbfit$param[3],digits=3), round(wbfit$param[4],digits=3), sep=", ")
res.den2 <- density(wbfit$res, adj=1)
r <- wbfit$param[1]
theo.den2 <- dweibull(xdata, shape = r, scale = 1 / gamma(1 + 1 / r))
emp.q2 <- quantile(wbfit$res, probs = seq(0.01, 0.99, 0.01))
theo.q2 <- qweibull(p=seq(0.01, 0.99, 0.01), shape = r, scale = 1 / gamma(1 + 1 / r), lower.tail = TRUE)
wbfit2 <- wbfit

wbfit <- wbfit1202
wbpar.str3 <- paste(round(wbfit$param[1],digits=3),round(wbfit$param[2],digits=3), round(wbfit$param[3],digits=3), round(wbfit$param[4],digits=3), sep=", ")
res.den3 <- density(wbfit$res, adj=1)
r <- wbfit$param[1]
theo.den3 <- dweibull(xdata, shape = r, scale = 1 / gamma(1 + 1 / r))
emp.q3 <- quantile(wbfit$res, probs = seq(0.01, 0.99, 0.01))
theo.q3 <- qweibull(p=seq(0.01, 0.99, 0.01), shape = r, scale = 1 / gamma(1 + 1 / r), lower.tail = TRUE)
wbfit3 <- wbfit

wbfit <- wbfit1209
wbpar.str4 <- paste(round(wbfit$param[1],digits=3),round(wbfit$param[2],digits=3), round(wbfit$param[3],digits=3), round(wbfit$param[4],digits=3), sep=", ")
res.den4 <- density(wbfit$res, adj=1)
r <- wbfit$param[1]
theo.den4 <- dweibull(xdata, shape = r, scale = 1 / gamma(1 + 1 / r))
emp.q4 <- quantile(wbfit$res, probs = seq(0.01, 0.99, 0.01))
theo.q4 <- qweibull(p=seq(0.01, 0.99, 0.01), shape = r, scale = 1 / gamma(1 + 1 / r), lower.tail = TRUE)
wbfit4 <- wbfit

wbfit <- wbfit0106
wbpar.str5 <- paste(round(wbfit$param[1],digits=3),round(wbfit$param[2],digits=3), round(wbfit$param[3],digits=3), round(wbfit$param[4],digits=3), sep=", ")
res.den5 <- density(wbfit$res, adj=1)
r <- wbfit$param[1]
theo.den5 <- dweibull(xdata, shape = r, scale = 1 / gamma(1 + 1 / r))
emp.q5 <- quantile(wbfit$res, probs = seq(0.01, 0.99, 0.01))
theo.q5 <- qweibull(p=seq(0.01, 0.99, 0.01), shape = r, scale = 1 / gamma(1 + 1 / r), lower.tail = TRUE)
wbfit5 <- wbfit

wbfit <- wbfit0303
wbpar.str6 <- paste(round(wbfit$param[1],digits=3),round(wbfit$param[2],digits=3), round(wbfit$param[3],digits=3), round(wbfit$param[4],digits=3), sep=", ")
res.den6 <- density(wbfit$res, adj=1)
r <- wbfit$param[1]
theo.den6 <- dweibull(xdata, shape = r, scale = 1 / gamma(1 + 1 / r))
emp.q6 <- quantile(wbfit$res, probs = seq(0.01, 0.99, 0.01))
theo.q6 <- qweibull(p=seq(0.01, 0.99, 0.01), shape = r, scale = 1 / gamma(1 + 1 / r), lower.tail = TRUE)
wbfit6 <- wbfit

# pdf(file=paste('sixwb_plot.pdf', sep=""), width = 12, height = 8)
par(mfrow=c(2,3), mar=c(5,5,4,2))
plot(theo.q1, emp.q1, ylim=c(0,6), xlim=c(0,5), main="Week 1", xlab="TheoQ", ylab="EmpQ")
lines(0:20, 0:20)
plot(theo.q2, emp.q2, ylim=c(0,6), xlim=c(0,5), main="Week 2", xlab="TheoQ", ylab="EmpQ")
lines(0:20, 0:20)
plot(theo.q3, emp.q3, ylim=c(0,6), xlim=c(0,5), main="Week 3", xlab="TheoQ", ylab="EmpQ")
lines(0:20, 0:20)
plot(theo.q4, emp.q4, ylim=c(0,6), xlim=c(0,5), main="Week 4", xlab="TheoQ", ylab="EmpQ")
lines(0:20, 0:20)
plot(theo.q5, emp.q5, ylim=c(0,6), xlim=c(0,5), main="Week 5", xlab="TheoQ", ylab="EmpQ")
lines(0:20, 0:20)
plot(theo.q6, emp.q6, ylim=c(0,6), xlim=c(0,5), main="Week 6", xlab="TheoQ", ylab="EmpQ")
lines(0:20, 0:20)
# dev.off()



# ========================================================
pdf(file=paste('sixwb_details_plot.pdf', sep=""), width = 12, height = 12)
par(mfrow=c(3,3), mar=c(5,5,5,5))
plot(theo.den1~xdata, type="l", col="blue", ylim=c(0, 1.5), lty="dashed", main = paste("WACD(1,1):", wbpar.str1))
lines(res.den1)
plot(theo.den1~xdata, type="l", col="blue", xlim=c(5, 8), ylim=c(0, 0.02), lty="dashed", main = "")
lines(res.den1)
acf(wbfit1$res, ylim=c(-0.03,0.03))

plot(theo.den2~xdata, type="l", col="blue", ylim=c(0, 1.5), lty="dashed", main = paste("WACD(1,1):", wbpar.str2))
lines(res.den2)
plot(theo.den2~xdata, type="l", col="blue", xlim=c(5, 8), ylim=c(0, 0.02), lty="dashed", main = "")
lines(res.den2)
acf(wbfit2$res, ylim=c(-0.03,0.03))

plot(theo.den3~xdata, type="l", col="blue", ylim=c(0, 1.5), lty="dashed", main = paste("WACD(1,1):", wbpar.str3))
lines(res.den3)
plot(theo.den3~xdata, type="l", col="blue", xlim=c(5, 8), ylim=c(0, 0.02), lty="dashed", main = "")
lines(res.den3)
acf(wbfit3$res, ylim=c(-0.03,0.03))

plot(theo.den4~xdata, type="l", col="blue", ylim=c(0, 1.5), lty="dashed", main = paste("WACD(1,1):", wbpar.str4))
lines(res.den4)
plot(theo.den4~xdata, type="l", col="blue", xlim=c(5, 8), ylim=c(0, 0.02), lty="dashed", main = "")
lines(res.den4)
acf(wbfit4$res, ylim=c(-0.03,0.03))

plot(theo.den5~xdata, type="l", col="blue", ylim=c(0, 1.5), lty="dashed", main = paste("WACD(1,1):", wbpar.str5))
lines(res.den5)
plot(theo.den5~xdata, type="l", col="blue", xlim=c(5, 8), ylim=c(0, 0.02), lty="dashed", main = "")
lines(res.den5)
acf(wbfit5$res, ylim=c(-0.03,0.03))

plot(theo.den6~xdata, type="l", col="blue", ylim=c(0, 1.5), lty="dashed", main = paste("WACD(1,1):", wbpar.str6))
lines(res.den6)
plot(theo.den6~xdata, type="l", col="blue", xlim=c(5, 8), ylim=c(0, 0.02), lty="dashed", main = "")
lines(res.den6)
acf(wbfit6$res, ylim=c(-0.03,0.03))
dev.off()


