id.str <- "optimtest"
testfit <- FitACDgeneral_optim(dur.adj, id.str, p = 1, q = 1, maxit = 100, init.param = c(2,0.1,0.2,0.6), distrib = "frechet")
mean(testfit$res)

id.str <- "2ndroundnewton"
# 1st round by optim()
frfit <- FitFACD(dur.adj, id.str = id.str, init.param = testfit$param, portmanteau = FALSE)
# 1st round by FACD with fixed gamma = 2
frfit <- FitFACD(dur.adj, id.str = id.str, init.param = fr2fit_2$param, portmanteau = FALSE)
