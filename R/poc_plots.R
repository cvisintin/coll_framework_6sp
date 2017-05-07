calibplot <- function(pred, negrug, posrug, ideal, ylim=c(0,1), xlim=c(0,1), capuci=TRUE, xlabel = "Predicted probability of presence", filename=NULL, title="Calibration plot", ...) {
  if (!is.null(filename)) png(filename)
  ylow <- pred$y - 2 * pred$se
  ylow[ylow<0] <- 0
  yhigh <- pred$y + 2 * pred$se
  if (capuci) yhigh[yhigh>1] <- 1
  plot(pred$x, ylow, type="l", col="orange", ylim=ylim, xlim=xlim,
       xlab=xlabel, lwd=2, ...)
  lines(pred$x, yhigh, lwd=2, col="orange")
  lines(pred$x, sapply(pred$x, ideal), lty="dashed")
  points(pred$x, pred$y, col="deepskyblue")
  rug(negrug)
  rug(posrug, col = "orange")
  title(title)
  if (!is.null(filename)) dev.off()
}

smoothingdf <- 6
smoothdist <- function(pred, res) {
  require(splines)
  gam1 <- glm(res ~ ns(pred, df=smoothingdf), weights=rep(1, length(pred)), family=binomial)
  x <- seq(min(pred), max(pred), length = 512)
  y <- predict(gam1, newdata = data.frame(pred = x), se.fit = TRUE,
               type = "response")
  data.frame(x=x, y=y$fit, se=y$se.fit)
}

# presence-only smoothed calibration plot
pocplot <- function(pred, back, linearize=TRUE, ...) {
  ispresence <- c(rep(1,length(pred)), rep(0, length(back)))
  predd <- smoothdist(c(pred,back), ispresence)
  c <- mean(back)*length(back)/length(pred)
  if (linearize) {
    fun <- function(x,y) c*y / (1-y)
    predd$y <- mapply(fun, predd$x, predd$y)
    predd$se <- mapply(fun, predd$x, predd$se)
    ideal <- function(x) x
    ylab <- "Relative probability of presence" 
  } 
  else {
    ideal <- function(x) x / (x + c)
    ylab <- "Probability of presence"
  }
  calibplot(predd, negrug=back, posrug=pred, ideal=ideal, ylab=ylab,
            capuci = FALSE, ...)
  predd
}

# presence-absence smoothed calibration plot
pacplot <- function(pred, pa, ...) {
  predd <- smoothdist(pred, pa)
  calibplot(predd, negrug=pred[pa==0], posrug=pred[pa==1], ideal=function(x) x, ylab="Probability of presence", ...)
}

# binned calibration plot with equal width bins
ecalp <- function(preds, acts, bins=10, do.plot=TRUE, do.clear=TRUE, filename=NULL, title="Binned calibration plot", ...){
  g <- floor(preds*bins)
  b <- 0:(bins-1)
  p <- sapply(b, function(x) if (length(acts[g==x])==0) -1 else sum(acts[g==x]) / length(acts[g==x]))
  mx <- sapply(b, function(x,g) mean(preds[g==x]), g)
  if(do.plot) {
    if (!is.null(filename)) png(filename)
    if (do.clear) {
      plot(mx, p, xlim=c(0,1), ylim=c(0,1), ...)
    } else {
      points(mx, p, xlim=c(0,1), ylim=c(0,1), ...)
    }
    rug(preds[acts==0])
    rug(preds[acts==1], col = "orange")
    abline(0,1,lty="dashed")
    title(title)
    if (!is.null(filename)) dev.off()
  }
  return(p)
}


# Code to recalibrate a model, without changing model discrimination
percent <- 95
rescaleinit <- function(plot) {
  mx <<- max(sapply(plot$y, function(y) y/(1-y)))
  len <<- length(plot$x)
  plotx <<- plot$x
  incr <<- isoreg(plot$y)$yf
}
rescale <- function(pred) {
  i <- findInterval(pred, plotx)
  if (i==0) i=1
  f <- incr[i]
  return((percent * (f / (mx * (1-f))) + (100-percent) * pred)/100)
}


#################Example#########################

# This code makes figures 1, 2 and 3, and calculates associated AUC and COR values

require(ROCR)
set.seed(122333444)

# model functions 
cube <- function(x) x*x*x
model1 <- function(x,y) 0.25 + 0.5 * x
model2 <- function(x,y) cube((x+y)/2)

# grid of x and y values for Figure 1
xg <- seq(0, 1, length=20)
yg <- xg
pg <- outer(xg, yg, function(x,y) (x+y)/2)
mg1 <- outer(xg, yg, model1)
mg2 <- outer(xg, yg, model2)

# Figure 1 a. First needs a "figs" directory to exist within your current working directory.
png("figs/truth.png")
persp(xg, yg, pg,col="deepskyblue", theta=30, phi=20, r=10000, d=0.1, expand=0.5, ltheta=90, lphi=180, shade=0.75, ticktype="detailed", nticks=5, xlab="Temperature", ylab="Precipitation", zlab="Probability", zlim=c(0,1))
title("(a) Truth")
dev.off()

# Figure 1 c
png("figs/pred_m1.png")
persp(xg, yg, mg1, col="deepskyblue", theta=30, phi=20, r=10000, d=0.1, expand=0.5, ltheta=90, lphi=180, shade=0.75, ticktype="detailed", nticks=5, xlab="Temperature", ylab="Precipitation", zlab="Predicted probability", zlim=c(0,1))
title("(c) Model 1")
dev.off()

# Figure 1 d
png("figs/pred_m2.png")
persp(xg, yg, mg2,col="deepskyblue", theta=30, phi=20, r=10000, d=0.1, expand=0.5, ltheta=90, lphi=180, shade=0.75, ticktype="detailed", nticks=5, xlab="Temperature", ylab="Precipitation", zlab="Predicted probability", zlim=c(0,1))
title("(d) Model 2")
dev.off()


# number of samples in data sets
ns <- 1000

# test data, for calculating statistics
xt <- runif(ns, min=0, max=1)
yt <- runif(ns, min=0, max=1)

# true probability of presence
pt <- (xt+yt)/2

# observed presence / absence, randomly drawn according to pt
ot <- rbinom(ns, 1, pt)

# model values
mt1 <- mapply(model1, xt, yt)
mt2 <- mapply(model2, xt, yt)

# Figure 1 b
palette(c("deepskyblue","orange"))
png("figs/testdata.png")
plot(xt,yt,col=ot+1,xlab="Temperature",ylab="Precipitation")
title("(b)")
dev.off()

# Figure 2 a, b, c (binned calibration plots):
ecalp(pt, ot, title="(a) Truth")
ecalp(mt1, ot, title="(b) Model 1")
ecalp(mt2, ot, title="(c) Model 2")


# Figure 2 d, e, f (smoothed calibration plots):
pacplot(pt, ot, title="(d) Truth")
pacplot(mt1, ot, title="(e) Model 1")
pacplot(mt2, ot, title="(f) Model 2")


# Predictions at 10000 random background points:
xb <- runif(10000)
yb <- runif(10000)
pb <- (xb+yb)/2
mb1 <- mapply(model1, xb, yb)
mb2 <- mapply(model2, xb, yb)


# Figure 3 a, b, c, d:
pocplot(mt2[ot==1], mb2, title="(a) Model 2", linearize=FALSE, ylim=c(0,0.2))
pocplot(mt2[ot==1], mb2, title="(b) Model 2")
pocplot(mt1[ot==1], mb1, title="(c) Model 1")
pocplot(pt[ot==1], pb, title="(d) Truth")


# Calculation of AUC and COR for Model 1 and 2 and true prob. of presence
performance(prediction(mt1, ot), "auc")@y.values[[1]] 
performance(prediction(mt2, ot), "auc")@y.values[[1]] 
performance(prediction(pt, ot), "auc")@y.values[[1]] 
cor(mt1, ot)
cor(mt2, ot)
cor(pt, ot)

# Generate fresh data to use while calibrating model 2
xc <- runif(ns, min=0, max=1)
yc <- runif(ns, min=0, max=1)
pc <- (xc+yc)/2
oc <- rbinom(ns, 1, pc)
mc2 <- mapply(model2, xc, yc)

# Recalibrate model 2
plot <- pocplot(mc2[oc==1], mb2, title="", linearize=FALSE)
rescaleinit(plot)
model2b <- function(x,y) {
  m <- model2(x,y)
  return(rescale(m))
}
mt2b <- mapply(model2b, xt, yt)

# Calculate AUC and COR for recalibrated Model 2
performance(prediction(mt2b, ot), "auc")@y.values[[1]] 
cor(mt2b, ot)



# Change the rescale function to recalibrate a model, but without using "isoreg"
rescaleinit.noisoreg <- function(plot) {
  mx <<- max(sapply(plot$y, function(y) y/(1-y)))
  len <<- length(plot$x)
  plotx <<- plot$x
  ploty <<- plot$y
}
rescale.noisoreg <- function(pred) {
  i <- findInterval(pred, plotx)
  if (i==0) i=1
  f <- ploty[i]
  return(f / (mx * (1-f)))
}

# Recalibrate model 2 again, without isoreg (i.e. using changed rescale function), and recalculate AUC and COR 
plot <- pocplot(mc2[oc==1], mb2, title="", linearize=FALSE)
rescaleinit.noisoreg(plot)
model2b <- function(x,y) {
  m <- model2(x,y)
  return(rescale.noisoreg(m))
}
mt2b <- mapply(model2b, xt, yt)
performance(prediction(mt2b, ot), "auc")@y.values[[1]] 
cor(mt2b, ot)
