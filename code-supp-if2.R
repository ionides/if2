##
##  R code to generate Figure S1 of 
##  "A new iterated filtering algorithm"
##  by E. L. Ionides, D. Nguyen, Y. Atchade, S. Stoev and A. A. King
##

require(pomp)

CLUSTER <- FALSE

if(CLUSTER){
 CORES <- 100    ## number of parallel processes
 JOBS <- 100     ## number of estimation runs
 NPMCMC <- 20000     ## number of PMCMC iterations per estimation run
 n <- 500       # number of time points for computing variances
 NP <- 1500     ## number of particles in pfilter operations
} else {
 CORES <- 12    ## number of parallel processes
 JOBS <- 12   ## number of estimation runs
 NPMCMC <- 1000    ## number of PMCMC iterations per estimation run
 n <-  100      # number of time points for computing variances
 NP <- 150     ## number of particles in pfilter operations
}

########### PMCMC for cholera model #############33

pompExample(dacca)

param.tab <- as.matrix(read.table(row.names=1,header=TRUE,text="
                  mle1 box_min box_max
gamma      20.800000000   10.00   40.00
eps        19.100000000    0.20   30.00
rho         0.000000000    0.00    0.00
delta       0.020000000    0.02    0.02
deltaI      0.060000000    0.03    0.60
clin        1.000000000    1.00    1.00
alpha       1.000000000    1.00    1.00
beta.trend -0.004980000   -0.01    0.00
log.beta1   0.747000000   -4.00    4.00
log.beta2   6.380000000    0.00    8.00
log.beta3  -3.440000000   -4.00    4.00
log.beta4   4.230000000    0.00    8.00
log.beta5   3.330000000    0.00    8.00
log.beta6   4.550000000    0.00    8.00
log.omega1 -1.692819521  -10.00    0.00
log.omega2 -2.543383579  -10.00    0.00
log.omega3 -2.840439389  -10.00    0.00
log.omega4 -4.691817993  -10.00    0.00
log.omega5 -8.477972478  -10.00    0.00
log.omega6 -4.390058806  -10.00    0.00
sd.beta     3.130000000    1.00    5.00
tau         0.230000000    0.10    0.50
S.0         0.621000000    0.00    1.00
I.0         0.378000000    0.00    1.00
Rs.0        0.000000000    0.00    0.00
R1.0        0.000843000    0.00    1.00
R2.0        0.000972000    0.00    1.00
R3.0        0.000000116    0.00    1.00
nbasis      6.000000000    6.00    6.00
nrstage     3.000000000    3.00    3.00
"))

dacca.pars <- c("gamma","eps","deltaI","beta.trend","log.beta1","log.beta2", "log.beta3","log.beta4", "log.beta5", "log.beta6", "log.omega1", "log.omega2","log.omega3","log.omega4","log.omega5","log.omega6", "sd.beta",   "tau")
dacca.ivps <- c("S.0","I.0","R1.0","R2.0","R3.0")
dacca.params <- c(dacca.pars,dacca.ivps)
dacca.rw.sd.pmcmc <- (param.tab[,"box_max"]-param.tab[,"box_min"])/100

dacca.hyperparams <- list(min=param.tab[,"box_min"],
                          max=param.tab[,"box_max"])

dacca.rprior <- function(params, hyperparams, ...)
{
  r <- runif(
             n=length(hyperparams$min),
             min=hyperparams$min,
             max=hyperparams$max
             )
  names(r) <- names(hyperparams$min)
  return(r)
}

dacca.dprior <- function(params, hyperparams, ..., log = FALSE)
{
  estpars <- (hyperparams$max-hyperparams$min)>0
  d <- sum(
           dunif(
                 x=params[estpars],
                 min=hyperparams$min[estpars],
                 max=hyperparams$max[estpars],
                 log=TRUE
                 )
           )
  if (log) d else exp(d)
}

dacca <- pomp(dacca,
              rprior=dacca.rprior,
              dprior=dacca.dprior,
              hyperparams=dacca.hyperparams)




binary.file <- "pmcmc-cholera.rda"
if(file.exists(binary.file)) load(binary.file) else {

 if (CLUSTER) {
   require(doSNOW)
   nodefile <- Sys.getenv("PBS_NODEFILE")
   hostlist <- tail(scan(nodefile,what=character(0)),-1)
   cl <- makeCluster(hostlist,type='SOCK')
   registerDoSNOW(cl)
 } else {
   require(doParallel)
   registerDoParallel(CORES)
 }

 tic <- Sys.time()

 pmcmc.overdispersed <- foreach(
                i=1:JOBS,
                .packages='pomp',
                .inorder=FALSE
                ) %dopar% {
                  
 set.seed(7777+i)
 th.draw <- rprior(dacca,coef(dacca))[,1]
 m <- try(
          pmcmc(
              dacca,
              Nmcmc=NPMCMC,
              Np=NP,
              max.fail=1000,
              start=th.draw,
              pars=c(dacca.pars,dacca.ivps),
              rw.sd=dacca.rw.sd.pmcmc
              )
          )
 list(pomp=m,start=th.draw)
 }

 pmcmc.underdispersed <- foreach(
                i=1:JOBS,
                .packages='pomp',
                .inorder=FALSE
                ) %dopar% {
                  
 set.seed(7777+i)
 th.draw <- coef(dacca)
 m <- try(
          pmcmc(
              dacca,
              Nmcmc=NPMCMC,
              Np=NP,
              max.fail=1000,
              start=th.draw,
              pars=c(dacca.pars,dacca.ivps),
              rw.sd=dacca.rw.sd.pmcmc
              )
          )
 list(pomp=m,start=th.draw)
 }

 toc <- Sys.time()
 pmcmc.time <- toc-tic

 pmcmc.out.overdispersed <- lapply(pmcmc.overdispersed,function(x)conv.rec(x$pomp))
 pmcmc.out.underdispersed <- lapply(pmcmc.underdispersed,function(x)conv.rec(x$pomp))

 save(pmcmc.out.overdispersed, pmcmc.out.underdispersed, pmcmc.time, file=binary.file,compress='xz')

 if (CLUSTER) stopCluster(cl)

}


range <- param.tab[,"box_max"]-param.tab[,"box_min"]
rw.scale <- range[dacca.params]/100 # sd of random walk for pmcmc proposal

binary.file <- "traceArray.rda"
if(file.exists(binary.file)) load(binary.file) else {
 traces.od <- lapply(pmcmc.out.overdispersed,function(x)x[,dacca.params])
 traces.ud <- lapply(pmcmc.out.underdispersed,function(x)x[,dacca.params])

 dt <- floor(dim(traces.od[[1]])[1] / n)
 trace.times <- seq(from=1,to=dim(traces.od[[1]])[1],by=dt)
 n <- length(trace.times)
 np <- length(dacca.params)

 trace.array.od <- array(0,c(length(traces.od),n,np))
 for(i in 1:length(traces.od)){
  trace.array.od[i,,] <- traces.od[[i]][trace.times,]
 }

 trace.array.ud <- array(0,c(length(traces.ud),n,np))
 for(i in 1:length(traces.ud)){
  trace.array.ud[i,,] <- traces.ud[[i]][trace.times,]
 }

 
 o.var.mat <- apply(trace.array.od,c(2,3),var)
 o.var.vec <- apply(o.var.mat,1,function(x)sum(x/(rw.scale^2)))

 u.var.mat <- apply(trace.array.ud,c(2,3),var)
 u.var.vec <- apply(u.var.mat,1,function(x)sum(x/(rw.scale^2)))


 save(trace.array.od,trace.array.ud,trace.times, o.var.vec, u.var.vec,
       file=binary.file,compress='xz')
}

pdf(file="pmcmc-varplot.pdf",height=5,width=8)
 par(mfrow=c(1,2))

 plot(y=u.var.vec,x=trace.times,ty="l",xlab="pmcmc iteration, m",
    ylab=expression(paste("weighted sum of variances, ", V[m])))
 plot.window(c(0,1),c(0,1))
 text(0.4,0.92,"A",cex=1.8)
 plot(y=o.var.vec,x=trace.times,ty="l",xlab="pmcmc iteration, m",
    ylab=expression(paste("weighted sum of variances, ", V[m])),
    ylim=c(0,max(o.var.vec)))
 lines(y=u.var.vec,x=trace.times,lty="dotted",lwd=2)
 plot.window(c(0,1),c(0,1))
 text(0.4,0.92,"B",cex=1.8)

dev.off()


##### Liu-West algorithm applied to the toy model #########

##### uses bsmc2() from pomp_0.54-1 #######################

theta.true <- c(th1=1,th2=1,x1.0=0,x2.0=0)
sd.y <- c(10,1)

toy.proc.sim <- function (x, t, params, ...) {
 th1 <- params["th1"]
 th2 <- params["th2"]
 xx <- c(exp(th1),exp(th1)*th2)
 names(xx) <- c("x1","x2")
 return(xx)
}

toy.meas.sim <- function (x, t, params, ...) {
 yy<- rnorm(n=2,mean=x,sd=sd.y)
 c(y1=yy[1],y2=yy[2])
}

toy.meas.dens <- function (y, x, t, params, log, ...) {
f <- dnorm(x=y,mean=x,sd=sd.y,log=log)
if(log) sum(f) else prod(f)
}

Times=100
dat=matrix(1,nrow=2,ncol=Times,dimnames=list(c("y1","y2"),NULL))

toy <- pomp(data=dat,times=1:Times, t0=0,
  rprocess=discrete.time.sim(step.fun=toy.proc.sim,delta.t=1),
  rmeasure=toy.meas.sim,dmeasure=toy.meas.dens)


loglik.exact <- function(po,params){
 th1 <- params["th1"]
 th2 <- params["th2"]
 xx <-  c(x1=exp(th1),x2=exp(th1)*th2)
 sum(apply(obs(po),2,toy.meas.dens,x=xx,t=0,params=params,log=T))
}


set.seed(555)
po <- simulate(toy,params=theta.true) 
mean.dat <- apply(obs(po),1, mean)
mle.exact <- c(th1=unname(log(mean.dat[1])),th2=unname(mean.dat[2]/mean.dat[1]))
max.exact <- loglik.exact(po,mle.exact)

L<-loglik.exact(po,coef(po))
L1 <- dnorm(obs(po)[1,],mean=exp(coef(po)["th1"]),sd=sd.y[1],log=T)
L2 <- dnorm(obs(po)[2,],mean=coef(po)["th2"]*exp(coef(po)["th1"]),sd=sd.y[2],log=T)

N1 <- 200
N2 <- 200
th1.vec <- seq(from=-2,to=2,length=N1)
th2.vec <- seq(from=0,to=10,length=N2)
lik.array <- matrix(NA,nr=N1,nc=N2)
for(i1 in 1:N1){
 for(i2 in 1:N2){
   th <- c(th1=th1.vec[i1],th2=th2.vec[i2])
   lik.array[i1,i2] <- loglik.exact(po,th)
 }
}

require(doParallel)
registerDoParallel(cores=9)
load('saved-toy.rda')
 th.min <- c(th1=-2,th2=0)
 th.max <- c(th1=2,th2=10)
toy.hyperparams <- list(min=th.min,max=th.max)
toy.dprior <- function(params, hyperparams, ..., log)
{
 f <- sum(dunif(params,
 min=hyperparams$min,
 max=hyperparams$max,
 log=TRUE))
 if (log) f else exp(f)
}

toy.rprior <- function(params, hyperparams,...)
{
 r<-c(runif(2,min=hyperparams$min,max=hyperparams$max),0,0)
 names(r) <- c("th1","th2", "x1.0", "x2.0" )
 return(r)
}

po <- pomp(po,dprior=toy.dprior,rprior=toy.rprior,hyperparams=toy.hyperparams)

delta <- 0.9999
NP <- 10000
BW <- 0.1

lwfit <- function(delta,np){
 a <- (3*delta-1)/(2*delta)
 h <- sqrt(1-a^2)
 # h = 0.1004 for delta = 0.99
 mpar <- foreach(i=1:8,.packages='pomp') %dopar% {
  set.seed(567+i)
  m <- bsmc2(po,Np=NP,est=c("th1","th2"),smooth=h,max.fail=1000)
  km <- density(m@post[1,],bw=BW)
  km
 }
 mpar
}

lw1 <- lwfit(0.99,np=10^4)
lw2 <- lwfit(0.999,np=10^4)
lw3 <- lwfit(0.9999,np=10^4)

XLIM <- c(-2,2)
YLIM <- c(0,1)
LINE.XAXIS <- 2.75
X.LABEL <- 0.08
Y.LABEL <- 0.92	

CEX.LAB <- 1.8
CEX.AXIS <- 1
CEX.AXIS.NUMBERS <- 0.8

k1x <- sapply(lw1,getElement,"x")
k1y <- sapply(lw1,getElement,"y")
k2x <- sapply(lw2,getElement,"x")
k2y <- sapply(lw2,getElement,"y")
k3x <- sapply(lw3,getElement,"x")
k3y <- sapply(lw3,getElement,"y")

lik.dev <- lik.array - max(lik.array)
post.array <- exp(lik.dev)/sum(exp(lik.dev))
post.th1 <- apply(post.array,1,sum)
th1.post <- ksmooth(y=post.th1/(th1.vec[2]-th1.vec[1]),x=th1.vec,bandwidth=0.2)


pdf(file="lw.pdf",width=6.5,height=3.5)
op <- par(mfrow=c(1,3),mai=c(0,0.1,0.1,0),omi=c(0.7,0,0,0.1))
par(cex=CEX.AXIS.NUMBERS)

######## delta = 0.99 ###########
matplot(k1x[,1:8],k1y[,1:8],lty=1,type='l',main="",axes=F,xlab='',ylab='',xlim=XLIM,ylim=YLIM)
box()
axis(side=1,cex=CEX.AXIS.NUMBERS)
mtext(side=1,bquote(theta[1]),line=LINE.XAXIS,cex=CEX.AXIS)
lines(th1.post,lty="dotted",lwd=2)
abline(h=0)
plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"A",cex=CEX.LAB)

######## delta = 0.999 ###########
matplot(k2x[,1:8],k2y[,1:8],lty=1,type='l',main="",axes=F,xlab='',ylab='',xlim=XLIM,ylim=YLIM)
box()
axis(side=1,cex=CEX.AXIS.NUMBERS)
mtext(side=1,bquote(theta[1]),line=LINE.XAXIS,cex=CEX.AXIS)
lines(th1.post,lty="dotted",lwd=2)
abline(h=0)
plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"B",cex=CEX.LAB)

######## delta = 0.9999 ###########
matplot(k3x[,1:8],k3y[,1:8],lty=1,type='l',main="",axes=F,xlab='',ylab='',xlim=XLIM,ylim=YLIM)
box()
axis(side=1,cex=CEX.AXIS.NUMBERS)
mtext(side=1,bquote(theta[1]),line=LINE.XAXIS,cex=CEX.AXIS)
lines(th1.post,lty="dotted",lwd=2)
abline(h=0)
plot.window(c(0,1),c(0,1))
text(x=X.LABEL,y=Y.LABEL,"C",cex=CEX.LAB)

dev.off()




