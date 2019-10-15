Bayesian models for mediation analysis
================
JM GALHARRET
10/13/2019

# Introduction

This R Markdown document provides the code to compute \(NDE(x)\) and
\(NIE(x)\) (\(x \in \{0,1\}\)) with the three models define in the
article Priors comparison in
[hal-02070053](https://hal.archives-ouvertes.fr/hal-02070053).

Download in this page the following files containing the three models
defined in Stan \[model1.stan, model2.stan, model3.stan. Copy and paste
this code into your \*.R file

The causal\_effects function provide the MCMC drawing from the posterior
distribution of NDE(0),NDE(1),NIE(0),NIE(1).

``` r
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
expit<-function(x) 1/(1+exp(-x))
causal_effects<-function(model,data){
  fit <- rstan::stan(file = model, data = data,control = list(adapt_delta = 0.99))
  draws<-rstan::extract(fit)
  effects<-function(thetaM,sigM,thetaY){
    B<-length(sigM)
    f_m<-function(x,m,a,sig) dnorm(m,a[1]+a[2]*x,sig)
    P_y<-function(x,m,b) 1/(1+exp(-(b[1]+b[2]*x+b[3]*m)))
    f1<-function(m,x,a,sig,b){f_m(x,m,a,sig) *(P_y(1,m,b)-P_y(0,m,b))}
    f2<-function(m,x,a,sig,b){(f_m(1,m,a,sig)-f_m(0,m,a,sig))*P_y(x,m,b)}
    NDE0<-rep(0,B)
    NDE1<-rep(0,B)
    NIE0<-rep(0,B)
    NIE1<-rep(0,B)
    for(t in 1:B) {
      NDE0[t]<-integrate(f1,lower=-Inf,upper=Inf,x=0,a=thetaM[t,],sig=sigM[t],
                         b=thetaY[t,])$value
      NDE1[t]<-integrate(f1,lower=-Inf,upper=Inf,x=1,a=thetaM[t,],sig=sigM[t],
                         b=thetaY[t,])$value
      NIE0[t]<-integrate(f2,lower=-Inf,upper=Inf,x=0,a=thetaM[t,],sig=sigM[t],
                         b=thetaY[t,])$value
      NIE1[t]<-integrate(f2,lower=-Inf,upper=Inf,x=1,a=thetaM[t,],sig=sigM[t],
                         b=thetaY[t,])$value
    }
    return(cbind(NDE0,NDE1,NIE0,NIE1))
  }

return(data.frame(effects(thetaM=draws$thetaM,sigM = sqrt(draws$varM),thetaY=draws$thetaY)))
}

summary_effects<-function(res){
  mean<-apply(res,2,mean)
  sd<-apply(res,2,sd)
  ci_lower<-apply(res,2, quantile,probs=.025)
  ci_upper<-apply(res,2, quantile,probs=.975)
  return(cbind(mean,sd,ci_lower,ci_upper))
}
```

# Example with Model 3:

You can run this code on your computer.

Remark: setwd() defines the directory in which the files
model1.stan,model2.stan, model3.stan are located.

``` r
setwd("~/Models/")

N<-50
X<-rep(0:1,N/2)
M<-rnorm(N,1-2*X,.75)
Y<-rbinom(N,1,expit(-0.5+1*M+1.5*X))
res<-causal_effects('model3.stan',data=list(N=N,X=X,M=M,Y=Y))
summary_effects(res)
```

You can plot the results:

``` r
library(tidyr)
summary_effects(res)
df<-gather(res,key="effects",value="posterior")
ggplot(df, aes(x=posterior,fill=effects))+geom_density(alpha=0.4)+theme_minimal()
```

# Example with Model 1:

For Model 1, we need historical data on which we calculate m and SIG

``` r
Xh<-rep(0:1,N/2)
Mh<-rnorm(N,1-2*Xh,.75)
Yh<-rbinom(N,1,expit(-0.5+1*Mh+1.5*Xh))
fit <- stan(file = 'model3.stan', data = list(N=N,X=Xh,M=Mh,Y=Yh))
draws<-extract(fit)

mu<-apply(cbind(draws$thetaM,draws$thetaY),2,mean)
SIG<-cov(cbind(draws$thetaM,draws$thetaY))
```

These values allow us to compute the effects:

``` r
N<-50
X<-rep(0:1,N/2)
M<-rnorm(N,1-2*X,.75)
Y<-rbinom(N,1,expit(-0.5+1*M+1.5*X))
causal_effects('model1.stan',data=list(N=N,X=X,M=M,Y=Y,mu=mu,SIG=SIG))
```
