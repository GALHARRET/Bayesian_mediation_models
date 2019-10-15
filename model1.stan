data {
int<lower=1> N;
vector[N] M;
vector[N] X;
int<lower=0, upper=1> Y[N];
vector [5] mu; // mean of theta=(thetaM,thetaY) obtained by mcmc on the historical data
matrix[5,5] SIG; // variance of theta obtained by mcmc on the historical data
}

transformed data{
matrix [N,2] phiM;
matrix [N,3] phiY;
for (n in 1:N){
phiM[n,1]=1;
phiM[n,2]=X[n];
phiY[n,1]=1;
phiY[n,2]=X[n];
phiY[n,3]=M[n];}
}

parameters {
vector[2] thetaM;
vector[3] thetaY;
vector<lower=0>[5] k;
//real<lower=0> lambda;
real <lower=0> varM;
}



model {
vector[5] MU;
vector[5] theta;
theta[1:2]=thetaM;
theta[3:5]=thetaY;
for(i in 1:5) MU[i]=k[i]*mu[i];
M~normal(phiM*theta[1:2],sqrt(varM));
Y~bernoulli_logit(phiY*theta[3:5]);
target+=multi_normal_lpdf(theta|MU,SIG);
for(i in 1:5) k[i]~normal(1,5);

}
