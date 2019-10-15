data {
int<lower=1> N;
vector [N] X;
vector [N] M;
int <lower=0,upper=1> Y[N];
}

transformed data{
matrix [N,2] phiM;
matrix[2,2] matM;
vector[N] Xc;
vector[N] Mc;
for (n in 1:N){
phiM[n,1]=1;
phiM[n,2]=X[n];}
matM=inverse(transpose(phiM)*phiM);
Xc=X-mean(X);
Mc=(M-mean(M))/(2*sd(M));
}

parameters {
vector[2] thetaM;
vector[3] thetaY;
real<lower=0> varM;
}

transformed parameters{
real thetaYc[3];
thetaYc[1]=thetaY[1]+thetaY[2]*mean(X)+thetaY[3]*mean(M);
thetaYc[3]=thetaY[3]*2*sd(M);
thetaYc[2]=thetaY[2];
}


model {
vector[2] e;
for (j in 1:2) e[j]=0;
thetaM~multi_normal(e,N*matM*varM);
M~normal(phiM*thetaM,sqrt(varM));
thetaYc[1]~cauchy(0,10);
thetaYc[2]~cauchy(0,2.5);
thetaYc[3]~cauchy(0,2.5);
Y~bernoulli_logit(thetaYc[1]+thetaYc[2]*Xc+thetaYc[3]*Mc);
}
