data {
int<lower=1> N;
vector [N] X;
vector [N] M;
int <lower=0,upper=1> Y[N];
}

transformed data{
matrix [N,2] phiM;
matrix [N,3] phiY;
matrix[2,2] matM;
matrix[3,3] matY;
for (n in 1:N){
phiM[n,1]=1;
phiM[n,2]=X[n];
phiY[n,1]=1;
phiY[n,2]=X[n];
phiY[n,3]=M[n];}
matM=inverse(transpose(phiM)*phiM);
matY=inverse(transpose(phiY)*phiY);
}



parameters {
vector[2] thetaM;
vector[3] thetaY;
real<lower=0> varM;
}


model {
vector[3] e;
for (j in 1:3) e[j]=0;
thetaM~multi_normal(e[1:2],N*matM*varM);
M~normal(phiM*thetaM,sqrt(varM));
thetaY~multi_normal(e,4*N^2*matY);
Y~bernoulli_logit(phiY*thetaY);
}
