library(glmnet)
library(mvtnorm)
library(SLOPE)
library(lpSolve)

#######generacja macierzy X


n<-100
L<-200;# length of chromosome in cM;
d<-5;# distance between neigboring markers in cM
r<-0.5*(1-exp(-0.02*d));#probability of recombination between different markers
X<-rep(0,L*n);
X<-matrix(X,nrow=n);
P<-rbinom(n,1,0.5);#genotype in the first marker
R<-rbinom(n*(L-1),1,r);
R<-matrix(R,nrow=n,ncol=(L-1));#matrix of recombination probabilities between neigboring markers
T<-cbind2(P,R);
T<-apply(T,1,'cumsum');
T<-t(T);
X<-T%%2;#final matrix of genotypes on ith chromosome

X<-scale(X);
n=100;
p=200;

#######

Y<-X%*%beta+5*rnorm(n);
#####sygnaly wystepuja na pierwszych 30 pozycjach

######
n=100;
p=200;

load('intro_final_2.Rdata')
###########################cluster example

######IR condition

####cluster example h

q=0.95;
seq1<-seq(1:p);
lambda<-qnorm(1-q*seq1/2/p);
lambda[1]<-lambda[1]+1;
wynsl<-rep(0,50);
wynla<-rep(0,50);

for (k in 2:50)
{
ind<-seq(1:k);
Xind<-X[,ind];
tildeX<-apply(Xind,1,'sum');

tildel<-sum(lambda[1:k]);

mat<-t(X)%*%tildeX*tildel/sum(tildeX^2);
smat<-sort(abs(mat), decreasing=TRUE);
cmat<-cumsum(smat);
cl<-cumsum(lambda);
wynsl[k]=max(cmat/cl);

wynla[k]<-max(abs(t(X)%*%X[,ind]%*%solve(t(X[,ind])%*%X[,ind])%*%rep(1,k)));

}

which(wynsl<=1+1E-10)
which(wynla<=1+1E-10)

###########################


beta<-rep(0,p);
beta[1:30]<-40;

obj2<-cv.glmnet(X,Y);
hatb2<-coef(obj2, s='lambda.min')[2:(p+1)];
plot(hatb2~beta, main='LASSO', ylab='estimator');
sum((hatb2-beta)^2)


lambda<-0.1*obj2$lambda.min;
###1.44337
obj3<-glmnet(X,Y,lambda=lambda, intercept=FALSE, thresh=1E-20);
hatb3<-coef(obj3, s='lambda.min')[2:(p+1)];
plot(hatb3~beta, main='LASSO', ylab='estimator');
sum((hatb3-beta)^2)

lambda<-0.06*obj2$lambda.min;
###1.44337
obj4<-glmnet(X,Y,lambda=lambda, intercept=FALSE, thresh=1E-20);
hatb4<-coef(obj3, s='lambda.min')[2:(p+1)];
plot(hatb4~beta, main='LASSO', ylab='estimator');
sum((hatb4-beta)^2)


q=0.95;
seq1<-seq(1:p);
lambda<-qnorm(1-q*seq1/2/p);
lambda[1]<-lambda[1]+1;


B<-SLOPE(X,Y, lambda=lambda, alpha=10/sqrt(n),intercept=FALSE, scale='none')
betahat=coefficients(B);
plot(betahat~beta, main='SLOPE', ylab='estimator');
sum((betahat-beta)^2)

plot(hatb3~beta, main='Cluster', ylab='estimator', xlab=expression(paste(beta)), pch=19, xlim=c(35,45), ylim=c(35,45),cex=2);
points(betahat~beta, col='red', pch=19, cex=2.5);
abline(h=40,lwd=2)
legend(35,45,c('LASSO','SLOPE'), pch=c(19,19), col=c('black','red'), pt.cex=c(2,2.5))





####checking identifiability condition


success<-rep(0,100)
beta<-rep(0,p);

obj<-rep(1,2*p);
constmat1<-cbind(X,-X);
constmat2<-cbind(diag(p),matrix(rep(0,p*p),nrow=p));
constmat3<-cbind(matrix(rep(0,p*p),nrow=p),diag(p));
constmat<-rbind(constmat1,constmat2,constmat3);
consteq<-c(rep('=',n),rep('>=',2*p));



for (k in 2:100)
{
beta[1:k]<-40;
Yn<-X%*%beta;
constright<-c(Yn,rep(0,2*p));
wyn1=lp('min',obj,constmat,consteq,constright);
wyn=wyn1$solution[1:p]-wyn1$solution[(p+1):(2*p)];
success[k]<-(max(abs(wyn-beta))<1E-7)
}

which(success==1)
###1


beta<-rep(0,p);
beta[1:90]<-40;
Y<-X%*%beta;

obj2<-cv.glmnet(X,Y);
hatb2<-coef(obj2, s='lambda.min')[2:(p+1)];
plot(hatb2~beta, main='LASSO', ylab='estimator');
sum((hatb2-beta)^2)

lambda<-0.005*obj2$lambda.min;
###1.44337
obj3<-glmnet(X,Y,lambda=lambda, intercept=FALSE, thresh=1E-20);
hatb3<-coef(obj3, s='lambda.min')[2:(p+1)];
plot(hatb3~beta, main='LASSO', ylab='estimator');
sum((hatb3-beta)^2)

*************************

beta<-rep(0,p);
beta[1:91]<-40;
Y<-X%*%beta;

obj2<-cv.glmnet(X,Y);
hatb2<-coef(obj2, s='lambda.min')[2:(p+1)];
plot(hatb2~beta, main='LASSO', ylab='estimator');
sum((hatb2-beta)^2)

lambda<-0.005*obj2$lambda.min;
###1.44337
obj3<-glmnet(X,Y,lambda=lambda, intercept=FALSE, thresh=1E-20);
hatb3<-coef(obj3, s='lambda.min')[2:(p+1)];
plot(hatb3~beta, main='LASSO', ylab='estimator');
sum((hatb3-beta)^2)

**************************
q=0.95;
seq1<-seq(1:p);
lambda<-qnorm(1-q*seq1/2/p);
lambda[1]<-lambda[1]+1;


B<-SLOPE(X,Y, lambda=lambda, alpha=0.05/sqrt(n),intercept=FALSE, scale='none')
betahat=coefficients(B);
plot(betahat~beta, main='SLOPE', ylab='estimator');
sum((betahat-beta)^2)

plot(hatb3~beta, main='Cluster', ylab='estimator', xlab=expression(paste(beta)), pch=19, xlim=c(0,45), ylim=c(-5,100),cex=2);
points(betahat~beta, col='red', pch=19, cex=2.5);
abline(h=40,lwd=2)
legend(20,60,c('LASSO','SLOPE'), pch=c(19,19), col=c('black','red'), pt.cex=c(2,2.5))




