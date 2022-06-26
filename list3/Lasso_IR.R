n<-100;
p<-300;
k<-5;

library(glmnet);

beta<-rep(0,p);
beta[1:k]<-3.5; 
     
X<-rnorm(n*p,0,1/sqrt(n));
X<-matrix(X,nrow=n);
Xsup<-X[,1:k];

IR=max(t(X)%*%X[,1:k]%*%solve(t(Xsup)%*%Xsup)%*%rep(1,k));

Y=X%*%beta;

obj<-cv.glmnet(X,Y);
betacv<-coefficients(obj, s='lambda.min')[2:(p+1),1];
plot(betacv~beta);

Y=X%*%beta+0.1*rnorm(n);

obj<-cv.glmnet(X,Y);
betacv<-coefficients(obj, s='lambda.min')[2:(p+1),1];
plot(betacv~beta);

lambda<-5*obj$lambda.min;
obj<-glmnet(X,Y,lambda=lambda, intercept=FALSE);
betacv<-coefficients(obj)[2:(p+1)];
plot(betacv~beta);


k<-20;


beta<-rep(0,p);
beta[1:k]<-3.5; 
     
X<-rnorm(n*p,0,1/sqrt(n));
X<-matrix(X,nrow=n);
Xsup<-X[,1:k];

IR=max(t(X)%*%X[,1:k]%*%solve(t(Xsup)%*%Xsup)%*%rep(1,k));

Y=X%*%beta;

obj<-cv.glmnet(X,Y);
betacv<-coefficients(obj, s='lambda.min')[2:(p+1),1];
plot(betacv~beta);

lambda<-10*obj$lambda.min;
obj2<-glmnet(X,Y,lambda=lambda);
betahat<-coefficients(obj2)[2:(p+1)];
plot(betahat~beta)


library('ADMM')

betahat<-admm.bp(X,Y);

