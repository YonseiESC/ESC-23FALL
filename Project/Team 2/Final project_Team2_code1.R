###Generating Data p=100, 300 500, 1000 /n=300
pZ<-1000;     N=pZ             #dimension
n<-300;      T=n              #sample size 
ncycle<-1000                #circle times


## Model 3
beta<-matrix(rnorm(pZ*r,0,1),pZ,r); diagc=diag(rep(sqrt(1*r),pZ))*6

#####   Gaussian
fm<-matrix(rnorm(r*n,0,1),r,n); epsilon<-matrix(rnorm(n*pZ,0,1),pZ,n)

Data<-beta%*%fm+diagc%*%epsilon
X=Data-matrix(rep(t(apply(Data,1,mean)),n),ncol=n)
X <- t(X)
sn<-cor(X)
sn
N<-nrow(X)
P<-ncol(X)
R=99

#x<-sn
eigenvalues <- eigen(sn)$values
eigenvalues<-sort(eigenvalues, decreasing=TRUE)
nge <- rep(0, P) 
# Repeat R times
for (r in 1:R) {
  for (j in 2:P) {
    X[,j] <- X[sample(N),j] # Permute columns 2, ..., P of x
  }
  permut_sn<-cor(X)
  permuted_eigenvalues <- eigen(permut_sn)$values # Obtain the eigenvalues of the permuted data
  permuted_eigenvalues <- sort(permuted_eigenvalues, decreasing=TRUE)
  nge <- nge + (permuted_eigenvalues >= eigenvalues) # Increase the count if the absolute value of the i-th permuted eigenvalue is greater than or equal the absolute value of the observed i-th eigenvalue
}
# Finally
p_values <- (nge + 1) / (R + 1) # Obtain p-values for each eigenvalue
p_values

sig_level=0.05
count=0
for (i in 1:pZ){
  if (p_values[i]<sig_level){
    count=count+1
  }
}
count ; count/pZ

