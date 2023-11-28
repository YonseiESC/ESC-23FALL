set.seed(2023) 
#set.seed(2022)
#set.seed(2021)
p<-50            # number of variables
# p=100, 300, 500, 1000
n<-300           # number of observations
f=5             # number of latent factor

## Model 2
beta<-matrix(rnorm(p*f,0,1),p,f); c<-runif(p,0,20); diagc<-diag(sqrt(c))

# Gaussian: Data generation of factor model(y=Bf+e)
fm<-matrix(rnorm(f*n,0,1),f,n); epsilon<-matrix(rnorm(n*p,0,1),p,n)  
Data<-beta%*%fm+diagc%*%epsilon                                        
X=Data-matrix(rep(t(apply(Data,1,mean)),n),ncol=n)
X <- t(X)           

X_cor <- cor(X)
X_cor
N<-nrow(X); P<-ncol(X)
R=99 # number of permutation

obs_eigenvalues <- eigen(X_cor)$values
obs_eigenvalues <- sort(obs_eigenvalues, decreasing=TRUE)
nge <- rep(0, p) 

# Repeat R times
A<- matrix(0, nrow=R+1, ncol=p)
A[1,] <- obs_eigenvalues
for (r in 1:R) {
  for (j in 2:p) {
    X[,j] <- X[sample(N),j] # Permute columns 2, ..., P of x
  }
  permutX_cor<-cor(X)
  permuted_eigenvalues <- eigen(permutX_cor)$values # Obtain the eigenvalues of the permuted data
  permuted_eigenvalues <- sort(permuted_eigenvalues, decreasing=TRUE)
  A[r+1,] <- permuted_eigenvalues
  nge <- nge + (permuted_eigenvalues >= obs_eigenvalues) # Increase the count if the absolute value of the i-th permuted eigenvalue is greater than or equal the absolute value of the observed i-th eigenvalue
}
# Finally
p_values <- (nge + 1) / (R + 1) # Obtain p-values for each eigenvalue
p_values

sig_level=0.05
sig_factor=0
for (i in 1:p){
  if (p_values[i]<sig_level){
    sig_factor=sig_factor+1
  }
}
sig_factor
head(obs_eigenvalues)
permutation_quantiles<-apply(A, 2, function(x) quantile(x, p = c(0,0.25,0.5,0.75,0.9,0.95,0.99)))
result<-cbind(obs_eigenvalues, t(permutation_quantiles), p_values)
rownames(result) <- c(seq(1,p,1))
result

