summary(iris)
plot(iris)


library(mlbench)


obj <- mlbench.spirals(100,1,0.025)

my.data1 = iris[,1:4]
my.data2 = 4 * obj$x

plot(my.data1[,1:2], col = iris$Species)
plot(my.data2)

km1<-kmeans(my.data1, centers = 3, iter.max = 10000)
km2<-kmeans(my.data2, centers = 2, iter.max = 10000)

plot(my.data1[,1:2], col=km1$cluster)
plot(my.data2[,1:2], col=km2$cluster)



s <- function(x1, x2, alpha=1) {
  exp(- alpha * norm(as.matrix(x1-x2), type="F"))
}

make.similarity <- function(my.data, similarity) {
  N <- nrow(my.data)
  S <- matrix(rep(NA,N^2), ncol=N)
  for(i in 1:N) {
    for(j in 1:N) {
      S[i,j] <- similarity(my.data[i,], my.data[j,])
    }
  }
  S
}


S1 <- make.similarity(my.data1, s)
S2 <- make.similarity(my.data2, s)

make.affinity <- function(S, n.neighbor=5) {
  N <- length(S[,1])
  if (n.neighbor >= N) { # fully connected
    W <- S
  } else {
    W <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighbor]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        W[i,j] <- S[i,j]
        W[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      } }
  }
  return(W)
}
W1 <- make.affinity(S1, 10)
W2 <- make.affinity(S2, 3)
# Diagonal matrix
D1 <- diag(apply(W1, 1, sum)) # sum rows 
D2 <- diag(apply(W2, 1, sum)) # sum rows 


eigen_result1<-eigen(D1, symmetric = TRUE)
eigen_result2<-eigen(D2, symmetric = TRUE)

D_minus_half1 <- eigen_result1$vectors %*% diag(1/sqrt(eigen_result1$values)) %*% t(eigen_result1$vectors)
D_minus_half2 <- eigen_result2$vectors %*% diag(1/sqrt(eigen_result2$values)) %*% t(eigen_result2$vectors)

L1 <- D1 - W1
L2 <- D2 - W2


L_sym1 <- (D_minus_half1 %*% L1 %*% D_minus_half1)
L_sym2 <- (D_minus_half2 %*% L2 %*% D_minus_half2)





L_eigen1 <- eigen(L1, symmetric=TRUE)
plot(rev(L_eigen1$values)[1:5])

L_eigen2 <- eigen(L2, symmetric=TRUE)
plot(rev(L_eigen2$values)[1:4])

L_sym_eigen1 <- eigen(L_sym1, symmetric=TRUE)
plot(rev(L_sym_eigen1$values)[1:4])

L_sym_eigen2 <- eigen(L_sym2, symmetric=TRUE)
plot(rev(L_sym_eigen2$values)[1:4])


k1 = 3
k2 = 2

Z_unnormal1 <- L_eigen1$vectors[,(ncol(L_eigen1$vectors)-k1+1):ncol(L_eigen1$vectors)]
Z_unnormal2 <- L_eigen2$vectors[,(ncol(L_eigen2$vectors)-k2+1):ncol(L_eigen2$vectors)]

Z_sym1 <- L_sym_eigen1$vectors[,(ncol(L_sym_eigen1$vectors)-k1+1):ncol(L_sym_eigen1$vectors)]
Z_sym2 <- L_sym_eigen2$vectors[,(ncol(L_sym_eigen2$vectors)-k2+1):ncol(L_sym_eigen2$vectors)]

library(stats)
result_unnormal1 <- kmeans(Z_unnormal1, centers=k1, nstart=5)
result_unnormal2 <- kmeans(Z_unnormal2, centers=k2, nstart=5)
result_sym1 <- kmeans(Z_sym1, centers=k1, nstart=5)
result_sym2 <- kmeans(Z_sym2, centers=k2, nstart=5)


plot(my.data1[,1:2],col=iris$Species)
title(main = "true data")

plot(my.data1[,1:2], col=result_unnormal1$cluster)
title(main = "Clustering with Unormalized laplacian")

plot(my.data1[,1:2], col=result_sym1$cluster)
title(main = "Clustering with Symmetric laplacian")



plot(my.data2[,1:2], col=result_unnormal2$cluster)
title(main = "Clustering with Unormalized laplacian")

plot(my.data2[,1:2], col=result_sym2$cluster)
title(main = "Clustering with Symmetric laplacian")

