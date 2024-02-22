######## library(truncnorm)
generate_data <- function(arg){
  ######################  initialization  ######################
  theta <- c(-2, 1, -1, 1)  
  X <- matrix(rtruncnorm(arg$N*3, a=-3, b=3, mean=0, sd=1), nrow=arg$N, ncol=3)
  logitp <- as.numeric(theta[1] + X%*%theta[2:4])
  p <- 1/ (1 + exp(-1 * logitp))
  Y <- rbinom(arg$N, size=1, prob=p)
  data_full <- data.frame(Y, X)
  colnames(data_full) <- c("Y", "X1", "X2", "X3")
  data2 <- data_full; data2$X3 <- NA
  data3 <- data_full; data3$X2 <- NA
  data4 <- data_full; data4$X2 <- NA; data4$X3 <- NA
  arg_init <- list(vars=arg$vars, norder=5, nbasis=10, cont="additive", cate="tensor_product")
  basis_2 <- generate_basis(data=data2, arg=arg_init)$basis
  basis_3 <- generate_basis(data=data3, arg=arg_init)$basis
  basis_4 <- generate_basis(data=data4, arg=arg_init)$basis

  ######################  propensity score  ######################
  odds <- matrix(0, nrow=arg$N, ncol=4)
  if (arg$missmech=="CCMV0"){
    odds[, 2] <- X[,1] + X[,2] - (Y==1) - 0.5
    odds[, 3] <- X[,1]/2 + X[,3] - 0.5*(Y==1) - 0.3
    odds[, 4] <- 1.5*X[,1] - (Y==1) - 0.4
  }
  if (arg$missmech=="CCMV1"){
    odds[, 2] <- (X[,1]+3)*(X[,1]+1.5)*(X[,1]-3)/5 + (X[,2]+3)*(X[,2]+1)*(X[,2]-3)/5 + (X[,1]+2)*(X[,2]+2)*(X[,2]-1)/10 + 3 - 2*(Y==1)
    odds[, 3] <- -(X[,1]+3)*(X[,1]+1)*(X[,1]-3)/5 + (X[,3]+3)*(X[,3]+1.5)*(X[,3]-3)/5  - 2*(Y==1)
    odds[, 4] <- -(X[,1]+2)*(X[,1]+0.5)*(X[,1]-4)/5 - 1 - 2*(Y==1) 
  }
  if (arg$missmech=="CCMV2"){
    odds[, 2] <- (X[,1]+3)*(X[,1]+2)*X[,1]*(X[,1]-2)*(X[,1]-3)/10 + (X[,2]+3)*(X[,2]+1)*(X[,2]-3)/10 - (X[,1]+2)*(X[,2]+2)*(X[,2]-1)/4 - 2*(Y==1)
    odds[, 3] <- (X[,1]+3)*(X[,1]+1)*(X[,1]-3)/10 + (X[,3]+3)*(X[,3]+2)*X[,3]*(X[,3]-2)*(X[,3]-3)/10 + (Y==1)*((X[,1]+1)*(X[,3]+2)*(X[,3]-2)/3 - 2)
    odds[, 4] <- ((X[,1]+3)*(X[,1]+2)*X[,1]*(X[,1]-2)*(X[,1]-3)/5 - 1)*(Y==0) - ((X[,1]+3)*(X[,1]+2.5)*(X[,1]+0.5)*(X[,1]-2.5)*(X[,1]-3)/10)*(Y==1)
  }
  
  ######################  missing pattern  ######################
  odds <- exp(odds)
  prop <- odds/rowSums(odds)
  pattern <- rep(NA, arg$N)
  for (i in 1:arg$N)
    pattern[i] <- sample(x=1:4, size=1, prob=prop[i, ])
  
  ######################  missing data  ######################
  data <- data_full
  data$X3[pattern==2] <- NA #X3 is missing
  data$X2[pattern==3] <- NA #X2 is missing
  data$X2[pattern==4] <- NA; data$X3[pattern==4] <- NA #both X2 X3 are missing
  return(list(data_full=data_full, data=data, pattern=pattern, odds=odds))
}