###################### initialization ######################
path1 <- paste(getwd(),"/result",sep=""); dir.create(path=path1, showWarnings = FALSE)
packagelist <- c("fda","prodlim","caret","nloptr","doParallel","doRNG","parallel","glmnet","pracma","plyr","truncnorm"); lapply(packagelist, require, character.only = TRUE)
n_iter <- 1000
maxcore <- detectCores()-2; cl <- makeCluster(maxcore, outfile="log.txt"); registerDoParallel(cl)

foreach (iter = 1:n_iter, .packages = packagelist) %dorng%{
  ###################### initialization ###################### 
  source("0data.R", local=TRUE); source("1basis.R", local=TRUE); source("2utility.R", local=TRUE)
  path2 <- paste(path1,"/",iter,sep=""); dir.create(path=path2, showWarnings = FALSE)

  ###################### generate data and basis ######################
  set.seed(iter)
  N=10000; M=4; vars=c(1:4); norder=4; nbasis=4; cont="tensor_product"; cate="tensor_product" 
  missmech="CCMV0"
  #missmech="CCMV1"
  #missmech="CCMV2"

  arg <- list(N=N, M=M, vars=vars, missmech=missmech)
  tmp1 <- generate_data(arg=arg)
  data_full <- tmp1$data_full; data <- tmp1$data; pattern <- tmp1$pattern; odds <- tmp1$odds
  arg <- append(arg, list(norder=norder, nbasis=nbasis, cont=cont, cate=cate))
  tmp2 <- generate_basis_list(data=data, pattern=pattern, arg=arg)
  basis_list <- tmp2$basis_list; penalty_list <- tmp2$penalty_list
  
  ###################### 1 full data ######################
  w_full <- rep(1, N)
  lm <- nloptr(rep(0,4), eval_f = loglik, eval_grad_f = loglik_grad, data=data_full, w=w_full, opts=list(algorithm="NLOPT_LD_LBFGS", maxeval=1e4, xtol_abs=1e-8))
  theta_full <- as.numeric(lm$solution)
  repeat{if(system2(command="mkdir", arg="lockdir", stderr=NULL)==0){break}}
  cat(theta_full, "\n", file=paste(path1,"/theta_full.txt",sep=""), append=TRUE)
  cat(theta_full, "\n", file=paste(path2,"/theta_full.txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
  
  ###################### 2 complete case ######################
  w_cc <- rep(0, N); w_cc[pattern==1] <- 1
  lm <- nloptr(rep(0,4), eval_f = loglik, eval_grad_f = loglik_grad, data=data, w=w_cc, opts=list(algorithm="NLOPT_LD_LBFGS", maxeval=1e4, xtol_abs=1e-8))
  theta_cc <- as.numeric(lm$solution)
  repeat{if(system2(command="mkdir", arg="lockdir", stderr=NULL)==0){break}}
  cat(theta_cc, "\n", file=paste(path1,"/theta_cc.txt",sep=""), append=TRUE)
  cat(theta_cc, "\n", file=paste(path2,"/theta_cc.txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
  
  ###################### 3 true propensity ######################
  #### generate weights and estimate parameters of interest
  w_true <- rep(0, N)
  w_true[pattern==1] <- rowSums(odds[pattern==1, ])
  lm <- nloptr(rep(0,4), eval_f = loglik, eval_grad_f = loglik_grad, data=data, w=w_true, opts=list(algorithm="NLOPT_LD_LBFGS", maxeval=1e4, xtol_abs=1e-8))
  theta_trueweight <- as.numeric(lm$solution)
  repeat{if(system2(command="mkdir", arg="lockdir", stderr=NULL)==0){break}}
  cat(theta_trueweight, "\n", file=paste(path1,"/theta_trueweight.txt",sep=""), append=TRUE)
  cat(theta_trueweight, "\n", file=paste(path2,"/theta_trueweight.txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
  
  ###################### 4 tailored loss ######################
  arg_cv <- append(arg, list(maxit=2500, nfld=5, gammas=c(0,0.1,0.5,0.9,1), lams_init=c(1e4, 1e-10), nlam=15, path=path1))
  arg_cv$loss <- "tailored"
  crossvalidation(iter=iter, arg=arg_cv, data=data, pattern=pattern, basis_list=basis_list, penalty_list=penalty_list)
  
  ###################### 5 entropy loss ######################
  arg_cv$loss <- "entropy"
  crossvalidation(iter=iter, arg=arg_cv, data=data, pattern=pattern, basis_list=basis_list, penalty_list=penalty_list)
  
  ###################### 6 entropy loss (parametric1) ######################
  data <- as.matrix(data, ncol=M) # convert the data.frame to matrix
  odds_ipw1 <- matrix(0, nrow=N, ncol=M)
  odds_ipw1[, 1] <- 1
  for (j in 2:M) {
    ind1j <- c(which(pattern==1), which(pattern==j))
    data1j <- data[ind1j, ]
    if (j==2) basis1j <- cbind(rep(1,length(ind1j)), data1j[, c(1,2,3)])
    if (j==3) basis1j <- cbind(rep(1,length(ind1j)), data1j[, c(1,2,4)])
    if (j==4) basis1j <- cbind(rep(1,length(ind1j)), data1j[, c(1,2)])
    pattern1j <- rep(1, length(ind1j))
    pattern1j[1:sum(pattern==1)] <- 0
    fit_ipw1 <- glm(pattern1j ~ basis1j+0, family=binomial)
    p_ipw1 <- fit_ipw1$fitted.values
    odds_ipw1[ind1j, j] <- p_ipw1/(1-p_ipw1)
  }
  w_ipw1 <- rep(0, arg$N); w_ipw1[pattern==1] <- rowSums(odds_ipw1)[pattern==1]
  lm <- nloptr(rep(0,4), eval_f = loglik, eval_grad_f = loglik_grad, data=data, w=w_ipw1, opts=list(algorithm="NLOPT_LD_LBFGS", maxeval=1e4, xtol_abs=1e-8))
  theta_ipw1 <- as.numeric(lm$solution)
  repeat{if(system2(command="mkdir", arg="lockdir", stderr=NULL)==0){break}}
  cat(theta_ipw1, "\n", file=paste(path1,"/theta_ipw1.txt",sep=""), append=TRUE)
  cat(theta_ipw1, "\n", file=paste(path2,"/theta_ipw1.txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
  
  ###################### 7 entropy loss (parametric2) ######################
  odds_ipw2 <- matrix(0, nrow=N, ncol=M)
  odds_ipw2[, 1] <- 1
  for (j in 2:M) {
    ind1j <- c(which(pattern==1), which(pattern==j))
    data1j <- data[ind1j, ]
    if (j==2) basis1j <- cbind(rep(1,length(ind1j)), data1j[, c(1,2,3)], data1j[, c(2,3)]^2, data1j[, c(2,3)]^3, data1j[, 2]*data1j[, 3], data1j[, 2]*data1j[, 3]^2)
    if (j==3) basis1j <- cbind(rep(1,length(ind1j)), data1j[, c(1,2,4)], data1j[, c(2,4)]^2, data1j[, c(2,4)]^3)
    if (j==4) basis1j <- cbind(rep(1,length(ind1j)), data1j[, c(1,2)], data1j[, 2]^2, data1j[, 2]^3)
    pattern1j <- rep(1, length(ind1j))
    pattern1j[1:sum(pattern==1)] <- 0
    fit_ipw2 <- glm(pattern1j ~ basis1j+0, family=binomial)
    p_ipw2 <- fit_ipw2$fitted.values
    odds_ipw2[ind1j, j] <- p_ipw2/(1-p_ipw2) 
  }
  w_ipw2 <- rep(0, arg$N); w_ipw2[pattern==1] <- rowSums(odds_ipw2)[pattern==1]
  lm <- nloptr(rep(0,4), eval_f = loglik, eval_grad_f = loglik_grad, data=data, w=w_ipw2, opts=list(algorithm="NLOPT_LD_LBFGS", maxeval=1e4, xtol_abs=1e-8))
  theta_ipw2 <- as.numeric(lm$solution)
  repeat{if(system2(command="mkdir", arg="lockdir", stderr=NULL)==0){break}}
  cat(theta_ipw2, "\n", file=paste(path1,"/theta_ipw2.txt",sep=""), append=TRUE)
  cat(theta_ipw2, "\n", file=paste(path2,"/theta_ipw2.txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
  
  ###################### 8 entropy loss (parametric3) ######################
  odds_ipw3 <- matrix(0, nrow=N, ncol=M)
  odds_ipw3[, 1] <- 1
  for (j in 2:M) {
    ind1j <- c(which(pattern==1), which(pattern==j))
    data1j <- data[ind1j, ]
    if (j==2) basis1j <- cbind(rep(1,length(ind1j)), data1j[, 1], data1j[, 2], data1j[, 2]^2, data1j[, 2]^3, data1j[, 2]^4, data1j[, 2]^5,
                               data1j[, 3], data1j[, 3]^2, data1j[, 3]^3, data1j[, 2]*data1j[, 3], data1j[, 2]*data1j[, 3]^2)
    if (j==3) basis1j <- cbind(rep(1,length(ind1j)), data1j[, 2], data1j[, 2]^2, data1j[, 2]^3,
                               data1j[, 4], data1j[, 4]^2, data1j[, 4]^3, data1j[, 4]^4, data1j[, 4]^5,
                               data1j[, 1], data1j[, 1]*data1j[, 4], data1j[, 1]*data1j[, 4]^2,
                               data1j[, 1]*data1j[, 2], data1j[, 1]*data1j[, 2]*data1j[, 4], data1j[, 1]*data1j[, 2]*data1j[, 4]^2)
    if (j==4) basis1j <- cbind(rep(1,length(ind1j)), data1j[, 2], data1j[, 2]^2, data1j[, 2]^3, data1j[, 2]^4, data1j[, 2]^5,
                               data1j[, 1]*rep(1,length(ind1j)), data1j[, 1]*data1j[, 2], data1j[, 1]*data1j[, 2]^2, data1j[, 1]*data1j[, 2]^3, data1j[, 1]*data1j[, 2]^4, data1j[, 1]*data1j[, 2]^5)
    pattern1j <- rep(1, length(ind1j))
    pattern1j[1:sum(pattern==1)] <- 0
    fit_ipw3 <- glm(pattern1j ~ basis1j+0, family=binomial)
    p_ipw3 <- fit_ipw3$fitted.values
    odds_ipw3[ind1j, j] <- p_ipw3/(1-p_ipw3)
  }
  w_ipw3 <- rep(0, arg$N); w_ipw3[pattern==1] <- rowSums(odds_ipw3)[pattern==1]
  lm <- nloptr(rep(0,4), eval_f = loglik, eval_grad_f = loglik_grad, data=data, w=w_ipw3, opts=list(algorithm="NLOPT_LD_LBFGS", maxeval=1e4, xtol_abs=1e-8))
  theta_ipw3 <- as.numeric(lm$solution)
  repeat{if(system2(command="mkdir", arg="lockdir", stderr=NULL)==0){break}}
  cat(theta_ipw3, "\n", file=paste(path1,"/theta_ipw3.txt",sep=""), append=TRUE)
  cat(theta_ipw3, "\n", file=paste(path2,"/theta_ipw3.txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
  
  ###################### 9 mean score  ######################
  data_sub <- data
  data_sub[, -1] <- as.numeric(data_sub[, -1] >= 0)
  odds_ms <- odds_meanscore(N=N, M=M, pattern=pattern, data_sub=data_sub)
  w_ms <- rep(0, arg$N); w_ms[pattern==1] <- rowSums(odds_ms)[pattern==1]
  lm <- nloptr(rep(0,4), eval_f = loglik, eval_grad_f = loglik_grad, data=data, w=w_ms, opts=list(algorithm="NLOPT_LD_LBFGS", maxeval=1e4, xtol_abs=1e-8))
  theta_ms <- as.numeric(lm$solution)
  repeat{if(system2(command="mkdir", arg="lockdir", stderr=NULL)==0){break}}
  cat(theta_ms, "\n", file=paste(path1,"/theta_ms.txt",sep=""), append=TRUE)
  cat(theta_ms, "\n", file=paste(path2,"/theta_ms.txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
  
  ###################### 10 mean score (error-prone)  ######################
  data_sub <- data
  err <- matrix(rnorm(N*3, mean=0, sd=0.5), ncol=3)
  data_sub[, -1] <- data_sub[, -1] + err
  data_sub[, -1] <- as.numeric(data_sub[, -1] >= 0)
  odds_ms2 <- odds_meanscore(N=N, M=M, pattern=pattern, data_sub=data_sub)
  w_ms2 <- rep(0, arg$N); w_ms2[pattern==1] <- rowSums(odds_ms2)[pattern==1]
  lm <- nloptr(rep(0,4), eval_f = loglik, eval_grad_f = loglik_grad, data=data, w=w_ms2, opts=list(algorithm="NLOPT_LD_LBFGS", maxeval=1e4, xtol_abs=1e-8))
  theta_ms2 <- as.numeric(lm$solution)
  repeat{if(system2(command="mkdir", arg="lockdir", stderr=NULL)==0){break}}
  cat(theta_ms2, "\n", file=paste(path1,"/theta_ms2.txt",sep=""), append=TRUE)
  cat(theta_ms2, "\n", file=paste(path2,"/theta_ms2.txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
  
}
stopCluster(cl)

