######## library(caret); library(glmnet); library(nloptr); library(plyr); library(pracma); library(prodlim)
######## standardization function
standard <- function(x){
  a <- range(x, na.rm=T)[1]
  b <- range(x, na.rm=T)[2]
  if (a==b) return(x) else return((x-a)/(b-a))
}

####### negative loglikelihood function of X and Y (parameter of interest)
loglik <- function(theta, data, w){
  X <- as.matrix(data[, -1], ncol=length(theta[-1]))
  X_tr <- cbind(rep(1,dim(X)[1]), X)
  Y <- data[, 1]
  logitp <- as.numeric(X_tr %*% theta)
  logliks <- -(Y==1)*logitp + log(1+exp(logitp))
  loglik_weighted_mean <- sum(w[w!=0] * logliks[w!=0])/length(w)
  return(loglik_weighted_mean)
}

loglik_grad <- function(theta, data, w){
  X <- as.matrix(data[, -1], ncol=length(theta[-1]))
  X_tr <- cbind(rep(1,dim(X)[1]), X)
  Y <- data[, 1]
  logitp <- as.numeric(X_tr %*% theta)
  logliks_grad <- -(Y==1)*X_tr + 1/(1+exp(-logitp))*X_tr
  loglik_grad_weighted_mean <- colMeans(matrix(w[w!=0] * logliks_grad[w!=0, ], ncol=length(theta)))
  return(as.numeric(loglik_grad_weighted_mean))
}


####### scoring rule (proposed method and penalty) #########
tailored_loss <- function(alpha, pattern, basis, penalty, lambda, gamma, j){
  loss1 <- exp(basis[pattern==1, ] %*% alpha)
  lossj <- -basis[pattern==j, ] %*% alpha
  loss <- (sum(loss1) + sum(lossj))/sum(pattern==1 | pattern==j)
  penalty1 <- sum(abs(alpha) * sqrt(penalty))
  penalty2 <- sum(alpha^2 * penalty)
  loss <- loss + lambda * (gamma*penalty1 + (1-gamma)*penalty2)
  return(loss)
}

tailored_loss_grad <- function(alpha, pattern, basis, penalty, lambda, gamma, j){
  loss1_grad <- as.numeric(exp(basis[pattern==1, ] %*% alpha)) * basis[pattern==1, ]
  lossj_grad <- -basis[pattern==j, ]
  loss1_grad <- matrix(loss1_grad, ncol=length(alpha))
  lossj_grad <- matrix(lossj_grad, ncol=length(alpha))
  loss_grad <- (colSums(loss1_grad) + colSums(lossj_grad))/sum(pattern==1 | pattern==j)
  penalty1_grad <- sign(alpha) * sqrt(penalty) 
  penalty2_grad <- 2 * alpha * penalty
  loss_grad <- loss_grad + lambda * (gamma*penalty1_grad + (1-gamma)*penalty2_grad)
  return(loss_grad)
}

####### restricted entropy (pairwise logistic regression) #########
entropy_loss <- function(alpha, pattern, basis, penalty, lambda, gamma, j){
  loss1 <- log(1+exp(basis[pattern==1, ] %*% alpha))
  lossj <- log(1+exp(-basis[pattern==j, ] %*% alpha))
  loss <- (sum(loss1) + sum(lossj))/sum(pattern==1 | pattern==j)
  penalty1 <- sum(abs(alpha) * sqrt(penalty))
  penalty2 <- sum(alpha^2 * penalty)
  loss <- loss + lambda * (gamma*penalty1 + (1-gamma)*penalty2)
  return(loss)
}

entropy_loss_grad <- function(alpha, pattern, basis, penalty, lambda, gamma, j){
  loss1_grad <- 1/(1+as.numeric(exp(-basis[pattern==1, ] %*% alpha))) * basis[pattern==1, ]
  lossj_grad <- -1/(1+as.numeric(exp(basis[pattern==j, ] %*% alpha))) * basis[pattern==j, ]
  loss1_grad <- matrix(loss1_grad, ncol=length(alpha))
  lossj_grad <- matrix(lossj_grad, ncol=length(alpha))
  loss_grad <- (colSums(loss1_grad) + colSums(lossj_grad))/sum(pattern==1 | pattern==j)
  penalty1_grad <- sign(alpha) * sqrt(penalty) 
  penalty2_grad <- 2 * alpha * penalty
  loss_grad <- loss_grad + lambda * (gamma*penalty1_grad + (1-gamma)*penalty2_grad)
  return(loss_grad)
}

####### cross validation for (either tailored or entropy) #########
crossvalidation <- function(iter, arg, data, pattern, basis_list, penalty_list){
  if (is.null(arg$maxit)) maxit <- 2500 else maxit <- arg$maxit
  if (is.null(arg$nfld)) nfld <- 5 else nfld <- arg$nfld
  if (is.null(arg$gammas)) gammas <- c(0, 0.5, 1) else gammas <- arg$gammas
  if (is.null(arg$lams_init)) lams_init <- c(1e3, 1e-8) else lams_init <- arg$lams_init
  if (is.null(arg$nlam)) nlam <- 12 else nlam <- arg$nlam
  flds <- createFolds(y=c(1:arg$N), k=nfld)
  lams <- exp(seq(log(lams_init[1]), log(lams_init[2]), len=nlam))
  alpha_star <- list()
  loss_cv <- list()
  for (j in 2:arg$M){
    alpha_star[[j]] <- list()
    loss_cv[[j]] <- matrix(0, nrow=length(gammas), ncol=nlam)
    for (igamma in 1:length(gammas)){
      for (ifld in 1:nfld){
        # cross validation : training data split
        ind_valid <- flds[[ifld]]; ind_train <- setdiff(c(1:arg$N), ind_valid)
        pattern_valid <- pattern[ind_valid]; pattern_train <- pattern[ind_train]
        basis_valid <- basis_list[[j]][ind_valid, ]; basis_train <- basis_list[[j]][ind_train, ]
        alpha_start <- rep(0.1, ncol(basis_list[[j]]))
        for (ilam in 1:nlam){ 
          # solver
          result_train <- nloptr(x0=alpha_start, eval_f=get(paste(arg$loss,"_loss",sep="")), eval_grad_f=get(paste(arg$loss,"_loss_grad",sep="")),
                                 pattern=pattern_train, basis=basis_train, penalty=penalty_list[[j]], lambda=lams[ilam], gamma=gammas[igamma], j=j,
                                 opts=list(algorithm="NLOPT_LD_LBFGS", maxeval=maxit, ftol_abs=1e-8)) # faster
          alpha0 <- result_train$solution
          if (result_train$status <= 0){
            result_train <- nloptr(x0=alpha_start, eval_f=get(paste(arg$loss,"_loss",sep="")), 
                                   pattern=pattern_train, basis=basis_train, penalty=penalty_list[[j]], lambda=lams[ilam], gamma=gammas[igamma], j=j,
                                   opts=list(algorithm="NLOPT_LN_COBYLA", maxeval=maxit, ftol_abs=1e-8)) # slower
            alpha0 <- result_train$solution
          }
          if (result_train$status <= 0)
            alpha0 <- rep(0, ncol(basis_list[[j]]))
          # validation
          if (arg$loss=="tailored"){
            loss1 <- exp(basis_valid[pattern_valid==1, ] %*% alpha0)
            lossj <- -basis_valid[pattern_valid==j, ] %*% alpha0
            loss1 <- as.numeric(loss1)
            lossj <- as.numeric(lossj)
          }else{
            loss1 <- log(1+exp(basis_valid[pattern_valid==1, ] %*% alpha0))
            lossj <- log(1+exp(-basis_valid[pattern_valid==j, ] %*% alpha0))
            loss1 <- as.numeric(loss1)
            lossj <- as.numeric(lossj)
          }  
          loss_cv[[j]][igamma, ilam] <- mean(c(loss1, lossj), na.rm=T) + loss_cv[[j]][igamma, ilam]
        }
      }
      # re-training all data to get the estimates
      alpha_star[[j]][[igamma]] <- matrix(NA, nrow=nlam, ncol=ncol(basis_list[[j]]))
      alpha_start <- rep(0.1, ncol(basis_list[[j]]))
      for (ilam in 1:nlam){
        # solver
        result <- nloptr(x0=alpha_start, eval_f=get(paste(arg$loss,"_loss",sep="")), eval_grad_f=get(paste(arg$loss,"_loss_grad",sep="")),
                         pattern=pattern, basis=basis_list[[j]], penalty=penalty_list[[j]], lambda=lams[ilam], gamma=gammas[igamma], j=j,
                         opts=list(algorithm="NLOPT_LD_LBFGS", maxeval=maxit, ftol_abs=1e-8)) # faster
        alpha0 <- result$solution
        if (result$status <= 0){
          result <- nloptr(x0=alpha_start, eval_f=get(paste(arg$loss,"_loss",sep="")),
                           pattern=pattern, basis=basis_list[[j]], penalty=penalty_list[[j]], lambda=lams[ilam], gamma=gammas[igamma], j=j,
                           opts=list(algorithm="NLOPT_LN_COBYLA", maxeval=maxit, ftol_abs=1e-8)) # slower
          alpha0 <- result$solution
        }
        if (result$status <= 0)
          alpha0 <- rep(0, ncol(basis_list[[j]]))
        # record the estimates
        alpha_star[[j]][[igamma]][ilam, ] <- alpha0
      }
    }
  }
  theta_true <- c(-2, 1, -1, 1) 
  odds <- matrix(0, nrow=arg$N, ncol=arg$M)
  for (j in 2:arg$M){
    igamma <- (which.min(loss_cv[[j]])-1) %% length(gammas) + 1
    ilam <- (which.min(loss_cv[[j]])-1) %/% length(gammas) + 1
    odds[, j] <- basis_list[[j]] %*% alpha_star[[j]][[igamma]][ilam, ]
  }
  odds <- exp(odds); w <- rep(0, arg$N); w[pattern==1] <- rowSums(odds)[pattern==1]
  lm <- nloptr(rep(0,4), eval_f=loglik, eval_grad_f=loglik_grad, data=data, w=w, opts=list(algorithm="NLOPT_LD_LBFGS", maxeval=1e4, xtol_abs=1e-8))
  theta <- as.numeric(lm$solution)
  sd <- sd_score(N=arg$N, M=arg$M, data=data, basis_list=basis_list, penalty_list=penalty_list, pattern=pattern, theta=theta, odds=odds, w=w)
  CI_lower <- theta - qnorm(0.975)*sd; CI_upper <- theta + qnorm(0.975)*sd
  CI_cover <- as.numeric((CI_lower < theta_true) & (theta_true < CI_upper))
  
  # output with file lock
  repeat{if(system2(command="mkdir", arg="lockdir",stderr)==0){break}}
  cat(theta, "\n", file=paste(arg$path,"/theta_",arg$loss,".txt",sep=""), append=TRUE)
  cat(theta, "\n", file=paste(arg$path,"/",iter,"/theta_",arg$loss,".txt",sep=""), append=TRUE)
  cat(theta, "\n", file=paste(arg$path,"/",iter,"/theta_est_",arg$loss,".txt",sep=""), append=TRUE)
  cat(sd, "\n", file=paste(arg$path,"/sd_",arg$loss,".txt",sep=""), append=TRUE)
  cat(sd, "\n", file=paste(arg$path,"/",iter,"/sd_",arg$loss,".txt",sep=""), append=TRUE)
  cat(CI_cover, "\n", file=paste(arg$path,"/CI_cover_",arg$loss,".txt",sep=""), append=TRUE)
  system2(command="rmdir", arg="lockdir")
}

sd_score <- function(N=NULL, M=NULL, data=NULL, basis_list=NULL, penalty_list=NULL, pattern=NULL, theta=NULL, odds=NULL, w=NULL){
  # require package "pracma", "glmnet"
  if(!"package:pracma" %in% search())
    require("pracma")
  if(!"package:glmnet" %in% search())
    require("glmnet")
  # fully-observed pattern
  ind1 <- which(pattern==1)
  #  psi_theta is derivative of loglikelihood
  psi <- matrix(0, nrow=N, ncol=length(theta))
  for (i in ind1)
    psi[i, ] <- loglik_grad(theta, data[i, ], w=1)
  # conditional expectation
  u <- NULL
  for (j in 2:M){
    u[[j]] <- matrix(NA, nrow=N, ncol=length(theta))
    for (k in 1:length(theta)){
      fitjk <- cv.glmnet(y=psi[ind1, k], x=basis_list[[j]][ind1, ], family="gaussian", penalty.factor=penalty_list[[j]], intercept = F)
      u[[j]][, k] <- predict(fitjk, newx = basis_list[[j]], s = "lambda.min")
    }
  }
  # Dtheta is the second derivative of loglikelihood. (hessian matrix)
  Dtheta <- hessian(f=loglik, x0=theta, data=data, w=w)
  inv_Dtheta <- solve(Dtheta)
  # Var_psi is the asymptotic variance of estimator (weighted sum of psi) 
  obj <- matrix(rep(-colMeans(w[ind1]*psi[ind1, ]),N), nrow=N, ncol=length(theta), byrow=TRUE)
  obj[ind1, ] <- obj[ind1, ] + w[ind1]*psi[ind1, ] 
  for (j in 2:M){
    indj <- which(pattern==j)
    obj[ind1, ] <- obj[ind1, ] - odds[ind1, j] * u[[j]][ind1, ]
    obj[indj, ] <- obj[indj, ] + u[[j]][indj, ]
  }
  # sqrt(n) * (theta_hat - theta0) -> N(0, invDtheta %*% Var_psi %*% invDtheta)
  Var_psi <- cov(obj)*(N-1)/N #we want E[obj^2],but cov gives var/(n-1)
  Var_psi <- Var_psi + colMeans(obj) %*% t(colMeans(obj))
  Var_theta <- inv_Dtheta %*% Var_psi %*% inv_Dtheta
  sd_theta <- sqrt(diag(Var_theta)/N)
  return(sd_theta = sd_theta)
}

odds_meanscore <- function(N=NULL, M=NULL, pattern=NULL, data_sub=NULL){
  odds_ms <- matrix(1, nrow=N, ncol=M)
  ind1 <- which(pattern==1)
  # counting numbers by class
  vars12 <- c(1,2,3)
  N2 <- count(df=data_sub[pattern==2, ], vars=vars12)
  N1 <- count(df=data_sub[pattern==1, ], vars=vars12)
  ratio2 <- rep(0, nrow(N2))
  for (i in 1:nrow(N2)){
    row_in_N1 <- row.match(N2[i, names(N2)!="freq"], N1[, names(N1)!="freq"])
    ratio2[i] <- N2[i, "freq"]/N1[row_in_N1, "freq"]
  }
  vars13 <- c(1,2,4)
  N3 <- count(df=data_sub[pattern==3, ], vars=vars13)
  N1 <- count(df=data_sub[pattern==1, ], vars=vars13)
  ratio3 <- rep(0, nrow(N3))
  for (i in 1:nrow(N3)){
    row_in_N1 <- row.match(N3[i, names(N3)!="freq"], N1[, names(N1)!="freq"])
    ratio3[i] <- N3[i, "freq"]/N1[row_in_N1, "freq"]
  }
  vars14 <- c(1,2)
  N4 <- count(df=data_sub[pattern==4, ], vars=vars14)
  N1 <- count(df=data_sub[pattern==1, ], vars=vars14)
  ratio4 <- rep(0, nrow(N4))
  for (i in 1:nrow(N4)){
    row_in_N1 <- row.match(N4[i, names(N4)!="freq"], N1[, names(N1)!="freq"])
    ratio4[i] <- N4[i, "freq"]/N1[row_in_N1, "freq"]
  }
  #### generate weights by counting numbers
  for (i in ind1){
    ind_match12 <- row.match(data_sub[i, vars12], N2[, -ncol(N2)])
    ind_match13 <- row.match(data_sub[i, vars13], N3[, -ncol(N3)])
    ind_match14 <- row.match(data_sub[i, vars14], N4[, -ncol(N4)])
    if (!is.na(ind_match12))
      odds_ms[i, 2] <- ratio2[ind_match12]
    if (!is.na(ind_match13))
      odds_ms[i, 3] <- ratio3[ind_match13]
    if (!is.na(ind_match14))
      odds_ms[i, 4] <- ratio4[ind_match14]
  }
  return(odds_ms=odds_ms)
}
