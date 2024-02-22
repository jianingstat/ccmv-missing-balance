The codes are designed for the simulation study.

# Preliminary setup
Suppose that we are interested in the coefficients of logistic regression:
$$\mathrm{logit}\{P(Y=1\mid X)\}=\theta_1X_1+\theta_2X_2+\theta_3X_3+\theta_4$$ 
where the true coefficients $\theta_0=(1,-1,1,-2)$ 
and $X_j,\ j=1,2,3$, are generated independently from a truncated standard normal distribution with support $[-3,3]$.
The response variable $Y$ is always observed while $X_1, X_2, X_3$ could be missing. 
Four non-monotone missing patterns are considered: 
| $R=1111$ | $R=1110$ | $R=1101$ | $R=1100$ |    
|---|---|---|---|
| $L^{1111}=(Y,X_1,X_2,X_3)$ | $L^{1110}=(Y,X_1,X_2)$ | $L^{1101}=(Y,X_1,X_3)$ | $L^{1000}=(Y,X_1)$ |

The categorical variable $R\in\\{1111,1110,1101,1100\\}$ indicating the missing patterns is generated from a multinomial distribution with the following probabilities. 
$$P(R=1111\mid L)=\frac{1}{1+\mathrm{Odds}^{1110}(L^{1110})+\mathrm{Odds}^{1101}(L^{1101})+\mathrm{Odds}^{1100}(L^{1100})},$$
$$P(R=1110\mid L)=\frac{\mathrm{Odds}^{1110}(L^{1110})}{1+\mathrm{Odds}^{1110}(L^{1110})+\mathrm{Odds}^{1101}(L^{1101})+\mathrm{Odds}^{1100}(L^{1100})},$$
$$P(R=1101\mid L)=\frac{\mathrm{Odds}^{1101}(L^{1101})}{1+\mathrm{Odds}^{1110}(L^{1110})+\mathrm{Odds}^{1101}(L^{1101})+\mathrm{Odds}^{1100}(L^{1100})},$$
$$P(R=1100\mid L)=\frac{\mathrm{Odds}^{1100}(L^{1100})}{1+\mathrm{Odds}^{1110}(L^{1110})+\mathrm{Odds}^{1101}(L^{1101})+\mathrm{Odds}^{1100}(L^{1100})},$$
where $\mathrm{Odds}^r(L^r)=\frac{P(R=r\mid L^r)}{P(R=1111\mid L^r)}$ is the true propensity odds, which only depends on the shared variables $L^r$, for $r=1110,1101,1100$. 
Therefore, the probabilities satisfy the complete-case missing variable condition (CCMV), which states that the propensity odds only depend on the commonly observed variables. 
More precisely, we have the following settings such that the logarithms of true propensity odds are functions of $L^r$.

### Setting 1
$$\log\mathrm{Odds}^{1110}(L^{1110})=X_1+X_2-I_{Y=1}-0.5,$$

$$\log\mathrm{Odds}^{1101}(L^{1101})=\frac{1}{2}X_1+X_3+-0.5I_{Y=1}-0.3,$$

$$\log\mathrm{Odds}^{1100}(L^{1100})=\frac{3}{2}X_1-I_{Y=1}-0.4.$$

### Setting 2
$$\log\mathrm{Odds}^{1110}(L^{1110})=\frac{1}{5}(X_1^2-9)(X_1+1.5)+\frac{1}{5}(X_2^2-9)(X_2+1)+\frac{1}{10}(X_1+2)(X_2+2)(X_2-1)-2I_{Y=1}+3,$$

$$\log\mathrm{Odds}^{1101}(L^{1101})=-\frac{1}{5}(X_1^2-9)(X_1+1)+\frac{1}{5}(X_3^2-9)(X_3+1.5)-2I_{Y=1},$$

$$\log\mathrm{Odds}^{1100}(L^{1100})=-\frac{1}{5}(X_1+2)(X_1+0.5)(X_1-4)-2I_{Y=1}-1.$$

### Setting 3
$$\log\mathrm{Odds}^{1110}(L^{1110})=\frac{1}{10}(X_1^2-9)(X_1^2-4)X_1+\frac{1}{10}(X_2^2-9)(X_2+1)+\frac{1}{4}(X_1+2)(X_2+2)(X_2-1)-2I_{Y=1},$$

$$\log\mathrm{Odds}^{1101}(L^{1101})=\frac{1}{10}(X_1^2-9)(X_1+1)+\frac{1}{10}(X_3^2-9)(X_3^2-4)X_3+I_{Y=1}\\{(X_1+1)(X_3^2-4)-2\\},$$

$$\log\mathrm{Odds}^{1100}(L^{1100})=I_{Y=0}\\{\frac{1}{5}(X_1^2-9)(X_1^2-4)X_1-1\\}-I_{Y=1}\\{\frac{1}{10}(X_1^2-9)(X_1^2-2.5^2)(X_1+0.5)\\}.$$

# Code example
```
################
#### set up ####
################
source("0data.R"); source("1basis.R"); source("2utility.R")
path1 <- paste(getwd(),"/result",sep=""); dir.create(path=path1, showWarnings = FALSE)
packagelist <- c("fda","prodlim","caret","nloptr","doParallel","doRNG","parallel","glmnet","pracma","plyr","truncnorm")
lapply(packagelist, require, character.only = TRUE)

#######################
#### simulate data ####
#######################
set.seed(123)
N=10000; M=4; vars=c(1:4); norder=4; nbasis=4
cont="tensor_product"; cate="tensor_product"; missmech="CCMV0"
arg <- list(N=N, M=M, vars=vars, missmech=missmech)
tmp1 <- generate_data(arg=arg)
data_full <- tmp1$data_full; data <- tmp1$data; pattern <- tmp1$pattern; odds <- tmp1$odds
arg <- append(arg, list(norder=norder, nbasis=nbasis, cont=cont, cate=cate))
tmp2 <- generate_basis_list(data=data, pattern=pattern, arg=arg)
basis_list <- tmp2$basis_list; penalty_list <- tmp2$penalty_list

###########################
##### balancing method ####
###########################
arg_cv <- append(arg, list(maxit=2500, nfld=5, gammas=c(0,0.1,0.5,0.9,1), lams_init=c(1e4, 1e-10), nlam=15, path=path1))
arg_cv$loss <- "tailored"
crossvalidation(iter=1, arg=arg_cv, data=data, pattern=pattern, basis_list=basis_list, penalty_list=penalty_list)
```

# Reference
Jianing Dong and Raymond K. W. Wong and Kwun Chuen Gary Chan. (2024) "Balancing Method for Non-monotone Missing Data". [link](https://arxiv.org/abs/2402.08873).










