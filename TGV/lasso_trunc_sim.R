library(glmnet)
library(R.utils)
source("funs.fixed.R")

source("funs.inf.R")
source("lasso_trunc_gauss.R")
source("lasso_thresh_intervals.R")
source("myfuns.R")

load("sim_param_list.RData")
num_sim = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

pinv=solve

options(error=dump.frames)
set.seed(322)
nsim=10

  #type="full"
  type="partial"

  param_settings = sim_list[[num_sim]]
  
  n = param_settings$n
  p = param_settings$p
  
  lam = param_settings$lam
  lam_high = param_settings$lam_high
  
  btrue = param_settings$btrue
  
  signal_strength = param_settings$signal_strength

  verbose=F
  
  # arrays to store confidence intervals
  ci_lee=ci_naive=ci_bonf=ci_lassothresh=ci_lsthresh=ci_nosign=array(NA,c(p,2,nsim))
  high_val_lassothresh = high_val_lsthresh = matrix(F, p, nsim)
  is_finite_filter = matrix(F,p,nsim)
  
  betaall=matrix(NA,nsim,p)
  btruepart=matrix(NA,p,nsim)
  btrueridge = matrix(NA,p,nsim)
  bthresh = matrix(NA, p ,nsim)
  bthresh_lasso = matrix(NA, p, nsim)
  aall=acterr=rep(NA,nsim)

  # noise level
  sigma= param_settings$sigma

  # desired miscoverage
  alpha= 0.1
  zalpha=-qnorm(alpha/2)
  alpha_bonf=alpha/p
  zalpha_bonf=-qnorm(alpha_bonf/2)
  t_thresh = zalpha_bonf

time1 = proc.time()
for(ii in 1:nsim){
  
  beta = 0
  while(sum(abs(beta)) < 1e-10){
    x=matrix(rnorm(n*p),n,p)
    x=scale(x,T,T)
    
    mutrue=x%*%btrue
    y=mutrue+sigma*rnorm(n)
    y=y-mean(y)
    
    a=glmnet(x,y,standardize=F, intercept = F, thresh = 1e-12, lambda = lam)
    beta=as.numeric(a$beta)
      
    aa = glmnet(x,y,standardize=F, intercept = F, thresh = 1e-12, lambda = lam_high)
    beta_high = as.numeric(aa$beta)
  }
  
  cat(ii)       
  betaall[ii,]=beta
  
  # least squares on selected variables
  a=lsfit(x[,beta!=0],y)
  aa=ls.diag(a)
  bhat=a$coef[-1]
  bhat0=a$coef[1]
  act=which(beta!=0)
  se=aa$std.err[-1]

  # true partial coefficients
  btruepart[,ii]=0
  btruepart[act,ii]=lsfit(x[,act,drop=F],mutrue)$coef[-1]

  # compute population ls_thresh coefficients
  bthresh[,ii] = 0
  tstat = compute_zstat(y, x, which(beta!=0), sigma)
  eta = form_ls_thresh_target(x, which(beta!=0), tstat, t_thresh)
  bthresh[which(beta!=0),ii] = eta%*%mutrue

  # compute population lasso thresh coefficients
  bthresh_lasso[,ii] = 0
  b_tmp = correct_lasso_estimate(x, beta, n*lam)
  b_high_tmp = correct_lasso_estimate(x, beta_high, n*lam_high)
  eta = form_lasso_thresh_target(x, b_tmp, b_high_tmp)
  bthresh_lasso[which(b_tmp!=0),ii] = eta%*%mutrue

  #naive intervals
  ci_naive[beta!=0,1,ii]=bhat-zalpha*se
  ci_naive[beta!=0,2,ii]=bhat+zalpha*se
  is_finite = which(ci_naive[,2,ii] - ci_naive[,1,ii] < Inf)
  
  #bonf-adj naive
  ci_bonf[beta!=0,1,ii]=bhat-zalpha_bonf*se
  ci_bonf[beta!=0,2,ii]=bhat+zalpha_bonf*se
  is_finite = intersect(is_finite, which(ci_bonf[,2,ii] - ci_bonf[,1,ii] < Inf))
  
  #lee et al intervals
  lee=get_lee_intervals(y,x,beta,sigma,lam*n, alpha, type = type)
  ci_lee[,,ii]=lee$ci
  is_finite = intersect(is_finite, which(ci_lee[,2,ii] - ci_lee[,1,ii] < Inf))
  
   # no conditioning on sign
   out = censored_intervals(y=y, x=x, beta = beta, sigma= sigma, lambda = n*lam,
                                     alpha = alpha, type = type, cluster_by_cor = F)
   ci_nosign[,,ii] = out$ci
   is_finite = intersect(is_finite, which(ci_nosign[,2,ii] - ci_nosign[,1,ii] < Inf))
   
   # lasso thresholding
   out = lasso_thresh_intervals(y, x, beta, beta_high, sigma, n*lam, n*lam_high, alpha)
   ci_lassothresh[,,ii] = out$ci
   high_val_lassothresh[out$high_value,ii] = T
   is_finite = intersect(is_finite, which(ci_lassothresh[,2,ii] - ci_lassothresh[,1,ii] < Inf))
   
   # t thresholding      
   out = ls_thresh_intervals(y=y, x=x, beta = beta, sigma= sigma, lambda = n*lam,
                                  alpha = alpha, type = type, t_thresh = t_thresh)
   ci_lsthresh[,,ii] = out$ci
   high_val_lsthresh[out$high_value,ii] = T
   is_finite = intersect(is_finite, which(ci_lsthresh[,2,ii] - ci_lsthresh[,1,ii] < Inf))
   
   is_finite_filter[is_finite,ii] = T
   
   param_settings = list(n=n, p = p, lam = lam, sigma = sigma, lam_high = lam_high,
                         alpha = alpha, type = type, t_thresh = t_thresh, nsim = ii, signal_strength = signal_strength)
   sim_name = paste("./data/n",n,"p",p,"lam",round(lam,2),"sig",
                   param_settings$signal_strength,"_type_",type,
                   ".RData", sep = "")
   save(btrue, btruepart, bthresh, bthresh_lasso, ci_naive,
        ci_bonf, ci_lee, ci_nosign, ci_lassothresh, ci_lsthresh,
        high_val_lassothresh, high_val_lsthresh, is_finite_filter, param_settings, file = sim_name)
}
time1 - proc.time()

