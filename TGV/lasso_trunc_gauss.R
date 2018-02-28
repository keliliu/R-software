require(intervals)

create_cor_x = function(n,p,p_active, rho, struct = "equi"){
  xactive = c()
  xnull = c()
  if(struct == "equi"){
    for(ii in 1:p_active){
      act_var = rnorm(n)
      xactive = cbind(xactive, act_var)
      
      num_spurious = floor((p - p_active)/p_active)
      xnull = cbind(xnull,rho * outer(act_var, rep(1, num_spurious)) + 
        sqrt(1-rho^2) * matrix(rnorm(n*num_spurious), nr = n))
    }
    x = cbind(xactive, xnull)
    if(dim(x)[2] < p) x = cbind(x, matrix(rnorm(n*(p-dim(x)[2])), nr = n))
  }else if(struct == "toeplitz"){
    cor_mat = toeplitz(rho^(0:(p-1)))
    x = matrix(rnorm(n*p), nr = n)%*%chol(cor_mat)
    x = x[,sample(p,p)]
  }else if(struct == "ind"){
    x = matrix(rnorm(n*p), nr = n) 
  }else{
    stop("Unknown covariance structure specified.")
  }
  return(x)
}

get_ci_length = function(ci_array, filter){
  nn = dim(filter)[2]
  length_vec = c()
  for(ii in 1:nn){
    tmp = ci_array[,2,ii] - ci_array[,1,ii]
    length_vec = c(length_vec,tmp[filter[,ii]])
  }
  return(length_vec)
}
  
choose_lam = function(n,p,sigma,btrue){
  lam_cv = quantile(find_lam_cv(n,p,sigma,btrue)$lam_min,0.5)
  lam_univ = quantile(find_lam_high(n,p,sigma), c(0.05,0.5,0.95))
  return(c(lam_cv,lam_univ))
}

# find distribution of lambda chosen by cv
find_lam_cv = function(n,p,sigma, btrue, simsize = 100, rho = NULL, 
                       struct = "ind"){
  set.seed(334)
  lam_cv = rep(NA, simsize)
  lam_1se = rep(NA,simsize)
  for(ii in 1:simsize){
    if(is.null(rho)){
      x=matrix(rnorm(n*p),n,p)
    }else{
      p_active = sum(btrue!=0)
      x = create_cor_x(n, p, p_active, rho, struct)
    }
    x=scale(x,T,T)
    
    y=x%*%btrue + sigma*rnorm(n)
    y=y-mean(y)
    
    fit = cv.glmnet(x,y)
    lam_cv[ii] = fit$lambda.min
    lam_1se[ii] = fit$lambda.1se
  }
  
  return(list(lam_min = lam_cv, lam_1se = lam_1se))
}

# find parameter for lam high
find_lam_high = function(n,p, sigma, simsize = 1000, rho = NULL,
                         struct = "ind", p_active = NULL){
  set.seed(334)
  max_inner = rep(NA, simsize)
  
  for(ii in 1:simsize){
    if(is.null(rho)){
      x=matrix(rnorm(n*p),n,p)
    }else{
      x = create_cor_x(n, p, p_active, rho, struct)
    }
    x=scale(x,T,T)
    
    y=sigma*rnorm(n)
    y=y-mean(y)
    
    max_inner[ii] = max(abs(t(y)%*%x))/n  
  }
  
  return(max_inner)
}

# returns A, b defining polyhedral constraints via Ay < b
compute_poly = function(x, active, s, lambda, active_only = F){
  # x is n x p matrix
  # active is n x 1 boolean vector indicating whether a variable is active
  # s is vector with values in {-1,1} with length equal to size of active sets
  # lambda is the lasso penalty parameter
  
  AA = c()
  bb = c()
  
  var_act = which(active)
  var_inact = which(!active)
  var_key = c(var_act, var_inact, var_inact)
  sign_key = c(rep(0,length(var_act)), rep(1, length(var_inact)),
               rep(-1, length(var_inact)))
  
  if(sum(active) > 0){
    xact = x[, active, drop = F]
    xxact_inv = solve(crossprod(xact))
    xact_pseudo = xxact_inv%*%t(xact)
    
    if(length(s)==1){
      dsign = matrix(s,1,1)
    }else{
      dsign = diag(s)
    }
    
    AA = -dsign %*%xact_pseudo
    bb = -lambda*dsign%*%xxact_inv%*%s
  }
  
  #if(active_only) return(list(A = AA, b = bb))
  
  if((dim(x)[2] - sum(active)) > 0){
    xinact = x[,!active, drop = F]
    
    if(sum(active) == 0){
      A0 = 1/lambda * t(xinact)
    }else{
      Proj_act = xact%*%xact_pseudo 
      A0 = 1/lambda * (t(xinact)%*%(diag(1,dim(x)[1]) - Proj_act))
    }
    
    AA = rbind(AA, A0, -A0)
    
    if(sum(active) == 0){
      b0 = rep(0, sum(1-active))
    }else{
      b0 = t(xinact)%*%t(xact_pseudo)%*%s
    }
    
    bb = c(bb, rep(1,length(b0)) - b0, rep(1,length(b0)) + b0)
  }
  
  return(list(A = AA, b = bb, var_key = var_key, sign_key = sign_key))
}

# given A,b for the constraint Ay < b, find interval endpoints for eta*y
find_endpoints = function(y, eta, A, b, var_key, sign_key, tol = 1e-12){
  # var_key gives key indicating which variable each row of A (a constraint)
  # corresponds to, needed to compute which variables become active/inactive
  
  eta = as.vector(eta)
  nn = length(y)
  
  cc = eta/sum(eta^2)
  z = y - cc*sum(eta*y)
  
  Ac = as.vector(A%*%cc)
  plus = which(Ac > (tol))
  var_key_plus = var_key[plus]
  sign_key_plus = sign_key[plus]
  minus = which(Ac < (-tol))
  var_key_minus = var_key[minus]
  sign_key_minus = sign_key[minus]
  zero = which(abs(Ac) <= tol)
  
  ee = b - as.vector(A%*%z)
  
  vplus = Inf
  ind_plus = NA # which variable becomes active/inactive at right boundary
  sign_plus = NA # sign of the new active/inactive variable
  if(length(plus) > 0){
    rr = ee[plus]/Ac[plus]
    ind_plus = which.min(rr) 
    vplus = rr[ind_plus]
    sign_plus = sign_key_plus[ind_plus]
    ind_plus = var_key_plus[ind_plus]
  }
  
  vminus = -Inf
  ind_minus = NA # which variable becomes active/inactive at left boundary
  sign_minus = NA
  if(length(minus) > 0){
    rr = ee[minus]/Ac[minus]
    ind_minus = which.max(rr) 
    vminus = rr[ind_minus]
    sign_minus = sign_key_minus[ind_minus]
    ind_minus = var_key_minus[ind_minus]
  }
  
  v0 = 1
  if(length(zero) > 0) v0 = min(ee[zero])
  
  if(v0 > (-tol)){
    return(list(int = c(vminus, vplus), var_left = ind_minus, var_right = ind_plus,
                sign_left = sign_minus, sign_right = sign_plus))
  }else{
    return(list(int = c(NA,NA), var_left = NA, var_right = NA, sign_left = NA,
                sign_right = NA))
  }
}

find_endpoints_defunct = function(y, eta, A, b, Sigmat = NULL, tol = 1e-12){
  eta = as.vector(eta)
  nn = length(y)
  if(is.null(Sigmat)) Sigmat = diag(nn)
  
  cc = as.vector(Sigmat%*%eta/as.numeric(t(eta)%*%Sigmat%*%eta)) 
  z = y - cc*sum(eta*y)
  
  Ac = as.vector(A%*%cc)
  plus = which(Ac > (tol))
  minus = which(Ac < (-tol))
  zero = which(abs(Ac) <= tol)
  
  ee = b - as.vector(A%*%z)
  
  vplus = Inf
  if(length(plus) > 0) vplus = min(ee[plus]/Ac[plus])
  
  vminus = -Inf
  if(length(minus) > 0) vminus = max(ee[minus]/Ac[minus])
  
  v0 = 1
  if(length(zero) > 0) v0 = min(ee[zero])
  
  if(v0 > (-tol)){
    return(c(vminus, vplus))
  }else{
    return(c(NA,NA))
  }
}

# checks KKT condition for the lasso
correct_lasso_estimate = function(x, beta, lambda, tol = 1e-12){
  active = (beta!=0)
  sob = sign(beta[active])
  poly_constraints = compute_poly(x = x, active = active, s = sob, lambda = lambda, active_only = F)
  
  constraints = poly_constraints$A%*%y - poly_constraints$b
  
  incorrect = poly_constraints$var_key[which(constraints > tol)]
  correction = poly_constraints$sign_key[which(constraints > tol)]
  
  beta[incorrect] = correction
  return(beta)
}

# returns P(Z > z | z in union of intervals)
tnorm.union.surv = function(z, mean, sd, intervals, bits = NULL){
  # intervals is a I x 2 matrix of disjoint intervals where the first column contains the lower endpoint
  
  pval = matrix(NA, nr = dim(intervals)[1], nc = length(mean))
  for(jj in 1:dim(intervals)[1]){
    if(z <= intervals[jj,1]){
      pval[jj,] = 1
    }else if(z >= intervals[jj,2]){
      pval[jj,] = 0
    }else{
      pval[jj,] = tnorm.surv(z, mean, sd, intervals[jj,1], intervals[jj,2], bits = bits)
    }
  }
  
  ww = matrix(NA, nr = dim(intervals)[1], nc = length(mean))
  for(jj in 1:dim(intervals)[1]){
    ww[jj,] = pnorm(intervals[jj,2], mean = mean, sd = sd) - pnorm(intervals[jj,1], mean = mean, sd = sd)
  }
  
  if(dim(intervals)[1]==1){
    ww = matrix(1,nr = 1, nc = length(mean))
  }else{
    ww = ww%*%diag(1/apply(ww,2,sum))
  }
  
  pval = apply(pval*ww,2,sum)
  
  return(as.numeric(pval))
}

# inverts truncated normal p-values to create an interval 
create_tnorm_interval = function(y, eta, sigma, alpha, intervals, gridrange=c(-100,100), gridpts = 100, 
                                 griddepth = 2, bits = NULL){
  z = sum(eta*y)
  sd = sqrt(sum(eta^2))*sigma
  
  grid = seq(gridrange[1]*sd,gridrange[2]*sd,length=gridpts)
  fun = function(x) { tnorm.union.surv(z, x, sd, intervals, bits) }
 
  int = grid.search(grid, fun, alpha/2, 1-alpha/2, gridpts, griddepth)
  return(int)
}

random_sign_intervals = function(y, x, beta, sigma, lambda, alpha, num_samples = 20, type = "full",
                                 active_only = F){
  
  active = (beta!=0)
  num_active = sum(active)
  if(num_active==0) stop("No non-zero coefficients")
  sob = sign(beta[active])
  
  if(type == "full"){
    eta = solve(crossprod(x))%*%t(x)
    eta = eta[active,,drop = F]
  }else if(type== "partial"){
    xa = x[, active, drop = F]
    eta = solve(crossprod(xa))%*%t(xa)
  }else{
    stop("Specified type does not exist.")
  }
  
  # normalize eta
  norm_eta = rep(NA, num_active)
  for(jj in 1:num_active){
    norm_eta[jj] = sqrt(sum(eta[jj,]^2))
    eta[jj,] = eta[jj,]/norm_eta[jj]
  }
  
  intervals_list = list()
  poly_constraints = compute_poly(x = x, active = active, s = sob, lambda = lambda, active_only = active_only)
  for(jj in 1:num_active){
    intervals_list[[jj]] = find_endpoints(y, eta[jj,], poly_constraints$A, poly_constraints$b,
                                          poly_constraints$var_key, poly_constraints$sign_key)$int
    if(is.na(intervals_list[[jj]][1]) | (intervals_list[[jj]][2] < intervals_list[[jj]][1])){
      stop(paste("Interval endpoints found for observed sign vector are not compatible.",
                 "Check to make sure intercept = F in glmnet", sep = ""))
    }
  }
  
  if(num_samples > 0){
    for(kk in 1:num_samples){
      poly_constraints = compute_poly(x = x, active = active, s = sob*sign(rnorm(num_active)), lambda = lambda,
                                      active_only = active_only)
      for(jj in 1:num_active){
        int = find_endpoints(y, eta[jj,], poly_constraints$A, poly_constraints$b,
                             poly_constraints$var_key, poly_constraints$sign_key)$int
        if(is.na(int[1]) | is.na(int[2])) next
        if(int[2] < int[1]) next
        intervals_list[[jj]] = rbind(intervals_list[[jj]], int)
      }
    }
  }
  
  for(jj in 1:num_active){
    intervals_list[[jj]] = Intervals(intervals_list[[jj]])
    intervals_list[[jj]] = interval_union(intervals_list[[jj]])
    intervals_list[[jj]] = as.matrix(intervals_list[[jj]])
  }
 
  ci = matrix(NA, nr = dim(x)[2], nc = 2)
  active = which(beta!=0)
  for(jj in 1:num_active){
    ci[active[jj],] = create_tnorm_interval(y, eta[jj,], sigma, alpha, intervals = intervals_list[[jj]])*
                      norm_eta[jj]
  }
 
  return(list(ci = ci, intervals_list = intervals_list))
}

draw_random_M_s = function(var, num_draw, kk, var_memo, sign_memo, max_attempt = 20){
  is_new = 0
  attempt = 0
  while(attempt < max_attempt & is_new == 0){
    is_new = 1
    candidate_var = sort(sample(var,num_draw))
    candidate_s = sign(rnorm(num_draw))
    for(ii in 1:kk){
      if(sum(candidate_var != var_memo[ii,]) == 0 & sum(candidate_s != sign_memo[ii,]) == 0){
        is_new = 0
        break
      }
    }
    attempt = attempt + 1
  }
  return(list(is_new = is_new, M = candidate_var, s = candidate_s))
}

random_model_intervals = function(y, x, beta, sigma, lambda, alpha, num_noise_var = 5, num_samples = 20, type = "full",
                                 active_only = F){
  
  active = (beta!=0)
  num_active = sum(active)
  if(num_active==0) stop("No non-zero coefficients")
  sob = sign(beta[active])
  
  # sample from inactive variables
  inactive = which(beta==0)
  if(length(inactive) > 0){
    num_inactive = length(inactive)
    noise_var = sample(inactive, min(num_noise_var, num_inactive))
  }else{
    noise_var = c()
  }
  
  var_of_interest = c(which(beta!=0), noise_var)
  num_interest = length(var_of_interest)
  
  if(type == "full"){
    eta = solve(crossprod(x))%*%t(x)
    eta = eta[var_of_interest,,drop = F]
  }else if(type== "partial"){
    xa = x[, var_of_interest, drop = F]
    eta = solve(crossprod(xa))%*%t(xa)
  }else{
    stop("Specified type does not exist.")
  }
  
  # normalize eta
  norm_eta = rep(NA, num_interest)
  for(jj in 1:num_interest){
    norm_eta[jj] = sqrt(sum(eta[jj,]^2))
    eta[jj,] = eta[jj,]/norm_eta[jj]
  }
  
  intervals_list = list()
  poly_constraints = compute_poly(x = x, active = active, s = sob, lambda = lambda, active_only = active_only)
  for(jj in 1:num_interest){
    intervals_list[[jj]] = matrix(find_endpoints(y, eta[jj,], poly_constraints$A, poly_constraints$b,
                                                 poly_constraints$var_key,
                                                 poly_constraints$sign_key)$int, nr=1,nc=2)
    if(is.na(intervals_list[[jj]][1,1]) | intervals_list[[jj]][1,2] < intervals_list[[jj]][1,1]){
      stop(paste("Interval endpoints found for observed sign vector are not compatible.",
                 "Check to make sure intercept = F in glmnet", sep = ""))
    }
  }
  
  if(num_samples > 0){
    # keep track of the (M,s) couples we've already sampled
    sampled_var_list = matrix(NA, nr = num_samples, nc = num_active)
    sampled_sign_list = matrix(NA, nr = num_samples, nc = num_active)
    sampled_var_list[1,] = sort(which(beta!=0))
    sampled_sign_list[1,] = sob[order(which(beta!=0))]
    
    for(kk in 1:num_samples){
      active = rep(F, length(beta))
      
      # draw a random (M,s)
      random_M_s = draw_random_M_s(var_of_interest, num_active, kk, sampled_var_list, sampled_sign_list)
      sampled_var_list[min(kk + 1, num_samples),] = random_M_s$M
      sampled_sign_list[min(kk + 1, num_samples),] = random_M_s$s
      if(random_M_s$is_new == 0) next
      
      active[random_M_s$M] = T
      poly_constraints = compute_poly(x = x, active = active, s = random_M_s$s, lambda = lambda,
                                      active_only = active_only)
      
      for(jj in 1:num_interest){
        int = find_endpoints(y, eta[jj,], poly_constraints$A, poly_constraints$b,
                             poly_constraints$var_key,
                             poly_constraints$sign_key)$int
        if(is.na(int[1]) | is.na(int[2])) next
        if(int[2] < int[1]) next
        intervals_list[[jj]] = rbind(intervals_list[[jj]], int)
      }
    }
  }
  
  ci = matrix(NA, nr = dim(x)[2], nc = 2)
  for(jj in 1:num_interest){
    ci[var_of_interest[jj],] = create_tnorm_interval(y, eta[jj,], sigma, alpha, intervals = intervals_list[[jj]])*
      norm_eta[jj]
  }
  
  return(list(ci = ci, intervals_list = intervals_list))
}

censored_intervals = function(y, x, beta, sigma, lambda, alpha, type = "full",
   Sigmat = NULL, cluster_by_cor = T, thresh = T, cor_thresh = 0.2, wcov = 0.9,
   trunc_type = "TG_V"){
    
    active = (beta!=0)
    num_active = sum(active)
    if(num_active==0) stop("No non-zero coefficients")
    sob = sign(beta[active])
    
    if(type == "full"){
      eta = solve(crossprod(x))%*%t(x)
      eta = eta[active,,drop = F]
    }else if(type== "partial"){
      if(cluster_by_cor){
        cluster_x = hclust(as.dist(1-abs(cor(x))))
        cluster_x = cutree(cluster_x, h = 1 - cor_thresh)
        active_cluster = unique(cluster_x[active])
        active_var = c()
        for(kk in 1:length(active_cluster)) active_var = c(active_var, which(cluster_x ==active_cluster[kk]))
        active_var = sort(active_var)
        do_inference = which(active[active_var])
        
        xa = x[, active_var, drop = F]
        #eta = solve(crossprod(xa))%*%t(xa)
        eta = solve(wcov*crossprod(xa) + (1-wcov)*diag(dim(xa)[2]))%*%t(xa)
        eta = eta[do_inference,,drop = F]
      }else{
        xa = x[, active, drop = F]
        eta = solve(crossprod(xa))%*%t(xa)  
      }
      
    }else{
      stop("Specified type does not exist.")
    }
    
    # normalize eta
    norm_eta = rep(NA, dim(eta)[1])
    for(jj in 1:dim(eta)[1]){
      norm_eta[jj] = sqrt(sum(eta[jj,]^2))
      eta[jj,] = eta[jj,]/norm_eta[jj]
    }
    
    
    
    ci = matrix(NA, nr = dim(x)[2], nc = 2)
    intervals_list = list()
    active = which(beta!=0)
    for(jj in 1:num_active){
      if(type == "full"){
        if(trunc_type == "TG_V"){
          marg_int = censored_marginal_intervals(y=y, x=x, eta = eta[jj,], jhat = active[jj],
                                                 beta = beta, sigma = sigma, lambda = lambda,
                                                 alpha = alpha)
        }else if(trunc_type == "TG_M"){
          marg_int = censored_marginal_intervals_defunct(y=y, x=x, eta = eta[jj,], jhat = active[jj],
                                                         beta = beta, sigma = sigma, lambda = lambda,
                                                         alpha = alpha)
        }else{
          stop("Specified trunc_type does not exist.")
        }
      }else{
        if(cluster_by_cor){
          marg_int = censored_marginal_intervals_cluster(y=y, x=x, eta = eta[jj,], jhat = active[jj],
                                                         beta = beta, sigma = sigma, lambda = lambda,
                                                         alpha = alpha, cor_thresh = cor_thresh)
        }else{
          marg_int = censored_marginal_intervals_defunct(y=y, x=x, eta = eta[jj,], jhat = active[jj],
                                                         beta = beta, sigma = sigma, lambda = lambda,
                                                         alpha = alpha)
        }
      }
      
      ci[active[jj],] = marg_int$ci*norm_eta[jj]
      intervals_list[[jj]] = marg_int$int
    }
    
    return(list(ci = ci, intervals_list = intervals_list))
}

get_lee_intervals = function(y, x, beta, sigma, lambda, alpha, type = "full"){
  
  active = (beta!=0)
  num_active = sum(active)
  if(num_active==0) stop("No non-zero coefficients")
  sob = sign(beta[active])
  
  if(type == "full"){
    eta = solve(crossprod(x))%*%t(x)
    eta = eta[active,,drop = F]
  }else if(type== "partial"){
      xa = x[, active, drop = F]
      eta = solve(crossprod(xa))%*%t(xa)  
  }else{
    stop("Specified type does not exist.")
  }
  
  # normalize eta
  norm_eta = rep(NA, dim(eta)[1])
  for(jj in 1:dim(eta)[1]){
    norm_eta[jj] = sqrt(sum(eta[jj,]^2))
    eta[jj,] = eta[jj,]/norm_eta[jj]
  }
  
  ci = matrix(NA, nr = dim(x)[2], nc = 2)
  intervals_list = list()
  active = which(beta!=0)
  for(jj in 1:num_active){
    marg_int = lee_marginal_intervals(y=y, x=x, eta = eta[jj,], jhat = active[jj],
                      beta = beta, sigma = sigma, lambda = lambda, alpha = alpha)
    ci[active[jj],] = marg_int$ci*norm_eta[jj]
    intervals_list[[jj]] = marg_int$int
  }
  
  return(list(ci = ci, intervals_list = intervals_list))
}


lee_marginal_intervals = function(y, x, eta, jhat, beta, sigma, lambda, alpha){
  eta = as.vector(eta)
  nn = length(y)
  
  cc = eta/sum(eta^2)
  z0 = sum(eta*y)
  nu = y - cc*z0
  
  active = (beta!=0)
  sob = sign(beta[active])
  
  poly_constraints = compute_poly(x = x, active = active, s = sob, lambda = lambda, active_only = active_only)
  
  endpoints = find_endpoints(y, eta, poly_constraints$A, poly_constraints$b,
                             poly_constraints$var_key, poly_constraints$sign_key)
  int = endpoints$int
  intervals_mat = matrix(int,1,2)
  
  ci= create_tnorm_interval(y, eta, sigma, alpha, intervals = intervals_mat)
  return(list(ci = ci, int = intervals_mat))
}

censored_marginal_intervals = function(y, x, eta, jhat, beta, sigma, lambda, alpha){
  
  eta = as.vector(eta)
  nn = length(y)
  
  cc = eta/sum(eta^2) 
  z0 = sum(eta*y)
  nu = y - cc*z0
  
  scale_factor = sum(x[,jhat]*cc)
  
  fit_lasso = glmnet(x = x[,-jhat] , y = nu, standardize=F, intercept = F, thresh = 1e-12,
                     lambda = lambda/nn)
  beta_star = as.vector(fit_lasso$beta)
  
  dd = sum(x[,jhat]*(x[,-jhat]%*%beta_star - nu))
  int = c(dd - lambda, dd + lambda)/scale_factor
  intervals_mat = rbind(c(-Inf,int[1]), c(int[2],Inf))
  
  ci= create_tnorm_interval(y, eta, sigma, alpha, intervals = intervals_mat)
  
  return(list(ci = ci, int = intervals_mat))
}

censored_marginal_intervals_defunct = function(y, x, eta, jhat, beta, sigma, lambda, alpha,
                                        zlow = -50, zhigh = 50, 
                                       digits = 5, active_only = F, max_attempt = 50){
  
  eta = as.vector(eta)
  nn = length(y)
  
  cc = eta/sum(eta^2)
  z0 = sum(eta*y)
  nu = y - cc*z0
  
  active = (beta!=0)
  sob = sign(beta[active])
  
  poly_constraints = compute_poly(x = x, active = active, s = sob, lambda = lambda, active_only = active_only)
  
  endpoints = find_endpoints(y, eta, poly_constraints$A, poly_constraints$b,
                             poly_constraints$var_key, poly_constraints$sign_key)
  int = endpoints$int
  if(is.na(int[1]) | int[2] < int[1]){
      stop(paste("Interval endpoints found for observed sign vector are not compatible.",
                 "Check to make sure intercept = F in glmnet or increase glmnet precision.", sep = ""))
  }
  
  intervals_mat = matrix(int,1,2)
  zlower = int[1]
  zupper = int[2]

  var_upper = endpoints$var_right
  var_lower = endpoints$var_left
  
  sign_upper = endpoints$sign_right
  sign_lower = endpoints$sign_left
  
  beta_star = beta
  attempt = 0
  while((zupper < zhigh) & (attempt < max_attempt)){
    # ystar = cc*(zupper + stepsize) + nu
    # fit_lasso = glmnet(x = x , y = ystar, standardize=F, intercept = F, thresh = 1e-12,
    #                    lambda = lambda/nn)
    # beta_star = as.vector(fit_lasso$beta)
    
    beta_star[var_upper] = sign_upper
    active_star = (beta_star!=0)
    s_star = sign(beta_star[active_star])
    
    poly_constraints = compute_poly(x = x, active = active_star, 
                                    s = s_star, lambda = lambda, active_only = active_only)
    endpoints = find_endpoints(y, eta, poly_constraints$A, poly_constraints$b,
                         poly_constraints$var_key, poly_constraints$sign_key)
    int = endpoints$int
    if(is.na(int[1]) | int[2] < int[1]){
      stop(paste("Interval endpoints found for observed sign vector are not compatible.",
                 "Check to make sure intercept = F in glmnet or increase glmnet precision.", sep = ""))
    }
    zupper = int[2]
    
    # the selected variable remains active in the new interval
    #if(active_star[jhat]) intervals_mat = rbind(intervals_mat,int)
    if(sum(active_star != active) == 0) intervals_mat = rbind(intervals_mat,int)
   
    var_upper = endpoints$var_right
    sign_upper = endpoints$sign_right
   
    attempt = attempt + 1
  }
 
  beta_star = beta
  attempt = 0
  while((zlower > zlow) & (attempt < max_attempt)){
    # ystar = cc*(zlower - stepsize) + nu
    # fit_lasso = glmnet(x = x , y = ystar, standardize=F, intercept = F, thresh = 1e-12,
    #                    lambda = lambda/nn)
    # beta_star = as.vector(fit_lasso$beta)
    
    beta_star[var_lower] = sign_lower
    active_star = (beta_star!=0)
    s_star = sign(beta_star[active_star])
    poly_constraints = compute_poly(x = x, active = active_star, 
                                    s = s_star, lambda = lambda, active_only = active_only)
    endpoints = find_endpoints(y, eta, poly_constraints$A, poly_constraints$b,
                         poly_constraints$var_key, poly_constraints$sign_key)
    int = endpoints$int
    if(is.na(int[1]) | int[2] < int[1]){
      stop(paste("Interval endpoints found for observed sign vector are not compatible.",
                 "Check to make sure intercept = F in glmnet or increase glmnet precision.", sep = ""))
    }
    zlower = int[1]
    
    # the selected variable remains active in the new interval
    #if(active_star[jhat]) intervals_mat = rbind(intervals_mat,int)
    if(sum(active_star != active) == 0) intervals_mat = rbind(intervals_mat,int)
    
    var_lower = endpoints$var_left
    sign_lower = endpoints$sign_left
    
    attempt = attempt + 1
  }
  
  # take union of intervals
  intervals_mat = round(intervals_mat, digits=  digits)
  intervals_mat = Intervals(intervals_mat)
  intervals_mat = interval_union(intervals_mat)
  intervals_mat = as.matrix(intervals_mat)
  
  ci= create_tnorm_interval(y, eta, sigma, alpha, intervals = intervals_mat)
  
  return(list(ci = ci, int = intervals_mat))
}

censored_marginal_intervals_cluster = function(y, x, eta, jhat, beta, sigma, lambda, alpha, cor_thresh,
                                               Sigmat = NULL, zlow = -100, zhigh = 100, stepsize = 0.5,
                                               digits = 5, active_only = F, max_attempt = 100){
  
  cluster_x = hclust(as.dist(1-abs(cor(x))))
  cluster_x = cutree(cluster_x, h = 1 - cor_thresh)
  active_cluster = rep(F, length(unique(cluster_x)))
  active_cluster[unique(cluster_x[beta!=0])] = T
  
  eta = as.vector(eta)
  nn = length(y)
  
  cc = eta/sum(eta^2) 
  z0 = sum(eta*y)
  nu = y - cc*z0
  
  active = (beta!=0)
  sob = sign(beta[active])
  
  poly_constraints = compute_poly(x = x, active = active, s = sob, lambda = lambda, active_only = active_only)
  
  int = find_endpoints(y, eta, poly_constraints$A, poly_constraints$b,
                       poly_constraints$var_key, poly_constraints$sign_key)$int
  if(is.na(int[1]) | int[2] < int[1]){
    stop(paste("Interval endpoints found for observed sign vector are not compatible.",
               "Check to make sure intercept = F in glmnet or increase glmnet precision.", sep = ""))
  }
  
  intervals_mat = matrix(int,1,2)
  zlower = int[1]
  zupper = int[2]
  
  attempt = 0
  while((zupper < zhigh) & attempt < max_attempt){
    ystar = cc*(zupper + stepsize) + nu
    fit_lasso = glmnet(x = x , y = ystar, standardize=F, intercept = F, thresh = 1e-12,
                       lambda = lambda/nn)
    beta_star = as.vector(fit_lasso$beta)
    active_star = (beta_star!=0)
    s_star = sign(beta_star[active_star])
    poly_constraints = compute_poly(x = x, active = active_star, 
                                    s = s_star, lambda = lambda, active_only = active_only)
    int = find_endpoints(y, eta, poly_constraints$A, poly_constraints$b,
                         poly_constraints$var_key, poly_constraints$sign_key)$int
    if(is.na(int[1]) | int[2] < int[1]){
      stop(paste("Interval endpoints found for observed sign vector are not compatible.",
                 "Check to make sure intercept = F in glmnet or increase glmnet precision.", sep = ""))
    }
    zupper = int[2]
    int
    # the selected variable remains active in the new interval
    active_cluster_star = rep(F, length(unique(cluster_x)))
    if(sum(active_star) > 0) active_cluster_star[unique(cluster_x[active_star])] = T
    if(sum(active_cluster_star != active_cluster) == 0) intervals_mat = rbind(intervals_mat,int)
    
    attempt = attempt + 1
  }
  
  attempt = 0
  while((zlower > zlow) & (attempt < max_attempt)){
    ystar = cc*(zlower - stepsize) + nu
    fit_lasso = glmnet(x = x , y = ystar, standardize=F, intercept = F, thresh = 1e-12,
                       lambda = lambda/nn)
    beta_star = as.vector(fit_lasso$beta)
    active_star = (beta_star!=0)
    s_star = sign(beta_star[active_star])
    poly_constraints = compute_poly(x = x, active = active_star, 
                                    s = s_star, lambda = lambda, active_only = active_only)
    int = find_endpoints(y, eta, poly_constraints$A, poly_constraints$b,
                         poly_constraints$var_key, poly_constraints$sign_key)$int
    if(is.na(int[1]) | int[2] < int[1]){
      stop(paste("Interval endpoints found for observed sign vector are not compatible.",
                 "Check to make sure intercept = F in glmnet or increase glmnet precision.", sep = ""))
    }
    zlower = int[1]
    
    # the selected variable remains active in the new interval
    active_cluster_star = rep(F, length(unique(cluster_x)))
    if(sum(active_star) > 0) active_cluster_star[unique(cluster_x[active_star])] = T
    if(sum(active_cluster_star != active_cluster) == 0) intervals_mat = rbind(intervals_mat,int)
    
    attempt = attempt + 1
  }
  
  # take union of intervals
  intervals_mat = round(intervals_mat, digits=  digits)
  intervals_mat = Intervals(intervals_mat)
  intervals_mat = interval_union(intervals_mat)
  intervals_mat = as.matrix(intervals_mat)
  
  ci= create_tnorm_interval(y, eta, sigma, alpha, intervals = intervals_mat)
  
  return(list(ci = ci, int = intervals_mat))
}

thresh_intervals = function(y, x, beta, sigma, lambda, alpha, type = "full",
                              Sigmat = NULL, beta_thresh = 0.2){
  
  active = (beta!=0)
  num_active = sum(active)
  if(num_active==0) stop("No non-zero coefficients")
  sob = sign(beta[active])
  
  if(type == "full"){
    eta = solve(crossprod(x))%*%t(x)
    eta = eta[active,,drop = F]
  }else if(type== "partial"){
      act_thresh = which(abs(beta) > beta_thresh)
      if(length(act_thresh) > 0){
        eta = matrix(NA, num_active, length(y))
        xa = x[, act_thresh, drop = F]
        Pa = (diag(length(y)) - xa%*%solve(crossprod(xa))%*%t(xa))
        eta = t(Pa%*%x[,active, drop = F])
        vars_thresh = which(abs(beta[active]) > beta_thresh)
        eta[vars_thresh,] = solve(crossprod(xa))%*%t(xa)
      }else{
        eta = t(x[,active,drop = F]) 
      }
    
  }else{
    stop("Specified type does not exist.")
  }
  
  # normalize eta
  norm_eta = rep(NA, dim(eta)[1])
  for(jj in 1:dim(eta)[1]){
    norm_eta[jj] = sqrt(sum(eta[jj,]^2))
    eta[jj,] = eta[jj,]/norm_eta[jj]
  }
  
  
  
  ci = matrix(NA, nr = dim(x)[2], nc = 2)
  intervals_list = list()
  active = which(beta!=0)
  for(jj in 1:num_active){
    if(type == "full"){
      marg_int = censored_marginal_intervals(y=y, x=x, eta = eta[jj,], jhat = active[jj],
                                             beta = beta, sigma = sigma, lambda = lambda,
                                             alpha = alpha)
    }else{
      marg_int = thresh_marginal_intervals(y=y, x=x, eta = eta[jj,], jhat = active[jj],
                                                       beta = beta, sigma = sigma, lambda = lambda,
                                                       alpha = alpha, beta_thresh = beta_thresh)
    }
    
    ci[active[jj],] = marg_int$ci*norm_eta[jj]
    intervals_list[[jj]] = marg_int$int
  }
  
  return(list(ci = ci, intervals_list = intervals_list))
}

thresh_marginal_intervals = function(y, x, eta, jhat, beta, sigma, lambda, alpha, beta_thresh = 0.2,
                                               zlow = -100, zhigh = 100, stepsize = 0.5,
                                               digits = 5, active_only = F, max_attempt = 100){
  
  eta = as.vector(eta)
  nn = length(y)
  
  cc = eta/sum(eta^2)
  z0 = sum(eta*y)
  nu = y - cc*z0
  
  active = (beta!=0)
  active_thresh = (abs(beta) > beta_thresh)
  sob = sign(beta[active])
  
  poly_constraints = compute_poly(x = x, active = active, s = sob, lambda = lambda, active_only = active_only)
  
  int = find_endpoints(y, eta, poly_constraints$A, poly_constraints$b,
                       poly_constraints$var_key, poly_constraints$sign_key)$int
  if(is.na(int[1]) | int[2] < int[1]){
    stop(paste("Interval endpoints found for observed sign vector are not compatible.",
               "Check to make sure intercept = F in glmnet or increase glmnet precision.", sep = ""))
  }
  
  intervals_mat = matrix(int,1,2)
  zlower = int[1]
  zupper = int[2]
  
  attempt = 0
  while((zupper < zhigh) & (attempt < max_attempt)){
    ystar = cc*(zupper + stepsize) + nu
    fit_lasso = glmnet(x = x , y = ystar, standardize=F, intercept = F, thresh = 1e-12,
                       lambda = lambda/nn)
    beta_star = as.vector(fit_lasso$beta)
    active_star = (beta_star!=0)
    active_thresh_star = (abs(beta_star) > beta_thresh)
    s_star = sign(beta_star[active_star])
    poly_constraints = compute_poly(x = x, active = active_star, 
                                    s = s_star, lambda = lambda, active_only = active_only)
    int = find_endpoints(y, eta, poly_constraints$A, poly_constraints$b,
                         poly_constraints$var_key, poly_constraints$sign_key)$int
    if(is.na(int[1]) | int[2] < int[1]){
      stop(paste("Interval endpoints found for observed sign vector are not compatible.",
                 "Check to make sure intercept = F in glmnet or increase glmnet precision.", sep = ""))
    }
    zupper = int[2]
    
    # the selected variable remains active in the new interval
    #if(active_star[jhat]) intervals_mat = rbind(intervals_mat,int)
    if((sum(active_thresh_star != active_thresh) == 0) & active_star[jhat]) intervals_mat = rbind(intervals_mat,int)
    
    attempt = attempt + 1
  }
  
  attempt = 0
  while((zlower > zlow) & (attempt < max_attempt)){
    ystar = cc*(zlower - stepsize) + nu
    fit_lasso = glmnet(x = x , y = ystar, standardize=F, intercept = F, thresh = 1e-12,
                       lambda = lambda/nn)
    beta_star = as.vector(fit_lasso$beta)
    active_star = (beta_star!=0)
    active_thresh_star = (abs(beta_star) > beta_thresh)
    s_star = sign(beta_star[active_star])
    poly_constraints = compute_poly(x = x, active = active_star, 
                                    s = s_star, lambda = lambda, active_only = active_only)
    int = find_endpoints(y, eta, poly_constraints$A, poly_constraints$b,
                         poly_constraints$var_key, poly_constraints$sign_key)$int
    if(is.na(int[1]) | int[2] < int[1]){
      stop(paste("Interval endpoints found for observed sign vector are not compatible.",
                 "Check to make sure intercept = F in glmnet or increase glmnet precision.", sep = ""))
    }
    zlower = int[1]
    
    # the selected variable remains active in the new interval
    #if(active_star[jhat]) intervals_mat = rbind(intervals_mat,int)
    if((sum(active_thresh_star != active_thresh) == 0) & active_star[jhat]) intervals_mat = rbind(intervals_mat,int)
    
    attempt = attempt + 1
  }
  
  # take union of intervals
  intervals_mat = round(intervals_mat, digits=  digits)
  intervals_mat = Intervals(intervals_mat)
  intervals_mat = interval_union(intervals_mat)
  intervals_mat = as.matrix(intervals_mat)
  
  ci= create_tnorm_interval(y, eta, sigma, alpha, intervals = intervals_mat)
  
  return(list(ci = ci, int = intervals_mat))
}

compute_zstat = function(y, x, active, sigma){
  xa = x[,active,drop = F]
  eta = xa%*%solve(crossprod(xa))
  eta = scale(eta, F, scale = sqrt(apply(eta^2,2,sum)))
  tstat = abs(as.vector(t(eta)%*%y/sigma))
  return(tstat)
}

form_ls_thresh_target = function(x, ind_active, tstat, t_thresh = 1){

  ind_thresh = which(tstat > t_thresh)
  act_thresh = ind_active[ind_thresh]
  eta = matrix(NA, length(ind_active), dim(x)[1])
  
  if(length(act_thresh) > 0){
    xa = x[, act_thresh, drop = F]
    Pa = (diag(dim(x)[1]) - xa%*%solve(crossprod(xa))%*%t(xa))
    eta = Pa%*%x[, ind_active, drop = F]
    eta = scale(eta, F, scale = apply(eta^2,2,sum))
    eta = t(eta)
    eta[ind_thresh,] = solve(crossprod(xa))%*%t(xa)
  }else{
    eta = x[,ind_active, drop = F]
    eta = scale(eta, F, scale = apply(eta^2,2,sum))
    eta = t(eta) 
  }
  
  return(eta)
}

ls_thresh_intervals = function(y, x, beta, sigma, lambda, alpha, type = "full",
                            t_thresh = 1){
  
  active = (beta!=0)
  num_active = sum(active)
  if(num_active==0) stop("No non-zero coefficients")
  sob = sign(beta[active])
  
  if(type == "full"){
    eta = solve(crossprod(x))%*%t(x)
    eta = eta[active,,drop = F]
    tstat=NA
  }else if(type== "partial"){
    tstat = compute_zstat(y, x, active, sigma)
    eta = form_ls_thresh_target(x, ind_active = which(beta!=0), tstat = tstat, t_thresh = t_thresh)
  }else{
    stop("Specified type does not exist.")
  }
  
  # normalize eta
  norm_eta = rep(NA, dim(eta)[1])
  for(jj in 1:dim(eta)[1]){
    norm_eta[jj] = sqrt(sum(eta[jj,]^2))
    eta[jj,] = eta[jj,]/norm_eta[jj]
  }
  
  ci = matrix(NA, nr = dim(x)[2], nc = 2)
  intervals_list = list()
  active = which(beta!=0)
  for(jj in 1:num_active){
    if(type == "full"){
      marg_int = censored_marginal_intervals(y=y, x=x, eta = eta[jj,], jhat = active[jj],
                                             beta = beta, sigma = sigma, lambda = lambda,
                                             alpha = alpha)
    }else{
      active_thresh = which(beta!=0)[which(tstat > t_thresh)]
      marg_int = ls_thresh_marginal_intervals(y=y, x=x, eta = eta[jj,], jhat = active[jj],
                                           beta = beta, sigma = sigma, lambda = lambda,
                                           alpha = alpha, active_thresh = active_thresh, t_thresh = t_thresh)
    }
    
    ci[active[jj],] = marg_int$ci*norm_eta[jj]
    intervals_list[[jj]] = marg_int$int
  }
  
  return(list(ci = ci, intervals_list = intervals_list, high_value = active[which(tstat > t_thresh)]))
}

ls_thresh_marginal_intervals = function(y, x, eta, jhat, beta, sigma, lambda, alpha, active_thresh, t_thresh = 1,
                                     zlow = -100, zhigh = 100, stepsize = 0.01,
                                     digits = 5, active_only = F, max_attempt = 100){
  
  eta = as.vector(eta)
  nn = length(y)
  
  cc = eta/sum(eta^2)
  z0 = sum(eta*y)
  nu = y - cc*z0
  
  active = (beta!=0)
  sob = sign(beta[active])
  
  poly_constraints = compute_poly(x = x, active = active, s = sob, lambda = lambda, active_only = active_only)
  
  endpoints = find_endpoints(y, eta, poly_constraints$A, poly_constraints$b,
                       poly_constraints$var_key, poly_constraints$sign_key)
  int = endpoints$int
  if(is.na(int[1]) | int[2] < int[1]){
    stop(paste("Interval endpoints found for observed sign vector are not compatible.",
               "Check to make sure intercept = F in glmnet or increase glmnet precision.", sep = ""))
  }
  
  # find subinterval over which significant set is stable
  active = which(beta!=0)
  ls_int = find_ls_interval(nu, cc, x, active = active, significant = active_thresh, zstat_thresh = t_thresh, sigma = sigma)
  
  intervals_mat = interval_intersection(ls_int, Intervals(int))  
  
  zlower = int[1]
  zupper = int[2]
  
  var_upper = endpoints$var_right
  var_lower = endpoints$var_left
  
  sign_upper = endpoints$sign_right
  sign_lower = endpoints$sign_left
  
  beta_star = beta
  
  attempt = 0
  while((zupper < zhigh) & (attempt < max_attempt)){
    # ystar = cc*(zupper + stepsize) + nu
    # fit_lasso = glmnet(x = x , y = ystar, standardize=F, intercept = F, thresh = 1e-12,
    #                    lambda = lambda/nn)
    # beta_star = as.vector(fit_lasso$beta)
    # tmp = as.vector(fit_lasso$beta)
    
    beta_star[var_upper] = sign_upper
    active_star = (beta_star!=0)
    s_star = sign(beta_star[active_star])
    
    # if(sum(active_star != (tmp!=0)) > 0){
    #   stop("predicted active set different from that found by glmnet")
    # }
    # if(sum(s_star != sign(tmp[tmp!=0])) > 0){
    #   stop("predicted signs incorrect")
    # }
    
    poly_constraints = compute_poly(x = x, active = active_star, 
                                    s = s_star, lambda = lambda, active_only = active_only)
    endpoints = find_endpoints(y, eta, poly_constraints$A, poly_constraints$b,
                         poly_constraints$var_key, poly_constraints$sign_key)
    int = endpoints$int
    if(is.na(int[1]) | int[2] < int[1]){
      stop(paste("Interval endpoints found for observed sign vector are not compatible.",
                 "Check to make sure intercept = F in glmnet or increase glmnet precision.", sep = ""))
    }
    zupper = int[2]
    
    ls_int = find_ls_interval(nu, cc, x, active = which(beta_star!=0), significant = active_thresh, zstat_thresh = t_thresh, sigma = sigma)
   
    # the selected variable remains active in the new interval
    if(active_star[jhat]){
      intervals_mat = interval_union(intervals_mat, interval_intersection(Intervals(int),ls_int))
    }

    var_upper = endpoints$var_right
    sign_upper = endpoints$sign_right
    
    attempt = attempt + 1
  }
  
  beta_star = beta
  
  attempt = 0
  while((zlower > zlow) & (attempt < max_attempt)){
    # ystar = cc*(zlower - stepsize) + nu
    # fit_lasso = glmnet(x = x , y = ystar, standardize=F, intercept = F, thresh = 1e-12,
    #                    lambda = lambda/nn)
    # # beta_star = as.vector(fit_lasso$beta)
    # tmp = as.vector(fit_lasso$beta)
    
    beta_star[var_lower] = sign_lower
    active_star = (beta_star!=0)
    s_star = sign(beta_star[active_star])
    
    # if(sum(active_star != (tmp!=0)) > 0){
    #   stop("predicted active set different from that found by glmnet")
    # }
    # if(sum(s_star != sign(tmp[tmp!=0])) > 0){
    #   stop("predicted signs incorrect")
    # }
    
    poly_constraints = compute_poly(x = x, active = active_star, 
                                    s = s_star, lambda = lambda, active_only = active_only)
    endpoints = find_endpoints(y, eta, poly_constraints$A, poly_constraints$b,
                         poly_constraints$var_key, poly_constraints$sign_key)
    int = endpoints$int
    if(is.na(int[1]) | int[2] < int[1]){
      stop(paste("Interval endpoints found for observed sign vector are not compatible.",
                 "Check to make sure intercept = F in glmnet or increase glmnet precision.", sep = ""))
    }
    zlower = int[1]
    
    ls_int = find_ls_interval(nu, cc, x, active = which(beta_star!=0), significant = active_thresh, zstat_thresh = t_thresh, sigma = sigma)
    
    # the selected variable remains active in the new interval
    if(active_star[jhat]){
      intervals_mat = interval_union(intervals_mat, interval_intersection(Intervals(int),ls_int))
    }
    
    var_lower = endpoints$var_left
    sign_lower = endpoints$sign_left
    
    attempt = attempt + 1
  }
  
  # take union of intervals
  intervals_mat = round(intervals_mat, digits=  digits)
  intervals_mat = interval_union(intervals_mat)
  intervals_mat = as.matrix(intervals_mat)
  
  ci= create_tnorm_interval(y, eta, sigma, alpha, intervals = intervals_mat)
  
  return(list(ci = ci, int = intervals_mat))
}


find_ls_interval = function(nu, cc, x, active, significant, zstat_thresh, sigma){
  
  # if a significant variable has been excluded from the active set, interval is null
  if(length(setdiff(significant, active)) >0){
    return(Intervals())
  }
  
  if(length(active) == 0) return(Intervals())
  
  xa = x[,active,drop = F]
  eta = xa%*%solve(crossprod(xa))
  eta = scale(eta, F, scale = sqrt(apply(eta^2,2,sum)))
  
  aa = t(eta)%*%nu
  bb = t(eta)%*%cc
  
  sig_id = c()
  if(length(significant)==0){
    not_sig = seq(1:length(active))
  }else{
    for(kk in 1:length(significant)){
      sig_id = c(sig_id, which(active == significant[kk]))
    }
    not_sig = setdiff(seq(1:length(active)), sig_id)
  }
  
  zupper = (-sign(bb)*aa + zstat_thresh*sigma)/abs(bb)
  zlower = (-sign(bb)*aa - zstat_thresh*sigma)/abs(bb)
  
  # find bounds for non-significant variables to remain non-significant
  int_mat = Intervals(c(-Inf,Inf))
  if(length(not_sig) > 0){
    MM = max(zlower[not_sig]) 
    mm = min(zupper[not_sig])
    if(MM >= mm){
      return(Intervals())
    }else{
      int_mat = matrix(c(MM,mm) ,1,2)
      int_mat = Intervals(int_mat)  
    }
  }
  
  # find bounds for significant variables to remain significant
  if(length(sig_id) > 0){
    for(kk in 1:length(sig_id)){
      feasible = Intervals(rbind(c(-Inf, zlower[sig_id[kk]]), c(zupper[sig_id[kk]], Inf)))
      int_mat = interval_intersection(int_mat, feasible)
    }
  }
  
  return(int_mat)
}

# poly_constraints = compute_poly(x = x, active = beta != 0, s = c(1,-1,1), lambda = n*lam)
# find_endpoints(y, (solve(t(x[,beta!=0])%*%x[,beta!=0])%*%t(x[,beta!=0]))[1,], poly_constraints$A, 
#                poly_constraints$b)
# tmp = random_sign_intervals(y, x, beta, sigma, n*lam, alpha, num_samples = 20, type = "partial")
