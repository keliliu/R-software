
form_lasso_thresh_target = function(x, beta, beta_high){
  
  ind_active = which(beta != 0)
  act_thresh = which(beta_high != 0)
  
  eta = matrix(NA, length(ind_active), dim(x)[1])
  if(length(act_thresh) > 0){
    xa = x[, act_thresh, drop = F]
    Pa = (diag(dim(x)[1]) - xa%*%solve(crossprod(xa))%*%t(xa))
    eta = Pa%*%x[, ind_active, drop = F]
    eta = scale(eta, F, scale = apply(eta^2,2,sum))
    eta = t(eta)
  
    ind_thresh = which(beta_high[ind_active] != 0 )
    pass_low = which(beta[act_thresh]!=0)
    eta[ind_thresh,] = (solve(crossprod(xa))%*%t(xa))[pass_low,]
  }else{
    eta = x[, ind_active, drop = F] 
    eta = scale(eta, F, scale = apply(eta^2,2,sum))
    eta = t(eta)
  }
  
  return(eta)
}

lasso_thresh_intervals = function(y, x, beta, beta_high, sigma, lambda, lambda_high, 
                               alpha){
  
  beta = correct_lasso_estimate(x, beta, lambda)
  beta_high = correct_lasso_estimate(x, beta_high, lambda_high)
  
  active = (beta!=0)
  num_active = sum(active)
  if(num_active==0) stop("No non-zero coefficients")
  sob = sign(beta[active])
  
  eta = form_lasso_thresh_target(x, beta, beta_high)

  # normalize eta
  norm_eta = rep(NA, dim(eta)[1])
  for(jj in 1:dim(eta)[1]){
    norm_eta[jj] = sqrt(sum(eta[jj,]^2))
    eta[jj,] = eta[jj,]/norm_eta[jj]
  }
  
  ci = matrix(NA, nr = dim(x)[2], nc = 2)
  intervals_list = list()
  active = which(beta!=0)
  active_thresh = which(beta_high!=0)
  for(jj in 1:num_active){
    marg_int = lasso_thresh_marginal_intervals(y=y, x=x, eta = eta[jj,], jhat = active[jj],
                                              beta = beta, beta_high = beta_high, sigma = sigma, 
                                              lambda = lambda, lambda_high = lambda_high, alpha = alpha)
    ci[active[jj],] = marg_int$ci*norm_eta[jj]
    intervals_list[[jj]] = marg_int$int
  }
  
  return(list(ci = ci, intervals_list = intervals_list, high_value = active_thresh))
}

find_truncation = function(y, x, eta, beta, lambda, interval_test, test_param, zlow = -50, zhigh = 50, max_attempt = 50, 
                           digits = 5, active_only = F){
  
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
  
  intervals_mat = Intervals(int)
  
  zlower = int[1]
  zupper = int[2]
  
  var_upper = endpoints$var_right
  var_lower = endpoints$var_left
  
  sign_upper = endpoints$sign_right
  sign_lower = endpoints$sign_left
  
  beta_star = beta
  attempt = 0
  while((zupper < zhigh) & (attempt < max_attempt)){
    
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
    
    if(interval_test(beta_star, test_param)){
      intervals_mat = interval_union(intervals_mat, Intervals(int))
    }
    
    var_upper = endpoints$var_right
    sign_upper = endpoints$sign_right
    attempt = attempt + 1
  }
  
  beta_star = beta
  attempt = 0
  while((zlower > zlow) & (attempt < max_attempt)){
    
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
  
    if(interval_test(beta_star, test_param)){
      intervals_mat = interval_union(intervals_mat, Intervals(int))
    }
    
    var_lower = endpoints$var_left
    sign_lower = endpoints$sign_left
    attempt = attempt + 1
  }
  
  # take union of intervals
  intervals_mat = round(intervals_mat, digits=  digits)
  intervals_mat = interval_union(intervals_mat)
  
  return(intervals_mat)
}

test_jhat_active = function(beta, test_param){
  if(beta[test_param$jhat] != 0) return(TRUE)
  return(FALSE)
}

test_model_identity = function(beta, test_param){
  if(sum((beta !=0 ) != test_param$active) == 0) return(TRUE)
  return(FALSE)
}


lasso_thresh_marginal_intervals = function(y, x, eta, jhat, beta, beta_high, sigma = sigma, 
                                           lambda, lambda_high, alpha = alpha){
  
  eta = as.vector(eta)
  nn = length(y)
  
  cc = eta/sum(eta^2)
  z0 = sum(eta*y)
  nu = y - cc*z0
  
  # find interval over which jhat is active
  jhat_active_int = find_truncation(y, x, eta, beta, lambda, test_jhat_active, list(jhat = jhat))
  
  # find interval over which high value targets remains unchanged
  high_stable_int = find_truncation(y, x, eta, beta_high, lambda_high, test_model_identity, 
                                    list(active = (beta_high !=0)))
  
  intervals_mat = interval_intersection(jhat_active_int, high_stable_int)
  intervals_mat = as.matrix(intervals_mat)
  
  ci= create_tnorm_interval(y, eta, sigma, alpha, intervals = intervals_mat)
  
  return(list(ci = ci, int = intervals_mat))
}
