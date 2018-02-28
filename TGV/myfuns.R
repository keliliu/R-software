foo=function(x,y,beta,lam,alpha,omega=rep(0,length(y)),type="partial")  {
    p=ncol(x)
    ci=vpts=array(NA,c(p,2))
vars=beta!=0
k=sum(vars)
 gridrange=c(-100,100)
 bits=NULL
pinv=solve

out = fixedLasso.poly(x,y,beta,lam,vars)
  G = out$G
    u = out$u
    u=u-G%*%omega
 if(type=="full"){ 
M = pinv(crossprod(x)) %*% t(x)
    M = M[vars,,drop=F]
 }
    if(type=="partial"){
        xa=x[,vars]
M = pinv(crossprod(xa)) %*% t(xa)
   
 }
varsn=which(vars)
  
  ss = rep(NA,k)
  for (j in 1:k) {
    if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))

    vj = M[j,]
    mj = sqrt(sum(vj^2))
    vj = vj / mj        # Standardize (divide by norm of vj)
   ss[j] = sign(sum(vj*y))
    vj = ss[j] * vj
  #  a = poly.pval(y,G,u,vj,sigma,bits)

 a = mypoly.int(y,G,u,vj,sigma,alpha,gridrange=gridrange,
      flip=(ss[j]==-1),bits=bits)
    ci[varsn[j],] = a$int * mj # Unstandardize (mult by norm of vj)
    vpts[varsn[j],]=c(a$vlo,a$vup)
  }
return(list(ci=ci,vpts=vpts))
    }

mypoly.int=
function(y, G, u, v, sigma, alpha, gridrange=c(-100,100),
                     gridpts=100, griddepth=2, flip=FALSE, bits=NULL) {
  
  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  
  rho = G %*% v / vv
  vec = (u - G %*% y + rho*z) / rho
  vlo = suppressWarnings(max(vec[rho>0]))
  vup = suppressWarnings(min(vec[rho<0]))
  
  xg = seq(gridrange[1]*sd,gridrange[2]*sd,length=gridpts)
  fun = function(x) { tnorm.surv(z,x,sd,vlo,vup,bits) }

  int = grid.search(xg,fun,alpha/2,1-alpha/2,gridpts,griddepth)
  tailarea = c(fun(int[1]),1-fun(int[2]))

  if (flip) {
    int = -int[2:1]
    tailarea = tailarea[2:1]
  }
  
  return(list(int=int,tailarea=tailarea,vlo=vlo,vup=vup))
}

foo2=
    function(x,y,beta,lam,alpha,omega=rep(0,length(y)),type="partial",vlims=NULL,nreps=NULL)  {
        #get interavls from mix of TN
    p=ncol(x)
    ci=v=array(NA,c(p,2))
vars=beta!=0
k=sum(vars)
 gridrange=c(-100,100)
 bits=NULL
pinv=solve
     
 if(type=="full"){ 
M = pinv(crossprod(x)) %*% t(x)
    M = M[vars,,drop=F]
 }
if(type=="partial"){
        xa=x[,vars]
M = pinv(crossprod(xa)) %*% t(xa)
   
 }
varsn=which(vars)

  for (j in 1:k) {
    if (verbose) cat(sprintf("Inference for variable %i ...\\n",vars[j]))

    vj = M[j,]
    mj = sqrt(sum(vj^2))
    vj = vj / mj        # Standardize (divide by norm of vj)
    sign[j] = sign(sum(vj*y))
    vj = sign[j] * vj
  #  a = poly.pval(y,G,u,vj,sigma,bits)

 a = mypoly.int2(y,G,u,vj,sigma,alpha,gridrange=gridrange,
      flip=(sign[j]==-1),bits=bits,vlims=vlims[varsn[j],],nrep=nrep)
    ci[varsn[j],] = a$int * mj # Unstandardize (mult by norm of vj)
  
  }
return(ci)
    }

foo3=function(x,y,beta,lam,alpha,signs,type="partial")  {
    p=ncol(x)
    ci=vpts=array(NA,c(p,2))
vars=beta!=0
k=sum(vars)
 gridrange=c(-100,100)
 bits=NULL
    pinv=solve
   sign=rep(NA,length(vars))

out = fixedLasso.poly(x,y,beta,lam,vars,signs)
  G = out$G
      #   u = out$u*signs
    u=out$u
  
 if(type=="full"){ 
M = pinv(crossprod(x)) %*% t(x)
    M = M[vars,,drop=F]
 }
    if(type=="partial"){
        xa=x[,vars]
M = pinv(crossprod(xa)) %*% t(xa)
   
 }
varsn=which(vars)

  for (j in 1:k) {
    if (verbose) cat(sprintf("Inference for variable %i ...\n",vars[j]))

    vj = M[j,] # eta for the jth coefficient
    mj = sqrt(sum(vj^2))
    vj = vj / mj        # Standardize (divide by norm of vj)
    sign[j] = sign(sum(vj*y))
   
 vj = sign[j] * vj
  #  a = poly.pval(y,G,u,vj,sigma,bits)

 a = mypoly.int(y,G,u,vj,sigma,alpha,gridrange=gridrange,
                flip=(sign[j]==-1),bits=bits)
  
    ci[varsn[j],] = a$int * mj # Unstandardize (mult by norm of vj)
    vpts[varsn[j],]=c(a$vlo,a$vup)
  }
return(list(ci=ci,vpts=vpts))
}

mypoly.int2=
function(y, G, u, v, sigma, alpha, gridrange=c(-100,100),
                     gridpts=100, griddepth=2, flip=FALSE, bits=NULL,vlims=NULL,nrep=nrep) {
  
  z = sum(v*y)
  vv = sum(v^2)
  sd = sigma*sqrt(vv)
  
 # rho = G %*% v / vv
 # vec = (u - G %*% y + rho*z) / rho
#  vlo = suppressWarnings(max(vec[rho>0]))
#  vup = suppressWarnings(min(vec[rho<0]))
  
  xg = seq(gridrange[1]*sd,gridrange[2]*sd,length=gridpts)
  # fun = function(x) { tnorm.surv(z,x,sd,vlo,vup,bits) }
   # browser()
    fun=function(x){
        val=0
      
        for(ii in 1:nrep){
            vlo=vlims[1]+.1*sigma
             vup=vlims[2]+.1*sigma
            val=val+tnorm.surv(z,x,sd,vlo,vup,bits)
        }
      
        return(val/nrep)
        }
#
  int = grid.search(xg,fun,alpha/2,1-alpha/2,gridpts,griddepth)
  tailarea = c(fun(int[1]),1-fun(int[2]))

  if (flip) {
    int = -int[2:1]
    tailarea = tailarea[2:1]
  }
  
  return(list(int=int,tailarea=tailarea))
}


relax=function(x,y,lam){
    p=ncol(x)
    n=length(y)
    if(p>1){
    a=glmnet(x,y,standardize=F)
    bhat=coef(a,s=lam,x=y,x=y)[-1]
    }
    if(p==1){
    coef=lsfit(x,y)$coef[2]
    bhat=sign(coef)*(abs(coef)-n*lam)*(abs(coef)>n*lam)
    }
    act=which(bhat!=0)
    bhat0=mean(y)
    bhat2=rep(0,length(bhat))
    if(sum(bhat!=0)>0){
        junk=lsfit(x[,bhat!=0,drop=F],y)
        bhat0=junk$coef[1]
        bhat2=junk$coef[-1]
        se=ls.diag(junk)$std.err[-1]
    }
    return(list(bhat0=bhat0,bhat=bhat2,act=act,se=se))
    }

cifunc=function(bhat,se,alpha){
    zalpha=-qnorm(alpha/2)
    p1=bhat-se*zalpha
    p2=bhat+se*zalpha
    return(c(p1,p2))
    }



bootfunc=function(x,y,lam,sigma,act,bhat0,bhat,nboot){
    #bhat is either  relaxed lasso  or debiased lasso \
    p=ncol(x)
    betastar=sestar=array(NA,c(p,nboot))
  
     muhat=rep(bhat0,nrow(x))+x[,act,drop=F]%*%bhat

for(kk in 1:nboot){
    ystar=muhat+sigma*rnorm(n)
    ystar=ystar-mean(ystar)
    a4=relax(x,ystar,lam)
    #if(sum(a4$bhat!=0)>0) {
        betastar[a4$act,kk]=a4$bhat
        sestar[a4$act,kk]=a4$se
       #}
    }
    return(list(betastar=betastar,sestar=sestar))
}

estalpha=function(bhat,act,betastar,sestar,alpha,type="full"){
    p=dim(betastar)[1]
    ci2=array(NA,c(p,2,nboot))
       alphalist=c(.0005, seq(.001,.009,by=.001),seq(.01,.1,by=.01))
    bhatfull=rep(0,p)
    if(type=="full") bhatfull[act]=bhat
    muhat=x[,act,drop=F]%*%bhat
        na=length(alphalist)

        mc=matrix(NA,p,na)
        for(k in 1:na){
            alpha2=alphalist[k]
            zalpha2=-qnorm(alpha2/2)
            for(kk in 1:nboot){
                if(type=="partial"){
                    actstar=!is.na(betastar[,kk])
                    bhatfull=rep(0,p)
                    bhatfull[actstar]=lsfit(x[,actstar,drop=F],muhat)$coef[-1]
                    }
          ci2[,1,kk]=betastar[,kk]-zalpha2*sestar[,kk]
                ci2[,2,kk]=betastar[,kk]+zalpha2*sestar[,kk]
                actstar=which(betastar[,kk]!=0)
          for(j in 1:length(actstar)){
              jj=actstar[j]
                mc[jj,k]=mean( (ci2[jj,1,]>bhatfull[jj] | ci2[jj,2,]< bhatfull[jj]),na.rm=T)
          }}}
        mcm=colMeans(mc,na.rm=T)
        o=which.min(abs(mcm-alpha))
        
        alphahat=alphalist[o]
        return(list(alphahat=alphahat,actual=mcm[o]))
        }
error.bars <-function(x, upper, lower, width = 0.02, ...) {
  xlim <- range(x)
  barw <- diff(xlim) * width
  segments(x, upper, x, lower, ...)
  segments(x - barw, upper, x + barw, upper, ...)
  segments(x - barw, lower, x + barw, lower, ...)
  range(upper, lower)
}

