# YZ modified
#' P-value calculation using empirical saddlepoint approximation
#'
#' Test for association between a SNP and a survival object
#' @param object        Object returned from Get_Null_Model()
#' @param G             Genotype matrix (0,1,2)
#' @param N             Sample size in G
#' @param AF            Allele frequency of G
#' @param cutoff        Threshold of sd for utilizing empirical saddlepoint approximation. Default is 2
#' @param set_nonzero   Index of individuals with non-zero genotypes
#' @param r             Calibration factor, ratio of GPG and GWG 
#' @param loco_chr      Indicator for leave-one-chromosome-out analysis. Default is NULL
#' @return a vector with three elements: 
#'  * pval
#'  * pval_nospa 
#'  * qstat: variance-adjusted test statistics
#'  * score: score statistics 
#'  * var: variance of score
#'  * converged: indicator for convergence
#' @export
Get_Pvalues_espa<-function(object, G, N, AF, cutoff=2, set_nonzero, r=1, loco_chr = NULL){

  N_zero = N-length(set_nonzero)
  G_center = as.vector(G-2*AF) # centered geno
  G_nonzero = G[set_nonzero]
  
  if(class(object) %in% c("glmmkin", "glmmai")){ # if obj is GMMAT null model 
    score = sum(G_center * object$residuals)
    var_score = var(object$residuals)*2*AF*(1-AF)*N 
    q1 = score/sqrt(var_score)
  }else if(class(object)=="tape_null"){ # if obj is TAPE null model
    if( object$LOCO ){
      # if loco
      resid = object$info_LOCO$resid_LOCO[,loco_chr]
      var_resid = object$info_LOCO$var_resid_LOCO[loco_chr]
    }else{
      # if not loco
      resid = object$resid
      var_resid = object$var_resid
    } # end if loco
    
    score = sum(G_nonzero * resid[set_nonzero]) #original
    var_score = as.vector(var_resid)*2*AF*(1-AF)*N * object$var_ratio
    q1 = score/sqrt(var_score)
 
  } # end if class(object)
  
  pval_norm = pnorm(abs(q1), lower.tail =FALSE )*2
  
  if(abs(q1) < cutoff){
    re <- list(pval = pval_norm, 
             pval_nospa = pval_norm, 
             qstat = q1,
             score = score, 
             var_score = var_score,
             converged = NA)
    return(re)
  }
  
  var_score = as.vector(var_resid) * sum(G_center^2) * object$var_ratio
  # G_norm = G_center/sqrt(var_score) # centered and scaled geno
  G_norm = G_center/sqrt(sum(G_center^2)) 
  q2 = score/sqrt(var_score)
  
  Gnorm_nonzero = G_norm[set_nonzero]
  # Gnorm_zero = -AF*2/sqrt(var_score)
  Gnorm_zero = -AF*2/sqrt(sum(G_center^2))
  pval1 = GetProb_SPA_sparse(object, Gnorm_nonzero, Gnorm_zero, set_nonzero, N_zero, abs(q2), lower.tail=FALSE)
  pval2 = GetProb_SPA_sparse(object, Gnorm_nonzero, Gnorm_zero, set_nonzero, N_zero, -abs(q2), lower.tail=TRUE)
  # pval1<-GetProb_SPA(object, G_norm, abs(q1), lower.tail =FALSE)
  # pval2<-GetProb_SPA(object, G_norm, -abs(q1), lower.tail =TRUE)
  pval = pval1$pval+pval2$pval
  
  re<-list(pval = pval, 
           pval_nospa = pval_norm, 
           qstat = q1,
           score = score, 
           var_score = var_score,
           converged = NA)
  return(re)
}


# Shawn 12/15/2019
Get_Pvalue<-function(G, obj){
  
  #obj<-out
  if(class(obj)!= "glmmkin_new"){
    stop("obj should be the returned object from glmmkin.ai_new")
  }
  
  score<-crossprod(obj$fittedobj$PY, G)[1,1]
  
  Sigma_iG <- crossprod(obj$fittedobj$Sigma_i, G)
  XSiG = crossprod(obj$fittedobj$Sigma_iX, G)
  PG = Sigma_iG - obj$fittedobj$XSi_XSiX_inv %*% XSiG
  GPG = crossprod(G, PG)
  
  T = score^2/GPG[1,1]
  pval<-pchisq(T, df=1, lower.tail = FALSE)
  
  re<-list(pval = pval, T=T, score=score, var_score=GPG)
  return(re)
  
}


K_org<-function(t, g, object){
  
  n_t = length(t)
  out = rep(0,n_t)
  for(i in 1:n_t){
    t1 = t[i]
    t2 = t1*g
    out[i] = sum(object$K_org_emp(t2))
  }
  return(out)
}


K1_adj<-function(t, g, q, object)
{
  n_t = length(t)
  out = rep(0,n_t)
  
  for(i in 1:n_t){
    t1 = t[i]
    t2 = t1*g
    out[i] = sum(g*object$K_1_emp(t2)) -q
  }
  return(out)
}

K2<-function(t, g, object)
{
  n_t<-length(t)
  out<-rep(0,n_t)
  
  for(i in 1:n_t){
    t1<-t[i]
    t2<-t1*g
    out[i]<-sum(g^2*object$K_2_emp(t2))
  }
  return(out)
}



GetProb_SPA <- function(object, g1, q, lower.tail =TRUE){
  
  out = uniroot(K1_adj, c(-1000,1000), g=g1, q=q, object=object,extendInt = "upX")
  zeta = out$root
  
  k1 = K_org(zeta,  g1, object=object)
  k2 = K2(zeta,  g1, object=object)
  
  temp1 = zeta * q - k1
  
  
  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}
  
  pval = pnorm(w + 1/w * log(v/w), lower.tail =lower.tail )
  pval_norm = pnorm(q, lower.tail =lower.tail )
  
  re <- c(pval,pval_norm)
  return(re)
  
}


GetProb_SPA_sparse <- function(object, Gnorm_nonzero, Gnorm_zero, set_nonzero, N_zero, q2, lower.tail){
  
  out = uniroot_K1_sparse(Gnorm_nonzero=Gnorm_nonzero, Gnorm_zero=Gnorm_zero, set_nonzero=set_nonzero, N_zero=N_zero, q2=q2, object=object)
  zeta = out$root
  
  k0 = K_org_sparse(zeta,  Gnorm_nonzero=Gnorm_nonzero, Gnorm_zero=Gnorm_zero, set_nonzero=set_nonzero, N_zero=N_zero, object=object)
  k2 = K2_sparse(zeta,  Gnorm_nonzero=Gnorm_nonzero, Gnorm_zero=Gnorm_zero, set_nonzero=set_nonzero, N_zero=N_zero, object=object)
  
  w = sign(zeta) * sqrt(2 * (zeta * q2 - k0))
  v = zeta * sqrt(k2)
  
  pval = pnorm(w + 1/w * log(v/w), lower.tail =lower.tail )
  pval_norm = pnorm(q2, lower.tail =lower.tail )
  
  re = list(pval=pval,pval_norm=pval_norm)
  return(re)
  
}


uniroot_K1_sparse = function(object, Gnorm_nonzero, Gnorm_zero, set_nonzero, N_zero, q2,
                           tol = .Machine$double.eps, maxiter = 1000, init=0){
  t = init
  K1_eval <- K1_adj_sparse(t, Gnorm_nonzero, Gnorm_zero, set_nonzero, N_zero, q2, object)
  
  prevJump <- Inf
  rep <- 1
  repeat {
    K2_eval <- K2_sparse(t, Gnorm_nonzero, Gnorm_zero, set_nonzero, N_zero, object)
    tnew <- t - K1_eval/K2_eval
    if (abs(tnew - t) < tol) {
      conv <- TRUE
      break
    }
    if (is.na(tnew) | rep == maxiter) {
      conv = FALSE
      break
    }
    newK1 <- K1_adj_sparse(tnew, Gnorm_nonzero, Gnorm_zero, set_nonzero, N_zero, q2, object)
    if (sign(K1_eval) != sign(newK1)) {
      if (abs(tnew - t) > prevJump - tol) {
        tnew <- t + sign(newK1 - K1_eval) * prevJump/2
        newK1 <- K1_adj_sparse(tnew, Gnorm_nonzero, Gnorm_zero, set_nonzero, N_zero, q2, object)
        prevJump <- prevJump/2
      }
      else {
        prevJump <- abs(tnew - t)
      }
    }
    rep <- rep + 1
    t <- tnew
    K1_eval <- newK1
  }
  return(list(root = t, iter = rep, converged = conv))
}

K_org_sparse<-function(t, Gnorm_nonzero, Gnorm_zero, set_nonzero, N_zero, object){
  
  n_t = length(t)
  out = rep(0,n_t)
  for(i in 1:n_t){
    t1<-t[i]
    t2_zero = t1*Gnorm_zero
    t2_nonzero = t1*Gnorm_nonzero
    out[i] = N_zero*object$K_org_emp(t2_zero) + sum(object$K_org_emp(t2_nonzero))
  }
  return(out)
}

K1_adj_sparse<-function(t, Gnorm_nonzero, Gnorm_zero, set_nonzero, N_zero, q2, object){
  
  n_t<-length(t)
  out<-rep(0,n_t)
  
  for(i in 1:n_t){
    t1 = t[i]
    t2_zero = t1*Gnorm_zero
    t2_nonzero = t1*Gnorm_nonzero
    out[i] = N_zero*Gnorm_zero*object$K_1_emp(t2_zero) + sum(Gnorm_nonzero*object$K_1_emp(t2_nonzero)) - q2
  }
  return(out)
}

K2_sparse<-function(t, Gnorm_nonzero, Gnorm_zero, set_nonzero, N_zero, object){
  
  n_t = length(t)
  out = rep(0,n_t)
  
  for(i in 1:n_t){
    t1 = t[i]
    t2_zero = t1*Gnorm_zero
    t2_nonzero = t1*Gnorm_nonzero
    out[i] = N_zero*Gnorm_zero^2*object$K_2_emp(t2_zero) + sum(Gnorm_nonzero^2*object$K_2_emp(t2_nonzero))
  }
  return(out)
}
