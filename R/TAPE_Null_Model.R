#' Null Model for Test with Adjusted Phenotype using Empirical Saddlepoint Approximation
#'
#' @param formula     Formula for the regression model 
#' @param data        Data frame to be used by the formula. Default is NULL.
#' @param KgenFile    String. PLINK file to construct the GRM
#' @param idstoIncludeFile  String. File including individuals to include for GRM construction. A column of integer individual index in genotype file in covariate file order.  
#' @param K           Matrix. Sparse matrix for close relatedness
#' @param KmatFile    String. K matrix file location
#' @param length_out  Number of samples in resampling. Default is 9999
#' @param range       Numeric vector with 2 elements. Range for the resampling. Default is c(-100,100)
#' @param VRgenFile   String. PLINK file of genotypes for estimation of variance ratio. Default is NULL (variance ratio 1)
#' @param idx_ratio   Integer. Number of variants used for variance ratio estimation. Default is 30.
#' @param tau         Vector. Variance component estimates. Default is (1,0)
#' @param fixtau      Vector. Indicators for fixed value for variance component estimates. Default is (0,0)
#' @param maxiter     Integer. Max number of iterations to fit the LMM model. Default is 500
#' @param tol         Numeric. Tolerance for convergence of variance components estimation. Default is 0.02
#' @param tol_pcg     Numeric. Tolerance for the PCG algorithm. Default is 1e-5
#' @param maxiter_pcg Integer. Max number of iterations for the PCG algorithm. Default is 500
#' @param nrun_trace  Integer. Number of random vectors used for trace estimation. Default is 30. If equals NULL, Hutchinson’s randomized trace estimator will not be used.
#' @param cutoff_trace Numeric. Threshold for the coefficient of variation for trace estimation. Default is 0.0025
#' @param outFile     Output file prefix
#' @param inv_normalize Indicator to perform the inverse normalization for the phentoype. Default is FALSE
#' @param verbose     Indicator for message output. Default is FALSE
#' @param n_memchunk  Integer. Number of memory chunks for genotype. Default is 1
#' @param loco        Indicator for leave-one-chromosome-out analysis. Default is FALSE
#' @param tol_tau     Numeric. Minimum threshold for variance components. Default is 1e-5
#' @param tol_coef    Numeric. Tolerance for convergence of coefficients. Default is 0.1
#' @param tol_vr      Numeric. Tolerance for variance ratio calculation. Default is 0.001
#' @param thresh_maf_vr Numeric. Minimum MAF for markers in variance ratio calculation. Default is 0.01.
#' @return a null model object of class "tape_null", containing the following elements
#' * resid
#' * var_resid
#' * var_ratio
#' * K_org_emp
#' * K_1_emp
#' * K_2_emp
#' * method
#' * theta: variance component estimates
#' * scaled_resid
#' * X: covariate matrix
#' * cov: XSigma_inv_X_inv
#' * W
#' * XSiX_inv_SiXt
#' * score0_nok
#' * mu
#' * XVX
#' * converged
#' * LOCO
#' * info_LOCO
#' @export
TAPE_Null_Model <- function(formula, data=NULL,KgenFile="", idstoIncludeFile="",K=NULL, KmatFile="", length_out = 9999,range=c(-100,100),
                            VRgenFile=NULL,idx_ratio =30, 
                            tau=c(1,0), fixtau = c(0,0),
                            maxiter = 500, tol =0.02, tol_pcg = 1e-5, maxiter_pcg = 500,
                            nrun_trace = 30, cutoff_trace = 0.0025, inv_normalize=FALSE,
                            verbose=FALSE, outFile="",n_memchunk=1, loco=FALSE, tol_tau=1e-5,
                            tol_coef=0.1, tol_vr=0.001,thresh_maf_vr=0.01){
  pt_start <- proc.time()
  if(inv_normalize){
    phenoCol=as.character(as.formula(formula)[2])
    cat("Perform the inverse nomalization for ",phenoCol,"\n")
    invPheno = qnorm((rank(data[,which(colnames(data) == phenoCol)], na.last="keep")-0.5)/sum(!is.na(data[,which(colnames(data) == phenoCol)])))
    data[,which(colnames(data) == phenoCol)] = invPheno
  }
  
  # ==== 1. Fixed effect estimation ====
  fit0 <- glm(formula, data=data,family=gaussian)
  pt_glm <- proc.time()
  if(verbose){
    cat("Fixed effect estimation completed in", pt_glm - pt_start, "secs\n")
    print(fit0)
    cat(paste0(rep("-",16),collapse=""),"\n\n")
  }
  
  # ==== 2. Variance component estimation ====
  nullmodel_lmm<-glmmaiNullModel(fit0,KgenFile=KgenFile, idstoIncludeFile=idstoIncludeFile,
                                 K=K,KmatFile=KmatFile, tau=tau, fixtau=fixtau,
                                 maxiter=maxiter, tol=tol, tol_pcg=tol_pcg,maxiter_pcg=maxiter_pcg,
                                 cutoff_trace = cutoff_trace, 
                                 verbose=verbose, nrun_trace=nrun_trace,n_memchunk=n_memchunk,
                                 loco=loco,tol_tau=tol_tau,tol_coef=tol_coef)
  pt_lmm <- proc.time()
  if(verbose){
    cat("Variance component estimation completed in", pt_lmm - pt_glm, "secs\n")
    cat(paste0(rep("-",16),collapse=""),"\n\n")
  }
  
  # ==== 3. CGF approximation ====
  idx0 = qcauchy(1:length_out/(length_out+1))
  idx1 = idx0*max(range)/max(idx0)
  resid = as.vector(nullmodel_lmm$residuals)
  var_resid = var(resid)
  scaled_resid=resid/sqrt(var_resid)
  tau = nullmodel_lmm$theta

  cumul<-NULL
  for(id in idx1){
    t<-id
    # e_resid<-exp(resid*t)
    e_resid<-exp(scaled_resid*t)
    M0<-mean(e_resid)
    # M1<-mean(resid*e_resid)
    # M2<-mean(resid^2*e_resid)
    M1<-mean(scaled_resid*e_resid)
    M2<-mean(scaled_resid^2*e_resid)
    k0<-log(M0)
    k1<-M1/M0
    k2<-(M0*M2-M1^2)/M0^2
    cumul<-rbind(cumul, c(t, k0, k1, k2))
  }

  K_org_emp<-approxfun(cumul[,1], cumul[,2], rule=2)
  K_1_emp<-approxfun(cumul[,1], cumul[,3], rule=2)
  K_2_emp<-approxfun(cumul[,1], cumul[,4], rule=2)
  
  pt_cgf <- proc.time()
  if(verbose){
    cat("CGF approximation completed in", pt_cgf - pt_lmm, "secs\n")
    cat(paste0(rep("-",16),collapse=""),"\n\n")
  }
  
  # ==== 4. Collect model object components ====
  mu = as.vector(nullmodel_lmm$linear_predictors)
  X = nullmodel_lmm$X
  XVX <- t(X) %*% (X)
  XVX_inv_VXt=solve(XVX) %*% t(X)
  score0_nok=colSums(X*fit0$residuals)
  
  pt_collect <- proc.time()
  if(verbose){
    cat("Model object components collected in", pt_collect - pt_cgf, "secs\n")
    cat(paste0(rep("-",16),collapse=""),"\n\n")
  }
  
  # ==== 5. Variance ratio estimation ====
  if(is.null(VRgenFile)){
    cat("VRgenFile is NULL, variance ratio set to 1\n")
    var_ratio <- 1
  }else{
    var_ratio <- getVR(null_object=nullmodel_lmm, fit0=fit0, genfile=VRgenFile, K=K, 
                       outfile=outFile, n_snp_sample=idx_ratio, thresh_maf=thresh_maf_vr,
                       n_memchunk=n_memchunk, tol_vr=tol_vr, tol_pcg = tol_pcg, maxiter_pcg = maxiter_pcg)
  }
  pt_vr <- proc.time()
  if(verbose){
    cat("Variance ratio estimation completed in", pt_vr - pt_collect, "secs\n")
    cat(paste0(rep("-",16),collapse=""),"\n\n")
  }
  
  # loco
  if(!is.null(nullmodel_lmm$result_LOCO)){
    resid_LOCO <- c()
    scaled_resid_LOCO <- c()
    var_resid_LOCO <- c()
    K_org_emp_LOCO <- c()
    K_1_emp_LOCO <- c()
    K_2_emp_LOCO <- c()
    mu_LOCO <- c()
    
    for(chr in 1:22){
      result_LOCO = nullmodel_lmm$result_LOCO[[chr]]
      resid_LOCO <- cbind(resid_LOCO, result_LOCO$residuals) # N by 22
      
      cumul<-NULL
      for(id in idx1){
        t<-id
        e_resid<-exp(result_LOCO$residuals*t)
        M0<-mean(e_resid)
        M1<-mean(result_LOCO$residuals*e_resid)
        M2<-mean(result_LOCO$residuals^2*e_resid)
        k0<-log(M0)
        k1<-M1/M0
        k2<-(M0*M2-M1^2)/M0^2
        cumul<-rbind(cumul, c(t, k0, k1, k2))
      }
      
      K_org_emp_LOCO <- c(K_org_emp_LOCO,approxfun(cumul[,1], cumul[,2], rule=2))
      K_1_emp_LOCO <- c(K_1_emp_LOCO, approxfun(cumul[,1], cumul[,3], rule=2))
      K_2_emp_LOCO <- c(K_2_emp_LOCO, approxfun(cumul[,1], cumul[,4], rule=2))
      
      var_resid_LOCO = c(var_resid_LOCO, var(result_LOCO$residuals))
      scaled_resid_LOCO = cbind(scaled_resid_LOCO, result_LOCO$residuals/sd(result_LOCO$residuals)) #N by 22
      mu_LOCO = cbind(mu_LOCO, as.vector(result_LOCO$linear_predictors)) #N by 22
    } # end for chr
    
    info_LOCO <- list(resid_LOCO = resid_LOCO,
                      scaled_resid_LOCO = scaled_resid_LOCO,
                      var_resid_LOCO = var_resid_LOCO,
                      K_org_emp_LOCO = K_org_emp_LOCO,
                      K_1_emp_LOCO = K_1_emp_LOCO,
                      K_2_emp_LOCO = K_2_emp_LOCO,
                      mu_LOCO = mu_LOCO)
    
    pt_loco <- proc.time()
    if(verbose){
      cat("LOCO-related calculation completed in", pt_loco - pt_vr, "secs\n")
      cat(paste0(rep("-",16),collapse=""),"\n\n")
    }
    
  }else{
    # if not loco
    info_LOCO <- NULL
  } # end if loco
  
  obj_tape=list(resid=resid, var_resid=var_resid,var_ratio=var_ratio,
                K_org_emp=K_org_emp, K_1_emp=K_1_emp, K_2_emp=K_2_emp,
                method = nullmodel_lmm$method, theta=nullmodel_lmm$theta,
                scaled_resid=scaled_resid,
                X=X, cov=nullmodel_lmm$cov,W=nullmodel_lmm$W,
                XSiX_inv_SiXt=nullmodel_lmm$XSiX_inv_SiXt, score0_nok =score0_nok, mu=mu, XVX=XVX,
                XVX_inv_VXt=XVX_inv_VXt, fit0=fit0,Sigma_iX=nullmodel_lmm$Sigma_iX,
                converged=nullmodel_lmm$converged, 
                LOCO=loco, info_LOCO = info_LOCO) 
  
  class(obj_tape)<-"tape_null"
  return(obj_tape)
}