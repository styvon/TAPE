# from SPAtest ScoreTest_wSaddleApprox_Get_X1
getX1 <- function(X1)
{
  q1<-ncol(X1)
  if(q1>=2)
  {
    if(sum(abs(X1[,1]-X1[,2]))==0)
    {
      X1=X1[,-2]
      q1<-q1-1
    }
  }
  qr1<-qr(X1)
  if(qr1$rank < q1){
    
    X1.svd<-svd(X1)
    X1 = X1.svd$u[,1:qr1$rank]
  } 
  
  return(X1)
}

#' Null Linear Mixed Model using AI-REML
#'
#' @param fit0        GLM object. 
#' @param KgenFile    String. PLINK file to construct the GRM
#' @param idstoIncludeFile  String. File including a column of rsids to include for GRM construction
#' @param K           Matrix (sparse). Kinship or genetic relatedness matrix
#' @param KmatFile    Vector of strings. K matrix file location
#' @param tau         Vector. Variance component estimates. Default is (1,0)
#' @param fixtau      Vector. Indicators for fixed value for variance component estimates. Default is (0,0)
#' @param maxiter     Integer. Max number of iterations to fit the LMM model. Default is 500
#' @param tol         Numeric. Tolerance for convergence of variance components. Default is 0.02
#' @param tol_pcg     Numeric. Tolerance for the PCG algorithm. Default is 1e-5
#' @param maxiter_pcg Integer. Max number of iterations for the PCG algorithm. Default is 500
#' @param nrun_trace  Integer. Number of random vectors used for trace estimation. Default is 30
#' @param cutoff_trace Numeric. Threshold for the coefficient of variation for trace estimation. Default is 0.0025
#' @param verbose     Indicator for message output. Default is FALSE
#' @param n_memchunk  Integer. Number of memory chunks for genotype. Default is 1
#' @param loco        Indicator for leave-one-chromosome-out analysis. Default is FALSE
#' @param tol_tau     Numeric. Minimum threshold for variance components. Default is 1e-5
#' @param tol_coef    Numeric. Tolerance for convergence of coefficients. Default is 0.1
#' @return List from getTau() or its variants
#' @export
glmmaiNullModel <- function(fit0, KgenFile="", idstoIncludeFile="",K = NULL, KmatFile="", 
                            tau=c(1, 0), fixtau = c(0,0),
                            maxiter = 500, tol =0.02, tol_pcg = 1e-5, maxiter_pcg = 500,
                            nrun_trace = 30, cutoff_trace = 0.0025, 
                            verbose=FALSE, n_memchunk = 1, loco=FALSE, tol_tau=1e-5, tol_coef=0.1){
  #===0 Get fix effects alpha and working Y ====
  y = fit0$y
  n <- length(y)
  offset <- fit0$offset
  if(is.null(offset)){
    offset = 0
  }
  if(is.null(offset)) offset <- rep(0, n)
  family <- fit0$family
  psi <- summary(fit0)$dispersion
  eta <- fit0$linear.predictors
  mu <- fit0$fitted.values
  mu.eta <- family$mu.eta(eta)
  Y <- eta - offset + (y - mu)/mu.eta
  sqrtW <- mu.eta/sqrt(1/as.vector(weights(fit0))*fit0$family$variance(mu))
  X <- model.matrix(fit0)
  X = SPAtest:::ScoreTest_wSaddleApprox_Get_X1(X)
  
  alpha <- fit0$coef
  
  if( (!is.null(K)) & (KgenFile!="") ){
    #===1. Use sparse kinmat & grm from PLINK geno file [2+1 tau]====
    if(verbose){
      cat("\nVariance components: block diagonal K matrix and dense GRM constructed from Plink file input \n")
    }
    #===1.0 Initialize geno obj in c++ ====
    subSampleInGeno = as.numeric(scan(idstoIncludeFile, character(), quote = "")) # vec of int
    if(verbose){
      print("Reading PLINK file for GRM construction...")
    }
    
    re1 = system.time({n_nomissing=setgeno(KgenFile, subSampleInGeno, memoryChunk=n_memchunk, isDiagofKinSetAsOne=TRUE)})
    
    if(loco){
      # ==== if loco ====
      startid_chr_LOCO = NULL
      endid_chr_LOCO = NULL
      bimData = data.table:::fread(paste0(KgenFile,".bim"),  header=F)
      for(ichr in 1:22){
        
        if(length(which(bimData[,1] == ichr)) > 0){
          startid_chr_LOCO = c(startid_chr_LOCO, min(which(bimData[,1] == ichr))-1)
          endid_chr_LOCO = c(endid_chr_LOCO, max(which(bimData[,1] == ichr))-1)
        }else{ 
          startid_chr_LOCO = c(startid_chr_LOCO, NA)
          endid_chr_LOCO = c(endid_chr_LOCO, NA)
        } # end if ichr has snps
        
      } # end for ichr
      if(verbose){
        cat("startid_chr_LOCO: ", startid_chr_LOCO, "\n")
        cat("endid_chr_LOCO: ", endid_chr_LOCO, "\n")
      }
      if(sum(!is.na(startid_chr_LOCO)) <= 1 | sum(!is.na(startid_chr_LOCO)) <= 1){
        cat("Number of autosomal chromosomes < 2, leave-one-chromosome-out analysis disabled\n")
        loco=FALSE
      } # end loco check
    }else{ 
      # ==== if no loco ====
      startid_chr_LOCO = rep(NA, 22)
      endid_chr_LOCO = rep(NA, 22)
    } # end if loco
    
    #===1.1 Initialize tau ====
    q_k <- 1 # total number of sparse variance components: 1
    # q_k <- length(K) # total number of sparse variance components
    q_f <- 1 # total number of dense variance components
    q <- q_k+q_f # total number of variance components
    if(verbose){
      cat("Total number of variance components:",q+1,"\n")
    }
    
    if(length(tau)!=(q+1)){
      msg <- paste("Input tau length:",length(tau), ", not matching required (",(q+1),"), updated to match length.")
      warning(msg)
      tau <- c(1,rep(0,q))
      fixtau <- rep(0,q+1)
    }
    if(family$family %in% c("poisson", "binomial")) {
      tau[1] <- 1
      fixtau[1] <- 1
    }
    
    idxtau <- which(fixtau == 0)
    if( any(fixtau==0) ){ 
      tau[fixtau == 0] <- psi/(q+1) # psi = var(Y)
    }
    if(family$family == "gaussian" & fixtau[1]==0){
      tau[1] <- psi # psi = var(Y)
    }
    if(verbose){
      cat("Initial tau is: ", paste(tau, collapse =","), "\n")
    }
    
    #===1.2 Set Sigma and Update tau ====
    re <- getTau_plink(n_nomissing= n_nomissing, sparsekin=K, tau=tau, fixtau=fixtau,
                 Y = Y, y=y, X=X,W=sqrtW^2,alpha=alpha,
                 maxiter=maxiter, tol=tol, tol_pcg=tol_pcg,maxiter_pcg=maxiter_pcg,
                 nrun_trace=nrun_trace, cutoff_trace=cutoff_trace, verbose=verbose, 
                 loco=loco, startid_chr_LOCO=startid_chr_LOCO, endid_chr_LOCO=endid_chr_LOCO, tol_tau=tol_tau,
                 tol_coef=tol_coef)
    closegeno()
  }else if( (!is.null(K)) & (KmatFile!="") ){
    #===2. Use sparse kinmat & grm from matrix file input====
    # K: sparse matrix, Kmatfile: file path to dense GRM
    if(verbose){
      cat("\nVariance components:  block diagonal K matrix and dense GRM file\n")
    }
    
    #===2.1 Initialize tau ====
    q_k <- 1 # total number of sparse variance components: 1
    # q_k <- length(K) # total number of sparse variance components
    q_f <- length(KmatFile) # total number of dense variance components
    q <- q_k+q_f # total number of variance components
    if(verbose){
      cat("Total number of variance components:",q+1,"\n")
    }
    
    if(length(tau)!=(q+1)){
      msg <- paste("Input tau length:",length(tau), ", not matching required (",(q+1),"), updated to match length.")
      warning(msg)
      tau <- c(1,rep(0,q))
      fixtau <- rep(0,q+1)
    }
    if(family$family %in% c("poisson", "binomial")) {
      tau[1] <- 1
      fixtau[1] <- 1
    }
    
    idxtau <- which(fixtau == 0)
    if( any(fixtau==0) ){ 
      tau[fixtau == 0] <- var(Y)/(q+1)
    }
    if(family$family == "gaussian" & fixtau[1]==0){
      tau[1] <- psi # psi = var(Y)
    }
    if(verbose){
      cat("Initial tau is: ", paste(tau, collapse =","), "\n")
    }
    
    
    #===2.2 Set Sigma and Update tau ====
    re <- getTau(sparsekin=K,kmatfile=KmatFile,tau=tau, fixtau=fixtau,
                 Y = Y, y=y, X=X,W=sqrtW^2,alpha=alpha,
                 maxiter=maxiter, tol=tol, tol_pcg=tol_pcg,maxiter_pcg=tol_pcg,
                 nrun_trace=nrun_trace, cutoff_trace=cutoff_trace, verbose=verbose)
    
  } else if( (!is.null(K)) & (KmatFile=="") ){
    #===3. If((!is.null(K)) & (KmatFile=="")) [single kinship cov component]====
    if(verbose){
      cat("Variance components:  block diagonal Kinship matrix\n")
    }
    
    n_nomissing = length(Y)
    
    #===3.1 Initialize tau ====
    q_k <- 1 # total number of sparse variance components: 1
    q_f <- 0 # total number of dense variance components
    q <- q_k+q_f # total number of variance components
    if(verbose){
      cat("Total number of variance components:",q+1,"\n")
    }
    
    if(length(tau)!=(q+1)){
      msg <- paste("Input tau length:",length(tau), ", not matching required (",(q+1),"), updated to match length.")
      warning(msg)
      tau <- c(1,rep(0,q))
      fixtau <- rep(0,q+1)
    }
    if(family$family %in% c("poisson", "binomial")) {
      tau[1] <- 1
      fixtau[1] <- 1
    }
    
    idxtau <- which(fixtau == 0)
    if( any(fixtau==0) ){ 
      tau[fixtau == 0] <- var(Y)/(q+1)
    }
    if(family$family == "gaussian" & fixtau[1]==0){
      tau[1] <- psi # psi = var(Y)
    }
    if(verbose){
      cat("Initial tau is: ", paste(tau, collapse =","), "\n")
    }
    
    #===3.2 Set Sigma and Update tau ====
    re <- getTau(n_nomissing=n_nomissing, sparsekin=K,kmatfile=KmatFile,tau=tau, fixtau=fixtau,
                 Y = Y, y=y, X=X,W=sqrtW^2,alpha=alpha,
                 maxiter=maxiter, tol=tol, tol_pcg=tol_pcg,maxiter_pcg=maxiter_pcg,
                 nrun_trace=nrun_trace, cutoff_trace=cutoff_trace, verbose=verbose)
    
  }else if((is.null(K)) & (KgenFile!="")){
    #=== 4. Use only grm from plink file ====
    if(verbose){
      cat("\nVariance components: dense GRM constructed from Plink file input \n")
    }
    #===4.0 Initialize geno obj in c++ ====
    subSampleInGeno = as.numeric(scan(idstoIncludeFile, character(), quote = "")) # vec of int
    if(verbose){
      print("Reading PLINK file for GRM construction...")
    }
    re1 = system.time({n_nomissing=setgeno(KgenFile, subSampleInGeno, memoryChunk=1, isDiagofKinSetAsOne=TRUE)})
    
    #===4.1 Initialize tau ====
    q_k <- 0 # total number of sparse variance components: 0
    q_f <- 1 # total number of dense variance components
    q <- q_k+q_f # total number of variance components
    if(verbose){
      cat("Total number of variance components:",q+1,"\n")
    }
    
    if(length(tau)!=(q+1)){
      msg <- paste("Input tau length:",length(tau), ", not matching required (",(q+1),"), updated to match length.")
      warning(msg)
      tau <- c(1,rep(0,q))
      fixtau <- rep(0,q+1)
    }
    if(family$family %in% c("poisson", "binomial")) {
      tau[1] <- 1
      fixtau[1] <- 1
    }
    
    idxtau <- which(fixtau == 0)
    if( any(fixtau==0) ){ 
      tau[fixtau == 0] <- var(Y)/(q+1)
    }
    if(family$family == "gaussian" & fixtau[1]==0){
      tau[1] <- var(Y)
    }
    if(verbose){
      cat("Initial tau is: ", paste(tau, collapse =","), "\n")
    }
    
    #===4.2 Set Sigma and Update tau ====
    re <- getTau_plink_nok(n_nomissing= n_nomissing, tau=tau, fixtau=fixtau,
                       Y = Y, y=y, X=X,W=sqrtW^2,alpha=alpha,
                       maxiter=maxiter, tol=tol, tol_pcg=tol_pcg,maxiter_pcg=maxiter_pcg,
                       nrun_trace=nrun_trace, cutoff_trace=cutoff_trace, verbose=verbose, tol_tau=tol_tau,
                       tol_coef=tol_coef)
    closegeno()
  }# end if (!is.null(K)) & (KmatFile!="") 
  
  re$W=sqrtW^2
  return(re)
  
  
}