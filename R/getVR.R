# suppressPackageStartupMessages( require(seqminer) ) # for read in PLINK file

#' Calculating variance ratio using genotype subsamples
#'
#' @param null_object   Null GLMM model object
#' @param fit0          Null GLM model object
#' @param genfile       String. PLINK file to genotype samples
#' @param K             Matrix (sparse). Kinship or genetic relatedness matrix
#' @param outfile       String. Output file prefix
#' @param n_snp_sample  Integer. Number of variants used for variance ratio estimation. Default is 30.
#' @param thresh_maf    Numeric. MAF threshold. Default is 0.01
#' @param verbose       Bool. Default is TRUE
#' @param n_memchunk    Integer. Number of memory chunks. Default is 1.
#' @param tol_vr      Numeric. Tolerance for variance ratio calculation. Default is 0.001
#' @param tol_pcg     Numeric. Tolerance for the PCG algorithm. Default is 1e-5
#' @param maxiter_pcg Integer. Max number of iterations for the PCG algorithm. Default is 500
#' @return Value of variance ratio (with details written to paste0(outfile,"_variance_ratio.txt"))
#' @export
getVR <- function(null_object, fit0, genfile, K, outfile, n_snp_sample=30, thresh_maf = 0.01, verbose=TRUE,
                  n_memchunk=1, tol_vr=0.001, tol_pcg = 1e-5, maxiter_pcg = 500){
    
  #==== prep output file ====
  file_vr = paste0(outfile,"_variance_ratio.txt")
  if(file.exists(file_vr)){
    file.remove(file_vr)
  } # end if file.exists(file_vr)
  output_header <- c("CHR", "SNPID", "BP", "ALT", "REF", "MAF", "VAR_EXACT", "VAR_APPROX")
  
  # ==== validate genotype file =======
  if( any( !file.exists( paste0(genfile,c(".bim",".bed",".fam")) )  )  ){
    msg = paste0("Missing file: ",genfile," (.bim/.bed/.fam)")
    stop(msg)
  }
  
  # ==== match IIDs==== 
  file_idx_indv = paste0(outfile,".idorderinvrgeno")
  temp_iids_geno <- as.vector(data.table::fread(paste0(genfile,".fam"))$V2)
  temp_iids_cov <- as.vector(fit0$data$IID) # Goncalo IID coding
  temp_ids <- match(temp_iids_cov,temp_iids_geno)
  contents <- paste0(paste0(temp_ids,collapse = "\n"),"\n")
  cat(contents,file = file_idx_indv)
  rm(list=c("temp_iids_geno","temp_iids_cov","temp_ids"))
  
  #==== load genfile ======
  idx_indv = as.numeric(scan(file_idx_indv, character(), quote = "")) # vec of int
  if(verbose){
    print("Reading PLINK file for variance ratio...")
  }
  system.time({re_geno=setgeno(genfile, idx_indv, memoryChunk=n_memchunk,isDiagofKinSetAsOne=T)}) 
  N = re_geno
  
  maf_vec <- getMAFR() 
  maf_vec <- 0.5 - abs(maf_vec-0.5)
  idx_variants_sample = sort(sample(which(maf_vec>=thresh_maf),n_snp_sample))
  
  info_variants <- data.table::fread(paste0(genfile,".bim"))[idx_variants_sample,c(1,2,4:6)]
  names(info_variants) <- c("CHR","SNPID","BP","ALT","REF")
  
  #==== fetch variables =====
  X = null_object$X
  Y = null_object$Y
  family = fit0$family
  eta = null_object$linear_predictors
  mu = null_object$fitted_values
  W = null_object$W
  tau_vec = null_object$theta
  Sigma_iX = null_object$Sigma_iX
  
  #==== iter through sample snps =====
  if(verbose){
    print("Iterating through sample snps...")
  }
  cv_vr <- tol_vr + 0.1
  is_first_try <- TRUE
  output <- c()
  while(cv_vr > tol_vr){
    
    if(!is_first_try){
      n_snp_sample = n_snp_sample + 10
      cat("CV larger than threshold, increase number of markers to", n_snp_sample,"\n")
      
      idx_variants_sample = sort(sample(which(maf_vec>=thresh_maf),n_snp_sample))
      
      info_variants <- data.table::fread(paste0(genfile,".bim"))[idx_variants_sample,c(1,2,4:6)]
      names(info_variants) <- c("CHR","SNPID","BP","ALT","REF")
    }
    
    for(piter in 1:length(idx_variants_sample)){
      id_variant <- idx_variants_sample[piter]
      geno <- getOneSNPR_nostd(id_variant-1)
      
      idx_nomiss <- which(geno!=-9)
      idx_miss <- which(geno==-9)
      propmiss <- mean(geno==-9)
      
      info_onevariant <- info_variants[piter,]
      maf <- mean(geno[idx_nomiss])/2
      geno[idx_miss] <- maf*2
      if(maf>0.5){
        maf <- 1-maf
        geno <- 2-geno
        info_onevariant[,c("ALT","REF")] <- info_onevariant[,c("REF","ALT")]
      }
      mac <- sum(geno)
      
      Gtilde = geno - X %*% (null_object$XSiX_inv_SiXt %*% geno) # covariate-adjusted geno
      gtilde = Gtilde/sqrt(mac)
      q = innerProduct(gtilde, Y)
      mu_eta = family$mu.eta(eta)
      if(!is.null(K)){
        Sigma_iG = pcgsolve_plink(K, W, tau_vec, Gtilde, "Jacobi", tol_pcg, maxiter_pcg)
        Sigma_ig = pcgsolve_plink(K, W, tau_vec, gtilde, "Jacobi", tol_pcg, maxiter_pcg)
      }else{
        Sigma_iG = pcgsolve_plink_nok(W, tau_vec, Gtilde, "Jacobi", tol_pcg, maxiter_pcg)
        Sigma_ig = pcgsolve_plink_nok(W, tau_vec, gtilde, "Jacobi", tol_pcg, maxiter_pcg)
      }
      
      var_exact_pre = sum(Gtilde*Sigma_iG) - sum(Gtilde * (Sigma_iX %*% (null_object$XSiX_inv_SiXt %*% Gtilde)) )
      var_exact = var_exact_pre/mac
      
      var_approx = sum(gtilde*Sigma_ig)
      
      output <- rbind(output, c(info_onevariant$CHR, info_onevariant$SNPID, info_onevariant$BP,
                                info_onevariant$ALT, info_onevariant$REF,
                                maf, var_exact, var_approx))
      
    } # end for piter
    
    output <- as.data.frame(output)
    names(output) <- output_header
    
    vrs <- as.numeric(output$VAR_EXACT)/as.numeric(output$VAR_APPROX)
    cv_vr <- sd(vrs)/mean(vrs)
    
    if(is_first_try){
      is_first_try = FALSE
    }
    
  } # end while cv_vr
  closegeno()
  gc()
  
  out_vr <- mean(vrs)
  write.table(output, file=file_vr, sep="\t", quote=F, col.names=T, row.names=F)
  
  if(verbose){
    cat("Variance ratio =",out_vr,"\n")
    cat("Detailed results written to:", file_vr,"\n")
  }
  
  
  return(out_vr)
}