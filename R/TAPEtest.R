#' Generating test results from null model
#'
#' Supports null model objects from TAPE or GMMAT.
#' @param null_object   Null model object
#' @param genfile       Genotype file prefix
#' @param samplefile    String, file with a column of IIDs in genfile
#' @param outfile       Path to output result file
#' @param genfile_format (Optional) String, format of genotype file. Default is 'bfile'
#' @param snpbatch_size (Optional) Number of SNPs per batch of analysis to be written to outfile. Default is 20000.
#' @param thresh_maf    (Optional) Threshold of minimum MAF to be analyzed. Default is 0.
#' @param thresh_mac    (Optional) Threshold of minimum MAC to be analyzed. Default is 0.
#' @param espa_nt       (Optional) Number of points in empirical saddlepoint approximation, used only if null object is from GMMAT. Default is 9999.
#' @param espa_range    (Optional) Range of support in approximation, used only if null object is from GMMAT. Default is c(-20,20). 
#' @param bgi_file      (Optional) Path or Suffix type for bgi file. 1="genfile.bgen.bgi", 2="genfile.bgi". Default is 1.
#' @param snplist_file  (Optional) String, file of a subset of snp ids to be included in test . Default is ''
#' @param verbose       (Optional) Boolean, indicator for output runtime. Default is FALSE.
#' @return Number of tested variants, and test results at outfile 
#' @export
TAPEtest <- function(null_object, genfile, samplefile, outfile, genfile_format="bfile", thresh_maf=0, thresh_mac=0, snpbatch_size=20000, espa_nt=9999, espa_range=c(-20,20), 
                     bgi_file="1",snplist_file="", verbose=FALSE){
  cat("Starting TAPEtest...\n")
  # ==== validate null model ======
  if(!(class(null_object) %in% c("glmmkin", "glmmai","tape_null"))){
    stop("Object class not supported")
  }else if(class(null_object) %in% c("glmmkin", "glmmai")){
    idx0 = qcauchy(1:espa_nt/(espa_nt+1))
    idx1 = idx0*max(espa_range)/max(idx0)
    resid <- null_object$residuals
    
    cumul<-NULL
    for(id in idx1){
      t<-id
      e_resid<-exp(resid*t)
      M0<-mean(e_resid)
      M1<-mean(resid*e_resid)
      M2<-mean(resid^2*e_resid)
      k0<-log(M0)
      k1<-M1/M0
      k2<-(M0*M2-M1^2)/M0^2
      cumul<-rbind(cumul, c(t, k0, k1, k2))
    }
    
    K_org_emp<-approxfun(cumul[,1], cumul[,2], rule=2)
    K_1_emp<-approxfun(cumul[,1], cumul[,3], rule=2)
    K_2_emp<-approxfun(cumul[,1], cumul[,4], rule=2)
    
    null_object$K_org_emp=K_org_emp
    null_object$K_1_emp=K_1_emp
    null_object$K_2_emp=K_2_emp
  }
  
  # ===== result headers ======
  if(genfile_format=="bfile"){
    output_header <- c("CHR", "SNPID", "BP", "ALT", "REF", "PROPMISS", "MAF", "SCORE", "VAR", "PVAL", "PVAL_NORM")
  }else if(genfile_format=="bgen"){
    output_header <- c("CHR","BP","SNPID","ALT","REF","AC","AF", "MAF", "SCORE", "VAR", "PVAL", "PVAL_NORM")
  }
  
  cat("Genotype format:", genfile_format, "\n")
  if(genfile_format=="bfile"){
    # ==== 1. PLINK format =======
    # ==== * validate genotype file =======
    if( any( !file.exists( paste0(genfile,c(".bim",".bed",".fam")) )  )  ){
      msg = paste0("Missing file: ",genfile," (.bim/.bed/.fam)")
      stop(msg)
    }
    
    info_variants <- data.table::fread(paste0(genfile,".bim"))[,c(1,2,4:6)]
    names(info_variants) <- c("CHR","SNPID","BP","ALT","REF")
    n_variants <- dim(info_variants)[1]
    n_variants_tested <- 0
    
    # ==== * split into batches ==== 
    n_batches <- (n_variants-1) %/% snpbatch_size + 1
    
    # ==== * match IIDs==== 
    cat("====== Matching IIDs ======\n")
    idstoIncludeFile = paste0(outfile,".idorderingeno")
    temp_iids_geno <- as.vector(data.table::fread(paste0(genfile,".fam"))$V2)
    temp_iids_cov <- as.vector(null_object$fit0$data$IID) # Goncalo IID coding
    temp_ids <- match(temp_iids_cov,temp_iids_geno)
    contents <- paste0(paste0(temp_ids,collapse = "\n"),"\n")
    cat(contents,file = idstoIncludeFile)
    rm(list=c("temp_iids_geno","temp_iids_cov","temp_ids"))
    
    # ==== * load genfile==== 
    cat("====== Loading genotype ======\n")
    tic <- proc.time()
    subSampleInGeno = as.numeric(scan(idstoIncludeFile, character(), quote = "")) # vec of int
    N = length(null_object$fit0$data$IID)
    cat(N, "samples matched in plink file\n")
    re_geno=setgeno_fast(genfile, subSampleInGeno, n_memchunk=1,n_snp_sample=-1, thresh_mac=-1)
    n_nomissing=re_geno$n_nomissing
    toc <- proc.time()
    time_init <- toc - tic
    if(verbose){
      cat("Proc time init: ",time_init,"\n")
    }
    
    # ==== * get pval in batches ==== 
    cat("====== Testing ======\n")
    for(ibatch in 1:n_batches){ 
      cat("------ Batch", ibatch,"/", n_batches," ------\n")
      output <- c()
      
      if(ibatch == n_batches){
        snps_in_batch <- ((ibatch-1)*snpbatch_size+1):n_variants
      }else{
        snps_in_batch <- ((ibatch-1)*snpbatch_size+1):(ibatch*snpbatch_size)
      }
      
      for(piter in snps_in_batch){
        if(verbose){
          tic <- proc.time()
          geno <- getOneSNPR_nostd(piter-1)
          toc <- proc.time()
          time_load <- toc - tic
          cat("Proc time load: ",time_load,"\n")
        }else{
          geno <- getOneSNPR_nostd(piter-1)
        } # end if verbose
        
        
        idx_nomiss <- which(geno!=-9)
        idx_miss <- which(geno==-9)
        propmiss <- mean(geno==-9)
        
        # render minor allele representation
        maf <- mean(geno[idx_nomiss])/2
        mac <- sum(geno[idx_nomiss])
        
        # check data validity
        if(mac <= thresh_mac | min(maf,1-maf) <=thresh_maf ){
          next
        }
        
        n_variants_tested = n_variants_tested + 1
        info_onevariant <- info_variants[piter,]
        if(maf>0.5){
          maf <- 1-maf
          geno <- abs(2-geno)
          info_onevariant[,c("ALT","REF")] <- info_onevariant[,c("REF","ALT")]
        }
        
        # get TAPE test results
        if(verbose){
          tic <- proc.time()
          set_nonzero = which(geno!=0)
          geno_mat <- as.matrix(geno)
          toc <- proc.time()
          time_pretest <- toc - tic
          cat("Proc time pretest: ",time_pretest,"\n")
          
          tic <- proc.time()
          tape_results <- Get_Pvalues_espa(object = null_object, G= geno_mat, set_nonzero = set_nonzero,cutoff = 2)
          toc <- proc.time()
          time_test <- toc - tic
          cat("Proc time test: ",time_test,"\n")
        }else{
          set_nonzero = which(geno!=0)
          geno_mat <- as.matrix(geno)
          tape_results <- Get_Pvalues_espa(object = null_object, G= geno_mat, set_nonzero = set_nonzero,cutoff = 2)
          
        } # end if(verbose)
        
        # append results: CHR, SNPID, BP, ALT, REF, PROPMISS, MAF, SCORE, VAR, PVAL, PVAL_NORM
        output <- rbind(output, c(info_onevariant, propmiss, maf, tape_results$score, tape_results$var_score, tape_results$pval, tape_results$pval_nospa))
      } # end for piter in snps_in_batch
      
      if(ibatch==1){
        colnames(output) <- output_header
        write.table(output, outfile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", append=FALSE, na="NA")
      }else{
        write.table(output, outfile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE, na="NA")
      }
    } # end for ibatch
    closegeno()
    
  }else if(genfile_format=="bgen"){
    
    if(bgi_file=="1"){
      bgi_suffix = ".bgen.bgi"
    }else if(bgi_file=="2"){
      bgi_suffix = ".bgi"
    }else{
      bgi_suffix = ""
    }
    cat("bgi_suffix=",bgi_suffix,"\n")
    # ==== 2. bgen format =======
    # ==== * validate genotype file =======
    if(bgi_suffix !=""){
      if( any( !file.exists( paste0(genfile,c(".bgen",bgi_suffix)) )  )  ){
        msg = paste0("Missing file: ",genfile," (.bgen/.bgi)")
        stop(msg)
      }
      file_bgi = paste0(genfile, bgi_suffix)
    }else{
      if( (!file.exists( paste0(genfile,".bgen") )) | (!file.exists( bgi_file) ) )  {
        msg = paste0("Missing file: ",genfile," (.bgen/.bgi)")
        stop(msg)
      }
      file_bgi = bgi_file
    }
    
    file_bgen = paste0(genfile, ".bgen")
    
    # # ==== * match IIDs==== 
    cat("====== Matching IIDs ======\n")
    iid_in_geno = data.table:::fread(samplefile, header=F, stringsAsFactors=FALSE)$V1
    iid_in_model = null_object$fit0$data$IID
    mapid_geno_to_model=match(iid_in_geno,iid_in_model)
    N = sum(!is.na(mapid_geno_to_model))
    cat(N, "samples matched in bgen file\n")
    mapid_geno_to_model[is.na(mapid_geno_to_model)] = -10 
    mapid_geno_to_model=mapid_geno_to_model-1
    
    # ==== * load genfile==== 
    names_info_variants <- c("CHR","BP","SNPID","ALT","REF","AC","AF")
    n_variants_tested <- 0
    
    cat("====== Loading genotype ======\n")
    if(snplist_file==""){
      # ---- ** if snplist_file =="" ----
      ids_to_exclude = as.character(vector())
      ids_to_include = as.character(vector())
      ranges_to_exclude = data.frame(chromosome = NULL, start = NULL, end = NULL)
      ranges_to_include = data.frame(chromosome = NULL, start = NULL, end = NULL)
      
      has_variant=TRUE
      if(verbose){
        tic <- proc.time()
        n_variants = setgenoTest_bgenDosage(file_bgen,file_bgi, ranges_to_exclude = ranges_to_exclude, ranges_to_include = ranges_to_include, ids_to_exclude= ids_to_exclude, ids_to_include=ids_to_include)
        if(n_variants == 0){
          has_variant = FALSE
        }
        is_query = getQueryStatus()
        SetSampleIdx(mapid_geno_to_model, N)
        toc <- proc.time()
        time_init <- toc - tic
        cat("Proc time init: ",time_init,"\n")
      }else{
        n_variants = setgenoTest_bgenDosage(file_bgen,file_bgi, ranges_to_exclude = ranges_to_exclude, ranges_to_include = ranges_to_include, ids_to_exclude= ids_to_exclude, ids_to_include=ids_to_include)
        if(n_variants == 0){
          has_variant = FALSE
        }
        is_query = getQueryStatus()
        SetSampleIdx(mapid_geno_to_model, N)
      } # end if verbose
      
      cat("n_variants","variants read in\n")
      
      
      # ==== *** split into batches ==== 
      n_batches <- (n_variants-1) %/% snpbatch_size + 1
      
      # ==== *** get pval in batches ==== 
      cat("====== Testing ======\n")
      for(ibatch in 1:n_batches){ 
        cat("------ Batch", ibatch,"/", n_batches," ------\n")
        output <- c()
        
        if(ibatch == n_batches){
          snps_in_batch <- ((ibatch-1)*snpbatch_size+1):n_variants
        }else{
          snps_in_batch <- ((ibatch-1)*snpbatch_size+1):(ibatch*snpbatch_size)
        }
        
        for(piter in snps_in_batch){
          if(verbose){
            tic <- proc.time()
            Gx = getDosage_bgen_noquery()
            toc <- proc.time()
            time_load <- toc - tic
            cat("Proc time load: ",time_load,"\n")
          }else{
            Gx = getDosage_bgen_noquery()
          } # end if verbose
          
          info_onevariant <- as.vector(unlist(Gx$variants))[-4]
          AC = Gx$variants$AC
          AF = Gx$variants$AF
          maf <- min(AF, 1-AF)
          mac <- round(maf * N)
          
          # check data validity
          if(mac <= thresh_mac | maf<=thresh_maf ){
            next
          }
          n_variants_tested = n_variants_tested + 1
          geno = Gx$dosages
          
          # get TAPE test results
          if(verbose){
            tic <- proc.time()
            set_nonzero = which(geno!=0)
            geno_mat <- as.matrix(geno)
            toc <- proc.time()
            time_pretest <- toc - tic
            cat("Proc time pretest: ",time_pretest,"\n")
            
            tic <- proc.time()
            tape_results <- Get_Pvalues_espa(object = null_object, G= geno_mat, set_nonzero = set_nonzero,cutoff = 2)
            toc <- proc.time()
            time_test <- toc - tic
            cat("Proc time test: ",time_test,"\n")
          }else{
            set_nonzero = which(geno!=0)
            geno_mat <- as.matrix(geno)
            tape_results <- Get_Pvalues_espa(object = null_object, G= geno_mat, set_nonzero = set_nonzero,cutoff = 2)
          } # end if verbose
          
          
          # append results: "CHR","BP","SNPID","ALT","REF","AC","AF", MAF, SCORE, VAR, PVAL, PVAL_NORM
          output <- rbind(output, c(info_onevariant, maf, tape_results$score, tape_results$var_score, tape_results$pval, tape_results$pval_nospa))
        } # end for piter in snps_in_batch
        
        if(ibatch==1){
          colnames(output) <- output_header
          write.table(output, outfile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", append=FALSE, na="NA")
        }else{
          write.table(output, outfile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE, na="NA")
        }
        
      } # end for ibatch
      
      closetestGenoFile_bgenDosage()
      
    }else{
      # ---- ** if snplist_file specified (use query) ----
      ids_to_exclude = as.character(vector())
      temp <- as.data.frame(data.table:::fread(snplist_file, header=F, sep=" ", stringsAsFactors=FALSE, colClasses=c("character")))
      ids_to_include = as.character(
        temp$V1
      )
      rm(temp)
      ranges_to_exclude = data.frame(chromosome = NULL, start = NULL, end = NULL)
      ranges_to_include = data.frame(chromosome = NULL, start = NULL, end = NULL)
      
      has_variant=TRUE
      if(verbose){
        tic <- proc.time()
        
        n_variants = setgenoTest_bgenDosage(file_bgen,file_bgi, ranges_to_exclude = ranges_to_exclude, ranges_to_include = ranges_to_include, ids_to_exclude= ids_to_exclude, ids_to_include=ids_to_include)
        if(n_variants == 0){
          has_variant = FALSE
        }
        is_query = getQueryStatus()
        SetSampleIdx(mapid_geno_to_model, N)
        
        toc <- proc.time()
        time_init <- toc - tic
        cat("Proc time init: ",time_init,"\n")
      }else{
        n_variants = setgenoTest_bgenDosage(file_bgen,file_bgi, ranges_to_exclude = ranges_to_exclude, ranges_to_include = ranges_to_include, ids_to_exclude= ids_to_exclude, ids_to_include=ids_to_include)
        if(n_variants == 0){
          has_variant = FALSE
        }
        is_query = getQueryStatus()
        SetSampleIdx(mapid_geno_to_model, N)
      } # end if verbose
      
      
      # return if no matching markers
      if(!is_query){
        cat("0 markers matched with snplist_file\n")
        return(n_variants_tested)
      }
      
      output <- c()
      for(piter in 1:n_variants){
        if(piter%%5000 ==1){
          cat("Variant",piter,"\n")
        }
        
        if(verbose){
          tic <- proc.time()
          Gx = getDosage_bgen_withquery()
          toc <- proc.time()
          time_load <- toc - tic
          cat("Proc time load: ",time_load,"\n")
        }else{
          Gx = getDosage_bgen_withquery()
        } # end if verbose
        
        info_onevariant <- as.vector(unlist(Gx$variants))[-4]
        geno = Gx$dosages
        AC = Gx$variants$AC
        AF = Gx$variants$AF
        maf <- min(AF, 1-AF)
        mac <- round(maf * N)
        
        # check data validity
        if(mac <= thresh_mac | maf<=thresh_maf ){
          next
        }
        n_variants_tested = n_variants_tested + 1
        
        # get TAPE test results
        if(verbose){
          tic <- proc.time()
          set_nonzero = which(geno!=0)
          geno_mat <- as.matrix(geno)
          toc <- proc.time()
          time_pretest <- toc - tic
          cat("Proc time pretest: ",time_pretest,"\n")
          
          tic <- proc.time()
          tape_results <- Get_Pvalues_espa(object = null_object, G= geno_mat, set_nonzero = set_nonzero,cutoff = 2)
          toc <- proc.time()
          time_test <- toc - tic
          cat("Proc time test: ",time_test,"\n")
        }else{
          set_nonzero = which(geno!=0)
          geno_mat <- as.matrix(geno)
          tape_results <- Get_Pvalues_espa(object = null_object, G= geno_mat, set_nonzero = set_nonzero,cutoff = 2)
        } # end if verbose
        
        
        # append results: "CHR","BP","SNPID","ALT","REF","AC","AF", MAF, SCORE, VAR, PVAL, PVAL_NORM
        output <- rbind(output, c(info_onevariant, maf, tape_results$score, tape_results$var_score, tape_results$pval, tape_results$pval_nospa))
      } # end for piter
      
      colnames(output) <- output_header
      write.table(output, outfile, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", append=FALSE, na="NA")
      
      closetestGenoFile_bgenDosage()
    } # end if snplist_file
      
    
  } # end if genfile_format
  
  cat("Test results saved to:",outfile,"\n")
  return(n_variants_tested)
}