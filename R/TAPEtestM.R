#' Generating test results from null model
#'
#' Supports null model objects from TAPE or GMMAT, less R-cpp overhead
#' @param null_object   Null model object
#' @param genfile       Genotype file prefix
#' @param samplefile    String, file with a column of IIDs in genfile
#' @param outfile       Path to output result file
#' @param genfile_format (Optional) String, format of genotype file. Default is 'bfile'
#' @param snpbatch_size (Optional) Number of SNPs per batch of analysis to be written to outfile. Default is 20000.
#' @param thresh_maf    (Optional) Threshold of minimum MAF to be analyzed. Default is 0.0001.
#' @param thresh_mac    (Optional) Threshold of minimum MAC to be analyzed. Default is 0.5.
#' @param espa_nt       (Optional) Number of points in empirical saddlepoint approximation, used only if null object is from GMMAT. Default is 9999.
#' @param espa_range    (Optional) Range of support in approximation, used only if null object is from GMMAT. Default is c(-20,20). 
#' @param bgi_file      (Optional) Path or Suffix type for bgi file. 1="genfile.bgen.bgi", 2="genfile.bgi". Default is 1.
#' @param snplist_file  (Optional) String, file of a subset of snp ids to be included in test . Default is ''
#' @param verbose       (Optional) Boolean, indicator for output runtime. Default is FALSE.
#' @param n_memchunk    (Optional) Integer. Number of memory chunks for PLINK genotype. Default is 1
#' @return Number of tested variants, and test results at outfile 
#' @export
TAPEtestM <- function(null_object, genfile, samplefile, outfile, genfile_format="bfile", thresh_maf=0.0001, thresh_mac=0.5, snpbatch_size=20000, espa_nt=9999, espa_range=c(-20,20), 
                     bgi_file="1",snplist_file="", verbose=FALSE,n_memchunk=1){
  # GRAB counterpart: GRAB.Marker
  if(verbose){
    cat("Starting TAPEtest...\n")
  }
  
  
  # ==== validate null model ======
  if(!(class(null_object) %in% c("glmmkin", "glmmai","tape_null"))){
    stop("Object class not supported")
  # }else if(class(null_object) %in% c("glmmkin", "glmmai")){
  }else{
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
    null_object$cumul = cumul
    
    # K_org_emp<-approxfun(cumul[,1], cumul[,2], rule=2)
    # K_1_emp<-approxfun(cumul[,1], cumul[,3], rule=2)
    # K_2_emp<-approxfun(cumul[,1], cumul[,4], rule=2)
    # 
    # null_object$K_org_emp=K_org_emp
    # null_object$K_1_emp=K_1_emp
    # null_object$K_2_emp=K_2_emp
  } # end if(!(class(null_object) %in% c("glmmkin", "glmmai","tape_null")))
  
  # ===== result headers ======
  if(genfile_format=="bfile"){
    output_header <- c("CHR", "SNPID", "BP", "ALT", "REF", "PROPMISS", "MAF", "SCORE", "VAR", "PVAL", "PVAL_NORM")
  }else if(genfile_format=="bgen"){
    # output_header <- c("CHR","BP","SNPID","ALT","REF","AC","AF", "MAF", "SCORE", "VAR", "PVAL", "PVAL_NORM")
    output_header <- c("MARKER","INFO","MAF", "SCORE", "PVAL")
  }
  write(output_header,file = outfile, ncolumns = length(output_header),sep="\t")
  
  if(verbose){
    cat("Genotype format:", genfile_format, "\n")
  }
  
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
    re_geno=setgeno_fast(genfile, subSampleInGeno, memoryChunk=n_memchunk)
    n_nomissing=re_geno
    toc <- proc.time()
    time_init <- toc - tic
    if(verbose){
      cat("Proc time init: ",time_init,"\n")
    }
    
    # ==== * get pval ==== 
    cat("====== Testing ======\n")
    
    output <- c()
    
    for(piter in 1:n_variants){
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
      chr = as.numeric(info_onevariant$CHR)
      
      geno[idx_miss] <- maf*2
      if(maf>0.5){
        maf <- 1-maf
        geno <- 2-geno
        info_onevariant[,c("ALT","REF")] <- info_onevariant[,c("REF","ALT")]
      }
      
      
      # get TAPE test results
      # if(verbose){
      #   tic <- proc.time()
      #   set_nonzero = which(geno!=0)
      #   geno_mat <- as.matrix(geno)
      #   toc <- proc.time()
      #   time_pretest <- toc - tic
      #   cat("Proc time pretest: ",time_pretest,"\n")
      #   
      #   tic <- proc.time()
      #   tape_results <- Get_Pvalues_espa(object = null_object, G= geno_mat, set_nonzero = set_nonzero,cutoff = 2)
      #   toc <- proc.time()
      #   time_test <- toc - tic
      #   cat("Proc time test: ",time_test,"\n")
      # }else{
      #   set_nonzero = which(geno!=0)
      #   geno_mat <- as.matrix(geno)
      #   tape_results <- Get_Pvalues_espa(object = null_object, G= geno_mat, set_nonzero = set_nonzero,cutoff = 2)
      #   
      # } # end if(verbose)
      set_nonzero = which(geno!=0)
      geno_mat <- as.matrix(geno)
      tape_results <- Get_Pvalues_espa(object = null_object, G= geno_mat, set_nonzero = set_nonzero,cutoff = 2, N=N, AF=maf, loco_chr=chr)
      
      # append results: CHR, SNPID, BP, ALT, REF, PROPMISS, MAF, SCORE, VAR, PVAL, PVAL_NORM
      temp_row <- c(info_onevariant$CHR, info_onevariant$CHR, info_onevariant$SNPID, info_onevariant$BP,
                    info_onevariant$ALT, info_onevariant$REF,
                    propmiss, maf, tape_results$score, tape_results$var_score, tape_results$pval, tape_results$pval_nospa)
      output <- rbind(output, temp_row)
      
      if((n_variants_tested %% snpbatch_size == 0) | (piter==n_variants)){
        cat("----- At",piter," ------\n")
        cat("Number of variants tested: ", n_variants_tested, "\n")
        output = as.data.frame(output)
        write.table(output, outfile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE, na="NA")
        
        output <- c()
        gc(verbose=F)
      } # end if((n_variants_tested %% snpbatch_size == 0) | (piter==n_variants))
      
    } # end for piter in snps_in_batch
    
    
    
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
    
    # Update control parameters
    control = list(range = espa_range,
                   length.out = espa_nt,
                   # allele.order = "ref-first",
                   allele.order = "alt-first",
                   impute_method = "mean",  
                   missing_cutoff = 0.15,
                   min_maf_marker = thresh_maf,
                   min_mac_marker = thresh_mac,
                   nMarkersEachChunk = snpbatch_size,
                   SPA_Cutoff = 2,
                   omp_num_threads = 1)
    if(snplist_file!=""){
      if(!file.exists(snplist_file)){
        stop(paste0("Cannot find file: ",snplist_file,"..."))
      }
      control$IDsToIncludeFile = snplist_file
    }else{
      control$IDsToIncludeFile = NULL
    } # end if(snplist_file!="")
    AlleleOrder = control$allele.order
    genoType = "BGEN"
    nMarkersEachChunk = control$nMarkersEachChunk
    outIndex = 1
    
    # # ==== * match IIDs==== 
    cat("====== Matching IIDs ======\n")
    iid_in_geno = data.table:::fread(samplefile, header=F, stringsAsFactors=FALSE)$V1
    iid_in_model = null_object$fit0$data$IID
    # mapid_geno_to_model=match(iid_in_geno,iid_in_model)
    # N = sum(!is.na(mapid_geno_to_model))
    # cat(N, "samples matched in bgen file\n")
    # mapid_geno_to_model[is.na(mapid_geno_to_model)] = -10 
    # mapid_geno_to_model=mapid_geno_to_model-1
    
    # ==== * load genfile==== 
    # check the setting of control, if not specified, the default setting will be used
    # The following functions are in 'control.R'
    
    cat("====== Loading genotype ======\n")
    # set up an object for genotype
    ## in objGeno = setGenoInput(GenoFile, GenoFileIndex, SampleIDs, control)
    
    # load variant info db
    db_con <- RSQLite::dbConnect(RSQLite::SQLite(), file_bgi)
    on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)
    bgiData = dplyr::tbl(db_con, "Variant")
    bgiData = as.data.frame(bgiData)
    
    if(AlleleOrder == "alt-first"){
      markerInfo = bgiData[,c(1,2,3,6,5,7)]  # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
    }else if(AlleleOrder == "ref-first"){
      markerInfo = bgiData[,c(1,2,3,5,6,7)]  # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
    }

    colnames(markerInfo) = c("CHROM", "POS", "ID", "REF", "ALT","genoIndex")
    
    if(any(!iid_in_model %in% iid_in_geno)){
      stop("At least one sample from 'SampleIDs' are not in 'GenoFile' and 'GenoFileIndex'.")
    }
    SampleIDs = iid_in_model
    
    setBGENobjInCPP(file_bgen, file_bgi, as.character(iid_in_geno), as.character(iid_in_model), F, F, AlleleOrder) # in linkBgen.cpp
    
    # names_info_variants <- c("CHR","BP","SNPID","ALT","REF","AC","AF")
    # n_variants_tested <- 0
    # variant_index_vec <- c()
    # gc(verbose=F)
    
    anyInclude = FALSE
    markersInclude = c()

    if(!is.null(control$IDsToIncludeFile)){
      IDsToInclude = data.table::fread(control$IDsToIncludeFile, 
                                       header = F, colClasses = c("character"))
      if(ncol(IDsToInclude) != 1){
        stop("IDsToIncludeFile should include one column.")
      }
      IDsToInclude = IDsToInclude[,1]
      
      posRows = which(markerInfo$ID %in% IDsToInclude)
      if(length(posRows) != 0){
        markersInclude = c(markersInclude, markerInfo$ID[posRows])
      }
      anyInclude = TRUE
    } # end if(!is.null(control$IDsToIncludeFile))
    markersInclude = unique(markersInclude)
    
    if(anyInclude){
      markerInfo = subset(markerInfo, ID %in% markersInclude)
    }
    
    genoList = list(genoType = genoType, 
                    markerInfo = markerInfo, 
                    SampleIDs = SampleIDs, 
                    AlleleOrder = AlleleOrder, 
                    GenoFile = file_bgen, 
                    GenoFileIndex = file_bgi)
    
    ## end in objGeno = setGenoInput(GenoFile, GenoFileIndex, SampleIDs, control)
    
    n = length(SampleIDs)
    MarkerIDs = markerInfo$ID

    GenoMat = getGenoInCPP(genoType, markerInfo, n) # in linkBgen.cpp
    
    CHROM = markerInfo$CHROM
    genoIndex = markerInfo$genoIndex
    
    # all markers were split into multiple chunks, 
    genoIndexList = splitMarker(genoIndex, snpbatch_size, CHROM);
    nChunks = length(genoIndexList)
    
    cat("Number of all markers to test:\t", nrow(markerInfo), "\n")
    cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
    cat("Number of chunks for all markers:\t", nChunks, "\n")
    
    chrom = "InitialChunk"
    cat("====== Testing ======\n")
    for(i in outIndex:nChunks){
      tempList = genoIndexList[[i]]
      genoIndex = tempList$genoIndex
      tempChrom = tempList$chrom
      
      # set up objects that do not change for different variants
      if(tempChrom != chrom){
        # setMarker(NullModelClass, objNull, control, chrom)
        # in linkBgen.cpp

        setMarker_GlobalVarsInCPP(control$impute_method,
                                  control$missing_cutoff,
                                  control$min_maf_marker,
                                  control$min_mac_marker,
                                  control$omp_num_threads)
        obj.setMarker = setMarker.TAPE(null_object, control)
        chrom = tempChrom
      } # end if(tempChrom != chrom)
      

     # ==== * get pval ==== 
      # GRAB counterpart: mainmarker.TAPE
      print(paste0("(",Sys.time(),") ---- Analyzing Chunk ", i, "/", nChunks, ": chrom ", chrom," ---- "))
      output <- getPvalues_cpp(genoType, genoIndex)  # in linkBgen.cpp, GRAB counterpart: mainMarker.TAPE: mainMarkerInCPP
      
      markerVec = output$markerVec   # marker IDs
      infoVec = output$infoVec       # marker infomation: CHR:POS:REF:ALT
      altFreqVec = output$altFreqVec       # minor allele frequencies (freq of ALT if flip=F, freq of REF if flip=T)
      BetaVec = output$BetaVec       # beta for ALT if flip=F, beta for REF if flip=T
      seBetaVec = output$seBetaVec   # sebeta
      pvalVec = output$pvalVec;      # marker-level p-values
      zScoreVec = output$zScoreVec
      
      # resMarker = data.frame(Marker = markerVec,
      #                        Info = infoVec,
      #                        AltFreq = altFreqVec,
      #                        Beta = BetaVec,
      #                        seBeta = seBetaVec,
      #                        Pval = pvalVec)
      resMarker = data.frame(Marker = markerVec,
                             Info = infoVec,
                             AltFreq = altFreqVec,
                             Score = zScoreVec,
                             Pval = pvalVec)
      
      # write summary statistics to output file
      if(i == 1){
        data.table::fwrite(resMarker, outfile, quote = F, sep = "\t", append = F, col.names = T, na="NA")
        # write.table(matrix(c("GRAB.outIndex", "Please_do_not_modify_this_file.", "Marker", 
        #                      format(nMarkersEachChunk, scientific=F), 1), 
        #                    ncol = 1), 
        #             OutputFileIndex, col.names = F, row.names = F, quote = F, append = F)
      }else{
        data.table::fwrite(resMarker, outfile, quote = F, sep = "\t", append = T, col.names = F, na="NA")
      #   write.table(matrix(i, ncol = 1), 
      #               OutputFileIndex, col.names = F, row.names = F, quote = F, append = T)
      } # end if(i == 1)
    } # end for(i in outIndex:nChunks)
      
    
    
    ## original slow ver
    # output <- c()
    # for(piter in 1:n_variants){
    #   
    #   Gx = getDosage_bgen_noquery()
    #   
    #   info_onevariant <- as.vector(unlist(Gx$variants))[-4]
    #   chr = as.numeric(Gx$variants$CHR)
    #   AC = Gx$variants$AC
    #   AF = Gx$variants$AF
    #   maf <- min(AF, 1-AF)
    #   mac <- round(maf * N)
    #   
    #   # check data validity
    #   if(mac <= thresh_mac | maf<=thresh_maf ){
    #     next
    #   }
    #   n_variants_tested = n_variants_tested + 1
    #   geno = Gx$dosages
    #   
    #   set_nonzero = which(geno!=0)
    #   geno_mat <- as.matrix(geno)
    #   tape_results <- Get_Pvalues_espa(object = null_object, G= geno_mat, set_nonzero = set_nonzero,cutoff = 2,N=N, AF=AF, loco_chr=chr)
    #   
    #   
    #   # append results: "CHR","BP","SNPID","ALT","REF","AC","AF", MAF, SCORE, VAR, PVAL, PVAL_NORM
    #   temp_row <- c(info_onevariant, maf, tape_results$score, tape_results$var_score, tape_results$pval, tape_results$pval_nospa)
    #   output <- rbind(output, temp_row)
    #   
    #   if((n_variants_tested %% snpbatch_size == 0) | (piter==n_variants)){
    #     cat("----- At",piter," ------\n")
    #     cat("Number of variants tested: ", n_variants_tested, "\n")
    #     output = as.data.frame(output)
    #     write.table(output, outfile, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t", append=TRUE, na="NA")
    #     
    #     output <- c()
    #     gc(verbose=F)
    #   } # end if((n_variants_tested %% snpbatch_size == 0) | (piter==n_variants))
    # } # end for piter in snps_in_batch
    
    # closetestGenoFile_bgenDosage()
    
    
  } # end if genfile_format
  
  cat("Test results saved to:",outfile,"\n")
  # return(n_variants_tested)
  return(nrow(markerInfo))
}


setMarker.TAPE = function(objNull, control)
{
  cumul = objNull$cumul
  mresid = objNull$resid
  N = length(mresid)
  SPA_Cutoff = control$SPA_Cutoff
  var_ratio = objNull$var_ratio
  
  # The following function is in  linkBgen.cpp (TAPE) / Main.cpp (GRAB)
  setTAPEobjInCPP(cumul,
                  mresid,
                  N,
                  SPA_Cutoff,
                  var_ratio)
}

## split 'markerInfo' into multiple chunks, each of which includes no more than 'nMarkersEachChunk' markers
splitMarker = function(genoIndex, nMarkersEachChunk, CHROM)
{
  genoIndexList = list()
  iTot = 1;
  
  uCHROM = unique(CHROM)
  for(chrom in uCHROM){
    pos = which(CHROM == chrom)
    gIdx = genoIndex[pos]
    M = length(gIdx)
    
    idxStart = seq(1, M, nMarkersEachChunk)
    idxEnd = idxStart + nMarkersEachChunk - 1
    
    nChunks = length(idxStart)
    idxEnd[nChunks] = M
    
    for(i in 1:nChunks){
      idxMarker = idxStart[i]:idxEnd[i]
      genoIndexList[[iTot]] = list(chrom = chrom,
                                   genoIndex = gIdx[idxMarker])
      iTot = iTot + 1;
    }
  }
  
  return(genoIndexList)
}

