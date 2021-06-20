// C++ functions for score test with bgen file

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
// std::this_thread::sleep_for (std::chrono::seconds(1));

#include "linkBgen.h"
#include "BGEN.hpp"
#include "UTIL.hpp"
#include "TAPE.hpp"

// global objects for different genotype formats
static BGEN::BgenClass* ptr_gBGENobj = NULL;

// global objects for different analysis methods
static TAPE::TAPEClass* ptr_gTAPEobj = NULL;

// global variables for analysis
std::string g_impute_method;      // "mean", "minor", or "drop"
double g_missingRate_cutoff;
unsigned int g_omp_num_threads;
double g_marker_minMAF_cutoff;
double g_marker_minMAC_cutoff;
double g_region_maxMAF_cutoff;
unsigned int g_region_maxMarkers_cutoff;   // maximal number of markers in one chunk, only used for region-based analysis

// [[Rcpp::export]]
void setMarker_GlobalVarsInCPP(std::string t_impute_method,
                               double t_missing_cutoff,
                               double t_min_maf_marker,
                               double t_min_mac_marker,
                               unsigned int t_omp_num_threads)
{
  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_marker_minMAF_cutoff = t_min_maf_marker;
  g_marker_minMAC_cutoff = t_min_mac_marker;
  g_omp_num_threads = t_omp_num_threads;
}

// [[Rcpp::export]]
void setBGENobjInCPP(std::string t_bgenFileName,
                     std::string t_bgenFileIndex,
                     std::vector<std::string> t_SampleInBgen,
                     std::vector<std::string> t_SampleInModel,
                     bool t_isSparseDosageInBgen,
                     bool t_isDropmissingdosagesInBgen,
                     std::string t_AlleleOrder)
{
  ptr_gBGENobj = new BGEN::BgenClass(t_bgenFileName,
                                     t_bgenFileIndex,
                                     t_SampleInBgen,
                                     t_SampleInModel,
                                     t_isSparseDosageInBgen,
                                     t_isDropmissingdosagesInBgen,
                                     t_AlleleOrder);
  int n = ptr_gBGENobj->getN();
  std::cout << "Number of samples:\t" << n << std::endl;
}



// [[Rcpp::export]]
void setTAPEobjInCPP(arma::mat t_cumul,
                     arma::vec t_mresid,
                     int t_N,
                     double t_SPA_Cutoff,
                     double t_var_ratio)
{
  // in TAPE.cpp / .hpp
  ptr_gTAPEobj = new TAPE::TAPEClass(t_cumul,t_mresid,t_N,t_SPA_Cutoff,t_var_ratio);
}


// a unified function to get single marker from genotype file
arma::vec Unified_getOneMarker(std::string t_genoType,   // "PLINK", "BGEN"
                               uint64_t t_gIndex,        // different meanings for different genoType
                               std::string& t_ref,       // REF allele
                               std::string& t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                               std::string& t_marker,    // marker ID extracted from genotype file
                               uint32_t& t_pd,           // base position
                               std::string& t_chr,       // chromosome
                               double& t_altFreq,        // frequency of ALT allele
                               double& t_altCounts,      // counts of ALT allele
                               double& t_missingRate,    // missing rate
                               double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                               bool t_isOutputIndexForMissing,               // if true, output index of missing genotype data
                               std::vector<uint32_t>& t_indexForMissing,     // index of missing genotype data
                               bool t_isOnlyOutputNonZero,                   // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                               std::vector<uint32_t>& t_indexForNonZero)     // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
{
  arma::vec GVec;
  
  if(t_genoType == "BGEN"){
    bool isBoolRead;
    GVec = ptr_gBGENobj->getOneMarker(t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo, 
                                      t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero,
                                      isBoolRead);
  }
  
  return GVec;
}

// [[Rcpp::export]]
arma::mat getGenoInCPP(std::string t_genoType,
                       Rcpp::DataFrame t_markerInfo,
                       int n)
{
  int q = t_markerInfo.nrow();         // number of markers requested
  std::vector<uint64_t> gIndexVec = t_markerInfo["genoIndex"];
  arma::mat GMat(n, q);
  
  std::string ref, alt, marker, chr;
  uint32_t pd;
  double altFreq, altCounts, missingRate, imputeInfo;
  std::vector<uint32_t> indexForMissing, indexForNonZero;
  
  for(int i = 0; i < q; i++){
    uint64_t gIndex = gIndexVec.at(i);
    arma::vec GVec = Unified_getOneMarker(t_genoType,    // "PLINK", "BGEN"
                                          gIndex,        // different meanings for different genoType
                                          ref,           // REF allele
                                          alt,           // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                                          marker,        // marker ID extracted from genotype file
                                          pd,            // base position
                                          chr,           // chromosome
                                          altFreq,       // frequency of ALT allele
                                          altCounts,     // counts of ALT allele
                                          missingRate,   // missing rate
                                          imputeInfo,    // imputation information score, i.e., R2 (all 1 for PLINK)
                                          false,         // if true, output index of missing genotype data
                                          indexForMissing,     // index of missing genotype data
                                          false,               // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                          indexForNonZero);    // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
    GMat.col(i) = GVec;
  }
  
  return GMat;
}

// a unified function to get marker-level p-value
void Unified_getMarkerPval(arma::vec t_GVec,
                           bool t_isOnlyOutputNonZero,
                           std::vector<uint32_t> t_indexForNonZero,
                           double& t_Beta, 
                           double& t_seBeta, 
                           double& t_pval,
                           double& t_zScore,
                           double t_altFreq)
{
  std::string t_method = "TAPE";
  if(t_method == "TAPE"){
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using TAPE method to calculate marker-level p-values, 't_isOnlyOutputNonZero' shold be false.");
    
    // in TAPE.cpp / .hpp
    ptr_gTAPEobj->getMarkerPval(t_GVec, t_Beta, t_seBeta, t_pval, t_altFreq, t_zScore);
  }
  
}

// [[Rcpp::export]]
Rcpp::List getPvalues_cpp(std::string t_genoType,
                          std::vector<uint32_t> t_genoIndex)  
{
  std::string t_method = "TAPE";
  int q = t_genoIndex.size();  // number of markers
  
  // set up output
  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT
  std::vector<std::string> markerVec(q);  // marker SNPIDs
  std::vector<std::string> refVec(q);    // marker information: REF
  std::vector<std::string> altVec(q);    // marker information: ALT
  std::vector<double> altFreqVec(q);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> altCountsVec(q);    // allele counts of ALT allele.
  std::vector<double> missingRateVec(q);  // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> BetaVec(q);         // beta value for ALT allele
  std::vector<double> seBetaVec(q);       
  std::vector<double> pvalVec(q);
  std::vector<double> zScoreVec(q);
  
  for(int i = 0; i < q; i++){
    
    if((i % 1000 == 0) || (i==q-1)){
      std::cout << "Completed " << i << "/" << q << " markers in the chunk." << std::endl;
    }
    
    // information of marker
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing, indexForNonZero;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;
    
    uint32_t gIndex = t_genoIndex.at(i);
    
    arma::vec GVec = Unified_getOneMarker(t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          true, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          false, // bool t_isOnlyOutputNonZero,
                                          indexForNonZero);
    
    int n = GVec.size();
    std::string info = chr+":"+std::to_string(pd)+":"+ref+":"+alt;
    
    // record basic information for the marker
    markerVec.at(i) = marker;               // marker IDs
    infoVec.at(i) = info;    // marker information: CHR:POS:REF:ALT
    altFreqVec.at(i) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    altCountsVec.at(i) = altCounts;         // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = missingRate;
    
    // MAF and MAC are for Quality Control (QC)
    double MAF = std::min(altFreq, 1 - altFreq);
    double MAC = MAF * n * (1 - missingRate);
    
    // Quality Control (QC) based on missing rate, MAF, and MAC
    if((missingRate > g_missingRate_cutoff) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff)){
      BetaVec.at(i) = arma::datum::nan;         
      seBetaVec.at(i) = arma::datum::nan;       
      pvalVec.at(i) = arma::datum::nan;
      continue;
    }
    
    if(missingRate != 0){
      // Function imputeGenoAndFlip is in UTIL.cpp
      flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing, g_impute_method);  // in UTIL.cpp
    }
    
    // analysis results for single-marker
    double Beta, seBeta, pval, zScore;
    
    Unified_getMarkerPval(GVec, 
                          false,
                          indexForNonZero, Beta, seBeta, pval, zScore, altFreq);
    
    BetaVec.at(i) = Beta * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true       
    seBetaVec.at(i) = seBeta;       
    pvalVec.at(i) = pval;
    zScoreVec.at(i) = zScore;
  }
  // }
  
  Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("markerVec") = markerVec,
                                          Rcpp::Named("infoVec") = infoVec,
                                          Rcpp::Named("altFreqVec") = altFreqVec,
                                          Rcpp::Named("altCountsVec") = altCountsVec,
                                          Rcpp::Named("missingRateVec") = missingRateVec,
                                          Rcpp::Named("BetaVec") = BetaVec,
                                          Rcpp::Named("seBetaVec") = seBetaVec,
                                          Rcpp::Named("pvalVec") = pvalVec,
                                          Rcpp::Named("zScoreVec") = zScoreVec);
  
  return OutList;  
}