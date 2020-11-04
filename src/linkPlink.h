#ifndef PKG_LINKPLINK_H
#define PKG_LINKPLINK_H


//This is a class with attritbutes about the genotype informaiton 
class genoClass{
private:
  	//COPY from RVTEST:
  	// we reverse the two bits as defined in PLINK format
  	const static unsigned char HOM_REF = 0x0;  // 0b00 ;
    const static unsigned char HET = 0x2;      // 0b10 ;
    const static unsigned char HOM_ALT = 0x3;  // 0b11 ;
    const static unsigned char MISSING = 0x1;  // 0b01 ;
public:
  	//to chunk the geno vector to avoid large continuous memory usage 
  	int numMarkersofEachArray;
    int numofGenoArray;
    int numMarkersofLastArray;
    std::vector< std::vector<unsigned char>* > genoVecofPointers;
	 
  	size_t M;
  	size_t N;
    size_t Nnomissing;
  	arma::fvec invstdvVec;
    std::vector<int> ptrsubSampleInGeno;
	
  	arma::fvec 	alleleFreqVec;
  	arma::fvec	m_OneSNP_Geno;
  	arma::fvec	m_OneSNP_StdGeno;
  	arma::fvec	m_DiagStd; // diagonal of GRM
    arma::fvec  m_DiagStd_LOCO;

  	unsigned char m_genotype_buffer[4];
    int geno_idx;
    int m_size_of_esi;
    unsigned char m_bits_val[8];

    // new in latest ver
    bool setKinDiagtoOne;
    arma::ivec	MACVec; //for variance ratio based on different MAC categories
    float minMAFtoConstructGRM = 0; // all markers used for grm construction
    int MwithMAFge_minMAFtoConstructGRM = 0;
    //LOCO
    size_t M_LOCO;
    int startid_LOCO;
    int endid_LOCO;
    int MwithMAFge_minMAFtoConstructGRM_LOCO = 0;

    // for submarker ops
    arma::ivec  subMarkerIndex;

    // funcs
    void setBit(unsigned char & ch, int ii, int aVal, int bVal);
    void setGenotype(unsigned char* c, const int pos, const int geno);
    void Init_OneSNP_Geno();
    arma::fvec * Get_OneSNP_Geno(size_t SNPIdx);
    arma::fvec * Get_OneSNP_Geno_atBeginning(size_t SNPIdx, std::vector<int> & indexNA, std::vector<unsigned char> & genoVecOneMarkerOld);
    int Get_OneSNP_StdGeno(size_t SNPIdx, arma::fvec * out );
    void Init_Diagof_StdGeno();
    void Init_Diagof_StdGeno_LOCO();

    //Function to assign values to all attributes
  	//This function is used instead of using a constructor because using constructor can not take
  	//genofile as an argument from .R 
    //genofile is the predix for plink bim, bed, fam, files   
  	void setGenoObj(std::string genofile, std::vector<int> subSampleInGeno, float memoryChunk, bool  isDiagofKinSetAsOne);

  	// param output funcs
  	void printFromgenoVec(unsigned char genoBinary0);
  	int getM() const{ return(M); }
    int getM_LOCO() const{ return(M_LOCO); }
  	int getN() const{ return(N); }
  	int getNnomissing() const{ return(Nnomissing); }
  	float getAC(int m){ return(alleleFreqVec[m]*2*Nnomissing); }
  	float getMAC(int m);
    int getMwithMAFge_minMAFtoConstructGRM() const{ return(MwithMAFge_minMAFtoConstructGRM); }
    arma::ivec getSubMarkerIndex(){ return(subMarkerIndex); }
    int getStartid_LOCO(){ return(startid_LOCO);}
    int getEndid_LOCO(){ return(endid_LOCO);}

  	//print out the vector of genotypes (first 100 indvs)
  	void printGenoVec();
  	//print out the vector of allele frequency
  	void printAlleleFreqVec();

};



struct crossProdGRM;
struct crossProdGRM_LOCO;




int setgeno(std::string genofile, std::vector<int> & subSampleInGeno, float memoryChunk, bool isDiagofKinSetAsOne);
void closegeno();
void setLOCOid(int startid_LOCO, int endid_LOCO);
int getNnomissingR();
float getMACR(int m);
arma::fvec getMAFR();
arma::fvec getDiagGRM();
arma::fvec getDiagGRM_LOCO();

// ======== Rcpp functions ref to geno object: pcgsolve ==============
arma::fvec getDiagOfSigma(const arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau);
arma::fvec getDiagOfSigma_nok(const arma::fvec& W, arma::fvec& tau);
arma::fvec getDiagOfSigma_LOCO(arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau);


float getMeanDiagOfGRM();

arma::fvec parallelCrossProdGRM(arma::fcolvec & bVec);
arma::fvec parallelCrossProdGRM_LOCO(arma::fcolvec & bVec);

arma::fvec getCrossprodGRM(arma::fcolvec& bVec);
arma::fvec getCrossprodGRM_LOCO(arma::fcolvec& bVec);

arma::fcolvec getCrossprodSig(arma::fcolvec& bVec, const arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau);
arma::fcolvec getCrossprodSig_nok(arma::fcolvec& bVec, const arma::fvec& W, arma::fvec& tau);
arma::fcolvec getCrossprodSig_LOCO(arma::fcolvec& bVec, arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau);

arma::fvec pcgsolve_plink(const arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau, arma::fvec & b, std::string preconditioner = "Jacobi", float tol = 1e-6, int maxIter = 1000);
arma::fvec pcgsolve_plink_LOCO(arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau, arma::fvec & b, std::string preconditioner = "Jacobi", float tol = 1e-6, int maxIter = 1000);
arma::fvec pcgsolve_plink_nok(const arma::fvec& W, arma::fvec& tau, arma::fvec & b, std::string preconditioner = "Jacobi", float tol = 1e-6, int maxIter = 1000);


arma::fmat pcgsolveMat_plink(arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau, arma::fmat & bMat, std::string preconditioner = "Jacobi", float tol = 1e-6, int maxIter = 1000);
arma::fmat pcgsolveMat_plink_LOCO(arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau, arma::fmat & bMat, std::string preconditioner = "Jacobi", float tol = 1e-6, int maxIter = 1000);
arma::fmat pcgsolveMat_plink_nok(const arma::fvec& W, arma::fvec& tau, arma::fmat & bMat, std::string preconditioner = "Jacobi", float tol = 1e-6, int maxIter = 1000);

// ======== functions ref to geno object: getTrace ==============
arma::fvec getTrace_plink(int q, arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec tau, arma::fmat& Sigma_iX, arma::fmat& XSiX_inv, int nrun_trace, int maxiter_pcg, float tol_pcg, float cutoff_trace);
arma::fvec getTrace_plink_nok(int q, const arma::fvec& W, arma::fvec tau, arma::fmat& Sigma_iX, arma::fmat& XSiX_inv, int nrun_trace, int maxiter_pcg, float tol_pcg, float cutoff_trace);

Rcpp::List getScore_plink(int n_nomissing, int q, arma::fvec& Y, arma::fmat& X, arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec tau, arma::fvec fixtau, bool is_AI, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace); 
Rcpp::List getScore_plink_nok(int n_nomissing, int q, arma::fvec& Y, arma::fmat& X, const arma::fvec& W, arma::fvec tau, arma::fvec fixtau, bool is_AI, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace); 
Rcpp::List getScore_plink_LOCO(int n_nomissing, int q, arma::fvec& Y, arma::fmat& X, arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec tau, arma::fvec fixtau, bool is_AI, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace); 

arma::fvec getOneSNPR(int i);
arma::fvec getOneSNPR_nostd(int i);









#endif