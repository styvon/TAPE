#include <RcppArmadillo.h>
// [[Rcpp:depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppParallel.h> 
//[[Rcpp::depends(RcppParallel)]]
#include <list>
#include <iterator>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <Rcpp/Benchmark/Timer.h>

#include "linkPlink.h"
#include "pcgUtil.h"

using namespace std;
using namespace Rcpp;
using namespace RcppParallel;


//create a geno object as a global variable
genoClass geno;

//http://gallery.rcpp.org/articles/parallel-inner-product/
struct crossProdGRM : public RcppParallel::Worker
{   
    // source vectors
    arma::fcolvec & m_bVec;
    unsigned int m_N;
    unsigned int m_M;

    // product that I have accumulated
    arma::fvec m_bout;
    int Msub_mafge1perc;  
  
    // constructors
    crossProdGRM(arma::fcolvec & y): m_bVec(y){
        m_M = geno.getM();
        m_N = geno.getNnomissing();
        m_bout.zeros(m_N);
        Msub_mafge1perc=0;
    } 
    crossProdGRM(const crossProdGRM& crossprodgrm, RcppParallel::Split): m_bVec(crossprodgrm.m_bVec){
          m_N = crossprodgrm.m_N;
          m_M = crossprodgrm.m_M;
          m_bout.zeros(m_N);
          Msub_mafge1perc=0;
    }
    // process just the elements of the range I've been asked to
    void operator()(std::size_t begin, std::size_t end) {
        arma::fvec vec;
        for(unsigned int i = begin; i < end; i++){
            if((geno.alleleFreqVec[i] > geno.minMAFtoConstructGRM) && (geno.alleleFreqVec[i] < 1-geno.minMAFtoConstructGRM)){
                geno.Get_OneSNP_StdGeno(i, &vec);
                float val1 = arma::dot(vec,  m_bVec);
                m_bout += val1 * (vec) ;
                Msub_mafge1perc += 1;
            }
        }
    }
  
    // join my value with that of another InnerProduct
	void join(const crossProdGRM & rhs) { 
		m_bout += rhs.m_bout;
	    Msub_mafge1perc += rhs.Msub_mafge1perc; 
	}  
};

//http://gallery.rcpp.org/articles/parallel-inner-product/
struct crossProdGRM_LOCO : public RcppParallel::Worker
{   
    // source vectors
    arma::fcolvec & m_bVec;
    unsigned int m_N;
    unsigned int m_M;
    unsigned int m_M_LOCO;
    int startid;
	int endid;

    // product that I have accumulated
    arma::fvec m_bout;
    int Msub_mafge1perc;  
  
    // constructors
    crossProdGRM_LOCO(arma::fcolvec & y): m_bVec(y){
        m_M = geno.getM();
        m_N = geno.getNnomissing();
        m_M_LOCO = geno.getM_LOCO();
        startid=geno.getStartid_LOCO();
        endid=geno.getEndid_LOCO();

        m_bout.zeros(m_N);
        Msub_mafge1perc=0;
    } 
    crossProdGRM_LOCO(const crossProdGRM_LOCO& crossprodgrm, RcppParallel::Split): m_bVec(crossprodgrm.m_bVec){
          m_N = crossprodgrm.m_N;
          m_M = crossprodgrm.m_M;
          m_M_LOCO = geno.getM_LOCO();
          startid=geno.getStartid_LOCO();
          endid=geno.getEndid_LOCO();

          m_bout.zeros(m_N);
          Msub_mafge1perc=0;
    }
    // process just the elements of the range I've been asked to
    void operator()(std::size_t begin, std::size_t end) {
        arma::fvec vec;
        float val1;
        for(unsigned int i = begin; i < end; i++){
            if((geno.alleleFreqVec[i] > geno.minMAFtoConstructGRM) && (geno.alleleFreqVec[i] < 1-geno.minMAFtoConstructGRM)){
                geno.Get_OneSNP_StdGeno(i, &vec);
                if(i >= startid && i <= endid){
                	val1 = 0;
                }else{
                	val1 = arma::dot(vec,  m_bVec);
                	Msub_mafge1perc += 1;
                } // end if within loco range             
                m_bout += val1 * (vec) ;
                
            } // end if within maf range
        } // end for snp
    }
  
    // join my value with that of another InnerProduct
	void join(const crossProdGRM_LOCO & rhs) { 
		m_bout += rhs.m_bout;
	    Msub_mafge1perc += rhs.Msub_mafge1perc; 
	}  
};



//==== member functions for genoClass =====
// geno class member function specification
void genoClass::setBit(unsigned char & ch, int ii, int aVal, int bVal){
	if (bVal == 1 && aVal == 1){
		ch ^= char(1 << ((ii*2) + 1)); //set a to be 1
	}else if(bVal == 0){
		ch ^= char(1 << (ii*2)); //change b to 0
		if(aVal == 1){
			ch ^= char(1 << ((ii*2) + 1)); //change a to 1
		}
	}
}

void genoClass::setGenotype(unsigned char* c, const int pos, const int geno) {
	(*c) |= (geno << (pos << 1));
}

void genoClass::Init_OneSNP_Geno(){
	m_size_of_esi = (Nnomissing+3)/4;
	int k = 8;
	while (k > 0){
		-- k;
		m_bits_val[k] = 1 << k;
	}
}

arma::fvec * genoClass::Get_OneSNP_Geno(size_t SNPIdx){
	m_OneSNP_Geno.zeros(Nnomissing);

	//avoid large continuous memory usage
	int indexOfVectorPointer = SNPIdx/numMarkersofEachArray;
    int SNPIdxinVec = SNPIdx % numMarkersofEachArray;
	////////////////

    size_t Start_idx = m_size_of_esi * SNPIdxinVec;
    size_t ind= 0;
    unsigned char geno1;
    int bufferGeno;
    for(size_t i=Start_idx; i< Start_idx+m_size_of_esi; i++){
    	geno1 = genoVecofPointers[indexOfVectorPointer]->at(i); //avoid large continuous memory usage
    	for(int j=0; j<4; j++){
    		int b = geno1 & 1 ;
            geno1 = geno1 >> 1;
            int a = geno1 & 1 ;
			bufferGeno = 2-(a+b);
			m_OneSNP_Geno[ind] = bufferGeno;
            ind++;
            geno1 = geno1 >> 1;
            if(ind >= Nnomissing){
            	return & m_OneSNP_Geno;
            }
        }
    }
    return & m_OneSNP_Geno;
}

arma::fvec * genoClass::Get_OneSNP_Geno_atBeginning(size_t SNPIdx, std::vector<int> & indexNA, std::vector<unsigned char> & genoVecOneMarkerOld){
	arma::fvec m_OneSNP_GenoTemp;
	m_OneSNP_GenoTemp.zeros(N);
	m_OneSNP_Geno.zeros(Nnomissing);
	int m_size_of_esi_temp = (N+3)/4;
	size_t ind= 0;
	unsigned char geno1;
	int bufferGeno;
	for(int i=0; i< m_size_of_esi_temp; i++){
		geno1 = genoVecOneMarkerOld[i];
		for(int j=0; j<4; j++){
			int b = geno1 & 1 ;
			geno1 = geno1 >> 1;
			int a = geno1 & 1 ;

			if (b == 1 && a == 0){
				bufferGeno = 3;
            }else if(b == 0 && a == 0){
            	bufferGeno = 2;
            }else if(b == 0 && a == 1){
            	bufferGeno = 1;
            }else if(b == 1 && a == 1){
            	bufferGeno = 0;
            }else{
                cout << "Error: Get_OneSNP_Geno_atBeginning: Invalid geno coding\n";
                break;
            }

            m_OneSNP_GenoTemp[ind] = bufferGeno;
			ind++;
            geno1 = geno1 >> 1;

            if(ind >= N){
            	int indxInOut = 0;
            	for(size_t indx=0; indx < Nnomissing; indx++){
            		m_OneSNP_Geno[indxInOut] = m_OneSNP_GenoTemp[ptrsubSampleInGeno[indx] - 1];
            		if(m_OneSNP_Geno[indxInOut] == 3){
            			indexNA.push_back(indxInOut);
            		}
            		indxInOut = indxInOut + 1;
            	}
            	return & m_OneSNP_Geno;
            } // end if(ind >= N)
        } // end for j	
	} // end for i

	return & m_OneSNP_Geno;
}

int genoClass::Get_OneSNP_StdGeno(size_t SNPIdx, arma::fvec * out ){
	//avoid large continuous memory usage
	int indexOfVectorPointer = SNPIdx/numMarkersofEachArray;
	int SNPIdxinVec = SNPIdx % numMarkersofEachArray;
	////////////////

	out->zeros(Nnomissing);

	size_t Start_idx = m_size_of_esi * SNPIdxinVec;
	size_t ind= 0;
	unsigned char geno1;

	float freq = alleleFreqVec[SNPIdx];
	float invStd = invstdvVec[SNPIdx];
	for(size_t i=Start_idx; i< Start_idx+m_size_of_esi; i++){
		geno1 = genoVecofPointers[indexOfVectorPointer]->at(i);
		for(int j=0; j<4; j++){
			int b = geno1 & 1 ;
			geno1 = geno1 >> 1;
			int a = geno1 & 1 ;
			(*out)[ind] = ((2-(a+b)) - 2*freq)* invStd;
			ind++;
			geno1 = geno1 >> 1;

			if(ind >= Nnomissing){
				return 1;
			}
		} // end for j
	} // end for i
	return 1;
}



void genoClass::Init_Diagof_StdGeno(){

	arma::fvec * temp = &m_OneSNP_StdGeno;
	int MminMAF = geno.getMwithMAFge_minMAFtoConstructGRM();

	if(arma::size(m_DiagStd)[0] != Nnomissing){
		m_DiagStd.zeros(Nnomissing);
		for(size_t j=0; j< M; j++){
			Get_OneSNP_StdGeno(j, temp);
			// temp = Get_OneSNP_Geno(j);
			m_DiagStd = m_DiagStd + (*temp) % (*temp);
		} // end for j
		for(size_t i=0; i<Nnomissing; i++){
			// m_DiagStd(i) = m_DiagStd(i)/float(MminMAF);
			if(m_DiagStd(i) < 1e-4){
				m_DiagStd(i) = 1e-4 ;
			}
		} // end for i
	} // end if(size(m_DiagStd)[0] != Nnomissing)

	// return & m_DiagStd;
}

void genoClass::Init_Diagof_StdGeno_LOCO(){

	arma::fvec * temp = &m_OneSNP_StdGeno;
	int MminMAF = geno.getM_LOCO();

	if(arma::size(m_DiagStd_LOCO)[0] != Nnomissing){
		m_DiagStd_LOCO.zeros(Nnomissing);
		for(size_t j=0; j< M; j++){
			if(j < startid_LOCO || j > endid_LOCO){
				Get_OneSNP_StdGeno(j, temp);
				// temp = Get_OneSNP_Geno(j);
				m_DiagStd_LOCO = m_DiagStd_LOCO + (*temp) % (*temp);
			} // end if j out of left-out chr
		} // end for j
		for(size_t i=0; i<Nnomissing; i++){
			// m_DiagStd(i) = m_DiagStd(i)/float(MminMAF);
			if(m_DiagStd_LOCO(i) < 1e-4){
				m_DiagStd_LOCO(i) = 1e-4 ;
			}
		} // end for i
	} // end if(size(m_DiagStd)[0] != Nnomissing)
}

//Function to assign values to all attributes
//This function is used instead of using a constructor because using constructor can not take
//genofile as an argument from .R 
//genofile is the predix for plink bim, bed, fam, files   
void genoClass::setGenoObj(std::string genofile, std::vector<int> subSampleInGeno, float memoryChunk, bool isDiagofKinSetAsOne){

	setKinDiagtoOne = isDiagofKinSetAsOne;   
	ptrsubSampleInGeno = subSampleInGeno;
	Nnomissing = subSampleInGeno.size(); 
	
	// reset
	alleleFreqVec.clear();
	MACVec.clear();
	invstdvVec.clear();
	M=0;
	N=0;
   	
	std::string bedfile = genofile+".bed";
	std::string bimfile = genofile+".bim"; 
	std::string famfile = genofile+".fam"; 
	std::string junk;

	//count the number of individuals
	ifstream test_famfile;
	test_famfile.open(famfile.c_str());
	if (!test_famfile.is_open()){
		printf("Error: setGenoObj: cannot open fam file");
		return ;
	}
	int indexRow = 0;
	while (std::getline(test_famfile,junk)){
    	indexRow = indexRow + 1;
    	junk.clear();
    }
	N = indexRow;
	test_famfile.clear();


	//count the number of markers
	ifstream test_bimfile;
	test_bimfile.open(bimfile.c_str());
	if (!test_bimfile.is_open()){
		printf("Error: setGenoObj: cannot open bim file");
		return ;
	}
	indexRow = 0;
	while (std::getline(test_bimfile,junk)){
        	indexRow = indexRow + 1;
        	junk.clear();
	}
	M = indexRow;
	test_bimfile.clear(); 

	junk.clear();

	// Init OneSNP Geno
	Init_OneSNP_Geno();
	indexRow = 0;
	std::vector<unsigned char> genoVecOneMarkerOld;
	std::vector<unsigned char> genoVecOneMarkerNew;
	/////////////////////////////

	// Added for reserve for genoVec
	size_t nbyteOld = ceil(float(N)/4);
	size_t nbyteNew = ceil(float(Nnomissing)/4);
	size_t reserve = ceil(float(Nnomissing)/4) * M + M*2;
	cout << "nbyte: " << nbyteOld << endl;
	cout << "nbyte: " << nbyteNew << endl;		
	cout << "reserve: " << reserve << endl;		

	genoVecOneMarkerOld.reserve(nbyteOld);
	genoVecOneMarkerOld.resize(nbyteOld);

	ifstream test_bedfile;
	test_bedfile.open(bedfile.c_str(), ios::binary);
	if (!test_bedfile.is_open()){
		printf("Error: setGenoObj: cannot open  bed file");
		return;
	}
	printf("\nM: %zu, N: %zu\n", M, N);

	//set up the array of vectors for genotype
	numMarkersofEachArray = floor((memoryChunk*pow (10.0, 9.0))/(ceil(float(N)/4)));
	//cout << "numMarkersofEachArray: " << numMarkersofEachArray << endl;
	if(M % numMarkersofEachArray == 0){
		numofGenoArray = M / numMarkersofEachArray;
		genoVecofPointers.resize(numofGenoArray);
        //cout << "size of genoVecofPointers: " << genoVecofPointers.size() << endl;
        for (int i = 0; i < numofGenoArray ; i++){
        	genoVecofPointers[i] = new std::vector<unsigned char>;
        	genoVecofPointers[i]->reserve(numMarkersofEachArray*ceil(float(N)/4));
        }
    }else{
    	numofGenoArray = M/numMarkersofEachArray + 1;
    	genoVecofPointers.resize(numofGenoArray);
    	numMarkersofLastArray = M - (numofGenoArray-1)*numMarkersofEachArray;
    	// cout << "size of genoVecofPointers: " << genoVecofPointers.size() << endl;
    	try{
    		for (int i = 0; i < numofGenoArray-1; i++){
    			genoVecofPointers[i] = new std::vector<unsigned char>;
    			genoVecofPointers[i]->reserve(numMarkersofEachArray*ceil(float(N)/4));
				//cout <<((*genoVecofPointers[i]).capacity()==numMarkersofEachArray*ceil(float(N)/4))<< endl;
    		}
    		genoVecofPointers[numofGenoArray-1] = new std::vector<unsigned char>;
			genoVecofPointers[numofGenoArray-1]->reserve(numMarkersofLastArray*ceil(float(N)/4));
		}
		catch(std::bad_alloc& ba){
			std::cerr << "setGenoObj: bad_alloc caught1: " << ba.what() << '\n';
			exit(EXIT_FAILURE);
		} // end catch
	} // end if(M % numMarkersofEachArray == 0)

	cout << "setGenoObj 1: setting genotype array pointers...";
	alleleFreqVec.zeros(M);
	invstdvVec.zeros(M);
	MACVec.zeros(M);
	float freq, Std, invStd;
	std::vector<int> indexNA;
	int lengthIndexNA;
	int indexGeno;
	int fillinMissingGeno;
	int b2;
	int a2;

    unsigned char geno1 = 0;
    int bufferGeno;
    int u;
    cout <<  "complete" << endl;

    cout << "setGenoObj 2: scanning file & NA filling...";
    for(size_t i = 0; i < M; i++){
    	genoVecOneMarkerOld.clear();
    	genoVecOneMarkerOld.reserve(nbyteOld);
    	genoVecOneMarkerOld.resize(nbyteOld);

    	test_bedfile.seekg(3+nbyteOld*i);
    	test_bedfile.read((char*)(&genoVecOneMarkerOld[0]),nbyteOld);

      	indexNA.clear();
      	Get_OneSNP_Geno_atBeginning(i, indexNA, genoVecOneMarkerOld);

		geno1 = 0;

		for(unsigned int j=0; j< Nnomissing; j++){
			u = j & (4 - 1);
			bufferGeno = m_OneSNP_Geno[j];
			if(bufferGeno == 0){
				setGenotype(&geno1, u, HOM_ALT);
			}else if(bufferGeno == 1){
				setGenotype(&geno1, u, HET);
			}else if(bufferGeno == 2){	
				setGenotype(&geno1, u, HOM_REF);	
			}else{
				setGenotype(&geno1, u, MISSING);
				m_OneSNP_Geno[j] = 0;  //12-18-2017 	
			}

			if(u == 3){
				genoVecofPointers[i/numMarkersofEachArray]->push_back(geno1); //avoid large continuous memory usage
				geno1 = 0;
			}
		} // end for j
				
		if(Nnomissing%4 != 0){

			genoVecofPointers[i/numMarkersofEachArray]->push_back(geno1); //avoid large continuous memory usage
			geno1 = 0;
		}

		lengthIndexNA = indexNA.size();
		freq = sum(m_OneSNP_Geno)/(2*(Nnomissing-lengthIndexNA));

		if (lengthIndexNA > 0){
			fillinMissingGeno = int(round(2*freq));
			if(fillinMissingGeno == 0){
				b2 = 1;
				a2 = 1;
        	}else if(fillinMissingGeno == 1){
                b2 = 0;
                a2 = 1;
        	}else{
                b2 = 0;
                a2 = 0;
        	}

        	for (int k=0; k<lengthIndexNA; k++){
        		indexGeno = indexNA[k];
        		m_OneSNP_Geno[indexGeno] = fillinMissingGeno;
        		setBit(genoVecofPointers[i/numMarkersofEachArray]->at((i%numMarkersofEachArray)*nbyteNew+(indexGeno/4)),indexGeno%4, a2, b2);
    		} // end for k
		}// end if (lengthIndexNA > 0)

        freq = float(sum(m_OneSNP_Geno))/(2*Nnomissing);
		Std = std::sqrt(2*freq*(1-freq));
		if(Std == 0){
			invStd= 0;
		}else {
			invStd= 1.0/Std;
		}

		alleleFreqVec[i] = freq;
		if(minMAFtoConstructGRM > 0){
			if(freq > minMAFtoConstructGRM && freq < (1-minMAFtoConstructGRM)){
				MwithMAFge_minMAFtoConstructGRM = MwithMAFge_minMAFtoConstructGRM + 1;
			}
		}else{
			MwithMAFge_minMAFtoConstructGRM = M;
		}// end if(minMAFtoConstructGRM > 0)	
			
		if(freq > 0.5){
		  MACVec[i] = (2*Nnomissing) - sum(m_OneSNP_Geno);
		}else{
		  MACVec[i] = sum(m_OneSNP_Geno);
		}
		
		invstdvVec[i] = invStd;
		m_OneSNP_Geno.clear();

	}//end for(int i = 0; i < M; i++)
	cout << MwithMAFge_minMAFtoConstructGRM << " markers with MAF > " << minMAFtoConstructGRM << " are used for GRM." << endl;

	test_bedfile.close();
	cout << "setGenoObj: completed" << endl;
}//end func setGenoObj

void genoClass::printFromgenoVec(unsigned char genoBinary0){
	unsigned char genoBinary = genoBinary0;
	for(int j=0; j<4; j++){
		int b = genoBinary & 1 ;
		genoBinary = genoBinary >> 1;
		int a = genoBinary & 1 ;
		genoBinary = genoBinary >> 1;
		cout << 2-(a+b) << " " << endl;
	}
	cout << endl;
}

float genoClass::getMAC(int m){
	if(alleleFreqVec[m] > 0.5){
		return((1-alleleFreqVec[m])*2*Nnomissing);
	}else{
		return(alleleFreqVec[m]*2*Nnomissing);
	}
}

//print out the vector of genotypes
void genoClass::printGenoVec(){
	for(unsigned int i=0; i<2; ++i){
		Get_OneSNP_Geno(i);
		for(unsigned int j=0; j< 100; j++){
			cout << m_OneSNP_Geno[j] << ' ';
		}
		cout << endl;
	}
	cout << "M = " << M << endl;
	cout << "N = " << N << endl;
}

//print out the vector of allele frequency
void genoClass::printAlleleFreqVec(){
	//for(int i=0; i<alleleFreqVec.size(); ++i)
	for(size_t i=(M-100); i<M; ++i){
		cout << alleleFreqVec[i] << ' ';
	}
	cout << endl;
}
// ==== end member functions for genoClass ======




// [[Rcpp::export]]
int setgeno(std::string genofile, std::vector<int> & subSampleInGeno, float memoryChunk, bool isDiagofKinSetAsOne){
	// int start_s=clock();
	Rcpp::Timer timer;
	timer.step("start");
	geno.setGenoObj(genofile, subSampleInGeno, memoryChunk, isDiagofKinSetAsOne);
	timer.step("setgeno");
	geno.Init_Diagof_StdGeno();
	timer.step("init_diag");
	NumericVector re_time(timer), re_time2;
	re_time2 = diff(re_time);
	Rcpp::Rcout << "Runtime : " << "\n";
	Rcpp::Rcout << "setgeno\tinit_diag" << "\n";
	Rcpp::Rcout << re_time2 << "\n";
	Rcpp::Rcout << "----------------\n\n";
	//geno.printAlleleFreqVec();
	//geno.printGenoVec();
	// int stop_s=clock();
	// cout << "time: " << (stop_s-start_s)/float(CLOCKS_PER_SEC)*1000 << endl;
	int out = geno.getNnomissing();
	return(out);
}

// [[Rcpp::export]]
int setgeno_fast(std::string genofile, std::vector<int> & subSampleInGeno, float memoryChunk){
	geno.setGenoObj(genofile, subSampleInGeno, memoryChunk, TRUE);
	int out = geno.getNnomissing();
	return(out);
}

// [[Rcpp::export]]
void closegeno(){
	for (int i = 0; i < geno.numofGenoArray; i++){
		(*geno.genoVecofPointers[i]).clear();	
    		delete geno.genoVecofPointers[i];
  	}

  	geno.genoVecofPointers.clear();

  	geno.invstdvVec.clear();
  	geno.ptrsubSampleInGeno.clear();
  	geno.alleleFreqVec.clear();
  	geno.m_OneSNP_Geno.clear();
  	geno.m_OneSNP_StdGeno.clear();
  	geno.m_DiagStd.clear();
  	printf("PLINK file closed\n");
}


// [[Rcpp::export]]
void setLOCOid(int startid_LOCO, int endid_LOCO){
  geno.startid_LOCO = startid_LOCO;
  geno.endid_LOCO = endid_LOCO;
  geno.M_LOCO = 0;
  for(size_t i=0; i< geno.M; i++){
	if(i < startid_LOCO || i > endid_LOCO){
  		if(geno.alleleFreqVec[i] > geno.minMAFtoConstructGRM && geno.alleleFreqVec[i] < 1-geno.minMAFtoConstructGRM){     
			geno.M_LOCO = geno.M_LOCO + 1;
  		}
	}
  } // end for i
  geno.Init_Diagof_StdGeno_LOCO();
}

// [[Rcpp::export]]
int getNnomissingR(){
	int n = geno.getNnomissing();
	return(n);
}

// [[Rcpp::export]]
float getMACR(int m){
	float out = geno.getMAC(m);
	return(out);
}

// [[Rcpp::export]]
arma::fvec getMAFR(){
	arma::fvec out = geno.alleleFreqVec;
	return(out);
}

// [[Rcpp::export]]
arma::fvec getDiagGRM(){
	arma::fvec out = geno.m_DiagStd;
	return(out);
}

// [[Rcpp::export]]
arma::fvec getDiagGRM_LOCO(){
	arma::fvec out = geno.m_DiagStd_LOCO;
	return(out);
}

// ======== functions ref to geno object ==============

//' Get diagonal elements from Sigma
//' Sigma = tau[0] * diag(1/W) + tau[1] * sparsekin + tau[2] * grm 
//' 
//' @param sparsekin  matrix of kinship coefficients (sparse).
//' @param W        numeric vector.
//' @param tau      numeric vector, initial values for variance components.
//' @return A vector with diagonal elements of Sigma
// [[Rcpp::export]]
arma::fvec getDiagOfSigma(const arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau){
  
	size_t Nnomissing = geno.getNnomissing();
	int MminMAF = geno.getMwithMAFge_minMAFtoConstructGRM();

	arma::fvec diagVec(Nnomissing);
	arma::fvec oneVec(Nnomissing);
	oneVec.ones();

	
	if(!(geno.setKinDiagtoOne)){ 
		// conversion for ops with sp_mat term: tau(1)*sparsekin.diag()
		// diagVec = tau(0)/W + tau(1)*sparsekin.diag() + tau(2)* (*geno.getDiagGRM()) /MminMAF;
		// diagVec = tau(0)/W + tau(1)*sparsekin.diag() + tau(2)* (*geno.getDiagGRM());
		diagVec = tau(0)/W + tau(1)*oneVec + tau(2)*geno.m_DiagStd/float(MminMAF);


	}else{
		// conversion for ops with sp_mat term: tau(1)*sparsekin.diag()=tau(1)
		// diagVec = tau(0)/W+ tau(1)*oneVec + tau(2)/MminMAF ;
		// diagVec = tau(0)/W+ tau(1)*oneVec + tau(2) ;
		diagVec = tau(0)/W+ tau(1)*oneVec + tau(2)*geno.m_DiagStd/float(MminMAF);
		
	}
	
	for(size_t i=0; i< Nnomissing; i++){
		if(diagVec(i) < 1e-4){
			diagVec(i) = 1e-4 ;
		}
	}
	return(diagVec);
}

//' Get diagonal elements from Sigma (LOCO)
//' Sigma = tau[0] * diag(1/W) + tau[1] * sparsekin + tau[2] * grm 
//' 
//' @param sparsekin  matrix of kinship coefficients (sparse).
//' @param W        numeric vector.
//' @param tau      numeric vector, initial values for variance components.
//' @return A vector with diagonal elements of Sigma
// [[Rcpp::export]]
arma::fvec getDiagOfSigma_LOCO(arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau){
  
	size_t Nnomissing = geno.getNnomissing();
	int MminMAF = geno.getM_LOCO();

	arma::fvec diagVec(Nnomissing);
	arma::fvec oneVec(Nnomissing);
	oneVec.ones();

	
	if(!(geno.setKinDiagtoOne)){ 
		// conversion for ops with sp_mat term: tau(1)*sparsekin.diag()
		// diagVec = tau(0)/W + tau(1)*sparsekin.diag() + tau(2)* (*geno.getDiagGRM()) /MminMAF;
		// diagVec = tau(0)/W + tau(1)*sparsekin.diag() + tau(2)* (*geno.getDiagGRM());
		diagVec = tau(0)/W + tau(1)*oneVec + tau(2)*geno.m_DiagStd_LOCO/float(MminMAF);


	}else{
		// conversion for ops with sp_mat term: tau(1)*sparsekin.diag()=tau(1)
		// diagVec = tau(0)/W+ tau(1)*oneVec + tau(2)/MminMAF ;
		// diagVec = tau(0)/W+ tau(1)*oneVec + tau(2) ;
		diagVec = tau(0)/W+ tau(1)*oneVec + tau(2)*geno.m_DiagStd_LOCO/float(MminMAF);
		
	}
	
	for(size_t i=0; i< Nnomissing; i++){
		if(diagVec(i) < 1e-4){
			diagVec(i) = 1e-4 ;
		}
	}
	return(diagVec);
}



//' Get diagonal elements from Sigma
//' without sparse kinship matrix
//' Sigma = tau[0] * diag(1/W) + tau[2] * grm 
//' 
//' @param W        numeric vector.
//' @param tau      numeric vector, initial values for variance components.
//' @return A vector with diagonal elements of Sigma
// [[Rcpp::export]]
arma::fvec getDiagOfSigma_nok(const arma::fvec& W, arma::fvec& tau){
  
	size_t Nnomissing = geno.getNnomissing();
	int MminMAF = geno.getMwithMAFge_minMAFtoConstructGRM();

	arma::fvec diagVec(Nnomissing);
	arma::fvec oneVec(Nnomissing);
	oneVec.ones();

	diagVec = tau(0)/W + tau(1)* geno.m_DiagStd/float(MminMAF);
	for(size_t i=0; i< Nnomissing; i++){
		if(diagVec(i) < 1e-4){
			diagVec(i) = 1e-4 ;
		}
	}
	return(diagVec);
}

//' Get mean of diagonal elements from Genetic relationship matrix
// [[Rcpp::export]]
float getMeanDiagOfGRM(){
  
	size_t Nnomissing = geno.getNnomissing();
	// int MminMAF = geno.getMwithMAFge_minMAFtoConstructGRM();

	arma::fvec diagVec(Nnomissing);
	diagVec = getDiagGRM();

	for(size_t i=0; i< Nnomissing; i++){
		if(diagVec(i) < 1e-4){
			diagVec(i) = 1e-4 ;
		}
	}

	float diagMean = arma::mean(diagVec);	
	return(diagMean);
}

//' Get product of GRM and vector parallel worker
//' 
//' @param bVec       numeric vector, used for product
//' @return vector
// [[Rcpp::export]]
arma::fvec parallelCrossProdGRM(arma::fcolvec & bVec) {
  
    // declare the InnerProduct instance that takes a pointer to the vector data
    int M = geno.getM();
    crossProdGRM crossprodgrm(bVec);    
    // call paralleReduce to start the work
    RcppParallel::parallelReduce(0, M, crossprodgrm); 
    // return crossprodgrm.m_bout/(crossprodgrm.Msub_mafge1perc);
    return crossprodgrm.m_bout/M;


}

//' Get product of GRM and vector parallel worker (LOCO)
//' 
//' @param bVec       numeric vector, used for product
//' @return vector
// [[Rcpp::export]]
arma::fvec parallelCrossProdGRM_LOCO(arma::fcolvec & bVec) {
  
    // declare the InnerProduct instance that takes a pointer to the vector data
    int M_LOCO = geno.getM_LOCO();
    crossProdGRM_LOCO crossprodgrm(bVec);    
    // call paralleReduce to start the work
    RcppParallel::parallelReduce(0, M_LOCO, crossprodgrm); 
    // return crossprodgrm.m_bout/(crossprodgrm.Msub_mafge1perc);
    return crossprodgrm.m_bout/M_LOCO;


}

//' Get product of GRM and vector
//' 
//' @param bVec       numeric vector, used for product
//' @return vector
// [[Rcpp::export]]
arma::fvec getCrossprodGRM(arma::fcolvec& bVec){
    arma::fvec crossProdVec;
    crossProdVec = parallelCrossProdGRM(bVec) ;
    return(crossProdVec);
}

//' Get product of GRM and vector (LOCO)
//' 
//' @param bVec       numeric vector, used for product
//' @return vector
// [[Rcpp::export]]
arma::fvec getCrossprodGRM_LOCO(arma::fcolvec& bVec){
    arma::fvec crossProdVec;
    crossProdVec = parallelCrossProdGRM_LOCO(bVec) ;
    return(crossProdVec);
}

//' Get product of Sigma and vector
//' Sigma = tau[0] * diag(1/W) + tau[1] * sparsekin + tau[2] * grm 
//' 
//' @param bVec       numeric vector, used for product
//' @param sparsekin  matrix of kinship coefficients (sparse).
//' @param W          numeric vector.
//' @param tau        numeric vector, initial values for variance components.
//' @return vector
// [[Rcpp::export]]
arma::fcolvec getCrossprodSig(arma::fcolvec& bVec, const arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau){
	int n = bVec.size();
    arma::fcolvec crossProdVec(n), temp_vec(n);
    arma::sp_mat bMat(n,1), temp_mat(n,1);
    // convert bVec to sparse
    for(int i=0; i<n; i++){
		bMat(i,0) = bVec(i);
	}
	// convert sparse matmult back to fvec
	temp_mat = sparsekin * bMat;
	for(int i=0; i<n; i++){
		temp_vec(i) = temp_mat(i,0);
	}
	// arma::fcolvec crossProdVec;
 //    arma::vec temp_vec_double;
 //    temp_vec_double = 0.0+sparsekin * arma::conv_to<arma::vec>::from(bVec);

    if(tau(2) == 0){
		// conversion for ops with sp_mat term: tau(1)*(bVec % sparsekin.diag())= tau(1)*bVec
		// crossProdVec = tau(0)*(bVec % (1.0/W)) + tau(1)* arma::conv_to<arma::fvec>::from(temp_vec_double); 
		crossProdVec = tau(0)*(bVec % (1.0/W)) + tau(1)* temp_vec; 
        return(crossProdVec);
    }
    arma::fvec crossProd1  = getCrossprodGRM(bVec);
    // conversion for ops with sp_mat term: tau(1)*(bVec % sparsekin.diag())= tau(1)*bVec
    crossProdVec = tau(0)*(bVec % (1.0/W)) + tau(1)* temp_vec + tau(2)*crossProd1; 
    // crossProdVec = tau(0)*(bVec % (1.0/W)) + tau(1)* arma::conv_to<arma::fvec>::from(temp_vec_double) + tau(2)*crossProd1; 

    return(crossProdVec);
}

//' Get product of Sigma and vector
//' Sigma = tau[0] * diag(1/W) + tau[1] * sparsekin + tau[2] * grm 
//' 
//' @param bVec       numeric vector, used for product
//' @param sparsekin  matrix of kinship coefficients (sparse).
//' @param W          numeric vector.
//' @param tau        numeric vector, initial values for variance components.
//' @return vector
// [[Rcpp::export]]
arma::fcolvec getCrossprodSig_LOCO(arma::fcolvec& bVec, arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau){
	int n = bVec.size();
    arma::fcolvec crossProdVec(n), temp_vec(n);
    arma::sp_mat bMat(n,1), temp_mat(n,1);
    // convert bVec to sparse
    for(int i=0; i<n; i++){
		bMat(i,0) = bVec(i);
	}
	// convert sparse matmult back to fvec
	temp_mat = sparsekin * bMat;
	for(int i=0; i<n; i++){
		temp_vec(i) = temp_mat(i,0);
	}
	// arma::fcolvec crossProdVec;
 //    arma::vec temp_vec_double;
 //    temp_vec_double = 0.0+sparsekin * arma::conv_to<arma::vec>::from(bVec); 

    if(tau(2) == 0){
		// conversion for ops with sp_mat term: tau(1)*(bVec % sparsekin.diag())= tau(1)*bVec
		crossProdVec = tau(0)*(bVec % (1.0/W)) + tau(1)* temp_vec; 
		// crossProdVec = tau(0)*(bVec % (1.0/W)) + tau(1)* arma::conv_to<arma::fvec>::from(temp_vec_double);
        return(crossProdVec);
    }
    arma::fvec crossProd1  = getCrossprodGRM_LOCO(bVec);

    // conversion for ops with sp_mat term: tau(1)*(bVec % sparsekin.diag())= tau(1)*bVec
    crossProdVec = tau(0)*(bVec % (1.0/W)) + tau(1)* temp_vec + tau(2)*crossProd1; 
    // crossProdVec = tau(0)*(bVec % (1.0/W)) + tau(1)* arma::conv_to<arma::fvec>::from(temp_vec_double) + tau(2)*crossProd1;

    return(crossProdVec);
}

//' Get product of Sigma and vector
//' without sparse kinship matrix
//' Sigma = tau[0] * diag(1/W) + tau[1] * grm 
//' 
//' @param bVec       numeric vector, used for product
//' @param W          numeric vector.
//' @param tau        numeric vector, initial values for variance components.
//' @return vector
// [[Rcpp::export]]
arma::fcolvec getCrossprodSig_nok(arma::fcolvec& bVec, const arma::fvec& W, arma::fvec& tau){

    arma::fcolvec crossProdVec;

    if(tau(1) == 0){
		crossProdVec = tau(0)*(bVec % (1.0/W)) ; 
        return(crossProdVec);
    }
    arma::fvec crossProd1  = getCrossprodGRM(bVec);

    // conversion for ops with sp_mat term: tau(1)*(bVec % sparsekin.diag())= tau(1)*bVec
    // crossProdVec = tau(0)*(bVec % (1/W)) + tau(1)*bVec + tau(2)*crossProd1;
    crossProdVec = tau(0)*(bVec % (1.0/W)) + tau(1)*crossProd1; 


    return(crossProdVec);
}


//' Preconditioned conjugate gradient method (from cPCG)
//'
//' Sigma = tau[0] * diag(1/W) + tau[1] * sparsekin + tau[2] * grm 
//' @title Solve for x in Sigma x = b using preconditioned conjugate gradient method.
//' @param sparsekin        matrix of kinship coefficients (sparse).
//' @param W                numeric vector.
//' @param tau              numeric vector, initial values for variance components.
//' @param b                vector, with same dimension as number of rows of A.
//' @param preconditioner   string, method for preconditioning: \code{"Jacobi"} (default), \code{"SSOR"}, or \code{"ICC"}.
//' @param tol              numeric, threshold for convergence, default is \code{1e-6}.
//' @param maxIter          numeric, maximum iteration, default is \code{1000}.
//' @return A vector representing solution x.
// [[Rcpp::export]]
arma::fvec pcgsolve_plink(const arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau, arma::fvec & b, std::string preconditioner, float tol, int maxIter) {

    // get number of columns of A=Sigma
    int C = geno.getNnomissing(); // num cols of Sigma
    int R = geno.getNnomissing(); // num rows of Sigma
    
    // get preconditioner M
    arma::fvec Mvec;
    if (preconditioner == "Jacobi"){
        // Mmat = arma::diagmat(A);
        Mvec = getDiagOfSigma(sparsekin, W, tau);
    }
    // to implement 
    // else if(preconditioner == "SSOR"){
    //   arma::fmat D = arma::diagmat(A);
    //   arma::fmat L = arma::trimatl(A);
    //   Mmat = (D+L) * D.i() * (D+L).t();
    // } else if(preconditioner == "ICC"){
    //   Mmat = icc(A);
    // }

  
    // initiate solution x as zeros
    arma::fvec x(C) ;
    x.zeros() ; 
    
    arma::fvec oneVec(C);
    oneVec.ones() ;
    
    // arma::fvec r = b - A * x;
    arma::fvec r = b;
    // arma::fmat Minv = M.i();
    arma::fvec Minv = 1.0/Mvec;

    arma::fvec z = Minv % r;
    arma::fvec p = z;
    float rz_old = arma::sum(r % z);
    float rz_new=1.0;
    // arma::fvec rz_ratio(1);

    
    arma::fvec Ap(R);
    float alpha, beta;
    // vector version of alpha
    // arma::fvec alphaVec(1);
    
    for(int iter = 0; (iter < maxIter) && (rz_new > tol); iter++){
        // Ap = A * p;
        Ap = getCrossprodSig(p, sparsekin, W, tau);
        alpha = rz_old / arma::sum(p % Ap); 
        x += alpha * p;
        r -= alpha * Ap;
        z = Minv % r;
        rz_new = arma::sum( z % r );
        beta = rz_new / rz_old; 

        p = z + beta * p;
        rz_old = rz_new;
        if (iter >= maxIter){
            Rcout << "pcg did not converge." << std::endl;
        }
    }

    
    return arma::round(x*10000.0)/10000.0;
  
} 

//' Preconditioned conjugate gradient method LOCO (from cPCG)
//'
//' Sigma = tau[0] * diag(1/W) + tau[1] * sparsekin + tau[2] * grm 
//' @title Solve for x in Sigma x = b using preconditioned conjugate gradient method.
//' @param sparsekin        matrix of kinship coefficients (sparse).
//' @param W                numeric vector.
//' @param tau              numeric vector, initial values for variance components.
//' @param b                vector, with same dimension as number of rows of A.
//' @param preconditioner   string, method for preconditioning: \code{"Jacobi"} (default), \code{"SSOR"}, or \code{"ICC"}.
//' @param tol              numeric, threshold for convergence, default is \code{1e-6}.
//' @param maxIter          numeric, maximum iteration, default is \code{1000}.
//' @return A vector representing solution x.
// [[Rcpp::export]]
arma::fvec pcgsolve_plink_LOCO(arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau, arma::fvec & b, std::string preconditioner, float tol, int maxIter) {

    // get number of columns of A=Sigma
    int C = geno.getNnomissing(); // num cols of Sigma
    int R = geno.getNnomissing(); // num rows of Sigma
    
    // get preconditioner M
    arma::fvec Mvec;
    if (preconditioner == "Jacobi"){
        // Mmat = arma::diagmat(A);
        Mvec = getDiagOfSigma_LOCO(sparsekin, W, tau);
    }
    // to implement 
    // else if(preconditioner == "SSOR"){
    //   arma::fmat D = arma::diagmat(A);
    //   arma::fmat L = arma::trimatl(A);
    //   Mmat = (D+L) * D.i() * (D+L).t();
    // } else if(preconditioner == "ICC"){
    //   Mmat = icc(A);
    // }

  
    // initiate solution x as zeros
    arma::fvec x(C) ;
    x.zeros() ; 
    
    arma::fvec oneVec(C);
    oneVec.ones() ;
    
    // arma::fvec r = b - A * x;
    arma::fvec r = b;
    // arma::fmat Minv = M.i();
    arma::fvec Minv = 1.0/Mvec;
    arma::fvec z = Minv % r;
    arma::fvec p = z;
    float rz_old = arma::sum(r % z);
    float rz_new=1.0;
    // arma::fvec rz_ratio(1);

    arma::fvec Ap(R);
    float alpha, beta;
    // vector version of alpha
    // arma::fvec alphaVec(1);
    
    for(int iter = 0; (iter < maxIter) && (rz_new > tol); iter++){
        // Ap = A * p;
        Ap = getCrossprodSig_LOCO(p, sparsekin, W, tau);
        alpha = rz_old / arma::sum(p % Ap); 
        x += alpha * p;
        r -= alpha * Ap;
        z = Minv % r;
        rz_new = arma::sum( z % r );
        beta = rz_new / rz_old; 

        p = z + beta * p;
        rz_old = rz_new;
        if (iter >= maxIter){
            Rcout << "pcg did not converge." << std::endl;
        }
    }

    return arma::round(x*10000.0)/10000.0;
  
} 

//' Preconditioned conjugate gradient method (from cPCG)
//' without sparse kinship matrix
//'
//' Sigma = tau[0] * diag(1/W) + tau[1] * sparsekin + tau[2] * grm 
//' @title Solve for x in Sigma x = b using preconditioned conjugate gradient method.
//' @param W                numeric vector.
//' @param tau              numeric vector, initial values for variance components.
//' @param b                vector, with same dimension as number of rows of A.
//' @param preconditioner   string, method for preconditioning: \code{"Jacobi"} (default), \code{"SSOR"}, or \code{"ICC"}.
//' @param tol              numeric, threshold for convergence, default is \code{1e-6}.
//' @param maxIter          numeric, maximum iteration, default is \code{1000}.
//' @return A vector representing solution x.
// [[Rcpp::export]]
arma::fvec pcgsolve_plink_nok(const arma::fvec& W, arma::fvec& tau, arma::fvec & b, std::string preconditioner, float tol, int maxIter) {

    // get number of columns of A=Sigma
    int C = geno.getNnomissing(); // num cols of Sigma
    int R = geno.getNnomissing(); // num rows of Sigma
    
    // get preconditioner M
    arma::fvec Mvec;
    if (preconditioner == "Jacobi"){
        // Mmat = arma::diagmat(A);
        Mvec = getDiagOfSigma_nok(W, tau);
    }
    // to implement 
    // else if(preconditioner == "SSOR"){
    //   arma::fmat D = arma::diagmat(A);
    //   arma::fmat L = arma::trimatl(A);
    //   Mmat = (D+L) * D.i() * (D+L).t();
    // } else if(preconditioner == "ICC"){
    //   Mmat = icc(A);
    // }

  
    // initiate solution x as zeros
    arma::fvec x(C) ;
    x.zeros() ; 
    
    arma::fvec oneVec(C);
    oneVec.ones() ;
    
    // arma::fvec r = b - A * x;
    arma::fvec r = b;
    // arma::fmat Minv = M.i();
    arma::fvec Minv = 1.0/Mvec;

    arma::fvec z = Minv % r;
    arma::fvec p = z;
    float rz_old = arma::sum(r % z);
    float rz_new=1.0;
    // arma::fvec rz_ratio(1);

    
    arma::fvec Ap(R);
    float alpha, beta;
    // vector version of alpha
    // arma::fvec alphaVec(1);
    
    for(int iter = 0; (iter < maxIter) && (rz_new > tol); iter++){
        // Ap = A * p;
        Ap = getCrossprodSig_nok(p, W, tau);
        alpha = rz_old / arma::sum(p % Ap); 
        x += alpha * p;
        r -= alpha * Ap;
        z = Minv % r;
        rz_new = arma::sum( z % r );
        beta = rz_new / rz_old; 

        p = z + beta * p;
        rz_old = rz_new;
        if (iter >= maxIter){
            Rcout << "pcg did not converge." << std::endl;
        }
    }

    
    return arma::round(x*10000.0)/10000.0;
  
} 


//' Preconditioned conjugate gradient method (from cPCG)
//'
//' Sigma = tau[0] * diag(1/W) + tau[1] * sparsekin + tau[2] * grm 
//' @title Solve for x in Sigma x = B (matrix) using preconditioned conjugate gradient method.
//' @param sparsekin        matrix of kinship coefficients (sparse).
//' @param W                numeric vector.
//' @param tau              numeric vector, initial values for variance components.
//' @param bMat                vector, with same dimension as number of rows of A.
//' @param preconditioner   string, method for preconditioning: \code{"Jacobi"} (default), \code{"SSOR"}, or \code{"ICC"}.
//' @param tol              numeric, threshold for convergence, default is \code{1e-6}.
//' @param maxIter          numeric, maximum iteration, default is \code{1000}.
//' @return A matrix representing solution Sigma_iB.
// [[Rcpp::export]]
arma::fmat pcgsolveMat_plink(arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau, arma::fmat & bMat, std::string preconditioner, float tol, int maxIter) {

    // get number of columns of A=Sigma
    int C = geno.getNnomissing(); // num cols of Sigma
    int ncol_b = bMat.n_cols;

    // // get preconditioner M
    // arma::fvec Mvec(C);
    
    // if (preconditioner == "Jacobi"){
    //     // Mmat = arma::diagmat(A);
    //     Mvec = getDiagOfSigma(sparsekin, W, tau);
    // }else{
    // 	Mvec.ones();
    // }
    // // to implement 
    // // else if(preconditioner == "SSOR"){
    // //   arma::fmat D = arma::diagmat(A);
    // //   arma::fmat L = arma::trimatl(A);
    // //   Mmat = (D+L) * D.i() * (D+L).t();
    // // } else if(preconditioner == "ICC"){
    // //   Mmat = icc(A);
    // // }

  
    // initiate solution Sigma_iB as zeros
    arma::fmat Sigma_iB(C,ncol_b);
    Sigma_iB.zeros() ; 
    
    arma::fvec oneVec(C);
    oneVec.ones() ;
    arma::fvec b(C), r(C), Minv(C), z(C), p(C);
    arma::fvec Ap(C);
    float rz_old, rz_new=1.0, alpha, beta;
    arma::fvec x(C); // temp output column of Sigma_iB
    x.zeros() ;
    for(int i = 0; i < ncol_b; i++){
        b = bMat.col(i);
        // // arma::fvec r = b - A * x;
        // r = b;
        // // arma::fmat Minv = M.i();
        // Minv = 1.0/Mvec;
        // z = Minv % r;
        // p = z;
        // rz_old = arma::sum(r % z);
        // rz_new=1.0;
         
        // for(int iter = 0; (iter < maxIter) && (rz_new > tol); iter++){
        //     // Ap = A * p;
        //     Ap = getCrossprodSig(p, sparsekin, W, tau);
        //     alpha = rz_old / arma::sum(p % Ap);
        //     x += alpha * p;
        //     r -= alpha * Ap;
        //     z = Minv % r;
        //     rz_new = arma::sum( z % r );
        //     beta = rz_new / rz_old; 
        //     p = z + beta * p;
        //     rz_old = rz_new;
        //     if (iter >= maxIter){
        //         Rcout << "pcg did not converge for column" << i << std::endl;
        //     }
        // }
        Sigma_iB.col(i) = pcgsolve_plink(sparsekin, W, tau, b, preconditioner, tol, maxIter);
        // Sigma_iB.col(i) = arma::round(x*10000.0)/10000.0;
        // Sigma_iB.col(i) = x;
        // Sigma_iX.col(i) = arma::solve(Sigma, XmatVecTemp, arma::solve_opts::likely_sympd );
    }  
    
    return Sigma_iB;
  
} 

//' Preconditioned conjugate gradient method (from cPCG)
//' without sparse kinship matrix
//'
//' Sigma = tau[0] * diag(1/W) + tau[1] * sparsekin + tau[2] * grm 
//' @title Solve for x in Sigma x = B (matrix) using preconditioned conjugate gradient method.
//' @param W                numeric vector.
//' @param tau              numeric vector, initial values for variance components.
//' @param bMat                vector, with same dimension as number of rows of A.
//' @param preconditioner   string, method for preconditioning: \code{"Jacobi"} (default), \code{"SSOR"}, or \code{"ICC"}.
//' @param tol              numeric, threshold for convergence, default is \code{1e-6}.
//' @param maxIter          numeric, maximum iteration, default is \code{1000}.
//' @return A matrix representing solution Sigma_iB.
// [[Rcpp::export]]
arma::fmat pcgsolveMat_plink_nok(const arma::fvec& W, arma::fvec& tau, arma::fmat & bMat, std::string preconditioner, float tol, int maxIter) {

    // get number of columns of A=Sigma
    int C = geno.getNnomissing(); // num cols of Sigma
    int ncol_b = bMat.n_cols;

    // // get preconditioner M
    // arma::fvec Mvec(C);
    
    // if (preconditioner == "Jacobi"){
    //     // Mmat = arma::diagmat(A);
    //     Mvec = getDiagOfSigma_nok(W, tau);
    // }else{
    // 	Mvec.ones();
    // }
    // to implement 
    // else if(preconditioner == "SSOR"){
    //   arma::fmat D = arma::diagmat(A);
    //   arma::fmat L = arma::trimatl(A);
    //   Mmat = (D+L) * D.i() * (D+L).t();
    // } else if(preconditioner == "ICC"){
    //   Mmat = icc(A);
    // }

  
    // initiate solution Sigma_iB as zeros
    arma::fmat Sigma_iB(C,ncol_b);
    Sigma_iB.zeros() ; 
    
    arma::fvec oneVec(C);
    oneVec.ones() ;
    arma::fvec b(C), r(C), Minv(C), z(C), p(C);
    arma::fvec Ap(C);
    float rz_old, rz_new=1.0, alpha, beta;
    arma::fvec x(C); // temp output column of Sigma_iB
    x.zeros() ;
    for(int i = 0; i < ncol_b; i++){
        b = bMat.col(i);
        // // arma::fvec r = b - A * x;
        // r = b;
        // // arma::fmat Minv = M.i();
        // Minv = 1.0/Mvec;
        // z = Minv % r;
        // p = z;
        // rz_old = arma::sum(r % z);
        // rz_new=1.0;
         
        // for(int iter = 0; (iter < maxIter) && (rz_new > tol); iter++){
        //     // Ap = A * p;
        //     Ap = getCrossprodSig_nok(p, W, tau);
        //     alpha = rz_old / arma::sum(p % Ap);
        //     x += alpha * p;
        //     r -= alpha * Ap;
        //     z = Minv % r;
        //     rz_new = arma::sum( z % r );
        //     beta = rz_new / rz_old; 
        //     p = z + beta * p;
        //     rz_old = rz_new;
        //     if (iter >= maxIter){
        //         Rcout << "pcg did not converge for column" << i << std::endl;
        //     }
        // }
        Sigma_iB.col(i) = pcgsolve_plink_nok(W, tau, b, preconditioner, tol, maxIter);
        // Sigma_iB.col(i) = arma::round(x*10000.0)/10000.0;
        // Sigma_iB.col(i) = x;
        // Sigma_iX.col(i) = arma::solve(Sigma, XmatVecTemp, arma::solve_opts::likely_sympd );
    }  
    
    return Sigma_iB;
  
} 

//' Preconditioned conjugate gradient method LOCO (from cPCG)
//'
//' Sigma = tau[0] * diag(1/W) + tau[1] * sparsekin + tau[2] * grm 
//' @title Solve for x in Sigma x = B (matrix) using preconditioned conjugate gradient method.
//' @param sparsekin        matrix of kinship coefficients (sparse).
//' @param W                numeric vector.
//' @param tau              numeric vector, initial values for variance components.
//' @param bMat                vector, with same dimension as number of rows of A.
//' @param preconditioner   string, method for preconditioning: \code{"Jacobi"} (default), \code{"SSOR"}, or \code{"ICC"}.
//' @param tol              numeric, threshold for convergence, default is \code{1e-6}.
//' @param maxIter          numeric, maximum iteration, default is \code{1000}.
//' @return A matrix representing solution Sigma_iB.
// [[Rcpp::export]]
arma::fmat pcgsolveMat_plink_LOCO(arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec& tau, arma::fmat & bMat, std::string preconditioner, float tol, int maxIter) {

    // get number of columns of A=Sigma
    int C = geno.getNnomissing(); // num cols of Sigma
    int ncol_b = bMat.n_cols;

    // get preconditioner M
    arma::fvec Mvec(C);
    
    if (preconditioner == "Jacobi"){
        // Mmat = arma::diagmat(A);
        Mvec = getDiagOfSigma_LOCO(sparsekin, W, tau);
    }else{
    	Mvec.ones();
    }
    // to implement 
    // else if(preconditioner == "SSOR"){
    //   arma::fmat D = arma::diagmat(A);
    //   arma::fmat L = arma::trimatl(A);
    //   Mmat = (D+L) * D.i() * (D+L).t();
    // } else if(preconditioner == "ICC"){
    //   Mmat = icc(A);
    // }

  
    // initiate solution Sigma_iB as zeros
    arma::fmat Sigma_iB(C,ncol_b);
    Sigma_iB.zeros() ; 
    
    arma::fvec oneVec(C);
    oneVec.ones() ;
    arma::fvec b(C), r(C), Minv(C), z(C), p(C);
    arma::fvec Ap(C);
    float rz_old, rz_new=1.0, alpha, beta;
    arma::fvec x(C); // temp output column of Sigma_iB
    x.zeros() ;
    for(int i = 0; i < ncol_b; i++){
        b = bMat.col(i);
        // arma::fvec r = b - A * x;
        r = b;
        // arma::fmat Minv = M.i();
        Minv = 1.0/Mvec;
        z = Minv % r;
        p = z;
        rz_old = arma::sum(r % z);
        rz_new=1.0;
         
        for(int iter = 0; (iter < maxIter) && (rz_new > tol); iter++){
            // Ap = A * p;
            Ap = getCrossprodSig_LOCO(p, sparsekin, W, tau);
            alpha = rz_old / arma::sum(p % Ap);
            x += alpha * p;
            r -= alpha * Ap;
            z = Minv % r;
            rz_new = arma::sum( z % r );
            beta = rz_new / rz_old; 
            p = z + beta * p;
            rz_old = rz_new;
            if (iter >= maxIter){
                Rcout << "pcg did not converge for column" << i << std::endl;
            }
        }

        Sigma_iB.col(i) = arma::round(x*10000.0)/10000.0;
        // Sigma_iB.col(i) = x;
        // Sigma_iX.col(i) = arma::solve(Sigma, XmatVecTemp, arma::solve_opts::likely_sympd );
    }  
    
    return Sigma_iB;
  
} 



//' Get trace using Hutchinson’s randomized trace estimator
//' @param n_nomissing	integer
//' @param q			integer, total number of sparse and dense correlation matrices
//' @param sparsekin 	matrix of kinship coefficients (sparse).
//' @param W      numeric vector.
//' @param tau      numeric vector, initial values for variance components.
//' @param Sigma_iX 	matrix
//' @param XSiX_inv 	matrix, covariance
//' @param nrun_trace  	Integer. Number of random vectors used for trace estimation.
//' @param maxiter_pcg 	Integer. Max number of iterations for the PCG algorithm.
//' @param tol_pcg     	Numeric. Tolerance for the PCG algorithm.
//' @param cutoff_trace Numeric. Threshold for the coefficient of variation for trace estimation.
// [[Rcpp::export]]
arma::fvec getTrace_plink(int n_nomissing, int q, arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec tau, arma::fmat& Sigma_iX, arma::fmat& XSiX_inv, int nrun_trace, int maxiter_pcg, float tol_pcg, float cutoff_trace){
  	set_seed(200);
  	arma::fmat Sigma_iXt = Sigma_iX.t();
  	int Nnomissing = n_nomissing;
  	arma::fmat temp_mat(nrun_trace, q+1);
  	arma::fvec temp_vec(nrun_trace);
  	arma::vec temp_vec_double(Nnomissing);
  	temp_mat.zeros();

    arma::fvec Sigma_iu;
    arma::fvec Pu;
    arma::fmat Au_mat(Nnomissing, q+1);
    arma::fvec u_vec;
    NumericVector u_vec0;

    int nrun_trace_start = 0;
    int nrun_trace_end = nrun_trace;
    arma::fvec trace_cv(q+1);
    trace_cv.fill(cutoff_trace + 0.1);

    // while((trace_cv > cutoff_trace) | (trace_cv0 > cutoff_trace)){
	while( arma::any(trace_cv > cutoff_trace) ){

    	for(int i = nrun_trace_start; i < nrun_trace_end; i++){

      		u_vec0 = nb(Nnomissing);
      		u_vec = as<arma::fvec>(u_vec0);
      		u_vec = u_vec*2 - 1;
   
      		Sigma_iu = pcgsolve_plink(sparsekin, W, tau, u_vec, "Jacobi", tol_pcg, maxiter_pcg);
          // Sigma_iu = arma::solve(Sigma, u_vec, arma::solve_opts::likely_sympd);

      		Pu = Sigma_iu - Sigma_iX * (XSiX_inv *  (Sigma_iXt * u_vec));
			Au_mat.col(0) = u_vec;
			temp_mat(i,0) = dot(Au_mat.col(0), Pu);

			// conversion for ops with sp_mat
			temp_vec_double = 0.0+sparsekin * arma::conv_to<arma::vec>::from(u_vec);
			Au_mat.col(1) = arma::conv_to<arma::fvec>::from(temp_vec_double);
			// Au_mat.col(1) = sparsekin * u_vec;
			temp_mat(i,1) = dot(Au_mat.col(1), Pu);

      		for(int j=2; j<q+1;j++){
        			Au_mat.col(j) = getCrossprodSig(u_vec, sparsekin, W, tau);
        			temp_mat(i,j) = dot(Au_mat.col(j), Pu);
      		} // end for j in 2:q
    		
  		} // end for i


  		// update trace cv vector
  		for(int i=0; i<q+1; i++){
    			temp_vec = temp_mat.col(i);
    			trace_cv(i) = calCV( temp_vec );
  		} 
  		
  		// trace_cv.elem( arma::find_nonfinite(trace_cv) ).zeros();

  		// if not converge, increase nrun_trace and rerun
  		if( arma::any(trace_cv > cutoff_trace) ){
  		// if((trace_cv > cutoff_trace) | (trace_cv0 > cutoff_trace)){
  			nrun_trace_start = nrun_trace_end;
  			nrun_trace_end = nrun_trace_end + 10;
  			temp_mat.reshape(nrun_trace_end,q+1);
  			//std::cout << "arma::mean(temp_mat0): " << arma::mean(temp_mat0) << std::endl;	
  			Rcout << "CV for trace random estimator using "<< nrun_trace_start << " runs is " << trace_cv <<  "(> " << cutoff_trace << std::endl;
  			Rcout << "try " << nrun_trace_end << "runs" << std::endl;
      } // end if arma::any(trace_cv > cutoff_trace)

    } // end while  arma::any(trace_cv > cutoff_trace)
    Au_mat.clear();
    Pu.clear();
    Sigma_iu.clear();
    u_vec.clear();
    temp_vec.clear();

  	arma::fvec traVec(q+1);
  	for(int i=0; i<q+1; i++){
  		traVec(i) = arma::mean(temp_mat.col(i));
  	}
	temp_mat.clear();
	
  	return(traVec);
}

//' Get trace using Hutchinson’s randomized trace estimator
//' without sparse kinship matrix
//' @param n_nomissing	integer
//' @param q			integer, total number of sparse and dense correlation matrices
//' @param W      numeric vector.
//' @param tau      numeric vector, initial values for variance components.
//' @param Sigma_iX 	matrix
//' @param XSiX_inv 	matrix, covariance
//' @param nrun_trace  	Integer. Number of random vectors used for trace estimation.
//' @param maxiter_pcg 	Integer. Max number of iterations for the PCG algorithm.
//' @param tol_pcg     	Numeric. Tolerance for the PCG algorithm.
//' @param cutoff_trace Numeric. Threshold for the coefficient of variation for trace estimation.
// [[Rcpp::export]]
arma::fvec getTrace_plink_nok(int n_nomissing, int q, const arma::fvec& W, arma::fvec tau, arma::fmat& Sigma_iX, arma::fmat& XSiX_inv, int nrun_trace, int maxiter_pcg, float tol_pcg, float cutoff_trace){
  	set_seed(200);
  	arma::fmat Sigma_iXt = Sigma_iX.t();
  	int Nnomissing = n_nomissing;
  	arma::fmat temp_mat(nrun_trace, q+1);
  	arma::fvec temp_vec(nrun_trace);
  	arma::vec temp_vec_double(Nnomissing);
  	temp_mat.zeros();

    arma::fvec Sigma_iu;
    arma::fvec Pu;
    arma::fmat Au_mat(Nnomissing, q+1);
    arma::fvec u_vec;
    NumericVector u_vec0;

    int nrun_trace_start = 0;
    int nrun_trace_end = nrun_trace;
    arma::fvec trace_cv(q+1);
    trace_cv.fill(cutoff_trace + 0.1);

    // while((trace_cv > cutoff_trace) | (trace_cv0 > cutoff_trace)){
	while( arma::any(trace_cv > cutoff_trace) ){

    	for(int i = nrun_trace_start; i < nrun_trace_end; i++){

      		u_vec0 = nb(Nnomissing);
      		u_vec = as<arma::fvec>(u_vec0);
      		u_vec = u_vec*2 - 1;
   
      		Sigma_iu = pcgsolve_plink_nok(W, tau, u_vec, "Jacobi", tol_pcg, maxiter_pcg);
          // Sigma_iu = arma::solve(Sigma, u_vec, arma::solve_opts::likely_sympd);

      		Pu = Sigma_iu - Sigma_iX * (XSiX_inv *  (Sigma_iXt * u_vec));
			Au_mat.col(0) = u_vec;
			temp_mat(i,0) = dot(Au_mat.col(0), Pu);

      		for(int j=1; j<q+1;j++){ // j=1 only, as q=1 for nok
        			Au_mat.col(j) = getCrossprodSig_nok(u_vec, W, tau);
        			temp_mat(i,j) = dot(Au_mat.col(j), Pu);
      		} // end for j in 1:q
    		
  		} // end for i


  		// update trace cv vector
  		for(int i=0; i<q+1; i++){
    			temp_vec = temp_mat.col(i);
    			trace_cv(i) = calCV( temp_vec );
  		} 
  		
  		// trace_cv.elem( arma::find_nonfinite(trace_cv) ).zeros();

  		// if not converge, increase nrun_trace and rerun
  		if( arma::any(trace_cv > cutoff_trace) ){
  		// if((trace_cv > cutoff_trace) | (trace_cv0 > cutoff_trace)){
  			nrun_trace_start = nrun_trace_end;
  			nrun_trace_end = nrun_trace_end + 10;
  			temp_mat.reshape(nrun_trace_end,q+1);
  			//std::cout << "arma::mean(temp_mat0): " << arma::mean(temp_mat0) << std::endl;	
  			Rcout << "CV for trace random estimator using "<< nrun_trace_start << " runs is " << trace_cv <<  "(> " << cutoff_trace << std::endl;
  			Rcout << "try " << nrun_trace_end << "runs" << std::endl;
      } // end if arma::any(trace_cv > cutoff_trace)

    } // end while  arma::any(trace_cv > cutoff_trace)
    Au_mat.clear();
    Pu.clear();
    Sigma_iu.clear();
    u_vec.clear();
    temp_vec.clear();

  	arma::fvec traVec(q+1);
  	for(int i=0; i<q+1; i++){
  		traVec(i) = arma::mean(temp_mat.col(i));
  	}
	temp_mat.clear();
	
  	return(traVec);
}

//' Get score within a variance component estimation iteration
//'
//' @title Get score value in GLMM
//' @param n_nomissing	integer, number of non-missing individuals as returned from setgeno()
//' @param q			integer, total number of sparse and dense correlation matrices
//' @param Y			column vector, working y.
//' @param X   			data matrix.
//' @param sparsekin 	matrix of kinship coefficients (sparse).
//' @param W 			numeric vector.
//' @param tau 			numeric vector, initial values for variance components.
//' @param fixtau 		vector with elements 0/1, indicator for whether tau elements are fixed.
//' @param is_AI        boolean. Whether AI trace estimation will be used.
//' @param tol_pcg     	Numeric. Tolerance for the PCG algorithm.
//' @param maxiter_pcg 	Integer. Max number of iterations for the PCG algorithm.
//' @param nrun_trace  	Integer. Number of random vectors used for trace estimation.
//' @param cutoff_trace Numeric. Threshold for the coefficient of variation for trace estimation.
//' @return A list with following elements:
//' * YPAPY
//' * Trace
//' * PY
//' * AI
//' * cov = XSiX_inv
//' * tau
//' * alpha
//' * dtau
//' * eta
//' * XSiX_inv_SiXt
// [[Rcpp::export]]
List getScore_plink(int n_nomissing, int q, arma::fvec& Y, arma::fmat& X, arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec tau, arma::fvec fixtau, bool is_AI, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace){

	Rcpp::Timer timer;
	// ==== variable declaration ====
	int n = n_nomissing;
	int q2 = arma::sum(fixtau==0);
	arma::uvec idxtau = arma::find(fixtau==0);
	
	int ncol_x = X.n_cols;
	arma::fvec tau0;

	arma::fmat Sigma_iX(n,ncol_x), Sigma_iXt(ncol_x,n), XSiX_inv(n,ncol_x);
	arma::fvec Sigma_iY(n);
	arma::fvec XmatVecTemp(n);
	arma::fvec PY(n);
	// arma::vec PY_temp(n); // for sp mat format conversion
	arma::sp_mat PYMat(n,1), temp_mat(n,1);
	arma::fvec temp_vec(n);

	arma::fmat APYmat(n,q+1); // store col vecs of APY for each var component
	arma::fvec PAPY_1(n), PAPY(n);
	// arma::fvec YPAPY(q+1); // for storing results from q+1 components
	arma::fmat AI(q+1,q+1);
	AI.zeros();
	arma::fvec YPAPY(q+1);
	YPAPY.zeros();
	arma::fvec Trace(q+1);
	Trace.zeros();

	tau0 = tau;
    float diagMean = getMeanDiagOfGRM();
    tau0(2) = tau0(2) / diagMean; // adjust by mean(diag(grm))

    timer.step("start");
	// ==== get computing components, alpha, eta ====
	Sigma_iX = pcgsolveMat_plink(sparsekin, W, tau0, X, "Jacobi", tol_pcg, maxiter_pcg);
    Sigma_iXt = Sigma_iX.t();
    Sigma_iY = pcgsolve_plink(sparsekin, W, tau0, Y, "Jacobi", tol_pcg, maxiter_pcg);
    
	arma::fmat XSiX = arma::symmatu( Sigma_iXt * X) ;
    try {
    	XSiX_inv = arma::inv_sympd(  XSiX  ); //=cov
    }
    catch(const std::exception& e){
    	std::cout << "getScore_plink: XSiX not symmetric positive definite, using Cholesky decomp. inverse" << std::endl;
    	XSiX_inv = arma::pinv(   XSiX  ); 
    }
    
    arma::fmat XSiX_inv_SiXt = XSiX_inv * Sigma_iXt;
    arma::fvec alpha = XSiX_inv_SiXt * Y;
    PY = Sigma_iY - Sigma_iX * alpha;
    // convert PY to sparsemat
    for(int i=0; i<n; i++){
		PYMat(i,0) = PY(i);
	}

    // arma::fvec diagp = getDiagOfSigma(sparsekin,W,tau0);
	arma::fvec diagp = tau0(0)/W;
	arma::fvec eta = Y - diagp % PY; 	

	timer.step("computing_components");

    if(q2>0){
		// ==== fill APYmat ====
		// ---- i=0 ----
 		APYmat.col(0) = PY/W;
		// ---- i=1 (sparsekin) ----
		// conversion for ops with sp_mat
		// PY_temp = 0.0+sparsekin * arma::conv_to<arma::vec>::from(PY);
		temp_mat = sparsekin * PYMat;
		for(int j=0; j<n; j++){
			temp_vec(j) = temp_mat(j,0);
		}
		// APYmat.col(1)= arma::conv_to<arma::fvec>::from(PY_temp);
		APYmat.col(1)= temp_vec ; 
		// ---- i=2 (GRM) ----
		APYmat.col(2)= getCrossprodGRM(PY) ;
	    
	    // ==== get YPAPY ====
	    for(int i=0; i<q+1; i++){
	    	XmatVecTemp = APYmat.col(i);
	        PAPY_1 = pcgsolve_plink(sparsekin, W, tau0, XmatVecTemp, "Jacobi", tol_pcg, maxiter_pcg);
			PAPY = PAPY_1 - Sigma_iX * (XSiX_inv * (Sigma_iXt * XmatVecTemp));
					
			// fill AI
			for(int j=0; j<=i; j++){
				AI(i,j) = arma::dot(APYmat.col(j), PAPY);
				// AI(i,j) = arma::sum(APYmat.col(j) % PAPY);
				
				if(j != i){
					AI(j,i) = AI(i,j);
				}	
			}

			YPAPY(i) = arma::dot(PY, APYmat.col(i));

	    } // end for i
	    timer.step("get_YPAPY");
		
		//====Calculate trace: 1/2====
		// use Hutchinson’s randomized trace estimator
		// q fixed as 2
		Trace = getTrace_plink(n_nomissing, 2, sparsekin, W, tau0, Sigma_iX, XSiX_inv, nrun_trace, maxiter_pcg, tol_pcg, cutoff_trace);
		timer.step("trace");
		
		//==== update tau: 1/2 ====
		// tau0 = tau0 +  tau0 % tau0/float(n) % (YPAPY-Trace);
		for(int i=0; i<q2; i++){
			tau0( idxtau(i) ) = tau0( idxtau(i) ) +  tau0( idxtau(i) ) * tau0( idxtau(i) )/float(n) * (YPAPY( idxtau(i) )-Trace( idxtau(i) ));
		}// end for i in q2
		

		// // ==== get AI: 2/2 with updated tau0 ====
		// Sigma_iX = pcgsolveMat_plink(sparsekin, W, tau0, X, "Jacobi", tol_pcg, maxiter_pcg);
		// Sigma_iXt = Sigma_iX.t();
	 //    Sigma_iY = pcgsolve_plink(sparsekin, W, tau0, Y, "Jacobi", tol_pcg, maxiter_pcg);  

		// XSiX = arma::symmatu( Sigma_iXt * X) ;
	 //    try {
	 //    	XSiX_inv = arma::inv_sympd(  XSiX  ); //=cov
	 //    }
	 //    catch(const std::exception& e){
	 //    	std::cout << "getScore_plink: XSiX not symmetric positive definite, using Cholesky decomp. inverse" << std::endl;
	 //    	XSiX_inv = arma::pinv(   XSiX  ); 
	 //    }   
	 //    XSiX_inv_SiXt = XSiX_inv * Sigma_iXt;
	 //    alpha = XSiX_inv_SiXt * Y;
	 //    PY = Sigma_iY - Sigma_iX * alpha;
	 //    // convert PY to sparsemat
	 //    for(int i=0; i<n; i++){
		// 	PYMat(i,0) = PY(i);
		// }

	 //    for(int i=0; i<q+1; i++){
	 //    	// fill APYmat
	 //    	if(i==0){ // identity mat
	 //    		APYmat.col(0) = PY/W;
	 //    	}
	 //    	else {
	 //        	if(i==1){ // sparsekin
	 //        		// conversion for ops with sp_mat
	 //        		// PY_temp = 0.0+sparsekin * arma::conv_to<arma::vec>::from(PY);
	 //        		// APYmat.col(1)= arma::conv_to<arma::fvec>::from(PY_temp);
	 //        		temp_mat = sparsekin * PYMat;
	 //        		for(int j=0; j<n; j++){
	 //        			temp_vec(j) = temp_mat(j,0);
	 //        		}
	 //        		// APYmat.col(1)= arma::conv_to<arma::fvec>::from(PY_temp);
	 //        		APYmat.col(i)= temp_vec ; 
	 //        	}
	 //        	else { //GRM
	 //        		APYmat.col(i)= getCrossprodGRM(PY) ;
	 //        	} // end if i==1 	
	 //        } // end if i==0

	 //        XmatVecTemp = APYmat.col(i);
	 //        PAPY_1 = pcgsolve_plink(sparsekin, W, tau0, XmatVecTemp, "Jacobi", tol_pcg, maxiter_pcg);
		// 	PAPY = PAPY_1 - Sigma_iX * (XSiX_inv * (Sigma_iXt * XmatVecTemp));
					
		// 	// fill AI
		// 	for(int j=0; j<=i; j++){
		// 		AI(i,j) = arma::dot(APYmat.col(j), PAPY);			
		// 		if(j != i){
		// 			AI(j,i) = AI(i,j);
		// 		}	
		// 	}
		// 	YPAPY(i) = arma::dot(PY, APYmat.col(i));
	 //    } // end for i
		// timer.step("AI");

		// //====Calculate trace: 2/2====
		// // use Hutchinson’s randomized trace estimator
		// // q fixed as 2
		// Trace = getTrace_plink(n_nomissing, 2, sparsekin, W, tau0, Sigma_iX, XSiX_inv, nrun_trace, maxiter_pcg, tol_pcg, cutoff_trace);
		// timer.step("trace");

		// ====== get dtau =======
		arma::fvec score_vec_update(q2);
		arma::fvec dtau(q+1), dtau_pre(q2);
		arma::fmat AI_update(q2,q2);
		
		for(int i=0; i<q2; i++){
			// fill AI_update
			for(int j=0; j<=i;j++){
				AI_update(i,j) = AI(idxtau(i),idxtau(j));
				if(j != i){
					AI_update(j,i) = AI_update(i,j);
				}
			} // end for j	
			// fill score_vec_update
			score_vec_update(i) = YPAPY(idxtau(i)) - Trace(idxtau(i));	
		} // end for i in q2
		try{
			// arma::fvec dtau = arma::solve(AI_mat, score_vec, arma::solve_opts::likely_sympd);
			dtau_pre = arma::solve(AI_update, score_vec_update, arma::solve_opts::allow_ugly);	
		}
		catch(std::runtime_error){
			std::cout << "getScore_plink: arma::solve(AI_mat, score_vec): AI seems singular, using less variant components matrix is suggested." << std::endl;
			dtau_pre.zeros();
		}

		// fill dtau using dtau_pre, padding 0
		int i2 = 0;
		for(int i=0; i<q+1; i++){
			if(fixtau(i)<1){ // not fixed
				dtau(i) = dtau_pre(i2);
				i2++;
			} else { // fixed
				dtau(i) = 0.0;
			}
		} // end for i
		
		timer.step("get_dtau");

		NumericVector re_time(timer), re_time2;
		re_time2 = diff(re_time);
		Rcpp::Rcout << "Runtime : " << "\n";
		// Rcpp::Rcout << "computing_components\tget_YPAPY\ttrace\tAI\ttrace\tget_dtau" << "\n";
		Rcpp::Rcout << "computing_components\tget_YPAPY\ttrace\tget_dtau" << "\n";
		Rcpp::Rcout << re_time2 << "\n";
		Rcpp::Rcout << "----------------\n\n";

		//====Update tau: 2/2 in glmmaiUpdate ====
		// List out = List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY, Named("AI") = AI, Named("cov") = XSiX_inv, Named("alpha") =alpha, Named("eta") = eta, Named("tau") = tau0, Named("PAPY")=PAPY_out,Named("APYmat")=APYmat);
		List out = List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY, Named("AI") = AI, Named("cov") = XSiX_inv, Named("tau") = tau0, Named("alpha") =alpha, Named("dtau")=dtau, Named("eta")=eta, Named("XSiX_inv_SiXt")=XSiX_inv_SiXt);
		return(out);
	} // end if q2>0

	List out = List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY, Named("AI") = NULL, Named("cov") = XSiX_inv, Named("tau") = tau0, Named("alpha") =alpha, Named("dtau")=NULL, Named("eta")=eta, Named("XSiX_inv_SiXt")=XSiX_inv_SiXt);

	return(out);
}

//' Get score within a variance component estimation iteration
//' without sparse kinship matrix
//' 
//' @title Get score value in GLMM
//' @param n_nomissing	integer, number of non-missing individuals as returned from setgeno()
//' @param q			integer, total number of sparse and dense correlation matrices
//' @param Y			column vector, working y.
//' @param X   			data matrix.
//' @param W 			numeric vector.
//' @param tau 			numeric vector, initial values for variance components.
//' @param fixtau 		vector with elements 0/1, indicator for whether tau elements are fixed.
//' @param is_AI        boolean. Whether AI trace estimation will be used.
//' @param tol_pcg     	Numeric. Tolerance for the PCG algorithm.
//' @param maxiter_pcg 	Integer. Max number of iterations for the PCG algorithm.
//' @param nrun_trace  	Integer. Number of random vectors used for trace estimation.
//' @param cutoff_trace Numeric. Threshold for the coefficient of variation for trace estimation.
//' @return A list with following elements:
//' * YPAPY
//' * Trace
//' * PY
//' * AI
//' * cov = XSiX_inv
//' * tau
//' * alpha
//' * dtau
//' * eta
//' * XSiX_inv_SiXt
// [[Rcpp::export]]
List getScore_plink_nok(int n_nomissing, int q, arma::fvec& Y, arma::fmat& X, const arma::fvec& W, arma::fvec tau, arma::fvec fixtau, bool is_AI, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace){

	// ==== variable declaration ====
	int n = n_nomissing;
	int q2 = arma::sum(fixtau==0);
	arma::uvec idxtau = arma::find(fixtau==0);
	
	int ncol_x = X.n_cols;
	arma::fvec tau0;

	arma::fmat Sigma_iX(n,ncol_x), Sigma_iXt(ncol_x,n), XSiX_inv(n,ncol_x);
	arma::fvec Sigma_iY(n);
	arma::fvec XmatVecTemp(n);
	arma::fvec PY(n);
	arma::vec PY_temp(n); // for sp mat format conversion

	arma::fmat APYmat(n,q+1); // store col vecs of APY for each var component
	arma::fvec PAPY_1(n), PAPY(n);
	// arma::fvec YPAPY(q+1); // for storing results from q+1 components
	arma::fmat AI(q+1,q+1);
	AI.zeros();
	arma::fvec YPAPY(q+1);
	YPAPY.zeros();
	arma::fvec Trace(q+1);
	Trace.zeros();

	tau0 = tau;
    float diagMean = getMeanDiagOfGRM();
    tau0(1) = tau0(1) / diagMean; // adjust by mean(diag(grm))

	// ==== get computing components, alpha, eta ====
	Sigma_iX = pcgsolveMat_plink_nok(W, tau0, X, "Jacobi", tol_pcg, maxiter_pcg);
    Sigma_iXt = Sigma_iX.t();
    Sigma_iY = pcgsolve_plink_nok(W, tau0, Y, "Jacobi", tol_pcg, maxiter_pcg);
    
	arma::fmat XSiX = arma::symmatu( Sigma_iXt * X) ;
    
    try {
    	XSiX_inv = arma::inv_sympd(  XSiX  ); //=cov
    }
    catch(const std::exception& e){
    	std::cout << "getScore_plink_nok: XSiX not symmetric positive definite, using Cholesky decomp. inverse"<< std::endl;
    	XSiX_inv = arma::pinv(   XSiX  ); //=cov
    }

    arma::fmat XSiX_inv_SiXt = XSiX_inv * Sigma_iXt;
    arma::fvec alpha = XSiX_inv_SiXt * Y;
    PY = Sigma_iY - Sigma_iX * alpha;

    // arma::fvec diagp = getDiagOfSigma(sparsekin,W,tau0);
	arma::fvec diagp = tau0(0)/W;
	arma::fvec eta = Y - diagp % PY; 	

    if(q2>0){
		// ==== get AI: 1/2 ====
	    for(int i=0; i<q+1; i++){
	    	// fill APYmat
	    	if(i==0){ // identity mat
	    		APYmat.col(0) = PY/W;
	    	}
	    	else {
	    		//GRM
	        	APYmat.col(i)= getCrossprodGRM(PY) ;	
	        } // end if i==0

	        XmatVecTemp = APYmat.col(i);
	        PAPY_1 = pcgsolve_plink_nok(W, tau0, XmatVecTemp, "Jacobi", tol_pcg, maxiter_pcg);
			// PAPY_1 = arma::solve(Sigma, XmatVecTemp, arma::solve_opts::likely_sympd	);
			PAPY = PAPY_1 - Sigma_iX * (XSiX_inv * (Sigma_iXt * XmatVecTemp));
			// PAPY = PAPY_1 - Sigma_iX * (XSiX_inv * (Sigma_iXt * PAPY_1));
					
			// fill AI
			for(int j=0; j<=i; j++){
				AI(i,j) = arma::dot(APYmat.col(j), PAPY);
				// AI(i,j) = arma::sum(APYmat.col(j) % PAPY);
				
				if(j != i){
					AI(j,i) = AI(i,j);
				}	
			}

			YPAPY(i) = arma::dot(PY, APYmat.col(i));

	    } // end for i
		
		//====Calculate trace: 1/2====
		// use Hutchinson’s randomized trace estimator
		// q fixed as 1
		Trace = getTrace_plink_nok(n_nomissing, 1, W, tau0, Sigma_iX, XSiX_inv, nrun_trace, maxiter_pcg, tol_pcg, cutoff_trace);

		//==== update tau: 1/2 ====
		// tau0 = tau0 +  tau0 % tau0/float(n) % (YPAPY-Trace);
		for(int i=0; i<q2; i++){
			tau0( idxtau(i) ) = tau0( idxtau(i) ) +  tau0( idxtau(i) ) * tau0( idxtau(i) )/float(n) * (YPAPY( idxtau(i) )-Trace( idxtau(i) ));
		}// end for i in q2

		// ==== get AI: 2/2 ====
		Sigma_iX = pcgsolveMat_plink_nok(W, tau0, X, "Jacobi", tol_pcg, maxiter_pcg);
	    Sigma_iXt = Sigma_iX.t();
	    Sigma_iY = pcgsolve_plink_nok(W, tau0, Y, "Jacobi", tol_pcg, maxiter_pcg);    
		XSiX = arma::symmatu( Sigma_iXt * X) ;
	    try {
	    	XSiX_inv = arma::inv_sympd(  XSiX  ); //=cov
	    }
	    catch(const std::exception& e){
	    	std::cout << "getScore_plink_nok: XSiX not symmetric positive definite, using Cholesky decomp. inverse"<< std::endl;
	    	XSiX_inv = arma::pinv(   XSiX  ); //=cov
	    }
	    XSiX_inv_SiXt = XSiX_inv * Sigma_iXt;
	    alpha = XSiX_inv_SiXt * Y;
	    PY = Sigma_iY - Sigma_iX * alpha;
	    for(int i=0; i<q+1; i++){
	    	// fill APYmat
	    	if(i==0){ // identity mat
	    		APYmat.col(0) = PY/W;
	    	}
	    	else {
	    		//GRM
	        	APYmat.col(i)= getCrossprodGRM(PY) ;	
	        } // end if i==0

	        XmatVecTemp = APYmat.col(i);
	        PAPY_1 = pcgsolve_plink_nok(W, tau0, XmatVecTemp, "Jacobi", tol_pcg, maxiter_pcg);
			PAPY = PAPY_1 - Sigma_iX * (XSiX_inv * (Sigma_iXt * XmatVecTemp));
					
			// fill AI
			for(int j=0; j<=i; j++){
				AI(i,j) = arma::dot(APYmat.col(j), PAPY);			
				if(j != i){
					AI(j,i) = AI(i,j);
				}	
			}
			YPAPY(i) = arma::dot(PY, APYmat.col(i));
	    } // end for i
		
		//====Calculate trace: 2/2====
		// use Hutchinson’s randomized trace estimator
		// q fixed as 1
		Trace = getTrace_plink_nok(n_nomissing, 1, W, tau0, Sigma_iX, XSiX_inv, nrun_trace, maxiter_pcg, tol_pcg, cutoff_trace);

		// ====== get dtau =======
		arma::fvec score_vec_update(q2);
		arma::fvec dtau(q+1), dtau_pre(q2);
		arma::fmat AI_update(q2,q2);

		for(int i=0; i<q2; i++){
			// fill AI_update
			for(int j=0; j<=i;j++){
				AI_update(i,j) = AI(idxtau(i),idxtau(j));
				if(j != i){
					AI_update(j,i) = AI_update(i,j);
				}
			} // end for j	
			// fill score_vec_update
			score_vec_update(i) = YPAPY(idxtau(i)) - Trace(idxtau(i));	
		} // end for i in q2
		try{
			// arma::fvec dtau = arma::solve(AI_mat, score_vec, arma::solve_opts::likely_sympd);
			dtau_pre = arma::solve(AI_update, score_vec_update, arma::solve_opts::allow_ugly);	
		}
		catch(std::runtime_error){
			std::cout << "getScore_plink_nok: arma::solve(AI_mat, score_vec): AI seems singular, using less variant components matrix is suggested." << std::endl;
			dtau_pre.zeros();
		}

		// fill dtau using dtau_pre, padding 0
		int i2 = 0;
		for(int i=0; i<q+1; i++){
			if(fixtau(i)<1){ // not fixed
				dtau(i) = dtau_pre(i2);
				i2++;
			} else { // fixed
				dtau(i) = 0.0;
			}
		} // end for i

		//====Update tau: 2/2 in glmmaiUpdate_nok ====
		// List out = List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY, Named("AI") = AI, Named("cov") = XSiX_inv, Named("alpha") =alpha, Named("eta") = eta, Named("tau") = tau0, Named("PAPY")=PAPY_out,Named("APYmat")=APYmat);
		List out = List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY, Named("AI") = AI, Named("cov") = XSiX_inv, Named("tau") = tau0, Named("alpha") =alpha, Named("dtau")=dtau, Named("eta")=eta, Named("XSiX_inv_SiXt")=XSiX_inv_SiXt);
		return(out);
	} // end if q2>0

	List out = List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY, Named("AI") = NULL, Named("cov") = XSiX_inv, Named("tau") = tau0, Named("alpha") =alpha, Named("dtau")=NULL, Named("eta")=eta, Named("XSiX_inv_SiXt")=XSiX_inv_SiXt);
	return(out);
}

//' Get score within a variance component estimation iteration
//'
//' @title Get score value in GLMM
//' @param n_nomissing	integer, number of non-missing individuals as returned from setgeno()
//' @param q			integer, total number of sparse and dense correlation matrices
//' @param Y			column vector, working y.
//' @param X   			data matrix.
//' @param sparsekin 	matrix of kinship coefficients (sparse).
//' @param W 			numeric vector.
//' @param tau 			numeric vector, initial values for variance components.
//' @param fixtau 		vector with elements 0/1, indicator for whether tau elements are fixed.
//' @param is_AI        boolean. Whether AI trace estimation will be used.
//' @param tol_pcg     	Numeric. Tolerance for the PCG algorithm.
//' @param maxiter_pcg 	Integer. Max number of iterations for the PCG algorithm.
//' @param nrun_trace  	Integer. Number of random vectors used for trace estimation.
//' @param cutoff_trace Numeric. Threshold for the coefficient of variation for trace estimation.
//' @return A list with following elements:
//' * YPAPY
//' * Trace
//' * PY
//' * AI
//' * cov = XSiX_inv
//' * tau
//' * alpha
//' * dtau
//' * eta
//' * XSiX_inv_SiXt
// [[Rcpp::export]]
List getScore_plink_LOCO(int n_nomissing, int q, arma::fvec& Y, arma::fmat& X, arma::sp_mat& sparsekin, const arma::fvec& W, arma::fvec tau, arma::fvec fixtau, bool is_AI, float tol_pcg, int maxiter_pcg, int nrun_trace, float cutoff_trace){

	// ==== variable declaration ====
	int n = n_nomissing;
	int q2 = arma::sum(fixtau==0);
	
	int ncol_x = X.n_cols;
	arma::fvec tau0;

	arma::fmat Sigma_iX(n,ncol_x), Sigma_iXt(ncol_x,n), XSiX_inv(n,ncol_x);
	arma::fvec Sigma_iY(n);
	arma::fvec XmatVecTemp(n);
	arma::fvec PY(n);
	arma::vec PY_temp(n); // for sp mat format conversion

	arma::fmat APYmat(n,q+1); // store col vecs of APY for each var component
	arma::fvec PAPY_1(n), PAPY(n);
	// arma::fvec YPAPY(q+1); // for storing results from q+1 components
	arma::fmat AI(q+1,q+1);
	AI.zeros();
	arma::fvec YPAPY(q+1);
	YPAPY.zeros();
	arma::fvec Trace(q+1);
	Trace.zeros();

	tau0 = tau;
    float diagMean = getMeanDiagOfGRM();
    tau0(2) = tau0(2) / diagMean; // adjust by mean(diag(grm))
	// ==== get computing components, alpha, eta ====
	Sigma_iX = pcgsolveMat_plink_LOCO(sparsekin, W, tau0, X, "Jacobi", tol_pcg, maxiter_pcg);
    // Sigma_iX = pcgsolveMat_plink(sparsekin, W, tau0, X, "", tol_pcg, maxiter_pcg);
    Sigma_iXt = Sigma_iX.t();
    // Sigma_iY = pcgsolve_plink(sparsekin, W, tau, Y, "Jacobi", tol_pcg, maxiter_pcg);
    Sigma_iY = pcgsolve_plink_LOCO(sparsekin, W, tau0, Y, "Jacobi", tol_pcg, maxiter_pcg);

	arma::fmat XSiX = arma::symmatu( Sigma_iXt * X) ;
    try {
    	XSiX_inv = arma::inv_sympd(  XSiX  ); //=cov
    }
    catch(const std::exception& e){
    	std::cout << "getScore_plink: XSiX not symmetric positive definite, using Cholesky decomp. inverse" << std::endl;
    	XSiX_inv = arma::pinv(   XSiX  ); 
    }

    arma::fmat XSiX_inv_SiXt = XSiX_inv * Sigma_iXt;
    arma::fvec alpha = XSiX_inv_SiXt * Y;
    PY = Sigma_iY - Sigma_iX * alpha;

    // arma::fvec diagp = getDiagOfSigma_LOCO(sparsekin,W,tau0);
	arma::fvec diagp = tau0(0)/W;
	arma::fvec eta = Y - diagp % PY; 

	if(q2>0){

		// ==== get AI: 1/2 ====
	    for(int i=0; i<q+1; i++){
	    	// fill APYmat
	    	if(i==0){ // identity mat
	    		APYmat.col(0) = PY/W;
	    	}
	    	else {
	        	if(i==1){ // sparsekin
	        		// conversion for ops with sp_mat
	        		PY_temp = 0.0+sparsekin * arma::conv_to<arma::vec>::from(PY);
	        		APYmat.col(1)= arma::conv_to<arma::fvec>::from(PY_temp);
	        	}
	        	else { //densekin
	        		APYmat.col(i)= getCrossprodGRM_LOCO(PY) ;
	        	} // end if i==1 	
	        } // end if i==0

	        XmatVecTemp = APYmat.col(i);
	        PAPY_1 = pcgsolve_plink_LOCO(sparsekin, W, tau0, XmatVecTemp, "Jacobi", tol_pcg, maxiter_pcg);
			// PAPY_1 = arma::solve(Sigma, XmatVecTemp, arma::solve_opts::likely_sympd	);
			PAPY = PAPY_1 - Sigma_iX * (XSiX_inv * (Sigma_iXt * XmatVecTemp));
			// PAPY = PAPY_1 - Sigma_iX * (XSiX_inv * (Sigma_iXt * PAPY_1));
					
			// fill AI
			for(int j=0; j<=i; j++){
				AI(i,j) = arma::dot(APYmat.col(j), PAPY);
				// AI(i,j) = arma::sum(APYmat.col(j) % PAPY);
				
				if(j != i){
					AI(j,i) = AI(i,j);
				}	
			}

			YPAPY(i) = arma::dot(PY, APYmat.col(i));

	    } // end for i
		
		//====Calculate trace: 1/2====
		// use Hutchinson’s randomized trace estimator
		// q fixed as 2
		Trace = getTrace_plink(n_nomissing, 2, sparsekin, W, tau0, Sigma_iX, XSiX_inv, nrun_trace, maxiter_pcg, tol_pcg, cutoff_trace);

		//==== update tau: 1/2 ====
		tau0 = tau0 +  tau0 % tau0/float(n) % (YPAPY-Trace);

		// ==== get AI: 2/2 ====
		Sigma_iX = pcgsolveMat_plink_LOCO(sparsekin, W, tau0, X, "Jacobi", tol_pcg, maxiter_pcg);
	    Sigma_iXt = Sigma_iX.t();
	    Sigma_iY = pcgsolve_plink_LOCO(sparsekin, W, tau0, Y, "Jacobi", tol_pcg, maxiter_pcg);
	    XSiX = arma::symmatu( Sigma_iXt * X) ;
	    try {
	    	XSiX_inv = arma::inv_sympd(  XSiX  ); //=cov
	    }
	    catch(const std::exception& e){
	    	std::cout << "getScore_plink: XSiX not symmetric positive definite, using Cholesky decomp. inverse" << std::endl;
	    	XSiX_inv = arma::pinv(   XSiX  ); 
	    }
	    XSiX_inv_SiXt = XSiX_inv * Sigma_iXt;
	    alpha = XSiX_inv_SiXt * Y;
	    PY = Sigma_iY - Sigma_iX * alpha;
	    for(int i=0; i<q+1; i++){
	    	// fill APYmat
	    	if(i==0){ // identity mat
	    		APYmat.col(0) = PY/W;
	    	}
	    	else {
	        	if(i==1){ // sparsekin
	        		// conversion for ops with sp_mat
	        		PY_temp = 0.0+sparsekin * arma::conv_to<arma::vec>::from(PY);
	        		APYmat.col(1)= arma::conv_to<arma::fvec>::from(PY_temp);
	        	}
	        	else { //densekin
	        		APYmat.col(i)= getCrossprodGRM_LOCO(PY) ;
	        	} // end if i==1 	
	        } // end if i==0

	        XmatVecTemp = APYmat.col(i);
	        PAPY_1 = pcgsolve_plink_LOCO(sparsekin, W, tau0, XmatVecTemp, "Jacobi", tol_pcg, maxiter_pcg);
			PAPY = PAPY_1 - Sigma_iX * (XSiX_inv * (Sigma_iXt * XmatVecTemp));
					
			// fill AI
			for(int j=0; j<=i; j++){
				AI(i,j) = arma::dot(APYmat.col(j), PAPY);			
				if(j != i){
					AI(j,i) = AI(i,j);
				}	
			}

			YPAPY(i) = arma::dot(PY, APYmat.col(i));

	    } // end for i
		
		//====Calculate trace: 2/2====
		// use Hutchinson’s randomized trace estimator
		// q fixed as 2
		Trace = getTrace_plink(n_nomissing, 2, sparsekin, W, tau0, Sigma_iX, XSiX_inv, nrun_trace, maxiter_pcg, tol_pcg, cutoff_trace);

		// ====== get dtau =======
		arma::fvec score_vec = YPAPY - Trace;
		arma::fvec dtau;
		try{
			// arma::fvec dtau = arma::solve(AI_mat, score_vec, arma::solve_opts::likely_sympd);
			dtau = arma::solve(AI, score_vec, arma::solve_opts::allow_ugly);	
		}
		catch(std::runtime_error){
			std::cout << "getScore_plink_LOCO: arma::solve(AI_mat, score_vec): AI seems singular, using less variant components matrix is suggested." << std::endl;
			dtau.zeros();
		}

		//====Update tau: 2/2 in glmmaiUpdate ====
		// List out = List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY, Named("AI") = AI, Named("cov") = XSiX_inv, Named("alpha") =alpha, Named("eta") = eta, Named("tau") = tau0, Named("PAPY")=PAPY_out,Named("APYmat")=APYmat);
		List out = List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY, Named("AI") = AI, Named("cov") = XSiX_inv, Named("tau") = tau0, Named("alpha") =alpha, Named("dtau")=dtau, Named("eta")=eta, Named("XSiX_inv_SiXt")=XSiX_inv_SiXt);
		return(out);
	} // end if q2>0

	// List out = List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY, Named("AI") = AI, Named("cov") = XSiX_inv, Named("alpha") =alpha, Named("eta") = eta, Named("tau") = tau0, Named("PAPY")=PAPY_out,Named("APYmat")=APYmat);
	List out = List::create(Named("YPAPY") = YPAPY, Named("Trace") = Trace, Named("PY") = PY, Named("AI") = NULL, Named("cov") = XSiX_inv, Named("tau") = tau0, Named("alpha") =alpha, Named("dtau")=NULL, Named("eta")=eta, Named("XSiX_inv_SiXt")=XSiX_inv_SiXt);

	return(out);
}


//' Get standardized genotype for a given index
//'
//' @param i	integer, index of the genotype to be extracted
//' @return A vector of standardized genotypes for Nnomissing individuals
// [[Rcpp::export]]
arma::fvec getOneSNPR(int i){
	arma::fvec vec;
	geno.Get_OneSNP_StdGeno(i, &vec);
	return(vec);
}

//' Get unstandardized genotype for a given index
//'
//' @param i	integer, index of the genotype to be extracted
//' @return A vector of standardized genotypes for Nnomissing individuals
// [[Rcpp::export]]
arma::fvec getOneSNPR_nostd(int i){
	arma::fvec * vecp;
	vecp=geno.Get_OneSNP_Geno(i);
	return(*vecp);
}