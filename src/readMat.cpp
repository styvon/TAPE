#include <RcppArmadillo.h>
// [[Rcpp:depends(RcppArmadillo)]]
#include <math.h>
#include <iostream>
#include <fstream>
#include <algorithm>    // std::max
#include <string>    // std::stod


using namespace Rcpp; 
// [[Rcpp::export]]
arma::fmat readMat(const std::string &filename, const std::string &delimeter = "\t")
{
    std::ifstream infile(filename);
    std::vector<std::vector<float>> datas;

    for(std::string line; std::getline(infile, line); ) {

        std::vector<float> data;

        // split string by delimeter
        auto start = 0U;
        auto end = line.find(delimeter);
        while (end != std::string::npos) {
            data.push_back(std::stod(line.substr(start, end - start)));
            start = end + delimeter.length();
            end = line.find(delimeter, start);
        }
        data.push_back(std::stod(line.substr(start, end)));
        datas.push_back(data);
    }

    arma::fmat data_mat = arma::zeros<arma::fmat>(datas.size(), datas[0].size());

    for (int i=0; i<datas.size(); i++) {
        arma::fmat r(datas[i]);
        data_mat.row(i) = r.t();
    }

    return data_mat;
}