#ifndef _UTILS_EXPORTED_
#define _UTILS_EXPORTED_

#include <Rcpp.h>
#include <vector>
#include <string>

Rcpp::String all_checks_collapse(const Rcpp::StringVector &checks);

double position_icscoreC(std::vector<double> pos,
    std::vector<double> bkg, std::string type = "PPM", double pseudocount = 1.0,
    double nsites = 100.0, bool relative_entropy = false);

std::string get_consensusC(std::vector<double> pos, std::string alph = "DNA",
    std::string type = "PPM", double pseudocount = 1.0);

std::string get_consensusAAC(std::vector<double> pos, std::string type = "PPM",
    double pseudocount = 0.0);

#endif
