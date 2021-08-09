#ifndef _UTILS_EXPORTED_
#define _UTILS_EXPORTED_

#include <Rcpp.h>
#include <vector>
#include <string>

Rcpp::String all_checks_collapse(const Rcpp::StringVector &checks);

double position_icscoreC(std::vector<double> pos,
    std::vector<double> bkg, const std::string &type = "PPM",
    double pseudocount = 1.0, double nsites = 100.0, bool relative_entropy = false);

std::string get_consensusC(std::vector<double> pos,
    const std::string &alphahet = "DNA", const std::string &type = "PPM",
    double pseudocount = 1.0);

std::string get_consensusAAC(std::vector<double> pos,
    const std::string &type = "PPM", double pseudocount = 0.0);

Rcpp::IntegerVector order_char_cpp(const Rcpp::CharacterVector x);

#endif
