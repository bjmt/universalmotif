#ifndef _UTILS_
#define _UTILS_

#include <Rcpp.h>
using namespace Rcpp;

NumericVector pcm_to_ppmC(NumericVector position, double pseudocount=0);

NumericVector ppm_to_pcmC(NumericVector position, double nsites);

NumericVector ppm_to_pwmC(NumericVector position, NumericVector bkg=0,
    double pseudocount=0, NumericVector nsites=NumericVector::create());

NumericVector pwm_to_ppmC(NumericVector position, NumericVector bkg=0);

NumericVector ppm_to_icmC(NumericVector position, NumericVector bkg=0,
    bool relative_entropy=false);

double position_icscoreC(NumericVector position, NumericVector bkg=0,
    String type="PPM", double pseudocount=0.8, double nsites=100,
    bool relative_entropy=false);

NumericVector icm_to_ppmC(NumericVector position);

String get_consensusC(NumericVector position, String alphabet="DNA",
    String type="PPM", double pseudocount=0.8);

NumericVector consensus_to_ppmC(String letter);

NumericVector consensus_to_ppmAAC(String letter);

String get_consensusAAC(NumericVector position, String type="PPM",
    double pseudocount=0.0);

StringVector collapse_rows_mat(CharacterMatrix seqs_k);

StringVector collapse_rows_df(DataFrame seqs_k);

StringVector strsplit_cpp(std::string x);

String all_checks_collapse(StringVector checks);

#endif
