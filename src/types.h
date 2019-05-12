#ifndef _TYPES_
#define _TYPES_

#include <vector>
#include <string>

/* length one character vector */
typedef std::string str_t;

/* logical vector */
typedef std::vector<bool> vec_bool_t;

/* numeric vector */
typedef std::vector<double> vec_num_t;
typedef std::vector<long double> vec_lnum_t;

/* character vector */
typedef std::vector<std::string> vec_str_t;

/* integer vector */
typedef std::vector<int> vec_int_t;

/* c-type char vector */
typedef std::vector<char> vec_char_t;

/* list of integer vectors OR integer matrix (mat[i][j]; i = col, j = row) */
typedef std::vector<std::vector<int>> list_int_t;

/* list of numeric vectors OR numeric matrix (mat[i][j]; i = col, j = row) */
typedef std::vector<std::vector<double>> list_num_t;
typedef std::vector<std::vector<long double>> list_lnum_t;

/* list of c-type char vectors */
typedef std::vector<std::vector<char>> list_char_t;

/* list of lists of integer vectors OR list of matrices */
typedef std::vector<std::vector<std::vector<int>>> list_mat_t;

#endif
