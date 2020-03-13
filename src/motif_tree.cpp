#include <Rcpp.h>
#include "types.h"

void reorder_internal_phylo(Rcpp::IntegerVector &new_order, const int &node,
    const Rcpp::IntegerVector &e1, const Rcpp::IntegerVector &e2,
    const R_xlen_t &n_tips, const vec_int_t &xi, const vec_int_t &xj,
    const vec_int_t &L, R_xlen_t &iii) {

  R_xlen_t i = node - n_tips - 1, k;
  for (R_xlen_t j = 0; j < xj[i]; ++j) {
    k = L[xi[i] + j];
    new_order[iii++] = k + 1;
    if (e2[k] > n_tips)
      reorder_internal_phylo(new_order, e2[k], e1, e2, n_tips, xi, xj, L, iii);
  }

}

Rcpp::IntegerVector get_new_phylo_order(const Rcpp::IntegerMatrix &edge,
    const R_xlen_t &n_tips) {

  // Something isn't working here

  R_xlen_t root = n_tips + 1, k, j;

  Rcpp::IntegerVector e1 = edge(Rcpp::_, 0);
  Rcpp::IntegerVector e2 = edge(Rcpp::_, 1);

  R_xlen_t m = Rcpp::max(e1), nnode = m - n_tips, n = edge.nrow(), iii = 0;
  vec_int_t L(n), pos(nnode), xi(nnode), xj(nnode);
  Rcpp::IntegerVector new_order(n);

  for (R_xlen_t i = 0; i < n; ++i) {
    ++xj[e1[i] - n_tips - 1];
  }

  for (R_xlen_t i = 1; i < nnode; ++i) {
    xi[i] = xi[i - 1] + xj[i - 1];
  }

  for (R_xlen_t i = 0; i < n; ++i) {
    k = e1[i] - n_tips - 1;
    j = pos[k];
    L[xi[k] + j] = i;
    ++pos[k];
  }

  reorder_internal_phylo(new_order, root, e1, e2, n_tips, xi, xj, L, iii);

  return new_order;

}

Rcpp::IntegerMatrix reorder_phylo_edge(const Rcpp::IntegerMatrix &edge,
    const Rcpp::IntegerVector new_order) {
  Rcpp::IntegerMatrix new_edge(edge.nrow(), edge.ncol());
  for (R_xlen_t i = 0; i < new_order.size(); ++i) {
    new_edge(new_order[i] - 1, Rcpp::_) = edge(i, Rcpp::_);
  }
  return new_edge;
}

Rcpp::NumericVector reorder_phylo_edge_length(const Rcpp::NumericVector &edge_length,
    const Rcpp::IntegerVector new_order) {
  Rcpp::NumericVector new_edge_length(edge_length.size());
  for (R_xlen_t i = 0; i < edge_length.size(); ++i) {
    new_edge_length[new_order[i] - 1] = edge_length[i];
  }
  return new_edge_length;
}

/* C++ ENTRY ---------------------------------------------------------------- */

// [[Rcpp::export(rng = false)]]
Rcpp::List hclust_to_phylo_cpp(const Rcpp::IntegerMatrix x_merge,
    const Rcpp::NumericVector x_height, const Rcpp::RObject x_labels) {

  R_xlen_t n = x_merge.nrow();
  Rcpp::IntegerMatrix edge(n * 2, 2);
  Rcpp::NumericVector edge_length(2 * n);
  vec_int_t node(n);
  R_xlen_t cur_nod = n + 3, j = 0, k, y;

  node[n - 1] = n + 2;

  for (R_xlen_t i = n - 1; i >= 0; --i) {
    edge(j, 0) = node[i];
    edge(j + 1, 0) = node[i];
    for (R_xlen_t l = 0; l <= 1; ++l) {
      k = j + l;
      y = x_merge(i, l);
      if (y > 0) {
        node[y - 1] = cur_nod;
        edge(k, 1) = cur_nod;
        ++cur_nod;
        edge_length[k] = x_height[i] - x_height[y];
      } else {
        edge(k, 1) = -y;
        edge_length[k] = x_height[i];
      }
    }
    j += 2;
  }

  Rcpp::StringVector labels(n + 1);
  if (x_labels.isNULL()) {
    for (int i = 1; i <= n + 1; ++i) {
      labels[i - 1] = i;
    }
  } else {
    labels = Rcpp::as<Rcpp::StringVector>(x_labels);
  }

  // Getting new order is behaving differently from ape
  Rcpp::IntegerVector new_order = get_new_phylo_order(edge, n + 1);
  edge = reorder_phylo_edge(edge, new_order);
  edge_length = reorder_phylo_edge_length(edge_length / 2, new_order);

  return Rcpp::List::create(
    Rcpp::_["edge"] = edge,
    Rcpp::_["edge.length"] = edge_length,
    Rcpp::_["tip.label"] = labels,
    Rcpp::_["Nnode"] = n
  );

}
