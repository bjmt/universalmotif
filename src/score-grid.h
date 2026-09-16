#ifndef UNIVERSALMOTIF_SCORE_GRID_H
#define UNIVERSALMOTIF_SCORE_GRID_H

#include <algorithm>
#include <cmath>
#include <limits>

// Scores are sums of integerised PWM entries (entry * 1000, truncated).
// A query/threshold must be rounded UP, not truncated, to implement >=.
// Snap only within floating-point round-trip error of an integer, so values
// such as (1001 / 1000.0) * 1000 recover their original grid position.
inline double score_grid_ceiling(double score) {
  double scaled = score * 1000.0;
  if (!std::isfinite(scaled)) return scaled;
  double nearest = std::round(scaled);
  double tolerance = 4 * std::numeric_limits<double>::epsilon() *
    std::max(1.0, std::abs(scaled));
  if (std::abs(scaled - nearest) <= tolerance) scaled = nearest;
  return std::ceil(scaled);
}

// Bounds and offset remain in integer score units throughout inversion.
// Empty score bins are valid cutoffs; thresholds beyond the highest possible
// score are represented explicitly as +Inf rather than clamped to a hit.
inline double score_grid_threshold(const double *cdf, double offset,
                                   long minimum, long maximum, double pvalue) {
  if (pvalue >= 1.0) return minimum / 1000.0;
  for (long score = minimum; score <= maximum; ++score) {
    long index = static_cast<long>(score - offset);
    if (cdf[index] <= pvalue) return score / 1000.0;
  }
  return std::numeric_limits<double>::infinity();
}

#endif
