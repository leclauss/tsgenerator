///\file tpm.hpp
///
///\brief File contains a top pair motif discovery algorithm.
///
///This algorithm computes the top pair motif of a time series by evaluating
///iteratively the diagonals of the distance matrix of the time series.

#ifndef TPM_HPP
#define TPM_HPP

#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <tsgtypes.hpp>


namespace tsg {

  ///\brief The top pair motif dicovery procedure.
  ///
  ///\param [in] &timeSeries Hands over the time series.
  ///\param [in] &sums_in Hands over the running sum.
  ///\param [in] &sumSquares_in Hands over the running sum of squares.
  ///\param [out] &pos0_out Returns the position of the first subsequence.
  ///\param [out] &pos1_out Returns the position of the second subsequence.
  ///\param [in] window_in Hands over the window size.
  ///
  ///This is a top pair motif discovery algorithm. This algorithm computes the
  ///top pair motif of a time series by evaluating iteratively the diagonals of
  ///the distance matrix of the time series.
  double tpm(const rseq &timeSeries_in, const rseq &sums_in, const rseq
      &sumSquares_in, int &pos0_out, int &pos1_out, const int window_in);
}

#endif
