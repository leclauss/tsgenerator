///\file tsm.hpp
///
///\brief File contains a top set motif discovery algorithm.
///
///This algorithm computes the top set motif of a time series by evaluating
///iteratively the diagonals of the distance matrix of the time series.

#ifndef TSM_HPP
#define TSM_HPP

#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <map>
#include <tsgtypes.hpp>


namespace tsg {

  ///\brief A PAA implementation.
  ///
  ///\param [in] &timeSeries Hands over the time series.
  ///\param [in] &sums_in Hands over the running sum.
  ///\param [in] &sumSquares_in Hands over the running sum of squares.
  ///\param [out] &paa_out Returns the PAA representation.
  ///\param [in] pos_in Hands over the subsequence position in the time series.
  ///\param [in] window_in Hands over the window size.
  ///\param [in] paa_in Hands over the paa size.
  ///
  ///This is a PAA implementation. This algorithm computes the PAA
  ///representation of a z normalized subsequence.
  void zNormPAA(const rseq &timeSeries_in, const rseq &sums_in, const rseq
      &sumSquares_in, rseq &paa_out, const int pos_in, const int window_in,
      const int paa_in = 6);

  ///\brief An inverse normal CDF implementation.
  ///
  ///\param [in] p_in Hands over the CDF value.
  ///\param [in] prec_in Hands over the absolute err.
  ///
  ///\return Hands over the inverse normal CDF.
  ///
  ///This is an inverse normal CDF implementation. It uses the gauss error
  ///function and a simple divide and conquer technique.
  double invNormalCDF(const double p_in, const double prec_in = 0.00001);

  ///\brief An inverse normal CDF implementation.
  ///
  ///\param [in] alphabet_in Hands over the number of areas with equal size
  ///under the gaussian curve.
  ///\param [out] &p_out Hands over the breakpoints.
  ///\param [in] prec_in Hands over the absolute err.
  ///
  ///This is an inverse normal CDF implementation. It uses the gauss error
  ///function and a simple divide and conquer technique.
  void invNormalCDF(const int alphabet_in, rseq &p_out, const double prec_in
      = 0.00001);

  ///\brief A SAX implementation.
  ///
  ///\param [in] &timeSeries Hands over the time series.
  ///\param [in] &sums_in Hands over the running sum.
  ///\param [in] &sumSquares_in Hands over the running sum of squares.
  ///\param [out] &sax_out Returns the SAX representation.
  ///\param [in] pos_in Hands over the subsequence position in the time series.
  ///\param [in] window_in Hands over the window size.
  ///\param [in] paa_in Hands over the paa size.
  ///\param [in] &breakpoints_in Hands over the breakpoints.
  ///
  ///This is a SAX implementation. This algorithm computes the SAX
  ///representation of a z normalized subsequence.
  void zNormalSAX(const rseq &timeSeries_in, const rseq &sums_in, const rseq
      &sumSquares_in, word &sax_out, const int pos_in, const int window_in,
      const int paa_in, const rseq &breakpoints_in
      = { -std::numeric_limits<double>::infinity(), -0.967422, -0.430727, 0.0,
      0.430727, 0.967422, std::numeric_limits<double>::infinity() });

  ///\brief A top set motif dicovery procedure.
  ///
  ///\param [in] &timeSeries Hands over the time series.
  ///\param [in] &sums_in Hands over the running sum.
  ///\param [in] &sumSquares_in Hands over the running sum of squares.
  ///\param [out] &motif_out Returns the position of the motifs and non-self
  ///matching subsequences.
  ///\param [in] window_in Hands over the window size.
  ///\param [in] paa_in Hands over the paa size.
  ///\param [in] &breakpoints_in Hands over the breakpoints.
  ///
  ///This is a top set motif discovery algorithm. This algorithm computes the
  ///top set motif of a time series by evaluating iteratively the diagonals of
  ///the distance matrix of the time series.
  void tsm(const rseq &timeSeries_in, const rseq &sums_in, const rseq
      &sumSquares_in, iseqs &motif_out, const int window_in, const int paa_in
      = 6, const rseq &breakpoints_in
      = { -std::numeric_limits<double>::infinity(), -0.967422, -0.430727, 0.0,
      0.430727, 0.967422, std::numeric_limits<double>::infinity() });
}

#endif
