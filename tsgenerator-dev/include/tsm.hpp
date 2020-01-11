///\file tsm.hpp
///
///\brief File contains a top set motif discovery algorithm.
///
///This algorithm computes the top set motif of a time series by evaluating
///iteratively the diagonals of the distance matrix of the time series.

#ifndef TSM_HPP
#define TSM_HPP

#include <cmath>
#include <limits>
#include <algorithm>
#include <random>
#include <chrono>
#include <tsgtypes.hpp>


namespace tsg {

  ///\brief This class encapsulates a top set motif discovery algorithm.
  ///
  ///The TSM class contains all data to compute the top set motif in
  ///an objective sequence alias time series.
  class TSM {

    protected:
      ///\brief This variable stores the time series.
      ///
      ///This variable stores the objective sequence alias time series.
      const rseq timeSeries;

      ///\brief This variable stores the running sum.
      ///
      ///This variable stores the running sum of the objective sequence.
      const rseq sums;

      ///\brief This variable stores the running sum of squares.
      ///
      ///This variable stores the running sum of squares of the objective
      ///sequence.
      const rseq sumSquares;

      ///\brief This variable stores the ADM distance matrix.
      ///
      ///This variable stores the ADM distance matrix of a neighborhood.
      rseqs admDist;

    public:
      ///\brief The constructor sets up the top set motif object.
      ///
      ///\param [in] &timeSeries Hands over the time series.
      ///\param [in] &sums_in Hands over the running sum.
      ///\param [in] &sumSquares_in Hands over the running sum of squares.
      ///
      ///The constructor initializes the top set motif object including the
      ///time series, running sum and sum of squares.
      TSM(const rseq &timeSeries_in, const rseq &sums_in, const rseq
          &sumSquares_in);

      ///\brief A PAA implementation.
      ///
      ///\param [out] &paa_out Returns the PAA representation.
      ///\param [in] pos_in Hands over the subsequence position in the time series.
      ///\param [in] window_in Hands over the window size.
      ///\param [in] paa_in Hands over the paa size.
      ///
      ///This is a PAA implementation. This algorithm computes the PAA
      ///representation of a z normalized subsequence.
      void zNormalPAA(rseq &paa_out, const int pos_in, const int window_in, const
          int paa_in = 6);

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
      ///\param [out] &sax_out Returns the SAX representation.
      ///\param [in] pos_in Hands over the subsequence position in the time series.
      ///\param [in] window_in Hands over the window size.
      ///\param [in] paa_in Hands over the paa size.
      ///\param [in] &breakpoints_in Hands over the breakpoints.
      ///
      ///This is a SAX implementation. This algorithm computes the SAX
      ///representation of a z normalized subsequence.
      void zNormalSAX(word &sax_out, const int pos_in, const int window_in,
          const int paa_in, const rseq &breakpoints_in
          = { -std::numeric_limits<double>::infinity(), -0.967422, -0.430727, 0.0,
          0.430727, 0.967422, std::numeric_limits<double>::infinity() });

      ///\brief A SAX minimum distance implementation.
      ///
      ///\param [in] &word0_in Hands over the first SAX word.
      ///\param [in] &word1_in Hands over the second SAX word.
      ///\param [in] window_in Hands over the window size.
      ///\param [in] &breakpoints_in Hands over the SAX breakpoints.
      ///
      ///\return Returns the minimum distance.
      ///
      ///This is an implementation of a minimum distance of two SAX words.
      double saxDist(const word &word0_in, const word word1_in, const int
          window_in, const rseq &breakpoints_in);

      ///\brief A z-normalized Euclidean distance implementation.
      ///
      ///\param [in] &pos0_in Hands over the position of the first subsequence.
      ///\param [in] &pos1_in Hands over the position of the second
      ///subsequence.
      ///\param [in] window_in Hands over the window size.
      ///
      ///\return Returns the z-normalized Euclidean distance.
      ///
      ///This is an implementation of the z-normalized Euclidean distance.
      double dist(const int pos0_in, const int pos1_in, const int window_in);

      ///\brief ADM distance getter function.
      ///
      ///\param [in] &i_in Hands over the index of the first neighbor
      ///subsequence.
      ///\param [in] &j_in Hands over the index of the second neighbor
      ///subsequence.
      ///
      ///\return Returns ADM distance.
      ///
      ///This getter function returns the ADM distance from the ith to the jth
      ///neighbor.
      double distADM(const int i_in, const int j_in);

      ///\brief A modified ADM algorithm implementation.
      ///
      ///\param [in] &neighborhood_in Hands over the subsequence indices.
      ///\param [in] window_in Hands over the window size.
      ///\param [in] range_in Hands over the motif range.
      ///
      ///This is a modified ADM algorithm implementation which computes the
      ///distance matrix of a set of subsequences.
      void adm(const iseq &neighborhood_in, const int window_in, const double
          range_in);

      ///\brief A top set motif dicovery procedure.
      ///
      ///\param [out] &motif_out Returns the position of the motifs and non-self
      ///matching subsequences.
      ///\param [in] window_in Hands over the window size.
      ///\param [in] range_in Hands over the motif range.
      ///\param [in] paa_in Hands over the paa size.
      ///\param [in] &breakpoints_in Hands over the breakpoints.
      ///
      ///\return The optimal number of subsequences in a set motif.
      ///
      ///This is a top set motif discovery algorithm. This algorithm computes the
      ///top set motif of a time series by evaluating iteratively the diagonals of
      ///the distance matrix of the time series.
      int tsm(iseq &motif_out, const int window_in, const double range_in,
          const int paa_in = 6, const rseq &breakpoints_in
          = { -std::numeric_limits<double>::infinity(), -0.967422, -0.430727, 0.0,
          0.430727, 0.967422, std::numeric_limits<double>::infinity() });
  };
}

#endif
