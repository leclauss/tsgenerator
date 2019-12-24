///\file scrimpplusplus.hpp
///
///\brief File contains the Scrimp++ algorithm.
///
///This is SCRIMP++.
///
///Details of the SCRIMP++ algorithm can be found at:
///(author information ommited for ICDM review),
///"SCRIMP++: Motif Discovery at Interactive Speeds", submitted to ICDM 2018.
///
///edited by Rafael Moczalla

#ifndef SCRIMPPLUSPLUS_HPP
#define SCRIMPPLUSPLUS_HPP

#include <cstdlib>
#include <fftw3.h>
#include <cmath>
#include <iomanip>
#include <limits>
#include <vector>
#include <algorithm>
#include <string>
#include <chrono>
#include <random>
#include <tsgtypes.hpp>


namespace tsg {

  ///\brief The Scrimp++ procedure.
  ///
  ///\param [in] &timeSeries Hands over the time series.
  ///\param [in] &ASum_in Hands over the running sum.
  ///\param [in] &ASumSq_in Hands over the running sum of squares.
  ///\param [out] &pos0_out Returns the position of the first subsequence.
  ///\param [out] &pos1_out Returns the position of the second subsequence.
  ///\param [in] window_in Hands over the window size.
  ///\param [in] stepSize_in Hands over the step size for PreScrimp.
  ///
  ///This is the Scrimp++ implementation edited to return the most similar
  ///subsequence pair in a time series.
  void scrimpPP(const rseq &timeSeries_in, const rseq &ASum_in, const rseq
      &ASumSq_in, int &pos0_out, int &pos1_out, const int window_in, const
      double stepSize_in = -1.0);
}

#endif
