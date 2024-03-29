///\file tsgutils.hpp
///
///\brief File contains additional utilities.
///
///This is a set of additional utilities for time series modifications and
///analysis.

#ifndef TSGUTILS_HPP
#define TSGUTILS_HPP

#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <tsgenerator.hpp>
#include <tsgtypes.hpp>
#include <basets.hpp>


namespace tsg {

  ///\brief An algorithm that computes the minimum and maximum distance.
  ///
  ///\param [in] &timeSeries Hands over the time series.
  ///\param [in] &sums_in Hands over the running sum.
  ///\param [in] &sumSquares_in Hands over the running sum of squares.
  ///\param [in] window_in Hands over the window size.
  ///\param [out] &min_out Returns the minimum distance of two subsequences.
  ///\param [out] &max_out Returns the maximum distance of two subsequences.
  ///
  ///This algotihm computes the minimum and maximum z-normalized Euclidean
  ///distance of two non overlapping subsequences.
  void minMaxDist(const rseq &timeSeries_in, const rseq &sums_in, const rseq
      &sumSquares_in, const int window_in, double &min_out, double &max_out) {

    int window = window_in;
    double rWindow = 1.0 / window;

    int length = (int)(timeSeries_in.size());

    //precompute running mean and standard deviation
    rseq mean;
    for (int i = 0; i < length - window + 1; i++)
      mean.push_back(sums_in[i] * rWindow);
    rseq var;
    for (int i = 0; i < length - window + 1; i++)
      var.push_back(sumSquares_in[i] * rWindow - mean[i] * mean[i]);
    rseq sigma;
    for (int i = 0; i < length - window + 1; i++)
      sigma.push_back(var[i] > 1.0 ? sqrt(var[i]) : 1.0);

    double q = 0.0;
    double distance;
    min_out = std::numeric_limits<double>::infinity();
    max_out = 0.0;

    //for each diagonal of non overlapping subsequences
    for (int k = window - 1; k < length - window + 1; k++) {

      //the dot product
      q = 0.0;

      //iterate throw the diagonal
      for (int i = 0; i < length - window + 1 - k; i++) {

        if (i == 0)
          //compute the initial dot product
          for (int j = 0; j < window; j++)
            q += timeSeries_in[i + j] * timeSeries_in[k + j];
        else
          //compute dot product iteratively
          q += timeSeries_in[i + window - 1] * timeSeries_in[i + k + window
            - 1] - timeSeries_in[i - 1] * timeSeries_in[i + k - 1];

        //compute distance with dot product
        distance = 2 * (window - (q - window * mean[i] * mean[i + k])
            / (sigma[i] * sigma[i + k]));

        //save the best so far distances
        if (distance < min_out)
          min_out = distance;

        if (distance > max_out)
          max_out = distance;
      }
    }
  };

  ///\brief Computes the running sum and sum of squares
  ///
  ///\param [in] &sequence_in Hands over the objective sequence.
  ///\param [in] window_in Hands over the window size.
  ///\param [out] &sums_out Returns the rolling sum.
  ///\param [out] &sumSquares_out Returns the rolling sum of squares.
  ///
  ///This algorithm computes the rolling sum and sum of squares for a given
  ///objective sequence and window.
  void calcRunnings(const rseq &sequence_in, const int window_in, rseq
      &sums_out, rseq &sumSquares_out) {

    // obtain running sum and sum of square
    if (!sums_out.empty() || !sumSquares_out.empty()) {

      sums_out.clear();
      sumSquares_out.clear();
    }

    sums_out.resize(sequence_in.size() - window_in + 1);
    sumSquares_out.resize(sequence_in.size() - window_in + 1);

    sums_out[0] = sequence_in[0];
    sumSquares_out[0] = sequence_in[0] * sequence_in[0];

    for (int i = 1; i < window_in; i++) {

      sums_out[0] += sequence_in[i];
      sumSquares_out[0] += sequence_in[i] * sequence_in[i];
    }

    for (int i = 0; i < (int)sequence_in.size() - window_in; i++) {

      sums_out[i + 1] = sums_out[i] - sequence_in[i] +
        sequence_in[i + window_in];

      sumSquares_out[i + 1] = sumSquares_out[i] - sequence_in[i] *
        sequence_in[i] + sequence_in[i + window_in] *
        sequence_in[i + window_in];
    }
  };

  ///\brief Computes a base time series
  ///
  ///\param [out] &timeSeries_out Returns the generated time series.
  ///\param [in] method_in Hands over the generation method.
  ///\param [in] length_in Hands over the time series length.
  ///\param [in] delta_in Hands over the maximum absolut difference of two
  ///consecutive values in the time series.
  ///\param [in] step_in Hands over the maximum step size for the next value in
  ///x direction
  ///\param [in] times_in Hands over the maximum number of values per step.
  ///\param [in] maxi_in Hands over the maximum absolut value.
  ///\param [in] noise_in Hands over the maximum absolut difference a value is
  ///displaced to add noise.
  ///
  ///This algorithm computes a times series using one of the methods defined in
  ///the BaseTS class.
  void generateBaseTimeSeries(rseq &timeSeries_out, const word method_in, const
      int length_in, const double delta_in, const double step_in, const int
      times_in, const double maxi_in, const double noise_in) {

    if (!timeSeries_out.empty()) {

      timeSeries_out.clear();
      timeSeries_out.resize(0);
    }

    BaseTS baseTS;

    // get the method
    int method = (int)(std::distance(methods.begin(),
          std::find(methods.begin(), methods.end(), method_in)));

    // check if method exists
    if (method >= (int) methods.size()) {

      std::cerr << "ERROR: Unknown method: " << method_in << std::endl;
      throw(EXIT_FAILURE);
    }

    switch (method) {

      case 0:
        baseTS.simpleRandomWalk(timeSeries_out, length_in, delta_in, noise_in);
        break;
      case 1:
        baseTS.realRandomWalk(timeSeries_out, length_in, delta_in, noise_in);
        break;
      case 2:
        baseTS.normalRandomWalk(timeSeries_out, length_in, delta_in, noise_in);
        break;
      case 3:
        baseTS.linearRandomWalk(timeSeries_out, length_in, delta_in, step_in,
            noise_in);
        break;
      case 4:
        baseTS.simpleRandomWalk(timeSeries_out, length_in, delta_in, maxi_in,
            noise_in);
        break;
      case 5:
        baseTS.realRandomWalk(timeSeries_out, length_in, delta_in, maxi_in,
            noise_in);
        break;
      case 6:
        baseTS.normalRandomWalk(timeSeries_out, length_in, delta_in, maxi_in,
            noise_in);
        break;
      case 7:
        baseTS.linearRandomWalk(timeSeries_out, length_in, delta_in, step_in, maxi_in,
            noise_in);
        break;
      case 8:
        baseTS.uniformRandom(timeSeries_out, length_in, delta_in, noise_in);
        break;
      case 9:
        baseTS.normalRandom(timeSeries_out, length_in, delta_in, noise_in);
        break;
      case 10:
        baseTS.piecewiseLinearRandom(timeSeries_out, length_in, delta_in,
            noise_in);
        break;
      case 11:
        baseTS.splineRepeated(timeSeries_out, length_in, delta_in, (int)step_in, times_in,
            noise_in);
        break;
        default:
        std::cerr << "ERROR: Unknown method: " << method << std::endl;
        throw(EXIT_FAILURE);
    }
  };

  ///\brief Computes the distance of two subsequences in a sequence.
  ///
  ///\param [in] &sequence_in Hands over the sequence.
  ///\param [in] pos0_in Hands over the position of the first subsequence.
  ///\param [in] pos1_in Hands over the position of the second subsequence.
  ///\param [in] window_in Hands over the window size.
  ///\param [in] sums_in Hands over the running sum.
  ///\param [in] sumSquares_in Hands over the running sum of squares.
  ///\param [in] bestSoFar_in Hands over the best similarity so far.
  ///
  ///\return The distance of the two z-normalized subsequences.
  ///
  ///This function computes the distance of two subsequences. Therefore, the
  ///subsequences are first z-normalized and the Euclidean Distance is
  ///computed. The return value is the distance of the two z-normalized
  ///subsequences.
  double zNormEuclDist(const rseq &sequence_in, const int pos0_in, const int
      pos1_in, const int window_in, const tsg::rseq sums_in, const tsg::rseq
      sumSquares_in, const double bestSoFar_in) {

    tsg::rseq sums = sums_in;
    tsg::rseq sumSquares = sumSquares_in;
    int length = sequence_in.size();
    int window = window_in;
    if (sequence_in.empty()) {

      std::cerr << "ERROR: Time series is empty!" << std::endl;
      throw(EXIT_FAILURE);
    }

    if (pos0_in + window > length || pos0_in < 0) {

      std::cerr << "ERROR: Position of first subsequence " <<
        pos0_in << " in similarity function is wrong!" << std::endl;
      throw(EXIT_FAILURE);
    }

    if (pos1_in + window > length || pos1_in < 0) {

      std::cerr << "ERROR: Position of second subsequence " <<
        pos1_in << " in similarity function is wrong!" << std::endl;
      throw(EXIT_FAILURE);
    }

    double rWindow = 1.0 / window;

    double mean0 = sums[pos0_in] * rWindow;
    double stdDev0 = sumSquares[pos0_in] * rWindow  - mean0 * mean0;
    stdDev0 = stdDev0 < 1.0 ? 1.0 : sqrt(stdDev0);

    double mean1 = sums[pos1_in] * rWindow;
    double stdDev1 = sumSquares[pos1_in] * rWindow  - mean1 * mean1;
    stdDev1 = stdDev1 < 1.0 ? 1.0 : sqrt(stdDev1);

    //calculate the similarity
    double sumOfSquares = 0.0;
    double bestSoFar = bestSoFar_in * bestSoFar_in;
    double norm0;
    double norm1;
    double diff;

    for (int i = 0; i < window && sumOfSquares < bestSoFar; i++) {

      norm0 = (sequence_in[pos0_in + i] - mean0) / stdDev0;
      norm1 = (sequence_in[pos1_in + i] - mean1) / stdDev1;
      diff = norm0 - norm1;
      sumOfSquares += diff * diff;
    }

    return sqrt(sumOfSquares);
  };
}

#endif
