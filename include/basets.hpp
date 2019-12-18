///\file basets.hpp
///
///\brief File contains the BaseTS class declaration.
///
///This is the header file of the BaseTS. The BaseTS is a collection of random
///time series generators without injecting any motifs.

#ifndef BASETS_HPP
#define BASETS_HPP

#include <iostream>
#include <random>
#include <chrono>
#include <vector>


using namespace std;


//extern bool check_if_float(const string string_in);

// --------------------------------------------------------------------------
///\brief This class represents the BaseTS.
///
///The BaseTS is a collection of generation methods for random time series
///withouat injecting any motifs.
// --------------------------------------------------------------------------
class BaseTS {

protected:

  // --------------------------------------------------------------------------
  ///\brief This variable stores the random engine.
  ///
  ///The pseudo random engine generates random numbers with the Mersene Twister
  ///19937 generator.
  // --------------------------------------------------------------------------
  mt19937 randomEngine;

public:

  // --------------------------------------------------------------------------
  ///\brief The constructor initializes the BaseTS.
  ///
  ///The constructor checks whether a true random engine is available.
  // --------------------------------------------------------------------------
  BaseTS();

  // --------------------------------------------------------------------------
  ///\brief Frees the memory allocated by the BaseTS.
  ///
  ///The destructor does actually nothing.
  // --------------------------------------------------------------------------
  ~BaseTS();

  ///\brief The simple random walk method.
  ///
  ///\param [out] &timeSeries_out Hands over the random synthetic time series.
  ///\param [in] length_in Hands over the length of the time series.
  ///\param [in] start_in Hands over the start value of the time series.
  ///\param [in] delta_in Hands over the average maximal difference between two
  ///consecutive values in the time series.
  ///\param [in] noise_in Hands over the noise option.
  ///
  ///This function generates a random synthetic time series by performing the
  ///simple random walk. The noise is a value between -noise_in / 2 and
  ///noise_in \ 2. When noise_in is set to 0.0 then the i+1th value v_{i+1} of
  ///the time series is either v_i - delta_in or v_i + delta_in.
  // --------------------------------------------------------------------------
  void simpleRandomWalk(vector<double> &timeSeries_out, int length_in, double
      start_in, double delta_in, double noise_in);

  ///\brief Generates a random syntheteic time series.
  ///
  ///\param [out] &timeSeries_out Hands over the random synthetic time series.
  ///\param [in] length_in Hands over the length of the time series.
  ///\param [in] start_in Hands over the start value of the time series.
  ///\param [in] delta_in Hands over the average maximal difference between two
  ///consecutive values in the time series.
  ///\param [in] noise_in Hands over the noise option.
  ///
  ///This function generates a random synthetic time series by performing the
  ///extended random walk. The noise is a value between -noise_in and noise_in.
  ///The i+1th value v_{i+1} of the time series is an uniform distributed value
  ///between v_i - delta_in and v_i + delta_in with +- noise_in / 2 noise.
  // --------------------------------------------------------------------------
  void realRandomWalk(vector<double> &timeSeries_out, int length_in, double
      start_in, double delta_in, double noise_in);

  ///\brief Generates a random syntheteic time series.
  ///
  ///\param [out] &timeSeries_out Hands over the random synthetic time series.
  ///\param [in] length_in Hands over the length of the time series.
  ///\param [in] start_in Hands over the start value of the time series.
  ///\param [in] delta_in Hands over the average maximal difference between two
  ///consecutive values in the time series.
  ///\param [in] noise_in Hands over the noise option.
  ///
  ///This function generates a random synthetic time series by performing the
  ///extended random walk. The noise is a value between -noise_in / 2 and
  ///noise_in /2. The i+1th value v_{i+1} of the time series is an normal
  ///distributed value between v_i - delta_in and v_i + delta_in with +-
  ///noise_in / 2 noise.
  // --------------------------------------------------------------------------
  void normalRandomWalk(vector<double> &timeSeries_out, int length_in, double
      start_in, double delta_in, double noise_in);

  ///\brief Generates a random syntheteic time series.
  ///
  ///\param [out] &timeSeries_out Hands over the random synthetic time series.
  ///\param [in] length_in Hands over the length of the time series.
  ///\param [in] start_in Hands over the start value of the time series.
  ///\param [in] delta_in Hands over the average maximal difference between two
  ///consecutive values in the time series.
  ///\param [in] step_in Hands over the step size.
  ///\param [in] noise_in Hands over the noise option.
  ///
  ///This function generates a random synthetic time series by performing the
  ///extended random walk. The noise is a value between -noise_in / 2 and
  ///noise_in /2. The I(i+1)th value v_{I(i+1)} of the time series is an normal
  ///distributed value between v_{I(i)} - delta_in and v_{I(i)} + delta_in with
  ///+- noise_in / 2 noise. The step size tells the function to generate
  ///a random next value v_{I(i+1)} poisson distributed difference I(i+1)
  ///- I(i) with mean step_in. Values with indices j, I(i) < j < I(i+1), are
  ///linear approximated.
  // --------------------------------------------------------------------------
  void linearRandomWalk(vector<double> &timeSeries_out, int length_in, double
      start_in, double delta_in, double step_in, double noise_in);

  ///\brief The simple random walk method.
  ///
  ///\param [out] &timeSeries_out Hands over the random synthetic time series.
  ///\param [in] length_in Hands over the length of the time series.
  ///\param [in] start_in Hands over the start value of the time series.
  ///\param [in] delta_in Hands over the average maximal difference between two
  ///consecutive values in the time series.
  ///\param [in] maxi_in Hands over the maximum absolut value.
  ///\param [in] noise_in Hands over the noise option.
  ///
  ///This function generates a random synthetic time series by performing the
  ///simple random walk. The noise is a value between -noise_in / 2 and
  ///noise_in \ 2. When noise_in is set to 0.0 then the i+1th value v_{i+1} of
  ///the time series is either v_i - delta_in or v_i + delta_in. The maxi_in
  ///value makes sure that the absolute value of a times series value is
  ///at most maxi_in.
  // --------------------------------------------------------------------------
  void simpleRandomWalk(vector<double> &timeSeries_out, int length_in, double
      start_in, double delta_in, double maxi_in, double noise_in);

  ///\brief Generates a random syntheteic time series.
  ///
  ///\param [out] &timeSeries_out Hands over the random synthetic time series.
  ///\param [in] length_in Hands over the length of the time series.
  ///\param [in] start_in Hands over the start value of the time series.
  ///\param [in] delta_in Hands over the average maximal difference between two
  ///consecutive values in the time series.
  ///\param [in] maxi_in Hands over the maximum absolut value.
  ///\param [in] noise_in Hands over the noise option.
  ///
  ///This function generates a random synthetic time series by performing the
  ///extended random walk. The noise is a value between -noise_in and noise_in.
  ///The i+1th value v_{i+1} of the time series is an uniform distributed value
  ///between v_i - delta_in and v_i + delta_in with +- noise_in / 2 noise. The
  ///maxi_in value makes sure that the absolute value of a times series value
  ///is at most maxi_in.
  // --------------------------------------------------------------------------
  void realRandomWalk(vector<double> &timeSeries_out, int length_in, double
      start_in, double delta_in, double maxi_in, double noise_in);

  ///\brief Generates a random syntheteic time series.
  ///
  ///\param [out] &timeSeries_out Hands over the random synthetic time series.
  ///\param [in] length_in Hands over the length of the time series.
  ///\param [in] start_in Hands over the start value of the time series.
  ///\param [in] delta_in Hands over the average maximal difference between two
  ///consecutive values in the time series.
  ///\param [in] maxi_in Hands over the maximum absolut value.
  ///\param [in] noise_in Hands over the noise option.
  ///
  ///This function generates a random synthetic time series by performing the
  ///extended random walk. The noise is a value between -noise_in / 2 and
  ///noise_in /2. The i+1th value v_{i+1} of the time series is an normal
  ///distributed value between v_i - delta_in and v_i + delta_in with +-
  ///noise_in / 2 noise. The maxi_in option makes sure that the absolute value
  ///of a times series value is at most maxi_in.
  // --------------------------------------------------------------------------
  void normalRandomWalk(vector<double> &timeSeries_out, int length_in, double
      start_in, double delta_in, double maxi_in, double noise_in);

  ///\brief Generates a random syntheteic time series.
  ///
  ///\param [out] &timeSeries_out Hands over the random synthetic time series.
  ///\param [in] length_in Hands over the length of the time series.
  ///\param [in] start_in Hands over the start value of the time series.
  ///\param [in] delta_in Hands over the average maximal difference between two
  ///consecutive values in the time series.
  ///\param [in] maxi_in Hands over the maximum absolut value.
  ///\param [in] step_in Hands over the step size.
  ///\param [in] noise_in Hands over the noise option.
  ///
  ///This function generates a random synthetic time series by performing the
  ///extended random walk. The noise is a value between -noise_in / 2 and
  ///noise_in /2. The i+1th value v_{i+1} of the time series is an normal
  ///distributed value between v_i - delta_in and v_i + delta_in with +-
  ///noise_in / 2 noise. The step size tells the function to generate a random
  ///next value v_{I(i+1)} poisson distributed difference I(i+1) - I(i) with
  ///mean step_in. Values with indices j, I(i) < j < I(i+1), are linear
  ///approximated. The maxi_in option makes sure that the absolute value of
  ///a times series value without noise is at most maxi_in.
  // --------------------------------------------------------------------------
  void linearRandomWalk(vector<double> &timeSeries_out, int length_in, double
      start_in, double delta_in, double step_in, double maxi_in, double
      noise_in);

  ///\brief Generates a random syntheteic time series.
  ///
  ///\param [out] &timeSeries_out Hands over the random synthetic time series.
  ///\param [in] length_in Hands over the length of the time series.
  ///\param [in] start_in Hands over the start value of the time series.
  ///\param [in] delta_in Hands over the average maximal difference between two
  ///consecutive values in the time series.
  ///\param [in] noise_in Hands over the noise option.
  ///
  ///This function generates a random synthetic time series by computing random
  ///values with uniform distribution. The values v_i are within -delta_in /
  ///2.0 <= v_i < delta_in / 2.0.
  // --------------------------------------------------------------------------
  void uniformRandom(vector<double> &timeSeries_out, int length_in, double
      start_in, double delta_in, double noise_in);

  ///\brief Generates a random syntheteic time series.
  ///
  ///\param [out] &timeSeries_out Hands over the random synthetic time series.
  ///\param [in] length_in Hands over the length of the time series.
  ///\param [in] start_in Hands over the start value of the time series.
  ///\param [in] delta_in Hands over the average maximal difference between two
  ///consecutive values in the time series.
  ///\param [in] noise_in Hands over the noise option.
  ///
  ///This function generates a random synthetic time series by computing random
  ///values with normal distribution variance delta_in / 2.0.
  void normalRandom(vector<double> &timeSeries_out, int length_in, double
      start_in, double delta_in, double noise_in);

  ///\brief Generates a random syntheteic time series.
  ///
  ///\param [out] &timeSeries_out Hands over the random synthetic time series.
  ///\param [in] length_in Hands over the length of the time series.
  ///\param [in] start_in Hands over the start value of the time series.
  ///\param [in] delta_in Hands over the average maximal difference between two
  ///consecutive values in the time series.
  ///\param [in] noise_in Hands over the noise option.
  ///
  ///This function generates a random synthetic time series by computing random
  ///values with piecewise distribution. The pieces are -delta_in / 2.0,
  ///-delta_in / 4.0, 0.0, delta_in / 4.0 and delta_in / 2.0. The corresponding
  ///weights are 0.0, 10.0, 0.0, 10.0, 0.0
  void piecewiseLinearRandom(vector<double> &timeSeries_out, int length_in,
      double start_in, double delta_in, double noise_in);
};

#endif
