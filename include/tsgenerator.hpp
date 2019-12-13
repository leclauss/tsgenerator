///\file tsgenerator.hpp
///
///\brief File contains the TSGenerator class declaration.
///
///This is the header file of the TSGenerator. The TSGenerator consists of
///setter functions and a run function. After configuring the TSGenerator with
///the setter functions one starts the time series generation by calling the
///run function.

#ifndef TSGENERATOR_HPP
#define TSGENERATOR_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <string>
#include <random>
#include <chrono>
#include <cmath>
#include <cfloat>
#include <limits>
#include <iterator>
#include <vector>
#include <tuple>
#include <fftw3.h>
#include <motifsetcollection.hpp>
#include <freepositions.hpp>
#include <basets.hpp>
#include <global.hpp>


using namespace std;


//extern bool check_if_float(const string string_in);

// --------------------------------------------------------------------------
///\brief This class represents the TSGenerator.
///
///The TSGenerator consists of setter functions and a run function as well as
///various variables. After configuring the TSGenerator with the setter
///functions one starts the time series generation by calling the run function.
// --------------------------------------------------------------------------
class TSGenerator {

public:

  // --------------------------------------------------------------------------
  ///\brief This struct represents a time series subsequence.
  ///
  ///The time series subsequence consists of a postion and length value.
  // --------------------------------------------------------------------------
  struct subsequence {

    int position = 0;
    int length = 0;
  };

protected:

  // --------------------------------------------------------------------------
  ///\brief This variable contains the base value for the time series values.
  ///
  ///This variable stores the base value of the time series values. The time
  ///series is generated around this value.
  // --------------------------------------------------------------------------
  double baseValue = 0.0;

  // --------------------------------------------------------------------------
  ///\brief This variable contains the length of the time series.
  ///
  ///The variable stores the length of the time series.
  // --------------------------------------------------------------------------
  int length = -1;

  // --------------------------------------------------------------------------
  ///\brief This variable the window size.
  ///
  ///The window size is used to generate the motifs such that the window size
  ///is the motif length.
  // --------------------------------------------------------------------------
  int window = -1;

  // --------------------------------------------------------------------------
  ///\brief This variable contains the maximum difference of two consecutive
  ///values.
  ///
  ///This variable stores the maximum absolute diffrenece of two consecutive
  ///values in the time series excluding noise and motifs.
  // --------------------------------------------------------------------------
  double delta = 1.0;

  // --------------------------------------------------------------------------
  ///\brief This variable contains the noise option.
  ///
  ///This variable stores the noise option that will be used to add noise to
  ///the time series. The resulting noise is a value between -noise / 2 and
  ///noise / 2.
  // --------------------------------------------------------------------------
  double noise = 0.1;

  // --------------------------------------------------------------------------
  ///\brief This list contains all motif types.
  ///
  ///This variable stores the names of all motifs available in this
  ///implementation. To add a motif one has to declare the motif in the
  ///motifcollection.hpp file and define the motif in the motifcollection.cpp
  ///file as well. Also, one has to modify the addMotif() function in the file
  ///tsgenerator.cpp.
  // --------------------------------------------------------------------------
  vector<string> motifTypes{"box", "triangle", "semicircle", "trapezoid",
    "positiveflank", "negativeflank", "sine", "cosine"};

  // --------------------------------------------------------------------------
  ///\brief This variable contains the motif type.
  ///
  ///This variable stores the motif type, i.e., the shape of the motif.
  // --------------------------------------------------------------------------
  int type = 0;

  // --------------------------------------------------------------------------
  ///\brief This variable contains the motif size.
  ///
  ///This variable stores the motif size selected, i.e., the number of
  ///subsequences non-self matched by the motif.
  // --------------------------------------------------------------------------
  int size = 3;

  // --------------------------------------------------------------------------
  ///\brief This variable contains the motif height.
  ///
  ///This variable stores the motif height, i.e., the maximum difference
  ///between two motif values.
  // --------------------------------------------------------------------------
  double height = 10.0;

  // --------------------------------------------------------------------------
  ///\brief This variable contains the free positions.
  ///
  ///This variable stores the positions of free subsequences in the time series
  ///between the motif sets subsequences.
  // --------------------------------------------------------------------------
  FreePositions freePositions;

  // --------------------------------------------------------------------------
  ///\brief This variable stores the random engine.
  ///
  ///The pseudo random engine generates random numbers with the Mersene Twister
  ///19937 generator.
  // --------------------------------------------------------------------------
  mt19937 randomEngine;

  // --------------------------------------------------------------------------
  ///\brief This variable contains a base time series generator.
  ///
  ///This variable stores a collection of methods to generate a basic time
  ///series without injected motifs..
  // --------------------------------------------------------------------------
  BaseTS baseTS;

  // --------------------------------------------------------------------------
  ///\brief Backup of the time series.
  ///
  ///This variable stores a backup of the time series.
  // --------------------------------------------------------------------------
  vector<double> timeSeriesBackup;

  // --------------------------------------------------------------------------
  ///\brief Means and standard deviations.
  ///
  ///This variable stores the means and standard deviations of all time series
  ///database time series.
  // --------------------------------------------------------------------------
  vector<tuple<double, double>> meansStds;


  // --------------------------------------------------------------------------
  ///\brief Calculates the rolling mean and standard deviation of a time series.
  ///
  ///\param [in] &timeSeries_in Hands over the time series.
  ///
  ///\return The rolling mean and standard deviation.
  ///
  ///This function calculates the rolling mean and standard deviation of a time
  ///series. The output tuple contains the means and standard deviations.
  // --------------------------------------------------------------------------
  tuple<vector<double>, vector<double>> calcRollingMeanStdDev(const
      vector<double> &timeSeries_in);

  // --------------------------------------------------------------------------
  ///\brief Update the rolling mean and standard deviation of a time series.
  ///
  ///\param [in] &timeSeries_in Hands over the time series.
  ///\param [in] pos_in Hands over the position of the injected subsequence.
  ///\param [in, out] &mean Hands over the rolling mean.
  ///\param [in, out] &stdDev Hands over the rolling std deviation.
  ///
  ///This function updates the rolling mean and standard deviation of a time
  ///series at a specific location.
  // --------------------------------------------------------------------------
  void updateRollingMeanStdDev(const vector<double> &timeSeries_in, int pos_in,
      vector<double> &mean, vector<double> & stdDev);

  // --------------------------------------------------------------------------
  ///\brief Calculates the similarity of two time series.
  ///
  ///\param [in] &timeSeriesOne_in Hands over the time series.
  ///\param [in] subsequenceOnePos_in Hands over the position of the first
  ///subsequence in the time series.
  ///\param [in] meanOne_in Hands over the mean of timeSeriesOne_in.
  ///\param [in] stdDevOne_in Hands over the standard deviation of
  ///timeSeriesOne_in.
  ///\param [in] subsequenceTwoPos_in Hands over the position of the second
  ///subsequence in the time series.
  ///\param [in] bestSoFar_in Hands over the best similarity so far.
  ///
  ///\return The similarity of the two z-normalized time series.
  ///
  ///This function calculates the similarity of two time series. Therefore, the
  ///time series are first z-normalized and the Euclidean Distance is
  ///calculated. The return value is the similarity of the two z-normalized
  ///time series.
  // --------------------------------------------------------------------------
  double similarity(const vector<double> &timeSeries_in, const int
      subsequenceOnePos_in, const int subsequenceTwoPos_in, const double
      bestSoFar_in, const vector<double> &means, const vector<double> &stds);

  // --------------------------------------------------------------------------
  ///\brief Calculates the mean and standard deviation.
  ///
  ///\param [in] &timeSeries_in Time series to calculate the mean and standard
  ///deviation.
  ///\param [in] subsequencePos_in Hands over the position of the subsequence
  ///in the time series.
  ///\param [out] &mean_out The mean of the time series.
  ///\param [out] &stddev_out The standard deviation of the time series.
  ///
  ///This function calculates the mean and standard deviation.
  // --------------------------------------------------------------------------
  void meanStdDev(const vector<double> &timeSeries_in, const int
      subsequencePos_in, double &mean_out, double &stdDev_out);

  // --------------------------------------------------------------------------
  ///\brief Calculates the similarity of two time series.
  ///
  ///\param [in] &timeSeriesOne_in Hands over the time series.
  ///\param [in] &subsequenceOne_in Hands over the first subsequence in the
  ///time series.
  ///\param [in] meanOne_in Hands over the mean of timeSeriesOne_in.
  ///\param [in] stdDevOne_in Hands over the standard deviation of
  ///timeSeriesOne_in.
  ///\param [in] subsequenceTwoPos_in Hands over the position of the second
  ///subsequence in the time series.
  ///\param [in] bestSoFar_in Hands over the best similarity so far.
  ///
  ///\return The similarity of the two z-normalized time series.
  ///
  ///This function calculates the similarity of two time series. Therefore, the
  ///time series are first z-normalized and the Euclidean Distance is
  ///calculated. The return value is the similarity of the two z-normalized
  ///time series.
  // --------------------------------------------------------------------------
  double similarity(const vector<double> &timeSeries_in, const vector<double>
      &subsequenceOne_in, const double meanOne_in, const double stdDevOne_in,
      const int subsequenceTwoPos_in, const double bestSoFar_in);

  // --------------------------------------------------------------------------
  ///\brief Calculates a motif set subsequence.
  ///
  ///\param [out] subsequence_out Hands over the calculated raw subsequence.
  ///\param [in] type_in Hands over the motif set type.
  ///\param [in] height_in Hands over the height of the subsequence.
  ///
  ///This function calculates a motif set subsequence according to the motif
  ///set type, height and window size. Attention! The range is random and the
  ///subsequence is not added to the synthetic time series.
  // --------------------------------------------------------------------------
  void calculateSubsequence(vector<double> &subsequence_out, int type_in,
      double height_in);

  // --------------------------------------------------------------------------
  ///\brief Search for unintentional matches in the time series with new
  ///injected subsequence.
  ///
  ///\param [in] &timeSereis_in Hands over the time series.
  ///\param [in] &subsequencePositions_in Hands over the position of the new
  ///subsequence.
  ///\param [in] &means_in Hands over the rolling means.
  ///\param [in] &stds_in Hands over the rolling standard deviations.
  ///\param [in] similarity_in Hands over the similarity to break.
  ///
  ///\return true if there exists an unintentional subsequence match with the
  ///new subequence.
  ///
  ///This function searches for unintentinal matches of subsequences in the
  ///time series with the new injected subsequence.
  // --------------------------------------------------------------------------
  bool searchForUnintentionalMatches(const vector<double> &timeSeries_in, const
      vector<int> &motifPositions_in, const vector<double> &means_in, const
      vector<double> &stds_in, double similarity_in);

  // --------------------------------------------------------------------------
  ///\brief Checks if there is a larger motif set.
  ///
  ///\param [in] &timeSereis_in Hands over the time series.
  ///\param [in] &subsequencePositions_in Hands over the position of the new
  ///subsequence.
  ///\param [in] &means_in Hands over the rolling means.
  ///\param [in] &stds_in Hands over the rolling standard deviations.
  ///\param [in] range_in Hands over the motif set range.
  ///
  ///\return true if there exists an larger motif set.
  ///
  ///This function checks if there is a larger motif set in the set of motif
  ///set positions including all overlapping subsequences in the time series.
  // --------------------------------------------------------------------------
  bool checkIfThereIsALargerMotifSet(const vector<double> &timeSeries_in, const
      vector<int> &motifPositions_in, const vector<double> &means_in, const
      vector<double> &stds_in, double range_in);

public:

  // --------------------------------------------------------------------------
  ///\brief The constructor initializes the TSGenerator.
  ///
  ///\param [in] length_in Hands over the time series length.
  ///\param [in] window_in Hands over the window size.
  ///\param [in] delta_in Hands over the maximum difference between two
  ///consecutive values in the time series.
  ///\param [in] noise_in Hands over the noise option, i.e., +-noise_in 2 will
  ///be added to each value of the time series.
  ///\param [in] type_in Hands over the motif type, i.e., the motif shape.
  ///\param [in] size_in Hands over the number of subsequences non-self matched
  ///by the motif.
  ///\param [in] height_in Hands over the maximum difference between two values
  ///of the motif.
  ///
  ///The constructor checks whether a true random engine is available and
  ///stores the result in the trueRandomEngineAvailable variable.
  // --------------------------------------------------------------------------
  TSGenerator(int length_in, int window_in, double delta_in, double noise_in,
      int type_in, int size_in, double height_in);

  // --------------------------------------------------------------------------
  ///\brief The constructor initializes the TSGenerator.
  ///
  ///\param [in] length_in Hands over the time series length.
  ///\param [in] window_in Hands over the window size.
  ///\param [in] delta_in Hands over the maximum difference between two
  ///consecutive values in the time series.
  ///\param [in] noise_in Hands over the noise option, i.e., +-noise_in 2 will
  ///be added to each value of the time series.
  ///\param [in] type_in Hands over the motif type, i.e., the motif shape.
  ///\param [in] size_in Hands over the number of subsequences non-self matched
  ///by the motif.
  ///\param [in] height_in Hands over the maximum difference between two values
  ///of the motif.
  ///
  ///The constructor checks whether a true random engine is available and
  ///stores the result in the trueRandomEngineAvailable variable.
  // --------------------------------------------------------------------------
  TSGenerator(int length_in, int window_in, double delta_in, double noise_in,
      string type_in, int size_in, double height_in);

  // --------------------------------------------------------------------------
  ///\brief Frees the memory allocated by the TSGenerator.
  ///
  ///The destructor does actually nothing.
  // --------------------------------------------------------------------------
  ~TSGenerator();

  // --------------------------------------------------------------------------
  ///\brief Generates a time series with defined time series motif sets.
  ///
  ///\param [out] &timeSeries_in Hands over the time series.
  ///\param [out] &d_in Hands over the range of each the motif sets.
  ///\param [out] &window_out Hands over the window size of each motif set.
  ///\param [out] &motifPositions_in Hands over the positions of each motif
  ///set.
  ///
  ///First, this function checks all variables and pointers for validity. Some
  ///random time series values are forwarded. Afterwards, the time series is
  ///generated and the time series and the motif positions are written into
  ///a file. A not defined motif tag is treated as a random motif. Therfore,
  ///a ranodm motif type is chosen. The default motif type is the random motif.
  // --------------------------------------------------------------------------
  void run(vector<double> &timeSeries_out, vector<double> &d_out,
      vector<int> &window_out, vector<vector<int>> &motifPositions_out);
};

#endif
