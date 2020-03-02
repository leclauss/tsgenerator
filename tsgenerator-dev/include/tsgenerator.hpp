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

#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>
#include <chrono>
#include <cmath>
#include <cfloat>
#include <limits>
#include <iterator>
#include <tsgtypes.hpp>
#include <motifsetcollection.hpp>
#include <freepositions.hpp>
#include <basets.hpp>
#include <tpm.hpp>
#include <tsm.hpp>


namespace tsg
{

  ///\brief This list contains all methods for generating base time series.
  ///
  ///This variable stores the names of all mehtods available in this
  ///implementation. To add a method one has to declare the function in the
  ///basets.hpp file and define the function in the basets.cpp file as well.
  ///Also, one has to modify the base time series method chooser in the file
  ///tsgenerator.cpp.
  const par methods{
    "simpleRandomWalk",
    "realRandomWalk",
    "normalRandomWalk",
    "linearRandomWalk",
    "boundedSimpleRandomWalk",
    "boundedRealRandomWalk",
    "boundedNormalRandomWalk",
    "boundedLinearRandomWalk",
    "uniformRandom",
    "normalRandom",
    "piecewiseLinearRandom",
    "splineRepeated"
  };

  ///\brief This list contains all motif types.
  ///
  ///This variable stores the names of all motifs available in this
  ///implementation. To add a motif one has to declare the motif in the
  ///motifcollection.hpp file and define the motif in the motifcollection.cpp
  ///file as well. Also, one has to modify the addMotif() function in the file
  ///tsgenerator.cpp.
  const par motifTypes{
    "box",
    "triangle",
    "semicircle",
    "trapezoid",
    "positiveflank",
    "negativeflank",
    "sine",
    "cosine"
  };

  ///\brief This list contains all available motif generation types.
  ///
  ///This variable stores the names of all motif generation types currently
  ///available.
  const par gens{
    "pair motif",
    "set motif",
    "latent motif"
  };


  ///\brief This class represents the TSGenerator.
  ///
  ///The TSGenerator consists of setter functions and a run function as well as
  ///various variables. After configuring the TSGenerator with the setter
  ///functions one starts the time series generation by calling the run
  ///function.
  class TSGenerator {

  protected:

    ///\brief This variable contains the base value for the time series values.
    ///
    ///This variable stores the base value of the time series values. The time
    ///series is generated around this value.
    double baseValue = 0.0;

    ///\brief This variable contains the length of the time series.
    ///
    ///The variable stores the length of the time series.
    int length = -1;

    ///\brief This variable the window size.
    ///
    ///The window size is used to generate the motifs such that the window size
    ///is the motif length.
    int window = -1;

    ///\brief This variable contains the maximum difference of two consecutive
    ///values.
    ///
    ///This variable stores the maximum absolute diffrenece of two consecutive
    ///values in the time series excluding noise and motifs.
    double delta = 1.0;

    ///\brief This variable contains the noise option.
    ///
    ///This variable stores the noise option that will be used to add noise to
    ///the time series. The resulting noise is a value between -noise / 2 and
    ///noise / 2.
    double noise = 0.1;

    ///\brief This variable contains the motif type.
    ///
    ///This variable stores the motif type, i.e., the shape of the motif.
    int type = 0;

    ///\brief This variable contains the motif size.
    ///
    ///This variable stores the motif size selected, i.e., the number of
    ///subsequences non-self matched by the motif.
    int size = 3;

    ///\brief This variable contains the motif height.
    ///
    ///This variable stores the motif height, i.e., the maximum difference
    ///between two motif values.
    double height = 10.0;

    ///\brief This variable contains the step size value.
    ///
    ///This variable stores the step size of random walk generation.
    double step = 1.0;

    ///\brief This variable contains the number of steps.
    ///
    ///This variable stores the number of values in a repeating section of the
    ///time series.
    int times = 2;

    ///\brief This variable contains the base time series generation method.
    ///
    ///This variable stores the method for the generation of the base times
    ///series.
    int method = 5;

    ///\brief This variable contains the time series maximum absolute value.
    ///
    ///This variable stores the maximum absolute value a base time serise may
    ///have.
    double maxi = 20.0;

    ///\brief This variable contains the generator type.
    ///
    ///This variable stores the generator type. Currently there are the box
    ///motif, set motif and latent motif generator types available for time
    ///series generation.
    int gen = 2;

    ///\brief This variable contains the free positions.
    ///
    ///This variable stores the positions of free subsequences in the time
    ///series between the motif sets subsequences.
    FreePositions freePositions;

    ///\brief This variable stores the random engine.
    ///
    ///The pseudo random engine generates random numbers with the Mersene
    ///Twister 19937 generator.
    std::mt19937 randomEngine;

    ///\brief This variable contains a base time series generator.
    ///
    ///This variable stores a collection of methods to generate a basic time
    ///series without injected motifs..
    BaseTS baseTS;

    ///\brief Running sum.
    ///
    ///This variable stores the running sum of the time series.
    rseq sums;

    ///\brief Running sum of squares.
    ///
    ///This variable stores the running sum of squares of the time series.
    rseq sumSquares;


    ///\brief Calculates the running sum and sum of squares of a sequence.
    ///
    ///\param [in] &sequence_in Hands over the sequence.
    ///
    ///This function computes the running sum and sum of squares of a time
    ///series.
    void calcRunnings(const rseq &sequence_in);

    ///\brief Update the running sum and sum of squares of a sequence.
    ///
    ///\param [in] &sequence_in Hands over the sequence.
    ///\param [in] pos_in Hands over the position of the injected subsequence.
    ///
    ///This function updates the running sum and sum of squares of a sequence at
    ///a specific location.
    void updateRunnings(const rseq &timeSeries_in, const int pos_in);

    ///\brief Computes the similarity of two subsequences in a sequence.
    ///
    ///\param [in] &sequence_in Hands over the sequence.
    ///\param [in] pos0_in Hands over the position of the first subsequence.
    ///\param [in] pos1_in Hands over the position of the second subsequence.
    ///\param [in] bestSoFar_in Hands over the best similarity so far.
    ///
    ///\return The similarity of the two z-normalized subsequences.
    ///
    ///This function computes the similarity of two subsequences. Therefore, the
    ///subsequences are first z-normalized and the Euclidean Distance is
    ///computed. The return value is the similarity of the two z-normalized
    ///subsequences.
    double similarity(const rseq &sequence_in, const int pos0_in, const int
        pos1_in, const double bestSoFar_in);

    ///\brief Computes the mean and standard deviation of a sequence.
    ///
    ///\param [in] &sequence_in Hands over the sequence.
    ///\param [out] &mean_out The mean of the sequence.
    ///\param [out] &stdDev_out The standard deviation of the sequence.
    ///
    ///This function calculates the mean and standard deviation of a sequence.
    void meanStdDev(const rseq &sequence_in, double &mean_out, double
        &variance_out);

    ///\brief Computes the similarity of a sequence and a subsequence in the time
    ///series.
    ///
    ///\param [in] &timeSeries_in Hands over the time series.
    ///\param [in] &sequence_in Hands over the sequence.
    ///\param [in] mean_in Hands over the mean of sequence_in.
    ///\param [in] stdDev_in Hands over the standard deviation of sequence_in.
    ///\param [in] pos_in Hands over the position of the subsequence in the time
    ///series.
    ///\param [in] bestSoFar_in Hands over the best similarity so far.
    ///
    ///\return The similarity of the two z-normalized sequence and subsequence.
    ///
    ///This function computes the similarity of a sequence and a subsequence in
    ///the time series. Therefore, the sequence and the subsequence are first
    ///z-normalized and the Euclidean Distance is computed. The return value is
    ///the similarity of the z-normalized sequence and subsequence.
    double similarity(const rseq &timeSeries_in, const rseq &subsequence_in,
        const double mean_in, const double stdDev_in, const int pos_in, const
        double bestSoFar_in);

    ///\brief Computes the similarity of two sequences.
    ///
    ///\param [in] &sequence0_in Hands over the first sequence.
    ///\param [in] &sequence1_in Hands over the second sequence.
    ///\param [in] bestSoFar_in Hands over the best similarity so far.
    ///
    ///\return The similarity of the two z-normalized sequences.
    ///
    ///This function computes the similarity of two sequences. Therefore, the
    ///sequences are first z-normalized and the Euclidean Distance is computed.
    ///The return value is the similarity of the two z-normalized sequneces.
    double similarity(const rseq &sequence0_in, const rseq &sequence1_in, const
        double bestSoFar_in);

    ///\brief Computes a motif set subsequence.
    ///
    ///\param [out] subsequence_out Hands over the calculated raw subsequence.
    ///\param [in] type_in Hands over the motif set type.
    ///\param [in] height_in Hands over the height of the subsequence.
    ///
    ///This function computes a motif set subsequence according to the motif
    ///type, height and window size. Attention! The range is random and the
    ///subsequence is not added to the synthetic time series.
    void calculateSubsequence(rseq &subsequence_out, int type_in, double
        height_in);

    ///\brief Calculates a base time series.
    ///
    ///\param [out] timeSeries_out Hands over the computed time series.
    ///
    ///This function computes a base times series according to the length, delta,
    ///maxi and noise values.
    void generateBaseTimeSeries(rseq &timeSeries_out);

    ///\brief Check if there is a better pair motif in the time series.
    ///
    ///\param [in] &timeSereis_in Hands over the time series.
    ///\param [in] &subsequencePositions_in Hands over the pair motif
    ///subsequence positions.
    ///\param [in] similarity_in Hands over the similarity to break.
    ///
    ///\return true if there exists a subsequence pair in the time series with
    ///smaller distance.
    ///
    ///This function checks if the subsequences overlapping the second injected
    ///pair motif sequence have another subsequence within range similarity_in.
    bool smallerDistance(const rseq &timeSeries_in, const iseq
        &motifPositions_in, const double similarity_in);

    ///\brief Checks if there is a larger motif set.
    ///
    ///\param [in] &timeSereis_in Hands over the time series.
    ///\param [in] &subsequencePositions_in Hands over the position of the new
    ///subsequence.
    ///\param [in] range_in Hands over the motif set range.
    ///
    ///\return true if there exists an larger motif set.
    ///
    ///This function checks if there is a larger motif set with given position.
    bool largerMotifSet(const rseq &timeSeries_in, const iseq
        &motifPositions_in, double range_in);

    ///\brief Checks if there is a larger motif set.
    ///
    ///\param [in] &timeSereis_in Hands over the time series.
    ///\param [in] &subsequencePositions_in Hands over the position of the new
    ///subsequence.
    ///\param [in] size_in Hands over the size of the injected motif.
    ///\param [in] range_in Hands over the motif set range.
    ///
    ///\return Size of the first largest set motif for the subsequence
    ///motifPositions_in.back().
    ///
    ///This function computes the largest set motif size for a given
    ///subsequence.
    int largerMotifSet(const rseq &timeSeries_in, const iseq motifPositions_in,
        const int size_in, const double range_in);

    ///\brief Injects a pair motif into the time series.
    ///
    ///\param [out] &timeSeries_out Hands over the time series.
    ///\param [out] &motif_out Hands over the motif sequences.
    ///
    ///tbd
    void injectPairMotif(rseq &timeSeries_out, rseqs &motif_out);

    ///\brief Injects a set motif into the time series.
    ///
    ///\param [out] &timeSeries_out Hands over the time series.
    ///\param [out] &motif_out Hands over the motif sequences.
    ///\param [out] &d_out Hands over the range of each the motif sets.
    ///\param [out] &pos_out Hands over the positions of each motif
    ///set.
    ///
    ///tbd
    void injectSetMotif(rseq &timeSeries_out, rseqs &motif_out, rseq &d_out,
        iseqs &pos_out);

    ///\brief Injects a latent motif into the time series.
    ///
    ///\param [out] &timeSeries_out Hands over the time series.
    ///\param [out] &motif_out Hands over the motif sequences.
    ///\param [out] &d_out Hands over the range of each the motif sets.
    ///\param [out] &pos_out Hands over the positions of each motif set.
    ///
    ///First the distance d of a most similar subsequence pair is computed and
    ///a random motif center sequence with similarity 3 / 2 d or more to all
    ///other subsequences is generated. Than sequences are injected within
    ///range d / 2 to any other subsequence. Finally the function checks wether
    ///there is a larget latent motif.
    void injectLatentMotif(rseq &timeSeries_out, rseqs &motif_out, rseq &d_out,
        iseqs &pos_out);

  public:

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
    ///\param [in] times_in Hands over the times of steps.
    ///\param [in] method_in Hands over the method for base time series
    ///generation.
    ///\param [in] max_in Hands over the maximum absolute value of the time
    ///series.
    ///\param [in] gen_in Hands over the motif generation type.
    ///
    ///The constructor checks whether a true random engine is available and
    ///stores the result in the trueRandomEngineAvailable variable.
    TSGenerator(const int length_in, const int window_in, const double
        delta_in, const double noise_in, const int type_in, const int size_in,
        const double height_in, const double step_in = 1.0, const int times_in
        = 3, const int method_in = 5, const double maxi_in = 100.0, const int
        gen_in = 2);

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
    ///\param [in] times_in Hands over the times of steps.
    ///\param [in] method_in Hands over the method for base time series
    ///generation.
    ///\param [in] maxi_in Hands over the maximum absolute value of the time
    ///series.
    ///\param [in] gen_in Hands over the motif generation type.
    ///
    ///The constructor checks whether a true random engine is available and
    ///stores the result in the trueRandomEngineAvailable variable.
    TSGenerator(const int length_in, const int window_in, const double
        delta_in, const double noise_in, const word type_in, const int size_in,
        const double height_in, const double step_in = 1.0, const int times_in
        = 3, const word method_in = "boundedNormalRandomWalk", const double
        maxi_in = 100.0, const word gen_in = "latent motif");

    ///\brief Frees the memory allocated by the TSGenerator.
    ///
    ///The destructor does actually nothing.
    ~TSGenerator();

    ///\brief Generates a time series with defined time series motif sets.
    ///
    ///\param [out] &timeSeries_out Hands over the time series.
    ///\param [out] &motif_out Hands over the motif sequences.
    ///\param [out] &d_out Hands over the range of each the motif sets.
    ///\param [out] &pos_out Hands over the positions of each motif
    ///set.
    ///
    ///First, this function checks all variables and pointers for validity.
    ///Some random time series values are forwarded. Afterwards, the time
    ///series is generated and the time series and the motif positions are
    ///written into a file. A not defined motif tag is treated as a random
    ///motif. Therfore, a ranodm motif type is chosen. The default motif type
    ///is the random motif.
    void run(rseq &timeSeries_out, rseqs &motif_out, rseq &d_out, iseqs
        &pos_out);
  };
}

#endif
