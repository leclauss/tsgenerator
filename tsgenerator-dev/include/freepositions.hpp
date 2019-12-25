///\file freepositions.hpp
///
///\brief File contains the FreePositions class declaration.
///
///This is the header file of the FreePositions. The FreePositions is a data
///structure to store the remaining free positions in the time series for
///subsequences to be injected.

#ifndef FREEPOSITIONS_HPP
#define FREEPOSITIONS_HPP

#include <random>
#include <chrono>
#include <tsgtypes.hpp>


namespace tsg {

  ///\brief This class represents the free positions in a time series.
  ///
  ///The FreePositons consists of setter functions and a run function as well
  ///as various variables. After configuring the TSGenerator with the setter
  ///functions one starts the time series generation by calling the run
  ///function.
  class FreePositions {

  protected:

    ///\brief This variable stores the window size.
    ///
    ///The window size is used to generate the motifs such that the window size
    ///is the motif length.
    int window = -1;

    ///\brief This variable counts the free positions.
    ///
    ///The free count is the number of free positions in the time series.
    int freeCount = -1;

    ///\brief This variable contains the free positions.
    ///
    ///This variable stores the positions of free subsequences in the time
    ///series between the motif sets subsequences.
    intervals freePositions;

    ///\brief This variable stores the random engine.
    ///
    ///The pseudo random engine generates random numbers with the Mersene
    ///Twister 19937 generator.
    std::mt19937 randomEngine;

    ///\brief This variable stores the random generated position.
    ///
    ///This is always the most current random position calculated with
    ///calculateRandomPosition().
    int randomPosition = -1;


    ///\brief Calculates a random number.
    ///
    ///\return A random number.
    ///
    ///This function calculates a random number in [0, freeCount].
    int randomNumber();

  public:

    ///\brief The constructor initializes the FreePositions.
    ///
    ///\param [in] length_in Length of the time series.
    ///\param [in] window_in Window size of the subsequences to be injected.
    ///
    ///The constructor checks whether a true random engine is available,
    ///stores the result in the trueRandomEngineAvailable variable and sets the
    ///length of the time series and the window size of the subsequences to be
    ///injected.
    FreePositions(const int length_in, const int window_in);

    ///\brief Frees the memory allocated by FreePositions.
    ///
    ///The destructor does actually nothing.
    ~FreePositions();

    ///\brief Calculates a random free position in the time series.
    ///
    ///\return The random position in the time series.
    ///
    ///This function calculates a random free position in the time series
    ///according to the window size.
    int calculateRandomPosition();

    ///\brief Removes position from free ones.
    ///
    ///This function removes the position in randomPosition.
    void removePosition();
  };
}

#endif
