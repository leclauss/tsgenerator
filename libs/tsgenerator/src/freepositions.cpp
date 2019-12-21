///\file freepositions.cpp
///
///\brief File contains the FreePositions class definition.
///
///This is the source file of the FreePositions class. The FreePositions stores
///the free positions in a time series for a subsequence to be injected.


#include <freepositions.hpp>


FreePositions::FreePositions(int length_in, int window_in)
  : window(window_in), randomEngine(random_device().entropy()
      ? random_device()()
      : chrono::system_clock::now().time_since_epoch().count()) {

  //create first large interval
  interval first;
  first.start = 0;
  first.end = length_in - window;
  freePositions.push_back(first);

  //initialize the free positions count
  freeCount = first.end + 1;
}

FreePositions::~FreePositions() {}

int FreePositions::randomNumber() {

  //initialize new distribution
  uniform_int_distribution<int> distribution(0, freeCount - 1);

  return distribution(randomEngine);
}

int FreePositions::calculateRandomPosition() {

  if (freeCount < 1)
    throw(EXIT_FAILURE);

  //calculate the free random position
  randomPosition = randomNumber();

  //index of the interval containing the position
  int iP = 0;

  //offset the position to the first interval
  randomPosition += freePositions[iP].start;

  //decode the random number into a position
  while (randomPosition > freePositions[iP].end) {

    randomPosition += freePositions[iP + 1].start - freePositions[iP].end - 1;
    iP++;
  }

  return randomPosition;
}

void FreePositions::removePosition() {

  //index of the interval containing the position
  int iP = 0;

  //determine the index of the interval containing the position
  while (randomPosition > freePositions[iP].end)
    iP++;

  //split the interval
  interval first;
  first.start = freePositions[iP].start;
  first.end = randomPosition - 2 * window;

  interval second;
  second.start = randomPosition + 2 * window;
  second.end = freePositions[iP].end;

  //remove the old interval & update free count
  freePositions.erase(freePositions.begin() + iP);
  freeCount--;

  //if there is no position in the second interval
  if (second.end - second.start < 0) {

    //remove all positions right from removed position
    freeCount -= second.end - randomPosition;
  }
  else {

    //else insert new interval and update position count
    freeCount -= 2 * window - 1;
    freePositions.insert(freePositions.begin() + iP, second);
  }

  //analogous for left
  if (first.end - first.start < 0) {

    freeCount -= randomPosition - first.start;
  }
  else {

    freeCount -= 2 * window - 1;
    freePositions.insert(freePositions.begin() + iP, first);
  }
}
