///\file basets.cpp
///
///\brief File contains the BaseTS class definition.
///
///This is the source file of the BaseTS. The BaseTS is a collection of methods
///for the generation of random time series without injected motifs.


#include <basets.hpp>


BaseTS::BaseTS()
  : randomEngine(random_device().entropy()
    ? random_device()()
    : chrono::system_clock::now().time_since_epoch().count()) { }

BaseTS::~BaseTS() { }

void BaseTS::simpleRandomWalk(vector<double> &timeSeries_out, int length_in,
    double start_in, double delta_in, double noise_in) {

  if (!(timeSeries_out.empty())) {

    timeSeries_out.clear();
    timeSeries_out.resize(0);
  }

  if (length_in < 1)
    return;

  //initialize the distributions for the random walk
  normal_distribution<double> distributionNoise(0.0, noise_in / 2.0);
  uniform_int_distribution<int> distributionRandomWalk(0, 1);

  //add the first value
  double value = start_in;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int itr = 1; itr < length_in; itr++) {

    value += delta_in * distributionRandomWalk(randomEngine) - 2 * delta_in
      + distributionNoise(randomEngine);

    timeSeries_out.push_back(value);
  }
}

void BaseTS::realRandomWalk(vector<double> &timeSeries_out, int length_in,
    double start_in, double delta_in, double noise_in) {

  if (!(timeSeries_out.empty())) {

    timeSeries_out.clear();
    timeSeries_out.resize(0);
  }

  if (length_in < 1)
    return;

  //initialize the distributions for the random walk
  normal_distribution<double> distributionNoise(0.0, delta_in / 2.0);
  uniform_real_distribution<double> distributionRandomWalk(-delta_in,
      delta_in);

  //add the first value
  double value = start_in;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int itr = 1; itr < length_in; itr++) {

    value += distributionRandomWalk(randomEngine)
      + distributionNoise(randomEngine);

    timeSeries_out.push_back(value);
  }
}

void BaseTS::normalRandomWalk(vector<double> &timeSeries_out, int length_in,
    double start_in, double delta_in, double noise_in) {

  if (!(timeSeries_out.empty())) {

    timeSeries_out.clear();
    timeSeries_out.resize(0);
  }

  if (length_in < 1)
    return;

  //initialize the distributions for noise and the random walk
  normal_distribution<double> distributionNoise(0.0, delta_in / 2.0);
  normal_distribution<double> distributionRandomWalk(0.0, delta_in);

  //add the first value
  double value = start_in;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int itr = 1; itr < length_in; itr++) {

    value += distributionRandomWalk(randomEngine)
      + distributionNoise(randomEngine);

    timeSeries_out.push_back(value);
  }
}
