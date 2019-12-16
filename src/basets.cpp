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
  normal_distribution<double> distributionNoise(0.0, noise_in / 2.0 > 0.0
      ? noise_in / 2.0 : numeric_limits<double>::min());
  uniform_int_distribution<int> distribution(0, 1);

  //add the first value
  double value = start_in;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = delta_in * distribution(randomEngine) - 2 * delta_in;

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    timeSeries_out.push_back(timeSeries_out[i - 1] + value);
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
  normal_distribution<double> distributionNoise(0.0, noise_in / 2.0 > 0.0
      ? noise_in / 2.0 : numeric_limits<double>::min());
  uniform_real_distribution<double> distribution(-delta_in, delta_in);

  //add the first value
  double value = start_in;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = distribution(randomEngine);

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    timeSeries_out.push_back(timeSeries_out[i - 1] + value);
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
  normal_distribution<double> distributionNoise(0.0, noise_in / 2.0 > 0.0
      ? noise_in / 2.0 : numeric_limits<double>::min());
  normal_distribution<double> distribution(0.0, delta_in > 0.0
      ? delta_in : numeric_limits<double>::min());

  //add the first value
  double value = start_in;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = 0.0;

    if (delta_in > 0.0)
      value += distribution(randomEngine);

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    timeSeries_out.push_back(timeSeries_out[i - 1] + value);
  }
}

void BaseTS::simpleRandomWalk(vector<double> &timeSeries_out, int length_in,
    double start_in, double delta_in, double maxi_in, double noise_in) {

  if (maxi_in - delta_in < 0.0) {

    cerr << "ERROR: Border " << maxi_in << " is to tied for delta " <<
      delta_in << endl;
    throw(EXIT_FAILURE);
  }

  if (!(timeSeries_out.empty())) {

    timeSeries_out.clear();
    timeSeries_out.resize(0);
  }

  if (length_in < 1)
    return;

  //initialize the distributions for the random walk
  normal_distribution<double> distributionNoise(0.0, noise_in / 2.0 > 0.0
      ? noise_in / 2.0 : numeric_limits<double>::min());
  uniform_int_distribution<int> distribution(0, 1);

  //add the first value
  double value = start_in;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = delta_in * distribution(randomEngine) - 2 * delta_in;

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    //check if border was crossed
    if (timeSeries_out[i - 1] + value < start_in - maxi_in ||
        timeSeries_out[i - 1] + value > start_in + maxi_in)
      value = -value;

    timeSeries_out.push_back(timeSeries_out[i - 1] + value);
  }
}

void BaseTS::realRandomWalk(vector<double> &timeSeries_out, int length_in,
    double start_in, double delta_in, double maxi_in, double noise_in) {

  if (maxi_in - delta_in < 0.0) {

    cerr << "ERROR: Border " << maxi_in << " is to tied for delta " <<
      delta_in << endl;
    throw(EXIT_FAILURE);
  }

  if (!(timeSeries_out.empty())) {

    timeSeries_out.clear();
    timeSeries_out.resize(0);
  }

  if (length_in < 1)
    return;

  //initialize the distributions for the random walk
  normal_distribution<double> distributionNoise(0.0, noise_in / 2.0 > 0.0
      ? noise_in / 2.0 : numeric_limits<double>::min());
  uniform_real_distribution<double> distribution(-delta_in, delta_in);

  //add the first value
  double value = start_in;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = distribution(randomEngine);

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    //check if border was crossed
    if (timeSeries_out[i - 1] + value < start_in - maxi_in ||
        timeSeries_out[i - 1] + value > start_in + maxi_in)
      value = -value;

    timeSeries_out.push_back(timeSeries_out[i - 1] + value);
  }
}

void BaseTS::normalRandomWalk(vector<double> &timeSeries_out, int length_in,
    double start_in, double delta_in, double maxi_in, double noise_in) {

  if (maxi_in - delta_in < 0.0) {

    cerr << "ERROR: Border " << maxi_in << " is to tied for delta " <<
      delta_in << endl;
    throw(EXIT_FAILURE);
  }

  if (!(timeSeries_out.empty())) {

    timeSeries_out.clear();
    timeSeries_out.resize(0);
  }

  if (length_in < 1)
    return;

  //initialize the distributions for noise and the random walk
  normal_distribution<double> distributionNoise(0.0, noise_in / 2.0 > 0.0
      ? noise_in / 2.0 : numeric_limits<double>::min());
  normal_distribution<double> distribution(0.0, delta_in > 0.0
      ? delta_in : numeric_limits<double>::min());

  //add the first value
  double value = start_in;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = 0.0;

    if (delta_in > 0.0)
      value += distribution(randomEngine);

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    //check if border was crossed
    if (timeSeries_out[i - 1] + value < start_in - maxi_in ||
        timeSeries_out[i - 1] + value > start_in + maxi_in)
      value = -value;

    timeSeries_out.push_back(timeSeries_out[i - 1] + value);
  }
}

void BaseTS::uniformRandom(vector<double> &timeSeries_out, int length_in,
    double start_in, double delta_in, double noise_in) {

  if (!(timeSeries_out.empty())) {

    timeSeries_out.clear();
    timeSeries_out.resize(0);
  }

  if (length_in < 1)
    return;

  //initialize the distributions for noise and the random walk
  normal_distribution<double> distributionNoise(0.0, noise_in / 2.0 > 0.0
      ? noise_in / 2.0 : numeric_limits<double>::min());
  uniform_real_distribution<double> distribution(-delta_in / 2.0, delta_in
      / 2.0);

  //add the first value
  double value = start_in;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = distribution(randomEngine);

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    timeSeries_out.push_back(value);
  }
}

void BaseTS::normalRandom(vector<double> &timeSeries_out, int length_in,
    double start_in, double delta_in, double noise_in) {

  if (!(timeSeries_out.empty())) {

    timeSeries_out.clear();
    timeSeries_out.resize(0);
  }

  if (length_in < 1)
    return;

  //initialize the distributions for noise and the random walk
  normal_distribution<double> distributionNoise(0.0, noise_in / 2.0 > 0.0
      ? noise_in / 2.0 : numeric_limits<double>::min());
  normal_distribution<double> distribution(0.0, delta_in / 2.0 > 0.0
      ? delta_in / 2.0 : numeric_limits<double>::min());

  //add the first value
  double value;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = start_in;

    if (delta_in / 2.0 > 0.0)
      value += distribution(randomEngine);

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    timeSeries_out.push_back(value);
  }
}

void BaseTS::piecewiseLinearRandom(vector<double> &timeSeries_out, int
    length_in, double start_in, double delta_in, double noise_in) {

  if (!(timeSeries_out.empty())) {

    timeSeries_out.clear();
    timeSeries_out.resize(0);
  }

  if (length_in < 1)
    return;

  //initialize the distributions for noise and the random walk
  normal_distribution<double> distributionNoise(0.0, noise_in / 2.0 > 0.0
      ? noise_in / 2.0 : numeric_limits<double>::min());

  vector<double> intervals {-delta_in / 2.0, -delta_in / 4.0, 0.0, delta_in
    / 4.0, delta_in / 2.0};
  vector<double> weights {0.0, 10.0, 0.0, 10.0, 0.0};
  piecewise_linear_distribution<double> distribution(intervals.begin(),
      intervals.end(), weights.begin());

  //add the first value
  double value = start_in;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = distribution(randomEngine);

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    timeSeries_out.push_back(value);
  }
}
