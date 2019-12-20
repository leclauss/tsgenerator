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

void BaseTS::cubicSpline(const vector<double> & x_in, const vector<double>
    &y_in, vector<double> &a_out, vector<double> &b_out, vector<double> &c_out,
    vector<double> &d_out) {

  if (x_in.size() != y_in.size()) {

    cerr << "ERROR: x and y values have different length!" << endl;
    throw(EXIT_FAILURE);
  }

  int length = x_in.size();

  if (!(a_out.empty())) {

    a_out.clear();
    a_out.resize(0);
  }

  for (int i = 0; i < length; i++)
    a_out.push_back(y_in[i]);

  if (!(b_out.empty()))
    b_out.clear();

  b_out.resize(length - 1);

  if (!(d_out.empty()))
    d_out.clear();

  d_out.resize(length - 1);

  vector<double> h;

  for (int i = 0; i < length - 1; i++)
    h.push_back(x_in[i + 1] - x_in[i]);

  vector<double> alpha;

  //dummy
  alpha.push_back(0);

  for (int i = 1; i < length - 1; i++)
    alpha.push_back(3.0 * (a_out[i + 1] - a_out[i]) / h[i] - 3.0 * (a_out[i]
          - a_out[i - 1]) / h[i - 1]);

  if (!(c_out.empty()))
    c_out.clear();

  c_out.resize(length);

  vector<double> l(length), mu(length), z(length);

  l[0] = 1.0;
  mu[0] = 0.0;
  z[0] = 0.0;

  for (int i = 1; i < length - 1; i++) {

    l[i] = 2.0 * (x_in[i + 1] - x_in[i - 1]) - h[i - 1] * mu[i - 1];
    mu[i] = h[i] / l[i];
    z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
  }

  l[length - 1] = 1.0;
  z[length - 1] = 0.0;
  c_out[length - 1] = 0.0;

  for (int i = length - 2; i >= 0; i--) {

    c_out[i] = z[i] - mu[i] * c_out[i + 1];
    b_out[i] = (a_out[i + 1] - a_out[i]) / h[i] - h[i] * (c_out[i + 1] + 2.0
        * c_out[i]) / 3.0;
    d_out[i] = (c_out[i + 1] - c_out[i]) / (3.0 * h[i]);
  }

  a_out.pop_back();
  c_out.pop_back();
}

void BaseTS::simpleRandomWalk(vector<double> &timeSeries_out, int length_in,
    double delta_in, double noise_in) {

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
  double value = 0.0;

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
    double delta_in, double noise_in) {

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
  double value = 0.0;

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
    double delta_in, double noise_in) {

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
  double value = 0.0;

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

void BaseTS::linearRandomWalk(vector<double> &timeSeries_out, int length_in,
    double delta_in, double step_in, double noise_in) {

  if (!(timeSeries_out.empty())) {

    timeSeries_out.clear();
    timeSeries_out.resize(0);
  }

  if (length_in < 1)
    return;

  int step = abs(step_in);

  //initialize the distributions for noise and the random walk
  normal_distribution<double> distributionNoise(0.0, noise_in / 2.0 > 0.0
      ? noise_in / 2.0 : numeric_limits<double>::min());
  normal_distribution<double> distribution(0.0, delta_in > 0.0
      ? delta_in : numeric_limits<double>::min());
  poisson_distribution<int> distributionStep((double)step);

  //add the first value
  double value = 0.0;

  timeSeries_out.push_back(value);

  //compute first step
  step = distributionStep(randomEngine);

  //store the last step and value
  double rStep = 1.0 / (double)step;
  double lValue = value;

  //get next value
  double nValue = lValue + distribution(randomEngine);

  //store the value difference
  double diff = nValue - lValue;

  //step iterator
  int is = 1;

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    //check if we need to linear approximate
    if (is < step) {

      value = lValue + diff * rStep * (double)is;

      //add noise
      if (noise_in / 2.0 > 0.0)
        value += distributionNoise(randomEngine);

      //update step iterator
      is++;
    }
    else { //check if we need a new step

      value = nValue;
      lValue = nValue;

      //compute next value
      if (delta_in > 0.0)
        nValue += distribution(randomEngine);

      //add noise
      if (noise_in / 2.0 > 0.0)
        value += distributionNoise(randomEngine);

      //get next step
      step = distributionStep(randomEngine);
      rStep = 1.0 / step;

      //update diff
      diff = nValue - lValue;

      //reset step iterator
      is = 1;
    }

    //add the value
    timeSeries_out.push_back(value);
  }
}

void BaseTS::simpleRandomWalk(vector<double> &timeSeries_out, int length_in,
    double delta_in, double maxi_in, double noise_in) {

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
  double value = 0.0;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = delta_in * distribution(randomEngine) - 2 * delta_in;

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    //check if border is crossed
    if (timeSeries_out[i - 1] + value < -maxi_in)
      value = abs(value);

    if (timeSeries_out[i - 1] + value > maxi_in)
      value = -abs(value);

    timeSeries_out.push_back(timeSeries_out[i - 1] + value);
  }
}

void BaseTS::realRandomWalk(vector<double> &timeSeries_out, int length_in,
    double delta_in, double maxi_in, double noise_in) {

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
  double value = 0.0;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = distribution(randomEngine);

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    //check if border is crossed
    if (timeSeries_out[i - 1] + value < -maxi_in)
      value = abs(value);

    if (timeSeries_out[i - 1] + value > maxi_in)
      value = -abs(value);

    timeSeries_out.push_back(timeSeries_out[i - 1] + value);
  }
}

void BaseTS::normalRandomWalk(vector<double> &timeSeries_out, int length_in,
    double delta_in, double maxi_in, double noise_in) {

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
  double value = 0.0;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = 0.0;

    if (delta_in > 0.0)
      value += distribution(randomEngine);

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    //check if border is crossed
    if (timeSeries_out[i - 1] + value < -maxi_in)
      value = abs(value);

    if (timeSeries_out[i - 1] + value > maxi_in)
      value = -abs(value);

    timeSeries_out.push_back(timeSeries_out[i - 1] + value);
  }
}

void BaseTS::linearRandomWalk(vector<double> &timeSeries_out, int length_in,
    double delta_in, double step_in, double maxi_in, double noise_in) {

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

  int step = abs(step_in);

  //initialize the distributions for noise and the random walk
  normal_distribution<double> distributionNoise(0.0, noise_in / 2.0 > 0.0
      ? noise_in / 2.0 : numeric_limits<double>::min());
  normal_distribution<double> distribution(0.0, delta_in > 0.0
      ? delta_in : numeric_limits<double>::min());
  poisson_distribution<int> distributionStep((double)step);

  //add the first value
  double value = 0.0;

  timeSeries_out.push_back(value);

  //compute first step
  step = distributionStep(randomEngine);

  //store the last step and value
  double rStep = 1.0 / (double)step;
  double lValue = value;

  //get next value
  value = distribution(randomEngine);

  //check if border is crossed
  if (lValue + value < -maxi_in)
    value = abs(value);

  if (lValue + value > maxi_in)
    value = -abs(value);

  double nValue = lValue + value;

  //store the value difference and create a noise variable
  double diff = nValue - lValue;

  //step iterator
  int is = 1;

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    //check if we need to linear approximate
    if (is < step) {

      value = lValue + diff * rStep * (double)is;

      //add noise
      if (noise_in / 2.0 > 0.0)
        value += distributionNoise(randomEngine);

      //update step iterator
      is++;
    }
    else { //check if we need a new step

      lValue = nValue;

      //compute next value
      if (delta_in > 0.0)
        value = distribution(randomEngine);

      //check if border is crossed
      if (lValue + value < -maxi_in)
        value = abs(value);

      if (lValue + value > maxi_in)
        value = -abs(value);

      nValue += value;
      value = lValue;

      //add noise
      if (noise_in / 2.0 > 0.0)
        value += distributionNoise(randomEngine);

      //get next step
      step = distributionStep(randomEngine);
      rStep = 1.0 / step;

      //update diff
      diff = nValue - lValue;

      //reset step iterator
      is = 1;
    }

    //add the value
    timeSeries_out.push_back(value);
  }
}

void BaseTS::uniformRandom(vector<double> &timeSeries_out, int length_in,
    double delta_in, double noise_in) {

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
  double value = 0.0;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = distribution(randomEngine);

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    timeSeries_out.push_back(value);
  }
}

void BaseTS::normalRandom(vector<double> &timeSeries_out, int length_in, double
    delta_in, double noise_in) {

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
  double value = 0.0;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = 0.0;

    if (delta_in / 2.0 > 0.0)
      value += distribution(randomEngine);

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    timeSeries_out.push_back(value);
  }
}

void BaseTS::piecewiseLinearRandom(vector<double> &timeSeries_out, int
    length_in, double delta_in, double noise_in) {

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
  double value = 0.0;

  timeSeries_out.push_back(value);

  //add the remaining values
  for (int i = 1; i < length_in; i++) {

    value = distribution(randomEngine);

    if (noise_in / 2.0 > 0.0)
      value += distributionNoise(randomEngine);

    timeSeries_out.push_back(value);
  }
}

void BaseTS::splineRepeated(vector<double> &timeSeries_out, int length_in,
    double delta_in, int step_in, int times_in, double noise_in) {

  if (!(timeSeries_out.empty())) {

    timeSeries_out.clear();
    timeSeries_out.resize(0);
  }

  if (length_in < 1)
    return;

  int step = abs(step_in);
  int times = abs(times_in);

  //initialize the distributions for noise and the random walk
  normal_distribution<double> distributionNoise(0.0, noise_in / 2.0 > 0.0
      ? noise_in / 2.0 : numeric_limits<double>::min());
  uniform_real_distribution<double> distribution(delta_in / 2.0 > 0.0
      ? -delta_in / 2.0 : -numeric_limits<double>::min(), delta_in / 2.0 > 0.0
      ? delta_in / 2.0 : numeric_limits<double>::min());
  uniform_int_distribution<int> distributionStep(1, step < 1 ? 1 : step);

  //generate the repeating sequence
  vector<double> x;
  vector<double> y;

  x.push_back(0);
  y.push_back(distribution(randomEngine));

  for (int i = 1; i < times; i++) {

    x.push_back(x[i - 1] + distributionStep(randomEngine));
    y.push_back(distribution(randomEngine));
  }

  x.push_back(x[times - 1] + distributionStep(randomEngine));
  y.push_back(y[0]);

  for (int i = 1; i < times; i++) {

    x.push_back(x[times] + x[i]);
    y.push_back(y[i]);
  }

  vector<double> a;
  vector<double> b;
  vector<double> c;
  vector<double> d;

  cubicSpline(x, y, a, b, c, d);

  int start = times / 2;
  int t = 0;

  //add spline parts until time series is full
  while (t < length_in) {

    //compute time series spline parts
    for (int i = start; i < start + times && t < length_in; i++) {

      for (int j = x[i]; j < x[i + 1] && t < length_in; j++) {

        //get spline value
        timeSeries_out.push_back(a[i] + b[i] * (j - x[i]) + c[i] * pow(j
              - x[i], 2.0) + d[i] * pow(j - x[i], 3.0));

        //add noise
        timeSeries_out[t] += distributionNoise(randomEngine);

        t++;
      }
    }
  }
}
