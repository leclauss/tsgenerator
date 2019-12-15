///\file tsgenerator.cpp
///
///\brief File contains the TSGenerator class definition.
///
///This is the source file of the TSGenerator. The TSGenerator consists of
///setter functions and the () function. After configuring the TSGenerator
///with the setter functions one starts the time series generation by calling
///the rue) function.


#include <tsgenerator.hpp>


TSGenerator::TSGenerator(int length_in, int window_in, double delta_in,
    double noise_in, int type_in, int size_in, double height_in, double
    start_in, int method_in, double maxi_in)
  : length(abs(length_in)), window(window_in), delta(delta_in),
  noise(noise_in), type(abs(type_in)), size(abs(size_in)), height(height_in),
  start(start_in), method(method_in), maxi(abs(maxi_in)), freePositions(length,
      window), randomEngine(random_device().entropy()
    ? random_device()()
    : chrono::system_clock::now().time_since_epoch().count()) {

  //check if type exists
  if (type >= (int) motifTypes.size()) {

    cerr << "ERROR: Wrong motif set type: " << type_in << endl;
    throw(EXIT_FAILURE);
  }

  //check if method exists
  if (method >= (int) methods.size()) {

    cerr << "ERROR: Unknown method: " << method_in << endl;
    throw(EXIT_FAILURE);
  }

  //check if the motif amount is set properly
  if (size < 3) {

    cerr << "ERROR: Wrong motif amount: " << size_in << endl;
    throw(EXIT_FAILURE);
  }

  //check if the time series length can handle all motif subsequences
  if (length < (2 * size - 1) * window) {

    cerr << "ERROR: The sum of the subsequences length is to" <<
      " large to fit into the time series!" << endl;
    throw(EXIT_FAILURE);
  }
}

TSGenerator::TSGenerator(int length_in, int window_in, double delta_in,
    double noise_in, string type_in, int size_in, double height_in, double
    start_in, string method_in, double maxi_in)
  : length(abs(length_in)), window(window_in), delta(delta_in),
  noise(noise_in), size(abs(size_in)), height(height_in), start(start_in),
  maxi(abs(maxi_in)), freePositions(length, window),
  randomEngine(random_device().entropy()
    ? random_device()()
    : chrono::system_clock::now().time_since_epoch().count()) {

  // get the type number
  type = distance(motifTypes.begin(), find(motifTypes.begin(),
        motifTypes.end(), type_in));

  // check if type exists
  if (type >= (int) motifTypes.size()) {

    cerr << "ERROR: Wrong motif set type: " << type_in << endl;
    throw(EXIT_FAILURE);
  }

  // get the method
  method = distance(methods.begin(), find(methods.begin(), methods.end(),
        method_in));

  // check if type exists
  if (method >= (int) methods.size()) {

    cerr << "ERROR: Unknown method: " << method_in << endl;
    throw(EXIT_FAILURE);
  }

  //check if the motif amount is set properly
  if (size < 3) {

    cerr << "ERROR: Wrong motif amount: " << size_in << endl;
    throw(EXIT_FAILURE);
  }

  //check if the time series length can handle all motif subsequences
  if (length < (2 * size - 1) * window) {

    cerr << "ERROR: The sum of the subsequences length is to" <<
      " large to fit into the time series!" << endl;
    throw(EXIT_FAILURE);
  }
}

TSGenerator::~TSGenerator() { }

void TSGenerator::calcRollingMeanStdDev(const vector<double> &timeSeries_in) {

  // obtain running means and stddevs
  if (!means.empty() || !stdDevs.empty()) {

    means.clear();
    stdDevs.clear();
  }

  means.resize(length - window + 1);
  stdDevs.resize(length - window + 1);

  double sum = 0;
  double squareSum = 0;

  // it is faster to multiply than to divide
  double rWindow = 1.0 / (double)window;

  for (int i = 0; i < length - window + 1; i++) {

    sum += timeSeries_in[i];
    squareSum += timeSeries_in[i] * timeSeries_in[i];
  }

  means[0] = sum * rWindow;
  double buf = squareSum * rWindow - means[0] * means[0];
  stdDevs[0] = buf > 0 ? sqrt(buf) : 0;

  for (int i = 1; i < length - window + 1; i++) {

    sum += timeSeries_in[i + window - 1] - timeSeries_in[i - 1];
    means[i] = sum * rWindow;

    squareSum += timeSeries_in[i + window - 1] * timeSeries_in[i
      + window - 1] - timeSeries_in[i - 1] * timeSeries_in[i - 1];
    buf = squareSum * rWindow - means[i] * means[i];
    stdDevs[i] = buf > 0 ? sqrt(buf) : 0;
  }
}

void TSGenerator::updateRollingMeanStdDev(const vector<double> &timeSeries_in,
    int pos_in) {

  int start = pos_in - window + 1;
  int end = pos_in + window - 1;

  // check if we need to update the very first mean and standard deviation
  if (start < 1) {
    means[0] = timeSeries_in[0];

    for (int i = 1; i < window; i++) {
      means[0] += timeSeries_in[i];
    }

    means[0] /= window;

    stdDevs[0] = timeSeries_in[0] * timeSeries_in[0];

    for (int i = 1; i < window; i++) {
      stdDevs[0] += timeSeries_in[i] * timeSeries_in[i];
    }

    stdDevs[0] /= window;

    stdDevs[0] -= means[0] * means[0];

    start = 1;
  }

  // check if we are at the end of the time series
  if (end > length - window)
    end = length - window;

  // update all means overlapping the new subsequence and all potentially
  // shifted once
  for (int i = start; i <= length - window; i++)
    means[i] = (length * means[i - 1] - timeSeries_in[i - 1] + timeSeries_in[i])
      / length;

  // update all changed variances
  for (int i = start; i <= end; i++)
    stdDevs[i] = stdDevs[i - 1] + means[i - 1] * means[i - 1] - means[i]
      * means[i] + (timeSeries_in[i] * timeSeries_in[i] - timeSeries_in[i - 1]
          * timeSeries_in[i - 1]) / window;
}

double TSGenerator::similarity(const vector<double> &timeSeries_in, const int
    subsequenceOnePos_in, const int subsequenceTwoPos_in, const double
    bestSoFar_in) {

  if (timeSeries_in.empty()) {

    cerr << "ERROR: Time series is empty!" << endl;
    throw(EXIT_FAILURE);
  }

  if (subsequenceOnePos_in + window > (int)timeSeries_in.size() ||
    subsequenceOnePos_in < 0) {

    cerr << "ERROR: Position of first subsequence " <<
      subsequenceOnePos_in << " in similarity function is wrong!" << endl;
    throw(EXIT_FAILURE);
  }

  if (subsequenceTwoPos_in + window > (int)timeSeries_in.size() ||
    subsequenceTwoPos_in < 0) {

    cerr << "ERROR: Position of second subsequence " <<
      subsequenceTwoPos_in << " in similarity function is wrong!" << endl;
    throw(EXIT_FAILURE);
  }

  double meanOne = means[subsequenceOnePos_in];
  double stdDevOne = stdDevs[subsequenceOnePos_in];

  double meanTwo = means[subsequenceTwoPos_in];
  double stdDevTwo = stdDevs[subsequenceTwoPos_in];

  //calculate similarity
  double sumOfSquares = 0.0;
  double bestSoFar = bestSoFar_in * bestSoFar_in;

  int itrOne = subsequenceOnePos_in;
  int itrTwo = subsequenceTwoPos_in;

        while (itrOne < subsequenceOnePos_in + window  &&
               itrTwo < subsequenceTwoPos_in + window  &&
               sumOfSquares < bestSoFar                        ) {

    double normed1 = (timeSeries_in[itrOne] - meanOne) / stdDevOne;
    double normed2 = (timeSeries_in[itrTwo] - meanTwo) / stdDevTwo;
    double diff = normed1 - normed2;
    sumOfSquares += diff * diff;

    itrOne++;
    itrTwo++;
  }

  return sqrt(sumOfSquares);
}

void TSGenerator::meanStdDev(const vector<double> &timeSeries_in, const int
    subsequencePos_in, double &mean_out, double &stdDev_out) {

  if (timeSeries_in.empty()) {

    cerr << "ERROR: Time series is empty!" << endl;
    throw(EXIT_FAILURE);
  }

  if (subsequencePos_in + window > (int)timeSeries_in.size() ||
    subsequencePos_in + window < window) {

    cerr << "ERROR: Position of subsequence " <<
      subsequencePos_in << " in meanStdDev function is wrong!" << endl;
    throw(EXIT_FAILURE);
  }


  //calculate mean
  mean_out = accumulate(timeSeries_in.begin() + subsequencePos_in,
      timeSeries_in.begin() + subsequencePos_in + window, 0.0)
    / window;


  //calculate standard deviations
  stdDev_out = 0.0;

  vector<double> diff(window);
  transform(timeSeries_in.begin() + subsequencePos_in, timeSeries_in.begin()
      + subsequencePos_in + window, diff.begin(), [mean_out](double x)
      { return x - mean_out; });
  stdDev_out = sqrt(inner_product(diff.begin(), diff.end(), diff.begin(), 0.0)
      / window);

  if (stdDev_out == 0)
    stdDev_out = 1; // Do not try to divide by a very small number
}

double TSGenerator::similarity(const vector<double> &timeSeries_in, const
    vector<double> &subsequenceOne_in, const double meanOne_in, const double
    stdDevOne_in, const int subsequenceTwoPos_in, const double bestSoFar_in) {

  if (timeSeries_in.empty()) {

    cerr << "ERROR: Time series is empty!" << endl;
    throw(EXIT_FAILURE);
  }

  if (subsequenceOne_in.empty()) {

    cerr << "ERROR: First subsequence is empty!" << endl;
    throw(EXIT_FAILURE);
  }

  if ((int)subsequenceOne_in.size() != window) {

    cerr << "ERROR: The first subsequence in similarity" <<
      " function has wrong size: " << subsequenceOne_in.size() << endl;
    throw(EXIT_FAILURE);
  }

  if (subsequenceTwoPos_in + window > (int)timeSeries_in.size() ||
    subsequenceTwoPos_in < 0) {

    cerr << "ERROR: Position of second subsequence " <<
      subsequenceTwoPos_in << " in similarity function is wrong!" << endl;
    throw(EXIT_FAILURE);
  }

  //calculate means and standard deviations
  double meanTwo = means[subsequenceTwoPos_in];
  double stdDevTwo = stdDevs[subsequenceTwoPos_in];

  //calculate similarity
  double sumOfSquares = 0.0;
  double bestSoFar = bestSoFar_in * bestSoFar_in;

  int itrOne = 0;
  int itrTwo = subsequenceTwoPos_in;

        while (itrOne < window && itrTwo < subsequenceTwoPos_in
            + window && sumOfSquares < bestSoFar) {

    double normed1 = (subsequenceOne_in[itrOne] - meanOne_in) / stdDevOne_in;
    double normed2 = (timeSeries_in[itrTwo] - meanTwo) / stdDevTwo;
    double diff = normed1 - normed2;
    sumOfSquares += diff * diff;

    itrOne++;
    itrTwo++;
  }

  return sqrt(sumOfSquares);
}

double TSGenerator::similarity(const vector<double> &sequence0_in, const
    vector<double> &sequence1_in, const double bestSoFar_in) {

  if (sequence0_in.empty()) {

    cerr << "ERROR: First sequence is empty!" << endl;
    throw(EXIT_FAILURE);
  }

  if (sequence1_in.empty()) {

    cerr << "ERROR: Second sequence is empty!" << endl;
    throw(EXIT_FAILURE);
  }

  //calculate means and standard deviations
  double mean0, stdDev0;
  meanStdDev(sequence0_in, 0, mean0, stdDev0);
  double mean1, stdDev1;
  meanStdDev(sequence1_in, 0, mean1, stdDev1);

  //calculate similarity
  double sumOfSquares = 0.0;
  double bestSoFar = bestSoFar_in * bestSoFar_in;

  for (int i = 0; i < window && sumOfSquares < bestSoFar; i++) {

    double normed0 = (sequence0_in[i] - mean0) / stdDev0;
    double normed1 = (sequence1_in[i] - mean1) / stdDev1;
    double diff = normed0 - normed1;
    sumOfSquares += diff * diff;
  }

  return sqrt(sumOfSquares);
}

void TSGenerator::calculateSubsequence(vector<double> &subsequence_out, int
    type_in, double height_in) {

  if (!subsequence_out.empty()) {

    subsequence_out.clear();
    subsequence_out.resize(0);
  }

  switch (type_in) {

    case 0:
      generateBoxMotif(subsequence_out, 1.0, 1.0, window, height_in);
      break;
    case 1:
      generateTriangleMotif(subsequence_out, 1.0, 1.0, window, height_in);
      break;
    case 2:
      generateSemicircleMotif(subsequence_out, 1.0, 1.0, window, height_in);
      break;
    case 3:
      generateTrapezoidMotif(subsequence_out, 1.0, 1.0, window, height_in);
      break;
    case 4:
      generatePositiveFlankMotif(subsequence_out, 1.0, 1.0, window, height_in);
      break;
    case 5:
      generateNegativeFlankMotif(subsequence_out, 1.0, 1.0, window, height_in);
      break;
    case 6:
      generateSineMotif(subsequence_out, 1.0, 1.0, window, height_in);
      break;
    case 7:
      generateCosineMotif(subsequence_out, 1.0, 1.0, window, height_in);
      break;
    default:
      cerr << "ERROR: Unknown motif type: " << type_in << endl;
      throw(EXIT_FAILURE);
  }
}

void TSGenerator::generateBaseTimeSeries(vector<double> &timeSeries_out) {

  if (!timeSeries_out.empty()) {

    timeSeries_out.clear();
    timeSeries_out.resize(0);
  }

  switch (method) {

    case 0:
      baseTS.simpleRandomWalk(timeSeries_out, length, start, delta, noise);
      break;
    case 1:
      baseTS.realRandomWalk(timeSeries_out, length, start, delta, noise);
      break;
    case 2:
      baseTS.normalRandomWalk(timeSeries_out, length, start, delta, noise);
      break;
    case 3:
      baseTS.simpleRandomWalk(timeSeries_out, length, start, delta, maxi,
          noise);
      break;
    case 4:
      baseTS.realRandomWalk(timeSeries_out, length, start, delta, maxi,
          noise);
      break;
    case 5:
      baseTS.normalRandomWalk(timeSeries_out, length, start, delta, maxi,
          noise);
      break;
    case 6:
      baseTS.uniformRandom(timeSeries_out, length, start, delta, noise);
      break;
    case 7:
      baseTS.normalRandom(timeSeries_out, length, start, delta, noise);
      break;
    case 8:
      baseTS.piecewiseLinearRandom(timeSeries_out, length, start, delta,
          noise);
      break;
    default:
      cerr << "ERROR: Unknown method: " << method << endl;
      throw(EXIT_FAILURE);
  }
}

bool TSGenerator::searchForUnintentionalMatches(const vector<double>
    &timeSeries_in, const vector<int> &motifPositions_in, double similarity_in)
{

  //lower and upper positions of the subsequences overlapping the new motif set
  int lowerBound = max(0, motifPositions_in.back() - window + 1);
  int upperBound = min(motifPositions_in.back() + window, length - window + 1);


  bool candidate;

  for (int itr = 0; itr < length - window + 1; itr++) {

    candidate = true;

    //check if itr is non overlapping to any motif set subsequence
    for (int position : motifPositions_in)
      if (abs(position - itr) < window) {

        candidate = false;
        break;
      }

    //if so check if itr is too similar to any subsequence overlapping with the
    //new added motif
    if (candidate)
      for (int newMotifItr = lowerBound; newMotifItr < upperBound;
          newMotifItr++)
            //non overlapping and ...
        if (abs(newMotifItr - itr) >= window &&
            //... matching
            similarity(timeSeries_in, itr, newMotifItr, similarity_in) <=
            similarity_in)
          return true;
  }

  return false;
}

bool TSGenerator::checkIfThereIsALargerMotifSet(const vector<double>
    &timeSeries_in, const vector<int> &motifPositions_in, double range_in) {

  vector<int> positions;
  int newMotif = motifPositions_in.back();

  for (int position : motifPositions_in)
    positions.push_back(position);

  sort(positions.begin(), positions.end());

  vector<subsequence> overlapping;


  //get all positions of time series subsequences overlapping the injected
  //motif positions
  for (int position : positions)
    if (position != newMotif) {

      //left of the injected motifs
      subsequence leftSubsequences;
      leftSubsequences.position = max(0, position - window + 1);
      leftSubsequences.length = position - leftSubsequences.position;
      overlapping.push_back(leftSubsequences);

      //right of the injected motifs
      subsequence rightSubsequences;
      rightSubsequences.position = position + 1;
      rightSubsequences.length   = min(position + window
          - rightSubsequences.position, length - window
          + 1 - rightSubsequences.position);
      overlapping.push_back(rightSubsequences);
    }


  int motifSize = 1;

  int lastPos = -window;

  //if the new motif has to much matches with left and right overlapping sets
  for (int itr = 0; itr < (int)overlapping.size(); itr++)
    for (int posItr = max(overlapping[itr].position, lastPos + window);
        posItr < overlapping[itr].position + overlapping[itr].length; posItr++)
      if (abs(motifPositions_in.back() - posItr) >= window &&
          similarity(timeSeries_in, motifPositions_in.back(), posItr,
            2 * range_in) <= 2 * range_in) {

        lastPos = posItr;
        motifSize++;
        break;
      }

  if (motifSize > (int)motifPositions_in.size())
    return true;


  subsequence leftNewSubsequence;
  leftNewSubsequence.position = max(0, motifPositions_in.back() - window
      + 1);
  leftNewSubsequence.length   = motifPositions_in.back()
    - leftNewSubsequence.position;

  subsequence rightNewSubsequence;
  rightNewSubsequence.position = motifPositions_in.back() + 1;
  rightNewSubsequence.length   = min(motifPositions_in.back() + window
      - rightNewSubsequence.position, length - window
      + 1 - rightNewSubsequence.position);


  //if one subsequence of the overlapping set of subsequences left of the new
  //motif has to much matches with the left and right overlapping sets
  for (int item = leftNewSubsequence.position; item
      < leftNewSubsequence.position + leftNewSubsequence.length; item++) {

    motifSize = 1;

    for (int itr = rightNewSubsequence.position; itr
        < rightNewSubsequence.position + rightNewSubsequence.length; itr++)
      if (abs(item - itr) >= window && similarity(timeSeries_in, item, itr,
            2 * range_in) <= 2 * range_in) {

        motifSize++;
        break;
      }

    lastPos = -window;

    for (int itr = 0; itr < (int)overlapping.size(); itr++)
      for (int posItr = max(overlapping[itr].position, lastPos + window);
          posItr < overlapping[itr].position + overlapping[itr].length;
          posItr++)
        if (abs(item - posItr) >= window &&
            similarity(timeSeries_in, item, posItr, 2 * range_in) <=
            2 * range_in) {

          lastPos = posItr;
          motifSize++;
          break;
        }

    if (motifSize > (int)motifPositions_in.size())
      return true;
  }


  //if one subsequence of the overlapping set of subsequences right of the new
  //motif has to much matches with the left and right overlapping sets
  for (int item = rightNewSubsequence.position; item
      < rightNewSubsequence.position + rightNewSubsequence.length; item++) {

    motifSize = 1;

    for (int itr = leftNewSubsequence.position; itr
        < leftNewSubsequence.position + leftNewSubsequence.length; itr++)
      if (abs(item - itr) >= window && similarity(timeSeries_in, item, itr,
            2 * range_in) <= 2 * range_in) {

        motifSize++;
        break;
      }

    lastPos = -window;

    for (int itr = 0; itr < (int)overlapping.size(); itr++)
      for (int posItr = max(overlapping[itr].position, lastPos + window);
          posItr < overlapping[itr].position + overlapping[itr].length;
          posItr++)
        if (abs(item - posItr) >= window && similarity(timeSeries_in, item,
              posItr, 2 * range_in) <= 2 * range_in) {

          lastPos = posItr;
          motifSize++;
          break;
        }

    if (motifSize > (int)motifPositions_in.size())
      return true;
  }

  return false;
}

void TSGenerator::run(vector<double> &timeSeries_out, vector<double>
    &d_out, vector<int> &window_out, vector<vector<int>>
    &motifPositions_out) {

  //clear the output buffers
  if (!d_out.empty()) {

    d_out.clear();
    d_out.resize(0);
  }

  if (!window_out.empty()) {

    window_out.clear();
    window_out.resize(0);
  }

  if (!motifPositions_out.empty()) {

    motifPositions_out.clear();
    motifPositions_out.resize(2);
  }
  else
    motifPositions_out.resize(2);

  //declaration stuff
  int positionOne = -1;
  int positionTwo = -1;
  vector<double> motifCenter;
  vector<double> subsequence(window, 0.0);
  bool repeatLoop = true;
  double mean;
  double stdDev;
  double d;
  normal_distribution<double> distribution(0.0, noise);

  do {

    repeatLoop = false;

    //generate a base time series
    generateBaseTimeSeries(timeSeries_out);

    //compute rolling mean and std dev
    calcRollingMeanStdDev(timeSeries_out);

    //determine similarity of the top motif pair in the random synthetic time
    //series
    scrimpPP(timeSeries_out, positionOne, positionTwo, window);
    d = similarity(timeSeries_out, positionOne, positionTwo,
        numeric_limits<double>::max());

    //calculate temporary motif set center subsequence in the window size
    //dimentional room of subsequence values
    calculateSubsequence(motifCenter, type, height);

    if (noise > 0.0) {

      for (auto& value : motifCenter)
        value += distribution(randomEngine);
    }

    mean = 0.0;
    stdDev = 0.0;
    meanStdDev(motifCenter, 0, mean, stdDev);

    //check if there is no other subsequence in the range 3.0 * d / 2.0 around
    //the center subsequence
    for (int itr = 0;
        itr < (int)timeSeries_out.size() - window + 1;
        itr++) {

      vector<double> tmpSubsequence(timeSeries_out.begin() + itr,
          timeSeries_out.begin() + itr + window);

      if (similarity(timeSeries_out, motifCenter, mean, stdDev, itr,
            numeric_limits<double>::max()) <= 3.0 * d / 2.0) {

        repeatLoop = true;
        break;
      }
    }
  } while(repeatLoop);


  //declaration stuff
  int position;
  int retryItr = 0;

  window_out.push_back(motifCenter.size());


  //calculate the subsequences of the motif set
  for (int motifItr = 0; motifItr < size; motifItr++) {

    //calculate the random position in the synthetic time series
    position = freePositions.calculateRandomPosition();

    motifPositions_out[0].push_back(position);

    //check if all motif set subsequences are in range smaller than r / 2.0
    while (retryItr < length + 100) {

      //copy the motif set center subsequence ...
      vector<double> newSubsequence(motifCenter);

      //... and add noise
      if (noise > 0.0)
        for (auto& value : newSubsequence)
          value += distribution(randomEngine);

      //check if the similarity of the subsequence and the center subsequence
      //is at most d / 2.0
      if (similarity(motifCenter, newSubsequence, d / 2.0) < d / 2.0) {

        //backup subsequence
        for (int i = 0; i < window; i++)
          subsequence[i] = timeSeries_out[position + i];

        //add subsequence to synthetic time series
        for (int i = 0; i < window; i++)
          timeSeries_out[i + position] = subsequence[0] + newSubsequence[i];

        //update the running mean and standard deviation
        updateRollingMeanStdDev(timeSeries_out, position);

        if (!searchForUnintentionalMatches(timeSeries_out,
              motifPositions_out[0], d) &&
            !checkIfThereIsALargerMotifSet(timeSeries_out,
              motifPositions_out[0], d / 2.0)) {

          break;
        }
        else {

          //restore old subsequence
          for (int i = 0; i < window; i++)
            timeSeries_out[position + i] = subsequence[i];

          //update the running mean and standard deviation
          updateRollingMeanStdDev(timeSeries_out, position);

          //get new random position in the synthetic time series
          position = freePositions.calculateRandomPosition();

          motifPositions_out[0].back() = position;
        }
      }
      else {

        noise -= noise / 100.0;

        if (noise < numeric_limits<double>::min())
          noise = 0.0;
      }

      if (retryItr == length + 100 - 1) {

        cerr << "ERROR: Cannot add motif set subsequence!" <<
          " Retry or change your settings!" << endl;
        throw(EXIT_FAILURE);
      }

      retryItr++;
    }

    //remove the position from available positions
    freePositions.removePosition();
  }

  d_out.push_back(d / 2.0);


cout << "in" << endl;
  //add top motif pair to output
  positionOne = 0;
  positionTwo = 0;

cout << "length: " << timeSeries_out.size() << endl;
  scrimpPP(timeSeries_out, positionOne, positionTwo, window);
cout << "scrimp pass" << endl;
  d = similarity(timeSeries_out, positionOne, positionTwo,
      numeric_limits<double>::max());
  d_out.push_back(d);

  motifPositions_out[1].push_back(positionOne);
  motifPositions_out[1].push_back(positionTwo);

  window_out.push_back(window);
cout << "out" << endl;
}

