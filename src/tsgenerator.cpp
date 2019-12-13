///\file tsgenerator.cpp
///
///\brief File contains the TSGenerator class definition.
///
///This is the source file of the TSGenerator. The TSGenerator consists of
///setter functions and the run() function. After configuring the TSGenerator
///with the setter functions one starts the time series generation by calling
///the run() function.


#include <tsgenerator.hpp>
#include <scrimpplusplus.hpp>


TSGenerator::TSGenerator(int length_in, int window_in, double delta_in,
    double noise_in, int type_in, int size_in, double height_in)
  : length(abs(length_in)), window(window_in), delta(delta_in),
  noise(noise_in), type(abs(type_in)), size(abs(size_in)), height(height_in),
  freePositions(length, window), randomEngine(random_device().entropy()
    ? random_device()()
    : chrono::system_clock::now().time_since_epoch().count()) {

  // check if type exists
  if (type >= (int) motifTypes.size()) {

    cerr << "ERROR: Wrong motif set type: " << type_in << endl;
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
    double noise_in, string type_in, int size_in, double height_in)
  : length(abs(length_in)), window(window_in), delta(delta_in),
  noise(noise_in), size(abs(size_in)), height(height_in), freePositions(length,
      window), randomEngine(random_device().entropy()
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

tuple<vector< double>, vector<double>> TSGenerator::calcRollingMeanStdDev(const
    vector<double> &timeSeries_in) {

  // obtain rolling/sliding means and stddevs
  vector<double> means(timeSeries_in.size() - window + 1);
  vector<double> stds(timeSeries_in.size() - window + 1);

  double sum = 0;
  double squareSum = 0;

  // it is faster to multiply than to divide
  double rWindowLength = 1.0 / (double)window;

  for (int ww = 0; ww < min(((int)timeSeries_in.size()), window); ww++) {

    sum += timeSeries_in[ww];
    squareSum += timeSeries_in[ww] * timeSeries_in[ww];
  }

  means[0] = sum * rWindowLength;
  double buf = squareSum * rWindowLength - means[0] * means[0];
  stds[0] = buf > 0 ? sqrt(buf) : 0;

  for (int w = 1, end = (int)timeSeries_in.size() - window + 1; w < end;
      w++) {

    sum += timeSeries_in[w + window - 1] - timeSeries_in[w - 1];
    means[w] = sum * rWindowLength;

    squareSum += timeSeries_in[w + window - 1] * timeSeries_in[w
      + window - 1] - timeSeries_in[w - 1] * timeSeries_in[w - 1];
    buf = squareSum * rWindowLength - means[w] * means[w];
    stds[w] = buf > 0 ? sqrt(buf) : 0;
  }

  return make_tuple(means, stds);
}

void TSGenerator::updateRollingMeanStdDev(const vector<double> &timeSeries_in,
    int pos_in, vector<double> &mean, vector<double> & stdDev) {

  int start = pos_in - window + 1;
  int end = pos_in + window - 1;

  // check if we need to update the very first mean and standard deviation
  if (start < 1) {
    mean[0] = timeSeries_in[0];

    for (int i = 1; i < window; i++) {
      mean[0] += timeSeries_in[i];
    }

    mean[0] /= window;

    stdDev[0] = timeSeries_in[0] * timeSeries_in[0];

    for (int i = 1; i < window; i++) {
      stdDev[0] += timeSeries_in[i] * timeSeries_in[i];
    }

    stdDev[0] /= window;

    stdDev[0] -= mean[0] * mean[0];

    start = 1;
  }

  // check if we are at the end of the time series
  if (end > (int) mean.size())
    end = mean.size();

  // update all means overlapping the new subsequence and all potentially
  // shifted once
  for (int i = start; i <= length - window; i++)
    mean[i] = (length * mean[i - 1] - timeSeries_in[i - 1] + timeSeries_in[i])
      / length;

  // update all changed variances
  for (int i = start; i <= end; i++)
    stdDev[i] = stdDev[i - 1] + mean[i - 1] * mean[i - 1] - mean[i] * mean[i]
      + (timeSeries_in[i] * timeSeries_in[i] - timeSeries_in[i - 1]
          * timeSeries_in[i - 1]) / window;
}

double TSGenerator::similarity(const vector<double> &timeSeries_in, const int
    subsequenceOnePos_in, const int subsequenceTwoPos_in, const double
    bestSoFar_in, const vector<double> &means_in, const vector<double>
    &stds_in) {

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

  double meanOne = means_in[subsequenceOnePos_in];
  double stdDevOne = stds_in[subsequenceOnePos_in];

  double meanTwo = means_in[subsequenceTwoPos_in];
  double stdDevTwo = stds_in[subsequenceTwoPos_in];

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
  double meanTwo = 0.0;
  double stdDevTwo = 0.0;

  meanStdDev(timeSeries_in, subsequenceTwoPos_in, meanTwo, stdDevTwo);


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

bool TSGenerator::searchForUnintentionalMatches(const vector<double>
    &timeSeries_in, const vector<int> &motifPositions_in, const
    vector<double> &means_in, const vector<double> &stds_in, double
    similarity_in) {

  //lower and upper positions of the subsequences overlapping the new motif set
  int lowerBound = max(0, motifPositions_in.back() - window + 1);
  int upperBound = min(motifPositions_in.back() + window, length
      - window + 1);


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
            similarity(timeSeries_in, itr, newMotifItr, similarity_in,
              means_in, stds_in) <= similarity_in)
          return true;
  }

  return false;
}

bool TSGenerator::checkIfThereIsALargerMotifSet(const vector<double>
    &timeSeries_in, const vector<int> &motifPositions_in, const
    vector<double> &means_in, const vector<double> &stds_in, double range_in) {

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
            2 * range_in, means_in, stds_in) <= 2 * range_in) {

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
            2 * range_in, means_in, stds_in) <= 2 * range_in) {

        motifSize++;
        break;
      }

    lastPos = -window;

    for (int itr = 0; itr < (int)overlapping.size(); itr++)
      for (int posItr = max(overlapping[itr].position, lastPos + window);
          posItr < overlapping[itr].position + overlapping[itr].length;
          posItr++)
        if (abs(item - posItr) >= window &&
            similarity(timeSeries_in, item, posItr, 2 * range_in, means_in,
              stds_in) <= 2 * range_in) {

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
            2 * range_in, means_in, stds_in) <= 2 * range_in) {

        motifSize++;
        break;
      }

    lastPos = -window;

    for (int itr = 0; itr < (int)overlapping.size(); itr++)
      for (int posItr = max(overlapping[itr].position, lastPos + window);
          posItr < overlapping[itr].position + overlapping[itr].length;
          posItr++)
        if (abs(item - posItr) >= window && similarity(timeSeries_in, item,
              posItr, 2 * range_in, means_in, stds_in) <= 2 * range_in) {

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

  //clear the buffers
  if (!timeSeries_out.empty()) {

    timeSeries_out.clear();
    timeSeries_out.resize(0);
  }

  if (!timeSeriesBackup.empty()) {

    timeSeriesBackup.clear();
    timeSeriesBackup.resize(length);
  }
  else {

    timeSeriesBackup.resize(length);
  }

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

  bool repeatLoop = true;

  double mean = 0.0;
  double stdDev = 0.0;

  vector<double> means;
  vector<double> stdDevs;

  vector<double> meansBackup;
  vector<double> stdDevsBackup;

  double d;

  normal_distribution<double> distribution(0.0, noise);


  do {

    repeatLoop = false;

    //generate a base time series
    baseTS.normalRandomWalk(timeSeries_out, length, 0.0, delta, noise);

    //backup time series
    for (int timeSeriesItr = 0; timeSeriesItr < length; timeSeriesItr++)
      timeSeriesBackup[timeSeriesItr] = timeSeries_out[timeSeriesItr];

    //compute rolling mean and std dev
    tie(means, stdDevs) = calcRollingMeanStdDev(timeSeries_out);

    //determine similarity of the top motif pair in the random synthetic time
    //series
    scrimpPP(timeSeries_out, positionOne, positionTwo, window);
    d = similarity(timeSeries_out, positionOne, positionTwo,
        numeric_limits<double>::max(), means, stdDevs);


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
  double difference;
  int retryItr = 0;

  window_out.push_back(motifCenter.size());


  //calculate the subsequences of the motif set
  for (int motifItr = 0; motifItr < size; motifItr++) {

    //calculate the random position in the synthetic time series
    position = freePositions.calculateRandomPosition();

    motifPositions_out[0].push_back(position);


    //calculate the difference between first and last value of the subsequence
    //in the time series
    difference = timeSeries_out[position] - timeSeries_out[position
      + window - 1];

    //update time series after subsequence
    for (int itr = position + window; itr < (int)timeSeries_out.size();
        itr++)
      timeSeries_out[itr] = timeSeries_out[itr] + difference;

    //store the first value of the subsequence
    double subsequenceFirstValue = timeSeries_out[position];

    //check if all motif set subsequences are in range smaller than r / 2.0
    while (retryItr < length + 100) {

      //copy the motif set center subsequence
      vector<double> newSubsequence(motifCenter);

      if (noise > 0.0) {

        for (auto& value : newSubsequence)
          value += distribution(randomEngine);
      }

      //add subsequence to synthetic time series
      for (int subsequenceItr = position; subsequenceItr < position
          + window; subsequenceItr++)
        timeSeries_out[subsequenceItr] = subsequenceFirstValue
          + newSubsequence[subsequenceItr - position];

      //check if the similarity of the subsequence and the center subsequence
      //is at most d / 2.0
      if (similarity(timeSeries_out, motifCenter, mean, stdDev, position,
            d / 2.0) < d / 2.0) {

        //calculate rolling means and standard deviations
        vector<double> means;
        vector<double> stds;
        tie(means, stds) = calcRollingMeanStdDev(timeSeries_out);


        if (!searchForUnintentionalMatches(timeSeries_out,
              motifPositions_out[0], means, stds, d) &&
            !checkIfThereIsALargerMotifSet(timeSeries_out,
              motifPositions_out[0], means, stds, d / 2.0)) {

          //update backup time series
          for (int timeSeriesItr = position; timeSeriesItr < length;
              timeSeriesItr++)
            timeSeriesBackup[timeSeriesItr] = timeSeries_out[timeSeriesItr];

          break;
        }
        else {

          //reset time series
          for (int timeSeriesItr = position; timeSeriesItr < length;
              timeSeriesItr++)
            timeSeries_out[timeSeriesItr] = timeSeriesBackup[timeSeriesItr];

          //calculate new random position in the synthetic time series
          position = freePositions.calculateRandomPosition();

          motifPositions_out[0].back() = position;


          //calculate the difference between first and last value of the
          //subsequence in the time series
          difference = timeSeries_out[position] - timeSeries_out[position
            + window - 1];

          //update time series after subsequence
          for (int itr = position + window; itr
              < (int)timeSeries_out.size(); itr++)
            timeSeries_out[itr] = timeSeries_out[itr] + difference;

          //store the first value of the subsequence
          subsequenceFirstValue = timeSeries_out[position];
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


  //add top motif pair to output
  positionOne = 0;
  positionTwo = 0;

  tie(means, stdDevs) = calcRollingMeanStdDev(timeSeries_out);

  scrimpPP(timeSeries_out, positionOne, positionTwo, window);
  d = similarity(timeSeries_out, positionOne, positionTwo,
      numeric_limits<double>::max(), means, stdDevs);
  d_out.push_back(d);
  motifPositions_out[1].push_back(positionOne);
  motifPositions_out[1].push_back(positionTwo);
  window_out.push_back(window);
}

