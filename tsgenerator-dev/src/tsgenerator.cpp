///\file tsgenerator.cpp
///
///\brief File contains the TSGenerator class definition.
///
///This is the source file of the TSGenerator. The TSGenerator consists of
///setter functions and the () function. After configuring the TSGenerator
///with the setter functions one starts the time series generation by calling
///the rue) function.


#include <tsgenerator.hpp>


namespace tsg {

  TSGenerator::TSGenerator(const int length_in, const int window_in, const
      double delta_in, const double noise_in, const int type_in, const int
      size_in, const double height_in, const double step_in, const int
      times_in, const int method_in, const double maxi_in, const int gen_in)
    : length(abs(length_in)), window(window_in), delta(delta_in),
    noise(noise_in), type(abs(type_in)), size(abs(size_in)), height(height_in),
    step(abs(step_in)), times(abs(times_in)), method(method_in),
    maxi(abs(maxi_in)), gen(abs(gen_in)), freePositions(length, window),
    randomEngine(std::random_device().entropy()
      ? std::random_device()()
      : (unsigned
        int)(std::chrono::system_clock::now().time_since_epoch().count())) {

    //check if type exists
    if (type >= (int) motifTypes.size()) {

      std::cerr << "ERROR: Wrong motif set type: " << type_in << std::endl;
      throw(EXIT_FAILURE);
    }

    //check if method exists
    if (method >= (int) methods.size()) {

      std::cerr << "ERROR: Unknown method: " << method_in << std::endl;
      throw(EXIT_FAILURE);
    }

    //check if the motif amount is set properly
    if (size < 3) {

      std::cerr << "ERROR: Wrong motif amount: " << size_in << std::endl;
      throw(EXIT_FAILURE);
    }

    //check if the generator type is set properly
    if (gen > 2) {

      std::cerr << "ERROR: Wrong generator type: " << gen_in << std::endl;
      throw(EXIT_FAILURE);
    }

    //check if the time series length can handle all motif subsequences
    if (length < (2 * size - 1) * window) {

      std::cerr << "ERROR: The sum of the subsequences length is to" <<
        " large to fit into the time series!" << std::endl;
      throw(EXIT_FAILURE);
    }
  }

  TSGenerator::TSGenerator(const int length_in, const int window_in, const
      double delta_in, const double noise_in, const word type_in, const int
      size_in, const double height_in, const double step_in, const int
      times_in, const word method_in, const double maxi_in, const word gen_in)
    : length(abs(length_in)), window(window_in), delta(delta_in),
    noise(noise_in), size(abs(size_in)), height(height_in), step(abs(step_in)),
    times(abs(times_in)), maxi(abs(maxi_in)), freePositions(length, window),
    randomEngine(std::random_device().entropy()
      ? std::random_device()()
      : (unsigned
        int)(std::chrono::system_clock::now().time_since_epoch().count())) {

    // get the type number
    type = (int)(std::distance(motifTypes.begin(),
          std::find(motifTypes.begin(), motifTypes.end(), type_in)));

    // check if type exists
    if (type >= (int) motifTypes.size()) {

      std::cerr << "ERROR: Wrong motif set type: " << type_in << std::endl;
      throw(EXIT_FAILURE);
    }

    // get the method
    method = (int)(std::distance(methods.begin(), std::find(methods.begin(),
            methods.end(), method_in)));

    // check if method exists
    if (method >= (int) methods.size()) {

      std::cerr << "ERROR: Unknown method: " << method_in << std::endl;
      throw(EXIT_FAILURE);
    }

    // get the generator
    gen = (int)(std::distance(gens.begin(), std::find(gens.begin(), gens.end(),
            gen_in)));

    // check if generator exists
    if (gen >= (int) gens.size()) {

      std::cerr << "ERROR: Unknown generator: " << gen_in << std::endl;
      throw(EXIT_FAILURE);
    }

    //check if the motif amount is set properly
    if (size < 3) {

      std::cerr << "ERROR: Wrong motif amount: " << size_in << std::endl;
      throw(EXIT_FAILURE);
    }

    //check if the time series length can handle all motif subsequences
    if (length < (2 * size - 1) * window) {

      std::cerr << "ERROR: The sum of the subsequences length is to" <<
        " large to fit into the time series!" << std::endl;
      throw(EXIT_FAILURE);
    }
  }

  TSGenerator::~TSGenerator() { }

  void TSGenerator::calcRunnings(const rseq &sequence_in) {

    // obtain running sum and sum of square
    if (!sums.empty() || !sumSquares.empty()) {

      sums.clear();
      sumSquares.clear();
    }

    sums.resize(sequence_in.size() - window + 1);
    sumSquares.resize(sequence_in.size() - window + 1);

    sums[0] = sequence_in[0];
    sumSquares[0] = sequence_in[0] * sequence_in[0];

    for (int i = 1; i < window; i++) {

      sums[0] += sequence_in[i];
      sumSquares[0] += sequence_in[i] * sequence_in[i];
    }

    for (int i = 0; i < (int)sequence_in.size() - window; i++) {

      sums[i + 1] = sums[i] - sequence_in[i] + sequence_in[i + window];

      sumSquares[i + 1] = sumSquares[i] - sequence_in[i] * sequence_in[i]
        + sequence_in[i + window] * sequence_in[i + window];
    }
  }

  void TSGenerator::updateRunnings(const rseq &sequence_in, const int pos_in) {

    int start = pos_in - window + 1;
    int end = pos_in + window - 1;

    // check if we need to update the very first sum and sum of squares
    if (start < 1) {

      sums[0] = sequence_in[0];
      sumSquares[0] = sequence_in[0] * sequence_in[0];

      for (int i = 1; i < window; i++) {

        sums[0] += sequence_in[i];
        sumSquares[0] += sequence_in[i] * sequence_in[i];
      }

      start = 1;
    }

    // check if we are at the end of the time series
    if (end > (int)sequence_in.size() - window)
      end = (int)(sequence_in.size()) - window;

    // update all changed runnings
    for (int i = start - 1; i < end; i++) {

      sums[i + 1] = sums[i] - sequence_in[i] + sequence_in[i + window];

      sumSquares[i + 1] = sumSquares[i] + sequence_in[i + window]
        * sequence_in[i + window] - sequence_in[i] * sequence_in[i];
    }
  }

  double TSGenerator::similarity(const rseq &timeSeries_in, const int pos0_in,
      const int pos1_in, const double bestSoFar_in) {

    if (timeSeries_in.empty()) {

      std::cerr << "ERROR: Time series is empty!" << std::endl;
      throw(EXIT_FAILURE);
    }

    if (pos0_in + window > length || pos0_in < 0) {

      std::cerr << "ERROR: Position of first subsequence " <<
        pos0_in << " in similarity function is wrong!" << std::endl;
      throw(EXIT_FAILURE);
    }

    if (pos1_in + window > length || pos1_in < 0) {

      std::cerr << "ERROR: Position of second subsequence " <<
        pos1_in << " in similarity function is wrong!" << std::endl;
      throw(EXIT_FAILURE);
    }

    double rWindow = 1.0 / window;

    double mean0 = sums[pos0_in] * rWindow;
    double stdDev0 = sumSquares[pos0_in] * rWindow  - mean0 * mean0;
    stdDev0 = stdDev0 < 1.0 ? 1.0 : sqrt(stdDev0);

    double mean1 = sums[pos1_in] * rWindow;
    double stdDev1 = sumSquares[pos1_in] * rWindow  - mean1 * mean1;
    stdDev1 = stdDev1 < 1.0 ? 1.0 : sqrt(stdDev1);

    //calculate the similarity
    double sumOfSquares = 0.0;
    double bestSoFar = bestSoFar_in * bestSoFar_in;
    double norm0;
    double norm1;
    double diff;

    for (int i = 0; i < window && sumOfSquares < bestSoFar; i++) {

      norm0 = (timeSeries_in[pos0_in + i] - mean0) / stdDev0;
      norm1 = (timeSeries_in[pos1_in + i] - mean1) / stdDev1;
      diff = norm0 - norm1;
      sumOfSquares += diff * diff;
    }

    return sqrt(sumOfSquares);
  }

  void TSGenerator::meanStdDev(const rseq &sequence_in, double &mean_out,
      double &stdDev_out) {

    if (sequence_in.empty()) {

      std::cerr << "ERROR: Time series is empty!" << std::endl;
      throw(EXIT_FAILURE);
    }

    double rWindow = 1.0 / window;

    //compute the mean and standard deviation
    double sum = sequence_in[0];
    double sumSquare = sequence_in[0] * sequence_in[0];

    for (int i = 1; i < (int)sequence_in.size(); i++) {

      sum += sequence_in[i];
      sumSquare += sequence_in[i] * sequence_in[i];
    }

    mean_out = sum * rWindow;
    stdDev_out = sumSquare * rWindow - mean_out * mean_out;
    stdDev_out = stdDev_out < 1.0 ? 1.0 : sqrt(stdDev_out);
  }

  double TSGenerator::similarity(const rseq &timeSeries_in, const rseq
      &sequence_in, const double mean_in, const double stdDev_in, const int
      pos_in, const double bestSoFar_in) {

    if (timeSeries_in.empty()) {

      std::cerr << "ERROR: Time series is empty!" << std::endl;
      throw(EXIT_FAILURE);
    }

    if (sequence_in.empty()) {

      std::cerr << "ERROR: First subsequence is empty!" << std::endl;
      throw(EXIT_FAILURE);
    }

    if ((int)sequence_in.size() != window) {

      std::cerr << "ERROR: The first subsequence in similarity" <<
        " function has wrong size: " << sequence_in.size() << std::endl;
      throw(EXIT_FAILURE);
    }

    if (pos_in + window > (int)timeSeries_in.size() ||
      pos_in < 0) {

      std::cerr << "ERROR: Position of second subsequence " <<
        pos_in << " in similarity function is wrong!" << std::endl;
      throw(EXIT_FAILURE);
    }

    double rWindow = 1.0 / window;

    //compute mean and standard deviation
    double mean1 = sums[pos_in] * rWindow;
    double stdDev1 = sumSquares[pos_in] * rWindow - mean1 * mean1;
    stdDev1 = stdDev1 < 1.0 ? 1.0 : sqrt(stdDev1);

    //compute the similarity
    double sumOfSquares = 0.0;
    double bestSoFar = bestSoFar_in * bestSoFar_in;

    double norm0;
    double norm1;
    double diff;

    for (int i = 0; i < window && sumOfSquares < bestSoFar; i++) {

      norm0 = (sequence_in[i] - mean_in) / stdDev_in;
      norm1 = (timeSeries_in[pos_in + i] - mean1) / stdDev1;
      diff = norm0 - norm1;
      sumOfSquares += diff * diff;
    }

    return sqrt(sumOfSquares);
  }

  double TSGenerator::similarity(const rseq &sequence0_in, const rseq
      &sequence1_in, const double bestSoFar_in) {

    if (sequence0_in.empty()) {

      std::cerr << "ERROR: First sequence is empty!" << std::endl;
      throw(EXIT_FAILURE);
    }

    if (sequence1_in.empty()) {

      std::cerr << "ERROR: Second sequence is empty!" << std::endl;
      throw(EXIT_FAILURE);
    }

    //compute the means and standard deviations
    double mean0, stdDev0;
    meanStdDev(sequence0_in, mean0, stdDev0);
    stdDev0 = sqrt(stdDev0);
    stdDev0 = stdDev0 < 1.0 ? 1.0 : stdDev0;

    double mean1, stdDev1;
    meanStdDev(sequence1_in, mean1, stdDev1);
    stdDev1 = sqrt(stdDev1);
    stdDev1 = stdDev1 < 1.0 ? 1.0 : stdDev1;

    //compute the similarity
    double sumOfSquares = 0.0;
    double bestSoFar = bestSoFar_in * bestSoFar_in;

    double norm0;
    double norm1;
    double diff;

    for (int i = 0; i < window && sumOfSquares < bestSoFar; i++) {

      norm0 = (sequence0_in[i] - mean0) / stdDev0;
      norm1 = (sequence1_in[i] - mean1) / stdDev1;
      diff = norm0 - norm1;
      sumOfSquares += diff * diff;
    }

    return sqrt(sumOfSquares);
  }

  void TSGenerator::calculateSubsequence(rseq &subsequence_out, int type_in,
      double height_in) {

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
        generatePositiveFlankMotif(subsequence_out, 1.0, 1.0, window,
            height_in);
        break;
      case 5:
        generateNegativeFlankMotif(subsequence_out, 1.0, 1.0, window,
            height_in);
        break;
      case 6:
        generateSineMotif(subsequence_out, 1.0, 1.0, window, height_in);
        break;
      case 7:
        generateCosineMotif(subsequence_out, 1.0, 1.0, window, height_in);
        break;
      default:
        std::cerr << "ERROR: Unknown motif type: " << type_in << std::endl;
        throw(EXIT_FAILURE);
    }
  }

  void TSGenerator::generateBaseTimeSeries(rseq &timeSeries_out) {

    if (!timeSeries_out.empty()) {

      timeSeries_out.clear();
      timeSeries_out.resize(0);
    }

    switch (method) {

      case 0:
        baseTS.simpleRandomWalk(timeSeries_out, length, delta, noise);
        break;
      case 1:
        baseTS.realRandomWalk(timeSeries_out, length, delta, noise);
        break;
      case 2:
        baseTS.normalRandomWalk(timeSeries_out, length, delta, noise);
        break;
      case 3:
        baseTS.linearRandomWalk(timeSeries_out, length, delta, step,
            noise);
        break;
      case 4:
        baseTS.simpleRandomWalk(timeSeries_out, length, delta, maxi,
            noise);
        break;
      case 5:
        baseTS.realRandomWalk(timeSeries_out, length, delta, maxi,
            noise);
        break;
      case 6:
        baseTS.normalRandomWalk(timeSeries_out, length, delta, maxi,
            noise);
        break;
      case 7:
        baseTS.linearRandomWalk(timeSeries_out, length, delta, step, maxi,
            noise);
        break;
      case 8:
        baseTS.uniformRandom(timeSeries_out, length, delta, noise);
        break;
      case 9:
        baseTS.normalRandom(timeSeries_out, length, delta, noise);
        break;
      case 10:
        baseTS.piecewiseLinearRandom(timeSeries_out, length, delta,
            noise);
        break;
      case 11:
        baseTS.splineRepeated(timeSeries_out, length, delta, (int)step, times,
            noise);
        break;
      default:
        std::cerr << "ERROR: Unknown method: " << method << std::endl;
        throw(EXIT_FAILURE);
    }
  }

  bool TSGenerator::searchForUnintentionalMatches(const rseq &timeSeries_in,
      const iseq &motifPositions_in, double similarity_in) {

    //lower and upper positions of the subsequences overlapping the new motif
    //set
    int lowerBound = std::max(0, motifPositions_in.back() - window + 1);
    int upperBound = std::min(motifPositions_in.back() + window, length
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

      //if so check if itr is too similar to any subsequence overlapping with
      //the new added motif
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

  bool TSGenerator::checkIfThereIsALargerMotifSet(const rseq &timeSeries_in,
      const iseq &motifPositions_in, double range_in) {

    iseq positions;
    int newMotif = motifPositions_in.back();

    for (int position : motifPositions_in)
      positions.push_back(position);

    sort(positions.begin(), positions.end());

    subsequences overlapping;


    //get all positions of time series subsequences overlapping the injected
    //motif positions
    for (int position : positions)
      if (position != newMotif) {

        //left of the injected motifs
        subsequence leftSubsequences;
        leftSubsequences.position = std::max(0, position - window + 1);
        leftSubsequences.length = position - leftSubsequences.position;
        overlapping.push_back(leftSubsequences);

        //right of the injected motifs
        subsequence rightSubsequences;
        rightSubsequences.position = position + 1;
        rightSubsequences.length   = std::min(position + window
            - rightSubsequences.position, length - window
            + 1 - rightSubsequences.position);
        overlapping.push_back(rightSubsequences);
      }


    int motifSize = 1;

    int lastPos = -window;

    //if the new motif has to much matches with left and right overlapping sets
    for (int itr = 0; itr < (int)overlapping.size(); itr++)
      for (int posItr = std::max(overlapping[itr].position, lastPos + window);
          posItr < overlapping[itr].position + overlapping[itr].length;
          posItr++)
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
    leftNewSubsequence.position = std::max(0, motifPositions_in.back() - window
        + 1);
    leftNewSubsequence.length   = motifPositions_in.back()
      - leftNewSubsequence.position;

    subsequence rightNewSubsequence;
    rightNewSubsequence.position = motifPositions_in.back() + 1;
    rightNewSubsequence.length   = std::min(motifPositions_in.back() + window
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
        for (int posItr = std::max(overlapping[itr].position, lastPos + window);
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


    //if one subsequence of the overlapping set of subsequences right of the
    //new motif has to much matches with the left and right overlapping sets
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
        for (int posItr = std::max(overlapping[itr].position, lastPos + window);
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

  void TSGenerator::injectPairMotif(rseq &timeSeries_out, rseq &d_out, iseq
      &window_out, iseqs &motifPositions_out) {

    std::cout << "pair motif not yet implemented" << std::endl;
  }

  void TSGenerator::injectSetMotif(rseq &timeSeries_out, rseq &d_out, iseq
      &window_out, iseqs &motifPositions_out) {

    std::cout << "set motif not yet implemented" << std::endl;
  }

  void TSGenerator::injectLatentMotif(rseq &timeSeries_out, rseq &d_out, iseq
      &window_out, iseqs &motifPositions_out) {

    //declaration stuff
    int positionOne = -1;
    int positionTwo = -1;
    rseq motifCenter;
    rseq subsequence(window, 0.0);
    bool repeatLoop = true;
    double mean;
    double stdDev;
    double d = std::numeric_limits<double>::max();
    std::normal_distribution<double> distribution(0.0, abs(noise) / 4.0 <= 0.0
        ? std::numeric_limits<double>::min() : noise / 4.0);

    for (int i = 0; i < 50 && repeatLoop; i++) {

      repeatLoop = false;

      //generate a base time series
      generateBaseTimeSeries(timeSeries_out);

      //compute running mean and std dev
      calcRunnings(timeSeries_out);

      //determine similarity of the top motif pair in the random synthetic time
      //series
      scrimpPP(timeSeries_out, sums, sumSquares, positionOne, positionTwo,
          window);
      d = similarity(timeSeries_out, positionOne, positionTwo,
          std::numeric_limits<double>::max());

      //calculate temporary motif set center subsequence in the window size
      //dimentional room of subsequence values
      calculateSubsequence(motifCenter, type, height);

      if (noise > 0.0)
        for (auto& value : motifCenter)
          value += distribution(randomEngine);

      mean = 0.0;
      stdDev = 0.0;
      meanStdDev(motifCenter, mean, stdDev);

      //check if there is no other subsequence in the range 3.0 * d / 2.0
      //around the center subsequence
      for (int itr = 0;
          itr < (int)timeSeries_out.size() - window + 1;
          itr++) {

        rseq tmpSubsequence(timeSeries_out.begin() + itr,
            timeSeries_out.begin() + itr + window);

        if (similarity(timeSeries_out, motifCenter, mean, stdDev, itr,
              std::numeric_limits<double>::max()) <= 3.0 * d / 2.0) {

          repeatLoop = true;
          break;
        }
      }
    }

    //declaration stuff
    int position;
    int retryItr = 0;
    double value;
    double min, max;
    double lth = noise / length;

    window_out.push_back((int)(motifCenter.size()));


    //inject sequences into the time series
    for (int motifItr = 0; motifItr < size; motifItr++) {

      //compute the random position for the subsequence
      position = freePositions.calculateRandomPosition();

      motifPositions_out[0].push_back(position);

      //try to inject another sequence
      while (retryItr < length + 100) {

        std::normal_distribution<double> distributionNoise(0.0, abs(noise)
            / 4.0 <= 0.0 ? std::numeric_limits<double>::min() : noise / 4.0);

        //copy another motif sequence ...
        rseq newSubsequence(motifCenter);

        //... with noise
        if (noise > 0.0)
          for (auto& item : newSubsequence)
            item += distributionNoise(randomEngine);

        //check if the sequence is within range d / 2.0 of the center
        if (similarity(motifCenter, newSubsequence, d / 2.0) < d / 2.0) {

          //backup subsequence at position
          for (int i = 0; i < window; i++)
            subsequence[i] = timeSeries_out[position + i];

          value = subsequence[0];

          //make sure we are in maxi when maxi can handle the motif height
          if (abs(height) <= 2 * maxi) {

            max = newSubsequence[0];
            min = max;

            for (int i = 0; i < window; i++) {

              if (newSubsequence[i] < min)
                min = newSubsequence[i];

              if (newSubsequence[i] > max)
                max = newSubsequence[i];
            }

            if (value + max > maxi)
              value = maxi - max;

            if (value - min < -maxi)
              value = -maxi - min;
          }

          //inject sequence into the time series
          for (int i = 0; i < window; i++)
            timeSeries_out[i + position] = value + newSubsequence[i];

          //update the running sum and sum of square
          updateRunnings(timeSeries_out, position);

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

            //update the running mean and variance
            updateRunnings(timeSeries_out, position);

            //get new random position in the synthetic time series
            position = freePositions.calculateRandomPosition();

            motifPositions_out[0].back() = position;

            noise -= lth;
          }
        }
        else {

          noise -= lth;
        }

        if (retryItr == length + 100 - 1) {

          std::cerr << "ERROR: Cannot add motif set subsequence!" <<
            " Retry or change your settings!" << std::endl;
          throw(EXIT_FAILURE);
        }

        retryItr++;
      }

      //remove the position from available positions
      freePositions.removePosition();
    }

    d_out.push_back(d / 2.0);
  }

  void TSGenerator::run(rseq &timeSeries_out, rseq &d_out, iseq &window_out,
      iseqs &motifPositions_out) {

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


    //choose correct generator
    switch (gen) {

      case 0:
        injectPairMotif(timeSeries_out, d_out, window_out, motifPositions_out);
        break;
      case 1:
        injectSetMotif(timeSeries_out, d_out, window_out, motifPositions_out);
        break;
      case 2:
        injectLatentMotif(timeSeries_out, d_out, window_out,
            motifPositions_out);
        break;
      default:
        std::cerr << "ERROR: Unknown motif generator: " << gen << std::endl;
        throw(EXIT_FAILURE);
    }


    //add top motif pair to output
    int positionOne = 0;
    int positionTwo = 0;

    scrimpPP(timeSeries_out, sums, sumSquares, positionOne, positionTwo,
        window);
    d_out.push_back(similarity(timeSeries_out, positionOne, positionTwo,
        std::numeric_limits<double>::max()));

    motifPositions_out[1].push_back(positionOne);
    motifPositions_out[1].push_back(positionTwo);

    window_out.push_back(window);
  }
}
