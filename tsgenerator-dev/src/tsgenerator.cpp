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
    step(abs(step_in)), times(abs(times_in)), method(abs(method_in)),
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
    method = (int)(std::distance(methods.begin(),
          std::find(methods.begin(), methods.end(), method_in)));

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

  bool TSGenerator::smallerDistance(const rseq &timeSeries_in, const iseq
      &motifPositions_in, const double similarity_in) {

    //lower and upper positions of the overlapping subsequences
    int lowerBound = std::max(0, motifPositions_in.back() - window + 1);
    int upperBound = std::min(motifPositions_in.back() + window, length
        - window + 1);


    //each overlapping subsequence
    for (int i = lowerBound; i < upperBound; i++) {

      //each subsequence of the time series
      for (int j = 0; j < length - window + 1; j++) {

        //the subsequences are not the injected sequences and non
        //overlapping
        if (!(j == motifPositions_in[0] && i == motifPositions_in[1]) &&
            abs(i - j) >= window) {

          //the similarity of both subsequences is smaller than the similarity
          //of both injected sequences
          if (similarity(timeSeries_in, i, j, similarity_in) <= similarity_in)
            return true;
        }
      }
    }

    return false;
  }

  bool TSGenerator::largerMotifSet(const rseq &timeSeries_in, const iseq
      &motifPositions_in, double range_in) {

    int motif = motifPositions_in.back();
    int start = motif - window + 1;
    int end = motif + window - 1;

    if (start < 0)
      start = 0;

    if (end > (int)timeSeries_in.size() - window)
      end = timeSeries_in.size() - window;

    //collect new matches
    iseq mats;

    //trivial matches
    for (int i = start; i <= end; i++)
      mats.push_back(i);

    //all other matches
    for (int i = start; i <= end; i++) {

      for (int j = 0; j <= (int)timeSeries_in.size() - window; j++) {

        if (similarity(timeSeries_in, i, j, range_in) <= range_in) {

          bool ins = false;

          for (int k = 0; k < (int)mats.size(); k++) {

            if (mats[k] == j) {

              ins = true;
              k = mats.size();
            }
            else if (j < mats[k]) {

              ins = true;
              mats.insert(mats.begin() + k, j);
              k = mats.size();
            }
          }

          if (!ins)
            mats.push_back(j);
        }
      }
    }

    //no matchings no larger set motif
    if ((int)mats.size() < 1)
      return false;

    //check matches for larger set motif
    for (int i = 0; i < (int)mats.size(); i++) {

      int size = 1;

      for (int j = 0; j <= (int)timeSeries_in.size() - window; j++) {

        //subseqeunce not overlapping the ith set motif
        if (!(abs(j - mats[i]) < window)) {

          //matching
          if (similarity(timeSeries_in, mats[i], j, range_in) <= range_in) {

            //filter overlaps
            j += window - 1;
            size++;

            //is the set motif larger?
            if (size > (int)motifPositions_in.size())
              return true;
          }
        }
      }
    }

    return false;
  }

  int TSGenerator::largerMotifSet(const rseq &timeSeries_in, const iseq
      motifPositions_in, const int size_in, const double range_in) {

    int motif = motifPositions_in.back();
    int start = motif - window + 1;
    int end = motif + window - 1;

    if (start < 0)
      start = 0;

    if (end > (int)timeSeries_in.size() - window)
      end = timeSeries_in.size() - window;

    //collect new matches
    iseq mats;

    //trivial matches
    for (int i = start; i <= end; i++)
      mats.push_back(i);

    //all other matches
    for (int i = start; i <= end; i++) {

      for (int j = 0; j <= (int)timeSeries_in.size() - window; j++) {

        if (similarity(timeSeries_in, i, j, range_in) <= range_in) {

          bool ins = false;

          for (int k = 0; k < (int)mats.size(); k++) {

            if (mats[k] == j) {

              ins = true;
              k = mats.size();
            }
            else if (j < mats[k]) {

              ins = true;
              mats.insert(mats.begin() + k, j);
              k = mats.size();
            }
          }

          if (!ins)
            mats.push_back(j);
        }
      }
    }

    //no matchings no larger set motif
    if ((int)mats.size() < 1)
      return 1;

    //check matches for larger set motif
    int largestSize = 1;

    for (int i = 0; i < (int)mats.size(); i++) {

      int size = 1;

      for (int j = 0; j <= (int)timeSeries_in.size() - window; j++) {

        //subsequences not overlapping the set motif i
        if (!(abs(j - mats[i]) < window)) {

          //matching
          if (similarity(timeSeries_in, mats[i], j, range_in) <= range_in) {

            //filter overlaps
            j += window - 1;
            size++;

            //is the set motif larger?
            if (size > size_in)
              return size;
          }
        }
      }

      if (size > largestSize)
        largestSize = size;
    }

    return largestSize;
  }

  void TSGenerator::injectPairMotif(rseq &timeSeries_out, rseqs &motif_out) {

    //empty motifs
    if (!motif_out.empty()) {

      for (auto &motif : motif_out)
        motif.clear();

      motif_out.clear();
      motif_out.resize(2);
      motif_out[0].resize(window);
      motif_out[1].resize(window);
    }
    else {

      motif_out.resize(2);
      motif_out[0].resize(window);
      motif_out[1].resize(window);
    }

    //declaration stuff
    int motifPos0, motifPos1;
    int pos0 = -1;
    int pos1 = -1;
    double mean, stdDev;
    rseq motif;
    rseq subsequence(window, 0.0);
    double d = std::numeric_limits<double>::max();
    std::normal_distribution<double> distribution(0.0, abs(noise) * 0.5 <= 0.0
        ? std::numeric_limits<double>::min() : noise * 0.5);
    double value;
    double min, max;
    std::uniform_real_distribution<double> distHeight(0.6, 1.4);
    double stretch = distHeight(randomEngine);

    //generate a base time series
    generateBaseTimeSeries(timeSeries_out);

    //compute first pair motif sequence
    calculateSubsequence(motif, type, height);
    motif_out[0] = motif;

    //stretch and add noise
    for (auto &item : motif) {

      item += distribution(randomEngine);
      item *= stretch;
    }

    //get new random position in the synthetic time series
    motifPos0 = freePositions.calculateRandomPosition();
    freePositions.removePosition();

    //inject sequence into the time series
    value = timeSeries_out[motifPos0];

    //make sure we are in maxi when maxi can handle the motif height
    if (method > 3
        && method < 8
        && abs(height) <= 2 * maxi) {

      max = motif[0];
      min = max;

      for (int i = 0; i < window; i++) {

        if (motif[i] < min)
          min = motif[i];

        if (motif[i] > max)
          max = motif[i];
      }

      if (value + max > maxi)
        value = maxi - max;

      if (value - min < -maxi)
        value = -maxi - min;
    }

    for (int i = 0; i < window; i++)
      timeSeries_out[i + motifPos0] = value + motif[i];

    //unstretch motif
    for (auto &item : motif)
      item /= stretch;

    //compute running mean and std dev
    calcRunnings(timeSeries_out);

    //determine simlarity of the top motif pair in the random synthetic time
    //series
    tpm(timeSeries_out, sums, sumSquares, pos0, pos1, window);
    d = similarity(timeSeries_out, pos0, pos1, d);

    //compute the second motif sequence
    double lth = noise / length;

    noise += lth;

    rseq second;
    second.resize(motif.size());

    do {

      stretch = distHeight(randomEngine);

      if (noise > 0.0)
        noise -= lth;
      else
        noise = 0.0;

      //copy motif
      for (int i = 0; i < window; i++)
        second[i] = motif[i];

      //stretch and add noise
      for (auto &item : second) {

        item += distribution(randomEngine);
        item *= stretch;
      }

      mean = 0.0;
      stdDev = 0.0;
      meanStdDev(second, mean, stdDev);

    } while(similarity(timeSeries_out, second, mean, stdDev, motifPos0, d) >=
        d);

    motif_out[1] = second;

    rseq tmp;
    tmp.resize(motif.size());

    //get a random position for the second motif sequence
    motifPos1 = freePositions.calculateRandomPosition();

    double dTmp;

    for(int i = 0; i < length + 1000; i++) {

      //inject the second motif sequence
      value = timeSeries_out[motifPos1];

      //make sure we are in maxi when maxi can handle the motif height
      if (method > 3
        && method < 8
        && abs(height) <= 2 * maxi) {

        max = second[0];
        min = max;

        for (int i = 0; i < window; i++) {

          if (second[i] < min)
            min = second[i];

          if (second[i] > max)
            max = second[i];
        }

        if (value + max > maxi)
          value = maxi - max;

        if (value - min < -maxi)
          value = -maxi - min;
      }

      for (int i = 0; i < window; i++) {

        tmp[i] = timeSeries_out[i + motifPos1];
        timeSeries_out[i + motifPos1] = value + second[i];
      }

      //update the runnings
      updateRunnings(timeSeries_out, motifPos1);

      dTmp = similarity(timeSeries_out, motifPos0, motifPos1, d);

      //check for success
      if(!smallerDistance(timeSeries_out, { motifPos0, motifPos1 },
            dTmp))
        break;

      if(i == length + 999) {

        std::cerr << "ERROR: Cannot inject second pair motif sequence!" <<
          std::endl;
        throw(EXIT_FAILURE);
      }

      //reset time series
      for (int i = 0; i < window; i++)
        timeSeries_out[i + motifPos1] = tmp[i];

      //update the runnings
      updateRunnings(timeSeries_out, motifPos1);

      //get new random position for the second motif sequence
      motifPos1 = freePositions.calculateRandomPosition();
    }
  }

  void TSGenerator::injectSetMotif(rseq &timeSeries_out, rseqs &motif_out, rseq
      &d_out, iseqs &pos_out) {

    //empty motifs
    if (!motif_out.empty()) {

      for (auto &motif : motif_out)
        motif.clear();

      motif_out.clear();
      motif_out.resize(1);
      motif_out[0].resize(window);
    }
    else {

      motif_out.resize(1);
      motif_out[0].resize(window);
    }

    //declaration stuff
    int positionOne = -1;
    int positionTwo = -1;
    rseq motif;
    rseq subsequence(window, 0.0);
    double d = std::numeric_limits<double>::max();
    std::normal_distribution<double> distNoise(0.0, abs(noise) * 0.5 <= 0.0
        ? std::numeric_limits<double>::min() : noise * 0.5);
    int position;
    double value;
    double min, max;
    std::uniform_real_distribution<double> distStretch(0.6, 1.4);
    double stretch = distStretch(randomEngine);

    //generate a base time series
    generateBaseTimeSeries(timeSeries_out);

    //compute the motif sequence
    calculateSubsequence(motif, type, height);

    motif_out[0] = motif;

    //stretch and add noise
    for (auto &item : motif) {

      item += distNoise(randomEngine);
      item *= stretch;
    }

    //get new random position in the synthetic time series
    position = freePositions.calculateRandomPosition();
    freePositions.removePosition();

    pos_out[0].push_back(position);

    //inject sequence into the time series
    value = timeSeries_out[position];

    //make sure we are in maxi when maxi can handle the motif height
    if (method > 3
        && method < 8
        && abs(height) <= 2 * maxi) {

      max = motif[0];
      min = max;

      for (int i = 0; i < window; i++) {

        if (motif[i] < min)
          min = motif[i];

        if (motif[i] > max)
          max = motif[i];
      }

      if (value + max > maxi)
        value = maxi - max;

      if (value - min < -maxi)
        value = -maxi - min;
    }

    for (int i = 0; i < window; i++)
      timeSeries_out[i + position] = value + motif[i];

    //unstretch
    for (auto &item : motif)
      item /= stretch;

    //compute running mean and std dev
    calcRunnings(timeSeries_out);

    //determine similarity of the top motif pair in the random synthetic time
    //series
    tpm(timeSeries_out, sums, sumSquares, positionOne, positionTwo, window);
    d = 0.9999999 * similarity(timeSeries_out, positionOne, positionTwo, d);

    d_out.push_back(d);

    //declaration stuff
    int retries = 20;
    double lth = noise / retries;
    double noiseBackup = noise;

    //inject sequences into the time series
    for (int motifItr = 1; motifItr < size; motifItr++) {

      //compute the random position for the subsequence
      position = freePositions.calculateRandomPosition();

      pos_out[0].push_back(position);

      //try to inject another sequence
      noise = noiseBackup;

      for (int retry = 0; retry <= retries; retry++) {

        if (retry == retries) {

          std::cerr << "ERROR: Cannot add another set motif subsequence!" <<
            " Retry or change your settings!" << std::endl;
          throw(EXIT_FAILURE);
        }

        //new random noise for retry
        std::normal_distribution<double> distNoiseRetry(0.0, abs(noise) * 0.5
            <= 0.0 ? std::numeric_limits<double>::min() : noise * 0.5);

        //copy another motif sequence ...
        rseq newSubsequence(motif);

        //... stretch and add noise
        stretch = distStretch(randomEngine);

        for (auto &item : newSubsequence) {

          item += distNoiseRetry(randomEngine);
          item *= stretch;
        }

        //check if the sequence is within range d of the motif
        if (similarity(motif, newSubsequence, d) < d) {

          //backup subsequence at position
          for (int i = 0; i < window; i++)
            subsequence[i] = timeSeries_out[position + i];

          value = subsequence[0];

          //make sure we are in maxi when maxi can handle the motif height
          if (method > 3
              && method < 8
              && abs(height) <= 2 * maxi) {

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

          if (!largerMotifSet(timeSeries_out, pos_out[0], d)) {

            break;
          }

          //restore old subsequence
          for (int i = 0; i < window; i++)
            timeSeries_out[position + i] = subsequence[i];

          //update the running mean and variance
          updateRunnings(timeSeries_out, position);

          //get new random position in the synthetic time series
          position = freePositions.calculateRandomPosition();

          pos_out[0].back() = position;
        }

        noise -= lth;
      }

      //remove the position from available positions
      freePositions.removePosition();
    }

    //inject smaller set motif to harden the algorithm
    position = freePositions.calculateRandomPosition();

    retries = size * 5;
    tsg::iseq secPos;
    tsg::rseq sec;
    secPos.push_back(position);

    for (int i = position; i < position + window; i++)
      sec.push_back(timeSeries_out[i] - timeSeries_out[position]);

    //remove the position from available positions
    freePositions.removePosition();

    int secSize = largerMotifSet(timeSeries_out, secPos, pos_out[0].size()
        - 1, d);

    if (secSize > (int)pos_out[0].size()) {

      std::cerr << "ERROR: Impossible error?" << std::endl;
      throw(EXIT_FAILURE);
    }

    //there is already a hidden motif with this subseqeunce
    if (secSize >= (int)pos_out[0].size() - 1)
      return;

    //harden the time series by injecting smaller motif
    for (int motifItr = 1; motifItr < size - 1; motifItr++) {

      //compute the random position for the subsequence
      position = freePositions.calculateRandomPosition();

      secPos.push_back(position);

      //try to inject another sequence
      for (int retry = 0; retry <= retries; retry++) {

        if (retry == retries) {

          std::cerr << "ERROR: Cannot add smaller motif set subsequence!" <<
            " Retry or change your settings!" << std::endl;
          throw(EXIT_FAILURE);
        }

        //backup subsequence at position
        for (int i = 0; i < window; i++)
          subsequence[i] = timeSeries_out[position + i];

        value = subsequence.back() - subsequence[0];
        value = subsequence[0] + value / 2.0;

        //inject sequence into the time series
        for (int i = 0; i < window; i++)
          timeSeries_out[i + position] = value + sec[i];

        //update the running sum and sum of square
        updateRunnings(timeSeries_out, position);

        secSize = largerMotifSet(timeSeries_out, secPos, pos_out[0].size()
            - 1, d);

        if (secSize <= (int)pos_out[0].size()) {

          if (secSize > (int)pos_out[0].size() - 1)
            return;

          break;
        }

        //restore old subsequence
        for (int i = 0; i < window; i++)
          timeSeries_out[position + i] = subsequence[i];

        //update the running mean and variance
        updateRunnings(timeSeries_out, position);

        //get new random position in the synthetic time series
        position = freePositions.calculateRandomPosition();

        pos_out[0].back() = position;
      }

      //remove the position from available positions
      freePositions.removePosition();
    }
  }

  void TSGenerator::injectLatentMotif(rseq &timeSeries_out, rseqs &motif_out,
      rseq &d_out, iseqs &pos_out) {

    //empty motifs
    if (!motif_out.empty()) {

      for (auto &motif : motif_out)
        motif.clear();

      motif_out.clear();
      motif_out.resize(1);
      motif_out[0].resize(window);
    }
    else {

      motif_out.resize(1);
      motif_out[0].resize(window);
    }

    //declaration stuff
    int positionOne = -1;
    int positionTwo = -1;
    rseq motifCenter;
    rseq subsequence(window, 0.0);
    double mean, stdDev;
    double d = std::numeric_limits<double>::max();
    std::normal_distribution<double> distNoise(0.0, abs(noise) * 0.5 <= 0.0
        ? std::numeric_limits<double>::min() : noise * 0.5);
    std::uniform_real_distribution<double> distStretch(0.6, 1.4);
    double stretch;

    //generate a base time series
    generateBaseTimeSeries(timeSeries_out);

    //compute running mean and std dev
    calcRunnings(timeSeries_out);

    //determine similarity of the top motif pair in the random synthetic time
    //series
    tpm(timeSeries_out, sums, sumSquares, positionOne, positionTwo, window);
    d = 0.49999999 * similarity(timeSeries_out, positionOne, positionTwo, d);

    d_out.push_back(d);

    //calculate motif set center subsequence in the window size dimentional
    //room of subsequence values
    calculateSubsequence(motifCenter, type, height);

    for (auto &item : motifCenter)
      item += distNoise(randomEngine);

    mean = 0.0;
    stdDev = 0.0;
    meanStdDev(motifCenter, mean, stdDev);

    motif_out[0] = motifCenter;

    //declaration stuff
    int position;
    double value;
    double min, max;
    int retries = 20;
    double lth = noise / retries;
    double noiseBackup = noise;


    //inject sequences into the time series
    for (int motifItr = 0; motifItr < size; motifItr++) {

      //compute the random position for the subsequence
      position = freePositions.calculateRandomPosition();

      pos_out[0].push_back(position);

      //try to inject another sequence
      noise = noiseBackup;

      for (int retry = 0; retry <= retries; retry++) {

        if (retry == retries) {

          std::cerr << "ERROR: Cannot add another latent motif subsequence!" <<
            " Retry or change your settings!" << std::endl;
          throw(EXIT_FAILURE);
        }

        //new random noise
        std::normal_distribution<double> distNoiseRetry(0.0, abs(noise) * 0.5
            <= 0.0 ? std::numeric_limits<double>::min() : noise * 0.5);

        //copy another motif sequence ...
        rseq newSubsequence(motifCenter);

        //... stretch and add noise
        stretch = distStretch(randomEngine);

        for (auto &item : newSubsequence) {

          item += distNoiseRetry(randomEngine);
          item *= stretch;
        }

        //check if the sequence is within range d of the center
        if (similarity(motifCenter, newSubsequence, d) < d) {

          //backup subsequence at position
          for (int i = 0; i < window; i++)
            subsequence[i] = timeSeries_out[position + i];

          value = subsequence[0];

          //make sure we are in maxi when maxi can handle the motif height
          if (method > 3
              && method < 8
              && abs(height) <= 2 * maxi) {

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

          if (!largerMotifSet(timeSeries_out, pos_out[0], 2.0 * d)) {

            break;
          }

          //restore old subsequence
          for (int i = 0; i < window; i++)
            timeSeries_out[position + i] = subsequence[i];

          //update the running mean and variance
          updateRunnings(timeSeries_out, position);

          //get new random position in the synthetic time series
          position = freePositions.calculateRandomPosition();

          pos_out[0].back() = position;
        }

        noise -= lth;
      }

      //remove the position from available positions
      freePositions.removePosition();
    }

    //inject smaller set motif to harden the algorithm
    position = freePositions.calculateRandomPosition();

    retries = size * 5;
    tsg::iseq secPos;
    tsg::rseq sec;
    secPos.push_back(position);

    for (int i = position; i < position + window; i++)
      sec.push_back(timeSeries_out[i] - timeSeries_out[position]);

    //remove the position from available positions
    freePositions.removePosition();

    int secSize = largerMotifSet(timeSeries_out, secPos, pos_out[0].size()
        - 1, 2.0 * d);

    if (secSize > (int)pos_out[0].size()) {

      std::cerr << "ERROR: Impossible error?" << std::endl;
      throw(EXIT_FAILURE);
    }

    //there is already a hidden motif with this subseqeunce
    if (secSize >= (int)pos_out[0].size() - 1)
      return;

    //harden the time series by injecting smaller motif
    for (int motifItr = 1; motifItr < size - 1; motifItr++) {

      //compute the random position for the subsequence
      position = freePositions.calculateRandomPosition();

      secPos.push_back(position);

      //try to inject another sequence
      for (int retry = 0; retry <= retries; retry++) {

        if (retry == retries) {

          std::cerr << "ERROR: Cannot add smaller motif set subsequence!" <<
            " Retry or change your settings!" << std::endl;
          throw(EXIT_FAILURE);
        }

        //backup subsequence at position
        for (int i = 0; i < window; i++)
          subsequence[i] = timeSeries_out[position + i];

        value = subsequence.back() - subsequence[0];
        value = subsequence[0] + value / 2.0;

        //inject sequence into the time series
        for (int i = 0; i < window; i++)
          timeSeries_out[i + position] = value + sec[i];

        //update the running sum and sum of square
        updateRunnings(timeSeries_out, position);

        secSize = largerMotifSet(timeSeries_out, secPos, pos_out[0].size()
            - 1, 2.0 * d);

        if (secSize <= (int)pos_out[0].size()) {

          if (secSize > (int)pos_out[0].size() - 1)
            return;

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

          pos_out[0].back() = position;
        }
      }

      //remove the position from available positions
      freePositions.removePosition();
    }
  }

  void TSGenerator::run(rseq &timeSeries_out, rseqs &motif_out, rseq &d_out,
      iseqs &pos_out) {

    //clear the output buffers
    if (!d_out.empty()) {

      d_out.clear();
      d_out.resize(0);
    }

    if (!pos_out.empty()) {

      pos_out.clear();
      pos_out.resize(gen ? 2 : 1);
    }
    else
      pos_out.resize(gen ? 2 : 1);


    //choose correct generator
    switch (gen) {

      case 0:
        injectPairMotif(timeSeries_out, motif_out);
        break;
      case 1:
        injectSetMotif(timeSeries_out, motif_out, d_out, pos_out);
        break;
      case 2:
        injectLatentMotif(timeSeries_out, motif_out, d_out, pos_out);
        break;
      default:
        std::cerr << "ERROR: Unknown motif generator: " << gen << std::endl;
        throw(EXIT_FAILURE);
    }


    //add top motif pair to output
    int positionOne = 0;
    int positionTwo = 0;

    tpm(timeSeries_out, sums, sumSquares, positionOne, positionTwo, window);
    d_out.push_back(similarity(timeSeries_out, positionOne, positionTwo,
        std::numeric_limits<double>::max()));

    pos_out[gen ? 1 : 0].push_back(positionOne);
    pos_out[gen ? 1 : 0].push_back(positionTwo);
  }
}
