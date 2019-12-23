///\file unittests.cpp
///
///\brief File contains all unit test.
///
///This is the source file to run all unit tests with stest.

#include <unittests.hpp>


void test_motifsetcollection() {

  TEST_GROUP_FUNCTION;

  vector<double> subsequence;

  generateBoxMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[5] == 5.0);
  TEST(subsequence[8] == 5.0);
  generateBoxMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[5] == -7.0);
  TEST(subsequence[8] == -7.0);
  generateTriangleMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 < subsequence[5]);
  TEST(subsequence[4] == subsequence[5]);
  TEST(subsequence[1] == subsequence[8]);
  TEST(subsequence[3] != subsequence[9]);
  generateTriangleMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[4] < 0.0);
  TEST(subsequence[4] == subsequence[5]);
  TEST(subsequence[1] == subsequence[8]);
  TEST(subsequence[3] != subsequence[9]);
  generateSemicircleMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 < subsequence[5]);
  TEST(subsequence[4] == subsequence[5]);
  TEST(subsequence[1] == subsequence[8]);
  TEST(subsequence[3] != subsequence[9]);
  generateSemicircleMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[4] < 0.0);
  TEST(subsequence[4] == subsequence[5]);
  TEST(subsequence[1] == subsequence[8]);
  TEST(subsequence[3] != subsequence[9]);
  generateTrapezoidMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 < subsequence[5]);
  TEST(subsequence[4] == subsequence[5]);
  TEST(subsequence[1] == subsequence[8]);
  TEST(subsequence[3] != subsequence[9]);
  generateTrapezoidMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[4] < 0.0);
  TEST(subsequence[4] == subsequence[5]);
  TEST(subsequence[1] == subsequence[8]);
  TEST(subsequence[3] != subsequence[9]);
  generatePositiveFlankMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 < subsequence[5]);
  TEST(subsequence[2] < subsequence[3]);
  TEST(subsequence[5] < subsequence[6]);
  TEST(subsequence[9] < subsequence[8]);
  generatePositiveFlankMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[4] < 0.0);
  TEST(subsequence[2] > subsequence[3]);
  TEST(subsequence[5] > subsequence[6]);
  TEST(subsequence[9] > subsequence[8]);
  generateNegativeFlankMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 < subsequence[5]);
  TEST(subsequence[2] > subsequence[3]);
  TEST(subsequence[5] > subsequence[6]);
  TEST(subsequence[0] < subsequence[1]);
  generateNegativeFlankMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[4] < 0.0);
  TEST(subsequence[2] < subsequence[3]);
  TEST(subsequence[5] < subsequence[6]);
  TEST(subsequence[0] > subsequence[1]);
  generateSineMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 < subsequence[3]);
  TEST(subsequence[8] < 0.0);
  TEST(subsequence[1] > subsequence[7]);
  TEST(subsequence[1] < subsequence[3]);
  TEST(abs(subsequence[1] + subsequence[8]) < 0.001);
  TEST(abs(subsequence[3] + subsequence[6]) < 0.001);
  generateSineMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 > subsequence[3]);
  TEST(subsequence[8] > 0.0);
  TEST(subsequence[1] < subsequence[7]);
  TEST(subsequence[1] > subsequence[3]);
  TEST(abs(subsequence[1] + subsequence[8]) < 0.001);
  TEST(abs(subsequence[3] + subsequence[6]) < 0.001);
  generateCosineMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[4] < 0.0);
  TEST(abs(subsequence[4] - subsequence[5]) < 0.001);
  TEST(abs(subsequence[3] - subsequence[6]) < 0.001);
  TEST(subsequence[3] != subsequence[9]);
  generateCosineMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 < subsequence[5]);
  TEST(abs(subsequence[4] - subsequence[5]) < 0.001);
  TEST(abs(subsequence[3] - subsequence[6]) < 0.001);
  TEST(subsequence[3] != subsequence[9]);
}

void test_freepositions() {

  TEST_GROUP_FUNCTION;

  { //test repeating without remove

    auto length = 10;
    auto window = 3;
    FreePositions freePos(length, window);

    for (int x = 0; x < 300; x++) {

      int pos = freePos.calculateRandomPosition();

      TEST(0 <= pos);
      TEST(pos <= length - window);
    }
  }

  //test normal behaviour
  for (int x = 0; x < 20; x++) {
    auto length = 44;
    auto window = 3;
    FreePositions freePos(length, window);
    vector<int> pos;

    for (int i = 0; i < 4; i++) {

      pos.push_back(freePos.calculateRandomPosition());
      freePos.removePosition();

      TEST(0 <= pos[i]);
      TEST(pos[i] <= length - window);

      for (int j = 0; j < (int) pos.size() - 1; j++) {

        TEST(abs(pos[j] - pos[i]) >= 2 * window);
      }
    }
  }

  { //test failure
    auto length = 3;
    auto window = 3;
    FreePositions freePos(length, window);

    auto pos = freePos.calculateRandomPosition();
    freePos.removePosition();

    TEST(0 <= pos);
    TEST(pos <= length - window);

    try {

      pos = freePos.calculateRandomPosition();
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }
  }
}

void test_basets() {

  TEST_GROUP_FUNCTION;

  BaseTS baseTS;

  vector<double> ts;

  //simple random walk
  baseTS.simpleRandomWalk(ts, 100, 1.0, 0.0);
  TEST(ts.size() == 100);

  baseTS.simpleRandomWalk(ts, 23, 1.0, 0.1);
  TEST(ts.size() == 23);

  baseTS.simpleRandomWalk(ts, 0, 1.0, 0.1);
  TEST(ts.size() == 0);

  //real random walk
  baseTS.realRandomWalk(ts, 100, 1.0, 0.0);
  TEST(ts.size() == 100);

  baseTS.realRandomWalk(ts, 23, 1.0, 0.1);
  TEST(ts.size() == 23);

  baseTS.realRandomWalk(ts, 0, 1.0, 0.1);
  TEST(ts.size() == 0);

  //normal random walk
  baseTS.normalRandomWalk(ts, 100, 1.0, 0.0);
  TEST(ts.size() == 100);

  baseTS.normalRandomWalk(ts, 23, 1.0, 0.1);
  TEST(ts.size() == 23);

  //linear random walk
  baseTS.linearRandomWalk(ts, 100, 1.0, 3.0, 0.0);
  TEST(ts.size() == 100);

  baseTS.linearRandomWalk(ts, 23, 1.0, 3.0, 0.1);
  TEST(ts.size() == 23);

  baseTS.linearRandomWalk(ts, 0, 1.0, 3.0, 0.1);
  TEST(ts.size() == 0);

  //bounded simple random walk
  baseTS.simpleRandomWalk(ts, 100, 1.0, 0.0);
  TEST(ts.size() == 100);

  baseTS.simpleRandomWalk(ts, 23, 1.0, 0.1);
  TEST(ts.size() == 23);

  baseTS.simpleRandomWalk(ts, 0, 1.0, 0.1);
  TEST(ts.size() == 0);

  //bounded real random walk
  baseTS.realRandomWalk(ts, 100, 1.0, 20.0, 0.0);
  TEST(ts.size() == 100);

  baseTS.realRandomWalk(ts, 23, 1.0, 20.0, 0.1);
  TEST(ts.size() == 23);

  baseTS.realRandomWalk(ts, 0, 1.0, 20.0, 0.1);
  TEST(ts.size() == 0);

  //bounded normal random walk
  baseTS.normalRandomWalk(ts, 100, 1.0, 20.0, 0.0);
  TEST(ts.size() == 100);

  baseTS.normalRandomWalk(ts, 23, 1.0, 20.0, 0.1);
  TEST(ts.size() == 23);

  baseTS.normalRandomWalk(ts, 0, 1.0, 20.0, 0.1);
  TEST(ts.size() == 0);

  //bounded linear random walk
  baseTS.linearRandomWalk(ts, 100, 1.0, 3.0, 20.0, 0.0);
  TEST(ts.size() == 100);

  baseTS.linearRandomWalk(ts, 23, 1.0, 3.0, 20.0, 0.1);
  TEST(ts.size() == 23);

  baseTS.linearRandomWalk(ts, 0, 1.0, 3.0, 20.0, 0.1);
  TEST(ts.size() == 0);

  //uniform random
  baseTS.uniformRandom(ts, 100, 1.0, 0.0);
  TEST(ts.size() == 100);

  baseTS.uniformRandom(ts, 23, 1.0, 0.1);
  TEST(ts.size() == 23);

  baseTS.uniformRandom(ts, 0, 1.0, 0.1);
  TEST(ts.size() == 0);

  //normal random
  baseTS.normalRandom(ts, 100, 1.0, 0.0);
  TEST(ts.size() == 100);

  baseTS.normalRandom(ts, 23, 1.0, 0.1);
  TEST(ts.size() == 23);

  baseTS.normalRandom(ts, 0, 1.0, 0.1);
  TEST(ts.size() == 0);

  //piecewise linear random
  baseTS.piecewiseLinearRandom(ts, 100, 1.0, 0.0);
  TEST(ts.size() == 100);

  baseTS.piecewiseLinearRandom(ts, 23, 1.0, 0.1);
  TEST(ts.size() == 23);

  baseTS.piecewiseLinearRandom(ts, 0, 1.0, 0.1);
  TEST(ts.size() == 0);

  //spline repeated
  baseTS.splineRepeated(ts, 1000, 20.0, 20.0, 10, 2.0);
  TEST(ts.size() == 1000);

  baseTS.splineRepeated(ts, 23, 20.0, 20.0, 10, 2.0);
  TEST(ts.size() == 23);

  baseTS.splineRepeated(ts, 0, 20.0, 20.0, 10, 2.0);
  TEST(ts.size() == 0);

}

void test_tsgenerator() {

  TEST_GROUP_FUNCTION;

  //redirect cerr to keep test output clean
  ostringstream local;

  auto cerr_buff = cerr.rdbuf(local.rdbuf());

  //tests for wrong input
  try {

    TestTSGenerator generator(300, 20, 20.0, 0.0, 2, 3, 50.0);
    TEST("Should work without throw!");
  }
  catch (...) {

    TEST(!"Should work without throw!");
  }

  try {

    TestTSGenerator generator(300, 20, 20.0, 0.0, 9, 3, 50.0);
    TEST(!"Has to throw an error!");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    TestTSGenerator generator(300, 20, 20.0, 0.0, 2, 2, 50.0);
    TEST(!"Has to throw an error!");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    TestTSGenerator generator(99, 20, 20, 0.0, 9, 3, 50.0);
    TEST(!"Has to throw an error!");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    TestTSGenerator generator(300, 500, 20, 0.0, 9, 3, 50.0);
    TEST(!"Has to throw an error!");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    TestTSGenerator generator(300, -1, 20, 0.0, 9, 3, 50.0);
    TEST(!"Has to throw an error!");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    TestTSGenerator generator(300, 20, 20.0, -1.0, 9, 3, 50.0);
    TEST(!"Has to throw an error!");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    TestTSGenerator generator(300, 20, 20.0, 0.01, 9, -1, 50.0);
    TEST(!"Has to throw an error!");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    TestTSGenerator generator(300, 20, 20.0, 0.0, 9, 20, 50.0);
    TEST(!"Has to throw an error!");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }


  {
    //test mean and standard deviation function
    TestTSGenerator generator(300, 20, 20.0, 0.0, 2, 3, 50.0);

    double mean;
    double stdDev;
    double rWindow = 0.1;

    generator.setWindow(10);

    try {

      generator.testMeanStdDev({}, mean, stdDev);
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    generator.testMeanStdDev({ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0 }, mean, stdDev);
    TEST(mean <= numeric_limits<double>::min() && mean >=
        -numeric_limits<double>::min());
    TEST(stdDev == 1.0);

    generator.testMeanStdDev({ 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
        3.0 }, mean, stdDev);
    TEST(mean <= 3.0 + numeric_limits<double>::min() && mean >= 3.0
        - numeric_limits<double>::min());
    TEST(stdDev == 1.0);

    generator.testMeanStdDev({ 3.0, -3.0, 1.75, 9.0, 33.0, 3.101, 0.03,
        3.99999, 4.23, -7.093 }, mean, stdDev);
    TEST(mean <= 4.801799 + 0.0000001 && mean >= 4.801799 - 0.0000001);
    TEST(stdDev <= 10.2679615 + 0.0000001 && stdDev >= 10.2679615 - 0.0000001);

    rWindow = 0.05;
    generator.setWindow(20);
    vector<double> subsequence(testTimeSeries.begin() + topMotifSetPos[0],
        testTimeSeries.begin() + topMotifSetPos[0] + 20);
    generator.testMeanStdDev(subsequence, mean, stdDev);
    TEST(mean <= topMotifSetMeans[0] + 0.0000001 && mean >= topMotifSetMeans[0]
        - 0.0000001);
    TEST(stdDev <= topMotifSetStdDevs[0] + 0.0000001 && stdDev >=
        topMotifSetStdDevs[0] - 0.0000001);

    vector<double> subsequence1(testTimeSeries.begin() + topMotifSetPos[1],
        testTimeSeries.begin() + topMotifSetPos[1] + 20);
    generator.testMeanStdDev(subsequence1, mean, stdDev);
    TEST(mean <= topMotifSetMeans[1] + 0.0000001 && mean >= topMotifSetMeans[1]
        - 0.0000001);
    TEST(stdDev <= topMotifSetStdDevs[1] + 0.0000001 && stdDev >=
        topMotifSetStdDevs[1] - 0.0000001);

    vector <double> subsequence2(testTimeSeries.begin() + topMotifSetPos[2],
        testTimeSeries.begin() + topMotifSetPos[2] + 20);
    generator.testMeanStdDev(subsequence2, mean, stdDev);
    TEST(mean <= topMotifSetMeans[2] + 0.0000001 && mean >= topMotifSetMeans[2]
        - 0.0000001);
    TEST(stdDev <= topMotifSetStdDevs[2] + 0.0000001 && stdDev >=
        topMotifSetStdDevs[2] - 0.0000001);


    //test running sum and sum of squares
    vector<double> sums;
    vector<double> sumSquares;
    double std;
    bool flagMeansStdsTest = true;

    generator.testCalcRunnings(testTimeSeries);
    sums = generator.getSums();
    sumSquares = generator.getSumSquares();

    if (TEST(testMeans.size() == sums.size())) {

      for (int itr = 0; itr < (int)sums.size(); itr++)
        if (!(testMeans[itr] - 0.000001 <= sums[itr] * rWindow  &&
              testMeans[itr] + 0.000001 >= sums[itr] * rWindow)) {

          flagMeansStdsTest = false;
          break;
        }

      TEST(flagMeansStdsTest);
    }

    flagMeansStdsTest = true;

    if (TEST(testStds.size() == sumSquares.size())) {

      for (int itr = 0; itr < (int)sumSquares.size(); itr++) {

        std = sums[itr] * rWindow;
        std = sqrt(sumSquares[itr] * rWindow - std * std);

        if (!(testStds[itr] - 0.000001 <= std &&
              testStds[itr] + 0.000001 >= std)) {

          flagMeansStdsTest = false;
          break;
        }
      }

      TEST(flagMeansStdsTest);
    }

    //test update running sum and sum of square
    {
      int pos = 59;
      vector<double> testSeries(testTimeSeries);

      //compute running sum and sum of square
      generator.testCalcRunnings(testSeries);

      //change the time series
      for (int i = pos; i < pos + 20; i++)
        testSeries[i] = 0;

      //update runnings
      generator.testUpdateRunnings(testSeries, pos);

      //update runnings again with original time series
      generator.testUpdateRunnings(testTimeSeries, pos);

      sums = generator.getSums();
      sumSquares = generator.getSumSquares();
      flagMeansStdsTest = true;

      if (TEST(testMeans.size() == sums.size())) {

        for (int itr = 0; itr < (int)sums.size(); itr++)
          if (!(testMeans[itr] - 0.000001 <= sums[itr] * rWindow  &&
                testMeans[itr] + 0.000001 >= sums[itr] * rWindow)) {

            flagMeansStdsTest = false;
            break;
          }

        TEST(flagMeansStdsTest);
      }

      flagMeansStdsTest = true;

      if (TEST(testStds.size() == sumSquares.size())) {

        for (int itr = 0; itr < (int)sumSquares.size(); itr++) {

          std = sums[itr] * rWindow;
          std = sqrt(sumSquares[itr] * rWindow - std * std);

          if (!(testStds[itr] - 0.000001 <= std &&
                testStds[itr] + 0.000001 >= std)) {

            flagMeansStdsTest = false;
            break;
          }
        }

        TEST(flagMeansStdsTest);
      }
    }


    //test the first similarity functions
    try {

      generator.testSimilarity({}, topMotifSetPos[0], topMotifSetPos[0],
          numeric_limits<double>::max());
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    try {

      generator.testSimilarity(testTimeSeries, 303, topMotifSetPos[0],
          numeric_limits<double>::max());
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    try {

      generator.testSimilarity(testTimeSeries, topMotifSetPos[0], 303,
          numeric_limits<double>::max());
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    generator.testCalcRunnings(testTimeSeries);

    double similarity = generator.testSimilarity(testTimeSeries,
        topMotifSetPos[0], topMotifSetPos[1], numeric_limits<double>::max());
    TEST(abs(similarity - testSequencesSimilarty) <= 0.0000001);

    similarity = generator.testSimilarity(testTimeSeries, topMotifSetPos[1],
        topMotifSetPos[2], numeric_limits<double>::max());
    TEST(abs(similarity) <= 0.0000001);

    similarity = generator.testSimilarity(testTimeSeries, topMotifSetPos[0],
        topMotifSetPos[2], 0.01);
    TEST(similarity >= 0.01);

    similarity = generator.testSimilarity(testTimeSeries, topMotifSetPos[1],
        topMotifSetPos[1], 0.01);
    TEST(abs(similarity) <= 0.0000001);


    //test the second similarity functions
    try {

      generator.testSimilarity({}, {}, topMotifSetMeans[0],
          topMotifSetStdDevs[0], topMotifSetPos[0],
          numeric_limits<double>::max());
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    try {

      generator.testSimilarity(testTimeSeries, vector<double>(),
          topMotifSetMeans[0], topMotifSetStdDevs[0], topMotifSetPos[0],
          numeric_limits<double>::max());
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    try {

      generator.testSimilarity(testTimeSeries, { 0.0, 0.0 },
          topMotifSetMeans[0], topMotifSetStdDevs[0], topMotifSetPos[0],
          numeric_limits<double>::max());
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    try {

      generator.testSimilarity(testTimeSeries, testSequenceOne,
          topMotifSetMeans[0], topMotifSetStdDevs[0], 303,
          numeric_limits<double>::max());
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    similarity = generator.testSimilarity(testTimeSeries,
        testMotifSetSubsequences[0], topMotifSetMeans[0],
        topMotifSetStdDevs[0], topMotifSetPos[1],
        numeric_limits<double>::max());
    TEST(abs(similarity - testSequencesSimilarty) <= 0.0000001);

    similarity = generator.testSimilarity(testTimeSeries,
        testMotifSetSubsequences[1], topMotifSetMeans[1],
        topMotifSetStdDevs[1], topMotifSetPos[2],
        numeric_limits<double>::max());
    TEST(abs(similarity) <= 0.0000001);

    similarity = generator.testSimilarity(testTimeSeries,
        testMotifSetSubsequences[0], topMotifSetMeans[0],
        topMotifSetStdDevs[0], topMotifSetPos[2], 0.01);
    TEST(similarity >= 0.01);

    similarity = generator.testSimilarity(testTimeSeries,
        testMotifSetSubsequences[1], topMotifSetMeans[1],
        topMotifSetStdDevs[1], topMotifSetPos[1], 0.01);
    TEST(abs(similarity) <= 0.0000001);


    //test calculate raw subsequence
    subsequence.clear();
    subsequence.resize(0);

    generator.testCalculateSubsequence(subsequence, 0, 5.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(subsequence[5] == 5.0);
    TEST(subsequence[18] == 5.0);
    generator.testCalculateSubsequence(subsequence, 0, -7.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(subsequence[5] == -7.0);
    TEST(subsequence[18] == -7.0);
    generator.testCalculateSubsequence(subsequence, 1, 5.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(0.0 < subsequence[5]);
    TEST(subsequence[4] == subsequence[15]);
    TEST(subsequence[1] == subsequence[18]);
    TEST(subsequence[3] != subsequence[19]);
    generator.testCalculateSubsequence(subsequence, 1, -7.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(subsequence[4] < 0.0);
    TEST(subsequence[4] == subsequence[15]);
    TEST(subsequence[1] == subsequence[18]);
    TEST(subsequence[3] != subsequence[19]);
    generator.testCalculateSubsequence(subsequence, 2, 5.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(0.0 < subsequence[5]);
    TEST(subsequence[4] == subsequence[15]);
    TEST(subsequence[1] == subsequence[18]);
    TEST(subsequence[3] != subsequence[19]);
    generator.testCalculateSubsequence(subsequence, 2, -7.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(subsequence[4] < 0.0);
    TEST(subsequence[4] == subsequence[15]);
    TEST(subsequence[1] == subsequence[18]);
    TEST(subsequence[3] != subsequence[19]);
    generator.testCalculateSubsequence(subsequence, 3, 5.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(0.0 < subsequence[5]);
    TEST(subsequence[4] == subsequence[15]);
    TEST(subsequence[1] == subsequence[18]);
    TEST(subsequence[3] != subsequence[19]);
    generator.testCalculateSubsequence(subsequence, 3, -7.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(subsequence[4] < 0.0);
    TEST(subsequence[4] == subsequence[15]);
    TEST(subsequence[1] == subsequence[18]);
    TEST(subsequence[3] != subsequence[19]);
    generator.testCalculateSubsequence(subsequence, 4, 5.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(0.0 < subsequence[15]);
    TEST(subsequence[2] < subsequence[13]);
    TEST(subsequence[5] < subsequence[16]);
    TEST(subsequence[19] < subsequence[18]);
    generator.testCalculateSubsequence(subsequence, 4, -7.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(subsequence[14] < 0.0);
    TEST(subsequence[2] > subsequence[13]);
    TEST(subsequence[5] > subsequence[16]);
    TEST(subsequence[19] > subsequence[18]);
    generator.testCalculateSubsequence(subsequence, 5, 5.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(0.0 < subsequence[15]);
    TEST(subsequence[2] > subsequence[13]);
    TEST(subsequence[5] > subsequence[16]);
    TEST(subsequence[0] < subsequence[1]);
    generator.testCalculateSubsequence(subsequence, 5, -7.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(subsequence[14] < 0.0);
    TEST(subsequence[2] < subsequence[13]);
    TEST(subsequence[5] < subsequence[16]);
    TEST(subsequence[0] > subsequence[1]);
    generator.testCalculateSubsequence(subsequence, 6, 5.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(0.0 < subsequence[3]);
    TEST(subsequence[18] < 0.0);
    TEST(subsequence[1] > subsequence[17]);
    TEST(subsequence[1] < subsequence[3]);
    TEST(abs(subsequence[1] + subsequence[18]) < 0.001);
    TEST(abs(subsequence[3] + subsequence[16]) < 0.001);
    generator.testCalculateSubsequence(subsequence, 6, -7.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(0.0 > subsequence[3]);
    TEST(subsequence[18] > 0.0);
    TEST(subsequence[1] < subsequence[17]);
    TEST(subsequence[1] > subsequence[3]);
    TEST(abs(subsequence[1] + subsequence[18]) < 0.001);
    TEST(abs(subsequence[3] + subsequence[16]) < 0.001);
    generator.testCalculateSubsequence(subsequence, 7, 5.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(subsequence[4] < 0.0);
    TEST(abs(subsequence[4] - subsequence[15]) < 0.001);
    TEST(abs(subsequence[3] - subsequence[16]) < 0.001);
    TEST(subsequence[3] != subsequence[19]);
    generator.testCalculateSubsequence(subsequence, 7, -7.0);
    TEST(subsequence.size() == 20);
    TEST(subsequence.front() == 0.0);
    TEST(subsequence.back() == 0.0);
    TEST(0.0 < subsequence[5]);
    TEST(abs(subsequence[4] - subsequence[15]) < 0.001);
    TEST(abs(subsequence[3] - subsequence[16]) < 0.001);
    TEST(subsequence[3] != subsequence[19]);


    generator.testCalculateSubsequence(subsequence, 1, 10.0);
    generator.testMeanStdDev(subsequence, mean, stdDev);
    vector<double> sequenceTwo;
    double meanTwo;
    double stdDevTwo;
    generator.testCalculateSubsequence(sequenceTwo, 7, 10.0);
    generator.testMeanStdDev(sequenceTwo, meanTwo, stdDevTwo);
    TEST(mean != meanTwo);


    //test the search for unintantional matches
    TEST(!generator.testSearchForUnintentionalMatches(testTimeSeries,
          topMotifSetPos, 2 * topMotifSetRange));
    TEST(generator.testSearchForUnintentionalMatches(testTimeSeries,
          { topMotifSetPos[0], topMotifSetPos[1] }, 2 * topMotifSetRange));


    //test the check for larger latent set motifs
    TEST(!generator.testCheckIfThereIsALargerMotifSet(testTimeSeries,
          topMotifSetPos,  topMotifSetRange));
    TEST(!generator.testCheckIfThereIsALargerMotifSet(testTimeSeries,
          { topMotifSetPos[0], topMotifSetPos[1], 249 }, topMotifSetRange));
    TEST(!generator.testCheckIfThereIsALargerMotifSet(testTimeSeries,
          { topMotifSetPos[0], topMotifSetPos[1], 249, 229 },
          topMotifSetRange));
    TEST(!generator.testCheckIfThereIsALargerMotifSet(testTimeSeries,
          { topMotifSetPos[0], topMotifSetPos[1], 249, 229, 100 },
          topMotifSetRange));
    generator.testCalcRunnings(testTimeSeriesCheckLargerMotifSet);
    TEST(generator.testCheckIfThereIsALargerMotifSet(
          testTimeSeriesCheckLargerMotifSet, { 34, 248, topMotifSetPos[1] },
          topMotifSetRange));
  }


  {
    //test synthetic time series motif set generation function
    TestTSGenerator simGenerator(3000, 100, 20, 0.0, 2, 3, 50.0);
    TestTSGenerator generator(3000, 100, 20, 0.0, 2, 3, 50.0);
    vector<double> timeSeries_out;
    vector<double> d_out;
    vector<int> windowSize_out;
    vector<vector<int>> positions_out;


    for (int itr = 0; itr < 3; itr++) {

      try {

        generator.run(timeSeries_out, d_out, windowSize_out,
            positions_out);
        TEST("Has to run without throwing an error!");

        simGenerator.testCalcRunnings(timeSeries_out);

        for (auto &pos0 : positions_out[0])
          for (auto &pos1 : positions_out[0])
            TEST(simGenerator.testSimilarity(timeSeries_out, pos0, pos1,
                  2 * d_out[0]) <= 2 * d_out[0]);

        TEST(d_out.size() == 2);
        TEST(windowSize_out.size() == 2);
        TEST(positions_out.size() == 2);
      }
      catch (...) {

        TEST(!"Has to run without throwing an error!");
      }
    }
  }


  //reset cerr
  cerr.rdbuf(cerr_buff);
}

int main() {

  try {
    TEST_WAIT(false);

    TEST_SECTION("motif collections");
    test_motifsetcollection();
    TEST_SECTION("free positions data type");
    test_freepositions();
    TEST_SECTION("base time series generation");
    test_basets();
    TEST_SECTION("time series generator");
    test_tsgenerator();

    TEST_SUMMARY;
  }
  catch (...) {
    TEST_EXCEPTION;
  }

  return EXIT_SUCCESS;
}
