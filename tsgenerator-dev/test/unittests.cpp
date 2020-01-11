///\file unittests.cpp
///
///\brief File contains all unit test.
///
///This is the source file to run all unit tests with stest.

#include <unittests.hpp>


void test_motifsetcollection() {

  TEST_GROUP_FUNCTION;

  tsg::rseq subsequence;

  tsg::generateBoxMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[5] == 5.0);
  TEST(subsequence[8] == 5.0);
  tsg::generateBoxMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[5] == -7.0);
  TEST(subsequence[8] == -7.0);
  tsg::generateTriangleMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 < subsequence[5]);
  TEST(subsequence[4] == subsequence[5]);
  TEST(subsequence[1] == subsequence[8]);
  TEST(subsequence[3] != subsequence[9]);
  tsg::generateTriangleMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[4] < 0.0);
  TEST(subsequence[4] == subsequence[5]);
  TEST(subsequence[1] == subsequence[8]);
  TEST(subsequence[3] != subsequence[9]);
  tsg::generateSemicircleMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 < subsequence[5]);
  TEST(subsequence[4] == subsequence[5]);
  TEST(subsequence[1] == subsequence[8]);
  TEST(subsequence[3] != subsequence[9]);
  tsg::generateSemicircleMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[4] < 0.0);
  TEST(subsequence[4] == subsequence[5]);
  TEST(subsequence[1] == subsequence[8]);
  TEST(subsequence[3] != subsequence[9]);
  tsg::generateTrapezoidMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 < subsequence[5]);
  TEST(subsequence[4] == subsequence[5]);
  TEST(subsequence[1] == subsequence[8]);
  TEST(subsequence[3] != subsequence[9]);
  tsg::generateTrapezoidMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[4] < 0.0);
  TEST(subsequence[4] == subsequence[5]);
  TEST(subsequence[1] == subsequence[8]);
  TEST(subsequence[3] != subsequence[9]);
  tsg::generatePositiveFlankMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 < subsequence[5]);
  TEST(subsequence[2] < subsequence[3]);
  TEST(subsequence[5] < subsequence[6]);
  TEST(subsequence[9] < subsequence[8]);
  tsg::generatePositiveFlankMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[4] < 0.0);
  TEST(subsequence[2] > subsequence[3]);
  TEST(subsequence[5] > subsequence[6]);
  TEST(subsequence[9] > subsequence[8]);
  tsg::generateNegativeFlankMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 < subsequence[5]);
  TEST(subsequence[2] > subsequence[3]);
  TEST(subsequence[5] > subsequence[6]);
  TEST(subsequence[0] < subsequence[1]);
  tsg::generateNegativeFlankMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[4] < 0.0);
  TEST(subsequence[2] < subsequence[3]);
  TEST(subsequence[5] < subsequence[6]);
  TEST(subsequence[0] > subsequence[1]);
  tsg::generateSineMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 < subsequence[3]);
  TEST(subsequence[8] < 0.0);
  TEST(subsequence[1] > subsequence[7]);
  TEST(subsequence[1] < subsequence[3]);
  TEST(abs(subsequence[1] + subsequence[8]) < 0.001);
  TEST(abs(subsequence[3] + subsequence[6]) < 0.001);
  tsg::generateSineMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(0.0 > subsequence[3]);
  TEST(subsequence[8] > 0.0);
  TEST(subsequence[1] < subsequence[7]);
  TEST(subsequence[1] > subsequence[3]);
  TEST(abs(subsequence[1] + subsequence[8]) < 0.001);
  TEST(abs(subsequence[3] + subsequence[6]) < 0.001);
  tsg::generateCosineMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST(subsequence.size() == 10);
  TEST(subsequence.front() == 0.0);
  TEST(subsequence.back() == 0.0);
  TEST(subsequence[4] < 0.0);
  TEST(abs(subsequence[4] - subsequence[5]) < 0.001);
  TEST(abs(subsequence[3] - subsequence[6]) < 0.001);
  TEST(subsequence[3] != subsequence[9]);
  tsg::generateCosineMotif(subsequence, 1.0, 1.0, 10, -7.0);
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
    tsg::FreePositions freePos(length, window);

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
    tsg::FreePositions freePos(length, window);
    tsg::iseq pos;

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
    tsg::FreePositions freePos(length, window);

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

  tsg::BaseTS baseTS;

  tsg::rseq ts;

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

void test_tsm() {

  TEST_GROUP_FUNCTION;

  tsg::rseq sums;
  tsg::rseq sumSquares;
  TestTSGenerator generator(300, 20, 20.0, 0.0, 2, 3, 50.0);

  generator.testCalcRunnings(testTimeSeries);
  sums = generator.getSums();
  sumSquares = generator.getSumSquares();

  tsg::TSM tsm(testTimeSeries, sums, sumSquares);

  //test the z-normalized PAA implementation
  tsg::rseq paa;
  tsg::rseq paaGT = { 0.1850, 0.1143, -0.1378, -0.0328, 0.1385, 0.0957 };

  try {

    tsm.zNormalPAA(paa, 56, 20, 6);

    for (int i = 0; i < (int)paaGT.size(); i++)
      TEST(paa[i] < paaGT[i] + 0.0001 && paa[i] > paaGT[i] - 0.0001);
  }
  catch (...) {

    TEST(!"Shouldn't throw an error!");
  }

  //test the invNormalCDF implementation
  tsg::rseq bGT = { -std::numeric_limits<double>::infinity(), -0.8416, -0.2533,
    0.2533, 0.8416, std::numeric_limits<double>::infinity() };

  double cdf = tsm.invNormalCDF(1.0 / 5.0, 0.000001);

  TEST(cdf < bGT[1] && cdf > bGT[1] - 0.0001);

  tsg::rseq b;

  try {

    tsm.invNormalCDF(5, b, 0.000001);

    TEST(b[0] == bGT[0]);
    TEST(b[5] == bGT[5]);

    for (int i = 1; i < (int)bGT.size() - 1; i++)
      TEST(b[i] < bGT[i] + 0.0001 && b[i] > bGT[i] - 0.0001);
  }
  catch (...) {

    TEST(!"Shouldn't throw an error!");
  }

  //test the z-normalized SAX implemenation
  tsg::word sax;

  try {

    tsm.zNormalSAX(sax, 148, 20, 6, bGT);
  }
  catch (...) {

    TEST(!"Shouldn't throw an error!");
  }

  //test the SAX minimum distance function
  double dist = -1.0;

  try {

    dist = tsm.saxDist({ 1, 1, 1, 1, 1 }, { 3, 3, 3, 3, 3}, 20, b);
    TEST(dist == sqrt(20.0 / 5.0 * 5.0 * (b[1] - b[2]) * (b[1] - b[2])));

    dist = tsm.saxDist({ 1, 1, 1, 1, 1 }, { 2, 2, 2, 2, 2}, 20, b);
    TEST(dist == 0.0);
  }
  catch (...) {

    TEST(!"Shouldn't throw an error!");
  }

  //test the z-normalized Euclidean distance function
  dist = tsm.dist(topMotifSetPos[0], topMotifSetPos[1], 20.0);
  TEST(abs(dist - testSequencesSimilarty) <= 0.0000001);

  dist = tsm.dist(topMotifSetPos[1], topMotifSetPos[2], 20.0);
  TEST(abs(dist) <= 0.0000001);

  //test the ADM implementation
  tsg::rseqs d;
  tsg::rseqs dGT = {
    { 0, 7.43966, 6.2613, 6.24954, 6.50969, 5.62853 },
    { 7.43966, 0, 4.00941, 5.40625, 4.59849, 7.11225 },
    { 6.2613, 4.00941, 0, 4.70289, 2.52727, 5.94591 },
    { 6.24954, 5.40625, 4.70289, 0, 4.4723, 6.72786 },
    { 6.50969, 4.59849, 2.52727, 4.4723, 0, 5.5463 },
    { 5.62853, 7.11225, 5.94591, 6.72786, 5.5463, 0 }
  };

  try {

    tsm.adm({ 4, 32, 86, 111, 140, 268 }, 20.0, 7.0, d);
    for (int i = 0; i < (int)dGT.size(); i++)
      for (int j = 0; j < (int)dGT.size(); j++)
        TEST(d[i][j] < dGT[i][j] + 0.0001 && d[i][j] > dGT[i][j] - 0.0001);
  }
  catch (...) {

    TEST(!"Shouldn't throw an error!");
  }

  //test the top set motif discovery implementation
  tsg::iseq motif;

  try {

    TEST(tsm.tsm(motif, 20, 0.5) == 3);
    TEST(tsm.tsm(motif, 20, 0.001) == 2);
    tsm.tsm(motif, 20, 6.0);
    for (int i = 0; i < (int)motif.size(); i++)
      for (int j = 0; j < (int)motif.size(); j++)
        if (i != j)
          TEST(abs(motif[i] - motif[j]) >= 20);
  }
  catch (...) {

    TEST(!"Shouldn't throw an error!");
  }
}

void test_tsgenerator() {

  TEST_GROUP_FUNCTION;

  //redirect cerr to keep test output clean
  std::ostringstream local;

  auto cerr_buff = std::cerr.rdbuf(local.rdbuf());

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
    TEST(mean <= std::numeric_limits<double>::min() && mean >=
        -std::numeric_limits<double>::min());
    TEST(stdDev == 1.0);

    generator.testMeanStdDev({ 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
        3.0 }, mean, stdDev);
    TEST(mean <= 3.0 + std::numeric_limits<double>::min() && mean >= 3.0
        - std::numeric_limits<double>::min());
    TEST(stdDev == 1.0);

    generator.testMeanStdDev({ 3.0, -3.0, 1.75, 9.0, 33.0, 3.101, 0.03,
        3.99999, 4.23, -7.093 }, mean, stdDev);
    TEST(mean <= 4.801799 + 0.0000001 && mean >= 4.801799 - 0.0000001);
    TEST(stdDev <= 10.2679615 + 0.0000001 && stdDev >= 10.2679615 - 0.0000001);

    rWindow = 0.05;
    generator.setWindow(20);
    tsg::rseq subsequence(testTimeSeries.begin() + topMotifSetPos[0],
        testTimeSeries.begin() + topMotifSetPos[0] + 20);
    generator.testMeanStdDev(subsequence, mean, stdDev);
    TEST(mean <= topMotifSetMeans[0] + 0.0000001 && mean >= topMotifSetMeans[0]
        - 0.0000001);
    TEST(stdDev <= topMotifSetStdDevs[0] + 0.0000001 && stdDev >=
        topMotifSetStdDevs[0] - 0.0000001);

    tsg::rseq subsequence1(testTimeSeries.begin() + topMotifSetPos[1],
        testTimeSeries.begin() + topMotifSetPos[1] + 20);
    generator.testMeanStdDev(subsequence1, mean, stdDev);
    TEST(mean <= topMotifSetMeans[1] + 0.0000001 && mean >= topMotifSetMeans[1]
        - 0.0000001);
    TEST(stdDev <= topMotifSetStdDevs[1] + 0.0000001 && stdDev >=
        topMotifSetStdDevs[1] - 0.0000001);

    tsg::rseq subsequence2(testTimeSeries.begin() + topMotifSetPos[2],
        testTimeSeries.begin() + topMotifSetPos[2] + 20);
    generator.testMeanStdDev(subsequence2, mean, stdDev);
    TEST(mean <= topMotifSetMeans[2] + 0.0000001 && mean >= topMotifSetMeans[2]
        - 0.0000001);
    TEST(stdDev <= topMotifSetStdDevs[2] + 0.0000001 && stdDev >=
        topMotifSetStdDevs[2] - 0.0000001);


    //test running sum and sum of squares
    tsg::rseq sums;
    tsg::rseq sumSquares;
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
      tsg::rseq testSeries(testTimeSeries);

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
          std::numeric_limits<double>::max());
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    try {

      generator.testSimilarity(testTimeSeries, 303, topMotifSetPos[0],
          std::numeric_limits<double>::max());
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    try {

      generator.testSimilarity(testTimeSeries, topMotifSetPos[0], 303,
          std::numeric_limits<double>::max());
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    generator.testCalcRunnings(testTimeSeries);

    double similarity = generator.testSimilarity(testTimeSeries,
        topMotifSetPos[0], topMotifSetPos[1], std::numeric_limits<double>::max());
    TEST(abs(similarity - testSequencesSimilarty) <= 0.0000001);

    similarity = generator.testSimilarity(testTimeSeries, topMotifSetPos[1],
        topMotifSetPos[2], std::numeric_limits<double>::max());
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
          std::numeric_limits<double>::max());
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    try {

      generator.testSimilarity(testTimeSeries, tsg::rseq(),
          topMotifSetMeans[0], topMotifSetStdDevs[0], topMotifSetPos[0],
          std::numeric_limits<double>::max());
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    try {

      generator.testSimilarity(testTimeSeries, { 0.0, 0.0 },
          topMotifSetMeans[0], topMotifSetStdDevs[0], topMotifSetPos[0],
          std::numeric_limits<double>::max());
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    try {

      generator.testSimilarity(testTimeSeries, testSequenceOne,
          topMotifSetMeans[0], topMotifSetStdDevs[0], 303,
          std::numeric_limits<double>::max());
      TEST(!"Has to throw an error!");
    }
    catch (int e) {

      TEST(e == EXIT_FAILURE);
    }

    similarity = generator.testSimilarity(testTimeSeries,
        testMotifSetSubsequences[0], topMotifSetMeans[0],
        topMotifSetStdDevs[0], topMotifSetPos[1],
        std::numeric_limits<double>::max());
    TEST(abs(similarity - testSequencesSimilarty) <= 0.0000001);

    similarity = generator.testSimilarity(testTimeSeries,
        testMotifSetSubsequences[1], topMotifSetMeans[1],
        topMotifSetStdDevs[1], topMotifSetPos[2],
        std::numeric_limits<double>::max());
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
    tsg::rseq sequenceTwo;
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
    tsg::rseq timeSeries_out;
    tsg::rseq d_out;
    tsg::iseqs positions_out;
    tsg::rseqs motif;


    for (int itr = 0; itr < 3; itr++) {

      try {

        generator.run(timeSeries_out, motif, d_out, positions_out);
        TEST("Has to run without throwing an error!");

        simGenerator.testCalcRunnings(timeSeries_out);

        for (auto &pos0 : positions_out[0])
          for (auto &pos1 : positions_out[0])
            TEST(simGenerator.testSimilarity(timeSeries_out, pos0, pos1,
                  2 * d_out[0]) <= 2 * d_out[0]);

        TEST(d_out.size() == 2);
        TEST(positions_out.size() == 2);
      }
      catch (...) {

        TEST(!"Has to run without throwing an error!");
      }
    }
  }


  //reset cerr
  std::cerr.rdbuf(cerr_buff);
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
    TEST_SECTION("top set motif");
    test_tsm();
    TEST_SECTION("time series generator");
    test_tsgenerator();

    TEST_SUMMARY;
  }
  catch (...) {
    TEST_EXCEPTION;
  }

  return EXIT_SUCCESS;
}
