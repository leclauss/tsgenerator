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
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(subsequence[5] == 5.0);
  TEST_R(subsequence[8] == 5.0);
  tsg::generateBoxMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(subsequence[5] == -7.0);
  TEST_R(subsequence[8] == -7.0);
  tsg::generateTriangleMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(0.0 < subsequence[5]);
  TEST_R(subsequence[4] == subsequence[5]);
  TEST_R(subsequence[1] == subsequence[8]);
  TEST_R(subsequence[3] != subsequence[9]);
  tsg::generateTriangleMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(subsequence[4] < 0.0);
  TEST_R(subsequence[4] == subsequence[5]);
  TEST_R(subsequence[1] == subsequence[8]);
  TEST_R(subsequence[3] != subsequence[9]);
  tsg::generateSemicircleMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(0.0 < subsequence[5]);
  TEST_R(subsequence[4] == subsequence[5]);
  TEST_R(subsequence[1] == subsequence[8]);
  TEST_R(subsequence[3] != subsequence[9]);
  tsg::generateSemicircleMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(subsequence[4] < 0.0);
  TEST_R(subsequence[4] == subsequence[5]);
  TEST_R(subsequence[1] == subsequence[8]);
  TEST_R(subsequence[3] != subsequence[9]);
  tsg::generateTrapezoidMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(0.0 < subsequence[5]);
  TEST_R(subsequence[4] == subsequence[5]);
  TEST_R(subsequence[1] == subsequence[8]);
  TEST_R(subsequence[3] != subsequence[9]);
  tsg::generateTrapezoidMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(subsequence[4] < 0.0);
  TEST_R(subsequence[4] == subsequence[5]);
  TEST_R(subsequence[1] == subsequence[8]);
  TEST_R(subsequence[3] != subsequence[9]);
  tsg::generatePositiveFlankMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(0.0 < subsequence[5]);
  TEST_R(subsequence[2] < subsequence[3]);
  TEST_R(subsequence[5] < subsequence[6]);
  TEST_R(subsequence[9] < subsequence[8]);
  tsg::generatePositiveFlankMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(subsequence[4] < 0.0);
  TEST_R(subsequence[2] > subsequence[3]);
  TEST_R(subsequence[5] > subsequence[6]);
  TEST_R(subsequence[9] > subsequence[8]);
  tsg::generateNegativeFlankMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(0.0 < subsequence[5]);
  TEST_R(subsequence[2] > subsequence[3]);
  TEST_R(subsequence[5] > subsequence[6]);
  TEST_R(subsequence[0] < subsequence[1]);
  tsg::generateNegativeFlankMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(subsequence[4] < 0.0);
  TEST_R(subsequence[2] < subsequence[3]);
  TEST_R(subsequence[5] < subsequence[6]);
  TEST_R(subsequence[0] > subsequence[1]);
  tsg::generateSineMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(0.0 < subsequence[3]);
  TEST_R(subsequence[8] < 0.0);
  TEST_R(subsequence[1] > subsequence[7]);
  TEST_R(subsequence[1] < subsequence[3]);
  TEST_R(abs(subsequence[1] + subsequence[8]) < 0.001);
  TEST_R(abs(subsequence[3] + subsequence[6]) < 0.001);
  tsg::generateSineMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(0.0 > subsequence[3]);
  TEST_R(subsequence[8] > 0.0);
  TEST_R(subsequence[1] < subsequence[7]);
  TEST_R(subsequence[1] > subsequence[3]);
  TEST_R(abs(subsequence[1] + subsequence[8]) < 0.001);
  TEST_R(abs(subsequence[3] + subsequence[6]) < 0.001);
  tsg::generateCosineMotif(subsequence, 1.0, 1.0, 10, 5.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(subsequence[4] < 0.0);
  TEST_R(abs(subsequence[4] - subsequence[5]) < 0.001);
  TEST_R(abs(subsequence[3] - subsequence[6]) < 0.001);
  TEST_R(subsequence[3] != subsequence[9]);
  tsg::generateCosineMotif(subsequence, 1.0, 1.0, 10, -7.0);
  TEST_R(subsequence.size() == 10);
  TEST_R(subsequence.front() == 0.0);
  TEST_R(subsequence.back() == 0.0);
  TEST_R(0.0 < subsequence[5]);
  TEST_R(abs(subsequence[4] - subsequence[5]) < 0.001);
  TEST_R(abs(subsequence[3] - subsequence[6]) < 0.001);
  TEST_R(subsequence[3] != subsequence[9]);
}

void test_freepositions() {

  TEST_GROUP_FUNCTION;

  { //test repeating without remove

    auto length = 10;
    auto window = 3;
    tsg::FreePositions freePos(length, window);

    for (int x = 0; x < 300; x++) {

      int pos = freePos.calculateRandomPosition();

      TEST_R(0 <= pos);
      TEST_R(pos <= length - window);
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

      TEST_R(0 <= pos[i]);
      TEST_R(pos[i] <= length - window);

      for (int j = 0; j < (int) pos.size() - 1; j++) {

        TEST_R(abs(pos[j] - pos[i]) >= 2 * window);
      }
    }
  }

  { //test failure
    auto length = 3;
    auto window = 3;
    tsg::FreePositions freePos(length, window);

    auto pos = freePos.calculateRandomPosition();
    freePos.removePosition();

    TEST_R(0 <= pos);
    TEST_R(pos <= length - window);

    try {

      pos = freePos.calculateRandomPosition();
      TEST_R(!"Has to throw an error!");
    }
    catch (int e) {

      TEST_R(e == EXIT_FAILURE);
    }
  }
}

void test_basets() {

  TEST_GROUP_FUNCTION;

  tsg::BaseTS baseTS;

  tsg::rseq ts;

  //simple random walk
  baseTS.simpleRandomWalk(ts, 100, 1.0, 0.0);
  TEST_R(ts.size() == 100);

  baseTS.simpleRandomWalk(ts, 23, 1.0, 0.1);
  TEST_R(ts.size() == 23);

  baseTS.simpleRandomWalk(ts, 0, 1.0, 0.1);
  TEST_R(ts.size() == 0);

  //real random walk
  baseTS.realRandomWalk(ts, 100, 1.0, 0.0);
  TEST_R(ts.size() == 100);

  baseTS.realRandomWalk(ts, 23, 1.0, 0.1);
  TEST_R(ts.size() == 23);

  baseTS.realRandomWalk(ts, 0, 1.0, 0.1);
  TEST_R(ts.size() == 0);

  //normal random walk
  baseTS.normalRandomWalk(ts, 100, 1.0, 0.0);
  TEST_R(ts.size() == 100);

  baseTS.normalRandomWalk(ts, 23, 1.0, 0.1);
  TEST_R(ts.size() == 23);

  //linear random walk
  baseTS.linearRandomWalk(ts, 100, 1.0, 3.0, 0.0);
  TEST_R(ts.size() == 100);

  baseTS.linearRandomWalk(ts, 23, 1.0, 3.0, 0.1);
  TEST_R(ts.size() == 23);

  baseTS.linearRandomWalk(ts, 0, 1.0, 3.0, 0.1);
  TEST_R(ts.size() == 0);

  //bounded simple random walk
  baseTS.simpleRandomWalk(ts, 100, 1.0, 0.0);
  TEST_R(ts.size() == 100);

  baseTS.simpleRandomWalk(ts, 23, 1.0, 0.1);
  TEST_R(ts.size() == 23);

  baseTS.simpleRandomWalk(ts, 0, 1.0, 0.1);
  TEST_R(ts.size() == 0);

  //bounded real random walk
  baseTS.realRandomWalk(ts, 100, 1.0, 20.0, 0.0);
  TEST_R(ts.size() == 100);

  baseTS.realRandomWalk(ts, 23, 1.0, 20.0, 0.1);
  TEST_R(ts.size() == 23);

  baseTS.realRandomWalk(ts, 0, 1.0, 20.0, 0.1);
  TEST_R(ts.size() == 0);

  //bounded normal random walk
  baseTS.normalRandomWalk(ts, 100, 1.0, 20.0, 0.0);
  TEST_R(ts.size() == 100);

  baseTS.normalRandomWalk(ts, 23, 1.0, 20.0, 0.1);
  TEST_R(ts.size() == 23);

  baseTS.normalRandomWalk(ts, 0, 1.0, 20.0, 0.1);
  TEST_R(ts.size() == 0);

  //bounded linear random walk
  baseTS.linearRandomWalk(ts, 100, 1.0, 3.0, 20.0, 0.0);
  TEST_R(ts.size() == 100);

  baseTS.linearRandomWalk(ts, 23, 1.0, 3.0, 20.0, 0.1);
  TEST_R(ts.size() == 23);

  baseTS.linearRandomWalk(ts, 0, 1.0, 3.0, 20.0, 0.1);
  TEST_R(ts.size() == 0);

  //uniform random
  baseTS.uniformRandom(ts, 100, 1.0, 0.0);
  TEST_R(ts.size() == 100);

  baseTS.uniformRandom(ts, 23, 1.0, 0.1);
  TEST_R(ts.size() == 23);

  baseTS.uniformRandom(ts, 0, 1.0, 0.1);
  TEST_R(ts.size() == 0);

  //normal random
  baseTS.normalRandom(ts, 100, 1.0, 0.0);
  TEST_R(ts.size() == 100);

  baseTS.normalRandom(ts, 23, 1.0, 0.1);
  TEST_R(ts.size() == 23);

  baseTS.normalRandom(ts, 0, 1.0, 0.1);
  TEST_R(ts.size() == 0);

  //piecewise linear random
  baseTS.piecewiseLinearRandom(ts, 100, 1.0, 0.0);
  TEST_R(ts.size() == 100);

  baseTS.piecewiseLinearRandom(ts, 23, 1.0, 0.1);
  TEST_R(ts.size() == 23);

  baseTS.piecewiseLinearRandom(ts, 0, 1.0, 0.1);
  TEST_R(ts.size() == 0);

  //spline repeated
  baseTS.splineRepeated(ts, 1000, 20.0, 20.0, 10, 2.0);
  TEST_R(ts.size() == 1000);

  baseTS.splineRepeated(ts, 23, 20.0, 20.0, 10, 2.0);
  TEST_R(ts.size() == 23);

  baseTS.splineRepeated(ts, 0, 20.0, 20.0, 10, 2.0);
  TEST_R(ts.size() == 0);

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
  tsg::rseq paaGT = { 0.2569, 0.1308, -0.5240, -0.3319, 0.0035, -0.2243 };

  try {

    tsm.zNormalPAA(paa, 56, 20, 6);

    for (int i = 0; i < (int)paaGT.size(); i++)
      TEST_R(paa[i] < paaGT[i] + 0.0001 && paa[i] > paaGT[i] - 0.0001);
  }
  catch (...) {

    TEST_R(!"Shouldn't throw an error!");
  }

  //test the invNormalCDF implementation
  tsg::rseq bGT = { -std::numeric_limits<double>::infinity(), -0.8416, -0.2533,
    0.2533, 0.8416, std::numeric_limits<double>::infinity() };

  double cdf = tsm.invNormalCDF(1.0 / 5.0, 0.000001);

  TEST_R(cdf < bGT[1] && cdf > bGT[1] - 0.0001);

  tsg::rseq b;

  try {

    tsm.invNormalCDF(5, b, 0.000001);

    TEST_R(b[0] == bGT[0]);
    TEST_R(b[5] == bGT[5]);

    for (int i = 1; i < (int)bGT.size() - 1; i++)
      TEST_R(b[i] < bGT[i] + 0.0001 && b[i] > bGT[i] - 0.0001);
  }
  catch (...) {

    TEST_R(!"Shouldn't throw an error!");
  }

  //test the z-normalized SAX implemenation
  tsg::word sax;

  try {

    tsm.zNormalSAX(sax, 148, 20, 6, bGT);
    TEST_R("Shouldn't throw an error!");
  }
  catch (...) {

    TEST_R(!"Shouldn't throw an error!");
  }

  //test the SAX minimum distance function
  double dist = -1.0;

  try {

    dist = tsm.saxDist({ 1, 1, 1, 1, 1 }, { 3, 3, 3, 3, 3}, 20, b);
    TEST_R(dist == sqrt(20.0 / 5.0 * 5.0 * (b[1] - b[2]) * (b[1] - b[2])));

    dist = tsm.saxDist({ 1, 1, 1, 1, 1 }, { 2, 2, 2, 2, 2}, 20, b);
    TEST_R(dist == 0.0);
  }
  catch (...) {

    TEST_R(!"Shouldn't throw an error!");
  }

  //test the z-normalized Euclidean distance function
  dist = tsm.dist(topMotifSetPos[0], topMotifSetPos[1], 20.0);
  TEST_R(abs(dist - testSequencesSimilarty) <= 0.0000001);

  dist = tsm.dist(topMotifSetPos[1], topMotifSetPos[2], 20.0);
  TEST_R(abs(dist) <= 0.0000001);

  //test the ADM implementation
  tsg::rseqs dGT = {
    { 0, 7.43966, 6.2613, 6.24954, 6.50969, 5.62853 },
    { 7.43966, 0, 4.00941, 5.40625, 4.59849, 7.11225 },
    { 6.2613, 4.00941, 0, 4.70289, 2.52727, 5.94591 },
    { 6.24954, 5.40625, 4.70289, 0, 4.4723, 6.72786 },
    { 6.50969, 4.59849, 2.52727, 4.4723, 0, 5.5463 },
    { 5.62853, 7.11225, 5.94591, 6.72786, 5.5463, 0 }
  };

  try {

    tsm.adm({4, 32, 86, 111, 140, 268 }, 20.0, 7.0);
    for (int i = 0; i < (int)dGT.size(); i++)
      for (int j = 0; j < (int)dGT.size(); j++)
        TEST_R(tsm.distADM(i, j) < dGT[i][j] + 0.0001 && tsm.distADM(i, j)
            > dGT[i][j] - 0.0001);
  }
  catch (...) {

    TEST_R(!"Shouldn't throw an error!");
  }

  //test the top set motif discovery implementation
  tsg::iseq motif;

  try {

    TEST_R(tsm.tsm(motif, 20, 0.5) == 3);
    TEST_R((int)motif.size() == 3);
    TEST_R(tsm.tsm(motif, 20, 0.001) == 2);
    TEST_R((int)motif.size() == 2);
    TEST_R(tsm.tsm(motif, 20, 6.0) == 13);
    TEST_R((int)motif.size() == 13);
    for (int i = 0; i < (int)motif.size(); i++)
      for (int j = 0; j < (int)motif.size(); j++)
        if (i != j)
          TEST_R(abs(motif[i] - motif[j]) >= 20);
  }
  catch (...) {

    TEST_R(!"Shouldn't throw an error!");
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
    TEST_R("Should work without throw!");
  }
  catch (...) {

    TEST_R(!"Should work without throw!");
  }

  try {

    TestTSGenerator generator(300, 20, 20.0, 0.0, 9, 3, 50.0);
    TEST_R(!"Has to throw an error!");
  }
  catch (int e) {

    TEST_R(e == EXIT_FAILURE);
  }

  try {

    TestTSGenerator generator(300, 20, 20.0, 0.0, 2, 2, 50.0);
    TEST_R(!"Has to throw an error!");
  }
  catch (int e) {

    TEST_R(e == EXIT_FAILURE);
  }

  try {

    TestTSGenerator generator(99, 20, 20, 0.0, 9, 3, 50.0);
    TEST_R(!"Has to throw an error!");
  }
  catch (int e) {

    TEST_R(e == EXIT_FAILURE);
  }

  try {

    TestTSGenerator generator(300, 500, 20, 0.0, 9, 3, 50.0);
    TEST_R(!"Has to throw an error!");
  }
  catch (int e) {

    TEST_R(e == EXIT_FAILURE);
  }

  try {

    TestTSGenerator generator(300, -1, 20, 0.0, 9, 3, 50.0);
    TEST_R(!"Has to throw an error!");
  }
  catch (int e) {

    TEST_R(e == EXIT_FAILURE);
  }

  try {

    TestTSGenerator generator(300, 20, 20.0, -1.0, 9, 3, 50.0);
    TEST_R(!"Has to throw an error!");
  }
  catch (int e) {

    TEST_R(e == EXIT_FAILURE);
  }

  try {

    TestTSGenerator generator(300, 20, 20.0, 0.01, 9, -1, 50.0);
    TEST_R(!"Has to throw an error!");
  }
  catch (int e) {

    TEST_R(e == EXIT_FAILURE);
  }

  try {

    TestTSGenerator generator(300, 20, 20.0, 0.0, 9, 20, 50.0);
    TEST_R(!"Has to throw an error!");
  }
  catch (int e) {

    TEST_R(e == EXIT_FAILURE);
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
      TEST_R(!"Has to throw an error!");
    }
    catch (int e) {

      TEST_R(e == EXIT_FAILURE);
    }

    generator.testMeanStdDev({ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0 }, mean, stdDev);
    TEST_R(mean <= std::numeric_limits<double>::min() && mean >=
        -std::numeric_limits<double>::min());
    TEST_R(stdDev == 1.0);

    generator.testMeanStdDev({ 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
        3.0 }, mean, stdDev);
    TEST_R(mean <= 3.0 + std::numeric_limits<double>::min() && mean >= 3.0
        - std::numeric_limits<double>::min());
    TEST_R(stdDev == 1.0);

    generator.testMeanStdDev({ 3.0, -3.0, 1.75, 9.0, 33.0, 3.101, 0.03,
        3.99999, 4.23, -7.093 }, mean, stdDev);
    TEST_R(mean <= 4.801799 + 0.0000001 && mean >= 4.801799 - 0.0000001);
    TEST_R(stdDev <= 10.2679615 + 0.0000001 && stdDev >= 10.2679615 - 0.0000001);

    rWindow = 0.05;
    generator.setWindow(20);
    tsg::rseq subsequence(testTimeSeries.begin() + topMotifSetPos[0],
        testTimeSeries.begin() + topMotifSetPos[0] + 20);
    generator.testMeanStdDev(subsequence, mean, stdDev);
    TEST_R(mean <= topMotifSetMeans[0] + 0.0000001 && mean >= topMotifSetMeans[0]
        - 0.0000001);
    TEST_R(stdDev <= topMotifSetStdDevs[0] + 0.0000001 && stdDev >=
        topMotifSetStdDevs[0] - 0.0000001);

    tsg::rseq subsequence1(testTimeSeries.begin() + topMotifSetPos[1],
        testTimeSeries.begin() + topMotifSetPos[1] + 20);
    generator.testMeanStdDev(subsequence1, mean, stdDev);
    TEST_R(mean <= topMotifSetMeans[1] + 0.0000001 && mean >= topMotifSetMeans[1]
        - 0.0000001);
    TEST_R(stdDev <= topMotifSetStdDevs[1] + 0.0000001 && stdDev >=
        topMotifSetStdDevs[1] - 0.0000001);

    tsg::rseq subsequence2(testTimeSeries.begin() + topMotifSetPos[2],
        testTimeSeries.begin() + topMotifSetPos[2] + 20);
    generator.testMeanStdDev(subsequence2, mean, stdDev);
    TEST_R(mean <= topMotifSetMeans[2] + 0.0000001 && mean >= topMotifSetMeans[2]
        - 0.0000001);
    TEST_R(stdDev <= topMotifSetStdDevs[2] + 0.0000001 && stdDev >=
        topMotifSetStdDevs[2] - 0.0000001);


    //test running sum and sum of squares
    tsg::rseq sums;
    tsg::rseq sumSquares;
    double std;
    bool flagMeansStdsTest = true;

    generator.testCalcRunnings(testTimeSeries);
    sums = generator.getSums();
    sumSquares = generator.getSumSquares();

    if (TEST_IF(testMeans.size() == sums.size())) {

      for (int itr = 0; itr < (int)sums.size(); itr++)
        if (!(testMeans[itr] - 0.000001 <= sums[itr] * rWindow  &&
              testMeans[itr] + 0.000001 >= sums[itr] * rWindow)) {

          flagMeansStdsTest = false;
          break;
        }

      TEST_R(flagMeansStdsTest);
    }

    flagMeansStdsTest = true;

    if (TEST_IF(testStds.size() == sumSquares.size())) {

      for (int itr = 0; itr < (int)sumSquares.size(); itr++) {

        std = sums[itr] * rWindow;
        std = sqrt(sumSquares[itr] * rWindow - std * std);

        if (!(testStds[itr] - 0.000001 <= std &&
              testStds[itr] + 0.000001 >= std)) {

          flagMeansStdsTest = false;
          break;
        }
      }

      TEST_R(flagMeansStdsTest);
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

      if (TEST_IF(testMeans.size() == sums.size())) {

        for (int itr = 0; itr < (int)sums.size(); itr++)
          if (!(testMeans[itr] - 0.000001 <= sums[itr] * rWindow  &&
                testMeans[itr] + 0.000001 >= sums[itr] * rWindow)) {

            flagMeansStdsTest = false;
            break;
          }

        TEST_R(flagMeansStdsTest);
      }

      flagMeansStdsTest = true;

      if (TEST_IF(testStds.size() == sumSquares.size())) {

        for (int itr = 0; itr < (int)sumSquares.size(); itr++) {

          std = sums[itr] * rWindow;
          std = sqrt(sumSquares[itr] * rWindow - std * std);

          if (!(testStds[itr] - 0.000001 <= std &&
                testStds[itr] + 0.000001 >= std)) {

            flagMeansStdsTest = false;
            break;
          }
        }

        TEST_R(flagMeansStdsTest);
      }
    }


    //test the first similarity functions
    try {

      generator.testSimilarity({}, topMotifSetPos[0], topMotifSetPos[0],
          std::numeric_limits<double>::max());
      TEST_R(!"Has to throw an error!");
    }
    catch (int e) {

      TEST_R(e == EXIT_FAILURE);
    }

    try {

      generator.testSimilarity(testTimeSeries, 303, topMotifSetPos[0],
          std::numeric_limits<double>::max());
      TEST_R(!"Has to throw an error!");
    }
    catch (int e) {

      TEST_R(e == EXIT_FAILURE);
    }

    try {

      generator.testSimilarity(testTimeSeries, topMotifSetPos[0], 303,
          std::numeric_limits<double>::max());
      TEST_R(!"Has to throw an error!");
    }
    catch (int e) {

      TEST_R(e == EXIT_FAILURE);
    }

    generator.testCalcRunnings(testTimeSeries);

    double similarity = generator.testSimilarity(testTimeSeries,
        topMotifSetPos[0], topMotifSetPos[1], std::numeric_limits<double>::max());
    TEST_R(abs(similarity - testSequencesSimilarty) <= 0.0000001);

    similarity = generator.testSimilarity(testTimeSeries, topMotifSetPos[1],
        topMotifSetPos[2], std::numeric_limits<double>::max());
    TEST_R(abs(similarity) <= 0.0000001);

    similarity = generator.testSimilarity(testTimeSeries, topMotifSetPos[0],
        topMotifSetPos[2], 0.01);
    TEST_R(similarity >= 0.01);

    similarity = generator.testSimilarity(testTimeSeries, topMotifSetPos[1],
        topMotifSetPos[1], 0.01);
    TEST_R(abs(similarity) <= 0.0000001);


    //test calculate raw subsequence
    subsequence.clear();
    subsequence.resize(0);

    generator.testCalculateSubsequence(subsequence);
    TEST_R(subsequence.size() == 20);
    TEST_R(subsequence.front() == 0.0);
    TEST_R(subsequence.back() == 0.0);
    TEST_R(0.0 < subsequence[5]);
    TEST_R(subsequence[4] == subsequence[15]);
    TEST_R(subsequence[1] == subsequence[18]);
    TEST_R(subsequence[3] != subsequence[19]);
  }

  {
    //test synthetic time series pair motif generation function
    TestTSGenerator simGenerator(1000, 20, 1.0, 2.0, 4, 3, 10.0, 1.0, 3, 5,
        20.0, 0);
    TestTSGenerator generator(1000, 20, 1.0, 2.0, 4, 3, 10.0, 1.0, 3, 5,
        20.0, 0);
    tsg::rseq timeSeries_out;
    tsg::rseq d_out;
    tsg::iseqs positions_out;
    tsg::rseqs motif;
    double d;


    for (int itr = 0; itr < 3; itr++) {

      try {

        generator.run(timeSeries_out, motif, d_out, positions_out);
        TEST_R("Has to run without throwing an error!");

        simGenerator.testCalcRunnings(timeSeries_out);
        d = simGenerator.testSimilarity(timeSeries_out, positions_out[0][0],
              positions_out[0][1], d_out[0]);
        TEST_R(d <= d_out[0] + 0.0000001 && d >= d_out[0] - 0.0000001);

        TEST_R(d_out.size() == 1);
        TEST_R(positions_out.size() == 1);
      }
      catch (...) {

        TEST_R(!"Has to run without throwing an error!");
      }
    }
  }

  {
    //test synthetic time series set motif generation function
    TestTSGenerator simGenerator(1000, 20, 1.0, 2.0, 4, 3, 10.0, 1.0, 3, 5,
        20.0, 1);
    TestTSGenerator generator(1000, 20, 1.0, 2.0, 4, 3, 10.0, 1.0, 3, 5,
        20.0, 1);
    tsg::rseq timeSeries_out;
    tsg::rseq d_out;
    tsg::iseqs positions_out;
    tsg::rseqs motif;


    for (int itr = 0; itr < 3; itr++) {

      try {

        generator.run(timeSeries_out, motif, d_out, positions_out);
        TEST_R("Has to run without throwing an error!");

        simGenerator.testCalcRunnings(timeSeries_out);

        for (auto &pos0 : positions_out[0])
          for (auto &pos1 : positions_out[0])
            TEST_R(simGenerator.testSimilarity(timeSeries_out, pos0, pos1,
                  2 * d_out[0]) <= 2 * d_out[0]);

        TEST_R(d_out.size() == 2);
        TEST_R(positions_out.size() == 2);
      }
      catch (...) {

        TEST_R(!"Has to run without throwing an error!");
      }
    }
  }

  {
    //test synthetic time series latent motif generation function
    TestTSGenerator simGenerator(1000, 20, 1.0, 2.0, 4, 3, 10.0, 1.0, 3, 5,
        20.0, 2);
    TestTSGenerator generator(1000, 20, 1.0, 2.0, 4, 3, 10.0, 1.0, 3, 5,
        20.0, 2);
    tsg::rseq timeSeries_out;
    tsg::rseq d_out;
    tsg::iseqs positions_out;
    tsg::rseqs motif;


    for (int itr = 0; itr < 3; itr++) {

      try {

        generator.run(timeSeries_out, motif, d_out, positions_out);
        TEST_R("Has to run without throwing an error!");

        simGenerator.testCalcRunnings(timeSeries_out);

        for (auto &pos0 : positions_out[0])
          for (auto &pos1 : positions_out[0])
            TEST_R(simGenerator.testSimilarity(timeSeries_out, pos0, pos1,
                  2 * d_out[0]) <= 2 * d_out[0]);

        TEST_R(d_out.size() == 2);
        TEST_R(positions_out.size() == 2);
      }
      catch (...) {

        TEST_R(!"Has to run without throwing an error!");
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
