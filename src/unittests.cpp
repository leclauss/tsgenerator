///\file unittests.cpp
///
///\brief File contains all unit test.
///
///This is the source file to run all unit tests with stest.

#include <unittests.hpp>


void test_check_if_float() {

  TEST_GROUP_FUNCTION;
  TEST(check_if_float("123456789"));
  TEST(check_if_float("-123456789"));
  TEST(check_if_float("123456789.123456789"));
  TEST(check_if_float("-123456789.123456789"));
  TEST(!check_if_float("abcd"));
  TEST(!check_if_float(""));
}

char **argv;
vector<string> argTokens;
vector<string> payload;

void test_parseArgs() {

  TEST_GROUP_FUNCTION;

  int argc = 10;
  char arg0[] = "--arg0";
  char arg1[] = "payload";
  char arg2[] = "-arg2";
  char arg3[] = "-1000.1";
  char arg4[] = "-empty";
  char arg5[] = "--multiplePayload";
  char arg6[] = "-100.000";
  char arg7[] = "string1 string2";
  char arg8[] = "";
  char arg9[] = "end";
  argv = new char*[argc];
  argv[0] = arg0;
  argv[1] = arg1;
  argv[2] = arg2;
  argv[3] = arg3;
  argv[4] = arg4;
  argv[5] = arg5;
  argv[6] = arg6;
  argv[7] = arg7;
  argv[8] = arg8;
  argv[9] = arg9;

  parseArgs(argc, argv, argTokens);

  TEST(argTokens[0] == "--arg0");
  TEST(argTokens[1] == "payload");
  TEST(argTokens[2] == "-arg2");
  TEST(argTokens[3] == "-1000.1");
  TEST(argTokens[4] == "-empty");
  TEST(argTokens[5] == "--multiplePayload");
  TEST(argTokens[6] == "-100.000");
  TEST(argTokens[7] == "string1 string2");
  TEST(argTokens[8] == "end");
}

void test_checkArg() {

  TEST_GROUP_FUNCTION;

  TEST(!checkArg(argTokens, "", payload));
  TEST(!checkArg(argTokens, "-1000.1", payload));
  TEST(checkArg(argTokens, "-empty", payload));
  TEST(payload.empty());
  TEST(checkArg(argTokens, "--arg0", payload));
  TEST(payload.size() == 1);
  TEST(payload[0] == "payload");
  TEST(checkArg(argTokens, "-arg2", payload));
  TEST(payload.size() == 1);
  TEST(payload[0] == "-1000.1");
  TEST(checkArg(argTokens, "--multiplePayload", payload));
  TEST(payload.size() == 3);
  TEST(payload[0] == "-100.000");
  TEST(payload[1] == "string1 string2");
  TEST(payload[2] == "end");
}

void test_print_version() {

  TEST_GROUP_FUNCTION;

  ostringstream local;
  ostringstream testStream;

  auto cout_buff = cout.rdbuf(local.rdbuf());

  print_version();

  testStream << PROGNAME << " " << VERSION << endl;

  testStream << "Copyright \302\251 2018 Rafael Moczalla" << endl;
  testStream << "Lizenz: Creative Commons - Attribution" <<
    " - Non-Commercial - Share Alike" << endl;
  testStream << "There are no guarantees as far as the law permits." <<
    endl << endl;
  testStream << "Written by Rafael Moczalla" << endl;

  cout.rdbuf(cout_buff);

  TEST(testStream.str() == local.str());
}

void test_print_help() {

  TEST_GROUP_FUNCTION;

  ostringstream local;
  ostringstream testStream;

  auto cout_buff = cout.rdbuf(local.rdbuf());

  print_help();

  testStream << "Call: " << endl;
  testStream << "    " << PROGNAME <<
    " -l INTEGER -w INTEGER -rd FLOAT [Options]" << endl;
  testStream << endl;
  testStream << "    -l,    --length INTEGER            " <<
    "                      Sets the length of the time series." << endl;
  testStream << "    -w,    --windowSize INTEGER        " <<
    "                      Sets the window size of the time series." << endl;
  testStream << "    -rd,   --randomness FLOAT          " <<
    "                      Sets the randomness of the time series." << endl;
  testStream << endl;
  testStream << "Options:" << endl;
  testStream << "    -lm,   --latentMotif (STRING INTEGER FLOAT)+" <<
    "             Sets the operating mode to latent motif and" << endl;
  testStream << "                                       " <<
    "                      the parameters for the latent motif namely" << endl;
  testStream << "                                       " <<
    "                      type, count and height." << endl;
  testStream << "    -o,    --out NAME                  " <<
    "                      Sets the output file name." << endl;
  testStream << "    -tsn,  --timeSeriesName NAME       " <<
    "                      Sets the time series name." << endl;
  testStream << "    -ho,   --horizontalOutput          " <<
    "                      Sets the output mode to horizontal." << endl;
  testStream << "    -r,    --range FLOAT FLOAT         " <<
    "                      Sets the range of the time series values." << endl;
  testStream << "    -h,    --help                      " <<
    "                      Prints these help text." << endl;
  testStream << "    -v,    --version                   " <<
    "                      Prints the version information." << endl;

  cout.rdbuf(cout_buff);

  TEST(testStream.str() == local.str());
}

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

void test_outputgenerator() {

  //redirect cerr to keep test output clean
  ostringstream local;

  auto cerr_buff = cerr.rdbuf(local.rdbuf());

  TEST_GROUP_FUNCTION;
  int folderNumber = 0;
  int secondFolderNumber;

  while(exists("time_series_" + to_string(folderNumber)))
    folderNumber++;

  create_directory("time_series_" + to_string(folderNumber));

  secondFolderNumber = folderNumber + 1;

  while(exists("time_series_" + to_string(secondFolderNumber)))
    secondFolderNumber++;


  TestOutputGenerator generator;

  generator.setTimeSeriesName("Test Name");
  generator.testOpenFile();
  try {

    generator.testOpenFile();
    TEST(!"Has to throw an error!");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  generator.closeOutputFileStream();
  generator.closeOutputFileMotifSetStream();
  generator.closeOutputFileGNUPlotScriptStream();

  TEST(generator.getOutputFolderNumber() == secondFolderNumber);
  TEST(generator.getOutputFileName() == ("./time_series_"
        + to_string(secondFolderNumber) + "/time_series"));
  TEST(generator.getOutputFileMotifSetsName() == ("./time_series_"
        + to_string(secondFolderNumber) + "/time_series_meta"));
  TEST(generator.getOutputFileGNUPlotScriptName() == ("./time_series_"
        + to_string(secondFolderNumber) + "/time_series_plot"));
  TEST(generator.getTimeSeriesName() == ("Test Name"));
  TEST(remove("time_series_" + to_string(secondFolderNumber) + "/time_series_"
        + to_string(secondFolderNumber) + ".csv"));
  TEST(remove("time_series_" + to_string(secondFolderNumber)
        + "/time_series_meta_" + to_string(secondFolderNumber) + ".csv"));
  TEST(remove("time_series_" + to_string(secondFolderNumber)
        + "/time_series_plot_" + to_string(secondFolderNumber) + ".plt"));
  TEST(remove("time_series_" + to_string(secondFolderNumber)));

  remove("time_series_" + to_string(folderNumber));


  TestOutputGenerator generatorHorizontal;

  TEST(generatorHorizontal.getOutputFolderNumber() == folderNumber);
  TEST(generatorHorizontal.getOutputFileName() == ("./time_series_"
        + to_string(folderNumber) + "/time_series"));
  TEST(generatorHorizontal.getOutputFileMotifSetsName() == ("./time_series_"
        + to_string(folderNumber) + "/time_series_meta"));
  TEST(generatorHorizontal.getOutputFileGNUPlotScriptName() ==
      ("./time_series_" + to_string(folderNumber) + "/time_series_plot"));
  TEST(generatorHorizontal.getTimeSeriesName() == ("Time Series"));

  //test if output generator throw exceptions correctly
  try {

    generatorHorizontal.printTimeSeriesHorizontal({}, { topMotifPairSimilarity
        }, { windowSize }, { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorHorizontal.printTimeSeriesHorizontal(testTimeSeries, {},
        { windowSize }, { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorHorizontal.printTimeSeriesHorizontal(testTimeSeries,
        { topMotifPairSimilarity, topMotifPairSimilarity }, { windowSize },
        { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorHorizontal.printTimeSeriesHorizontal(testTimeSeries,
        { topMotifPairSimilarity }, {}, { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorHorizontal.printTimeSeriesHorizontal(testTimeSeries,
        { topMotifPairSimilarity }, { windowSize, windowSize },
        { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorHorizontal.printTimeSeriesHorizontal(testTimeSeries,
        { topMotifPairSimilarity }, { windowSize }, { {} });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  //remove folder
  TEST(remove("time_series_" + to_string(folderNumber)));


  TestOutputGenerator generatorVertical;

  TEST(generatorVertical.getOutputFolderNumber() == folderNumber);
  TEST(generatorVertical.getOutputFileName() == ("./time_series_"
        + to_string(folderNumber) + "/time_series"));
  TEST(generatorVertical.getOutputFileMotifSetsName() == ("./time_series_"
        + to_string(folderNumber) + "/time_series_meta"));
  TEST(generatorVertical.getOutputFileGNUPlotScriptName() == ("./time_series_"
        + to_string(folderNumber) + "/time_series_plot"));
  TEST(generatorVertical.getTimeSeriesName() == ("Time Series"));

  //test if output generator throw exceptions correctly
  try {

    generatorVertical.printTimeSeriesVertical({}, { topMotifPairSimilarity },
      { windowSize }, { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorVertical.printTimeSeriesVertical(testTimeSeries, {},
      { windowSize }, { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorVertical.printTimeSeriesVertical(testTimeSeries,
        { topMotifPairSimilarity, topMotifPairSimilarity }, { windowSize },
        { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorVertical.printTimeSeriesVertical(testTimeSeries,
        { topMotifPairSimilarity }, {}, { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorVertical.printTimeSeriesVertical(testTimeSeries,
        { topMotifPairSimilarity }, { windowSize, windowSize },
        { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorVertical.printTimeSeriesVertical(testTimeSeries,
        { topMotifPairSimilarity }, { windowSize }, { {} });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  //remove folder
  TEST(remove("time_series_" + to_string(folderNumber)));


  TestOutputGenerator generatorHorizontalTwo;

  TEST(generatorHorizontalTwo.getOutputFolderNumber() == folderNumber);
  TEST(generatorHorizontalTwo.getOutputFileName() == ("./time_series_"
        + to_string(folderNumber) + "/time_series"));
  TEST(generatorHorizontalTwo.getOutputFileMotifSetsName() == ("./time_series_"
        + to_string(folderNumber) + "/time_series_meta"));
  TEST(generatorHorizontalTwo.getOutputFileGNUPlotScriptName() ==
      ("./time_series_" + to_string(folderNumber) + "/time_series_plot"));
  TEST(generatorHorizontalTwo.getTimeSeriesName() == ("Time Series"));

  generatorHorizontalTwo.printTimeSeriesHorizontal(testTimeSeries,
      { topMotifSetRange, topMotifPairSimilarity }, { windowSize, windowSize },
      { topMotifSetPos, topMotifPairPos });

  //test if time series file has correct content
  ostringstream correctStream;
  correctStream << fixed;
  correctStream.str("");
  correctStream.clear();

  bool start = true;

  for (double item : testTimeSeries)
    if (start) {

      correctStream << item;
      start = false;
    }
    else
      correctStream << ", " << item;

  ifstream iTestTimeSeries("time_series_" + to_string(folderNumber)
      + "/time_series_" + to_string(folderNumber) + ".csv");

  string line;

  if (TEST_IF(getline(iTestTimeSeries, line) ? true : false))
    TEST(line == correctStream.str());

  TEST(!getline(iTestTimeSeries, line));

  iTestTimeSeries.close();

  //test if motif positions file has correct content
  ifstream iTestPositions("time_series_" + to_string(folderNumber)
      + "/time_series_meta_" + to_string(folderNumber) + ".csv");

  correctStream << defaultfloat;
  correctStream.str("");
  correctStream.clear();

  correctStream << " , \"range/similarity\", \"window size\"," <<
    " \"position 0\", \"position 1\", \"position 2\"";

  if(TEST_IF(getline(iTestPositions, line) ? true : false))
    TEST(line == correctStream.str());

  correctStream.str("");
  correctStream.clear();

  correctStream << fixed << "\"Top Latent Motif Locations\", " <<
    0.240684 << ", " << windowSize << ", ";
  correctStream << defaultfloat << topMotifSetPos[0] << ", " <<
    topMotifSetPos[1] << ", " << topMotifSetPos[2];

  if(TEST_IF(getline(iTestPositions, line) ? true : false))
    TEST(line == correctStream.str());

  correctStream.str("");
  correctStream.clear();

  correctStream << fixed << "\"Top Pair Motif Locations\", " <<
    0.000000 << ", " << windowSize << ", ";
  correctStream << defaultfloat << topMotifPairPos[0] << ", " <<
    topMotifPairPos[1];

  if(TEST_IF(getline(iTestPositions, line) ? true : false))
    TEST(line == correctStream.str());

  TEST(!getline(iTestPositions, line));

  iTestPositions.close();

  ifstream iTestPlot("time_series_" + to_string(folderNumber)
      + "/time_series_plot_" + to_string(folderNumber) + ".plt");

  if(TEST_IF(getline(iTestPlot, line) ? true : false))
#ifdef _WIN32
    TEST(line == "set terminal wxt size 1600, 300");
#elif __APPLE__
    TEST(line == "set terminal aqua size 1600, 300");
#else
    TEST(line == "set terminal qt size 1600, 300 persist");
#endif
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set xrange[0:299]");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set yrange[-67.793538:110.480620]");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set grid");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set obj rect from 35, graph 0 to 54,"
        " graph 1 back fc rgb \"#c5c5c5\" fs border rgb \"#c5c5c5\"");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set obj rect from 150, graph 0 to 169,"
        " graph 1 back fc rgb \"#c5c5c5\" fs border rgb \"#c5c5c5\"");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set obj rect from 250, graph 0 to 269,"
        " graph 1 back fc rgb \"#c5c5c5\" fs border rgb \"#c5c5c5\"");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set obj rect from 150, graph 0 to 169,"
        " graph 1 back fc rgb \"#656565\" "
        "fs pattern 1 border rgb \"#656565\"");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set obj rect from 250, graph 0 to 269,"
        " graph 1 back fc rgb \"#656565\" "
        "fs pattern 1 border rgb \"#656565\"");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == ("plot \"time_series_" + to_string(folderNumber) + ".csv\""
        " matrix title \"Time Series " + to_string(folderNumber) + "\" with"
        " lines"));
#ifdef _WIN32
#elif __APPLE__
#else
  TEST_IF(getline(iTestPlot, line) ? true : false);
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "pause mouse close");
#endif

  TEST(!getline(iTestPlot, line));

  iTestPlot.close();

  //remove files and folder
  TEST(remove("time_series_" + to_string(folderNumber) + "/time_series_"
        + to_string(folderNumber) + ".csv"));
  TEST(remove("time_series_" + to_string(folderNumber) + "/time_series_meta_"
        + to_string(folderNumber) + ".csv"));
  TEST(remove("time_series_" + to_string(folderNumber) + "/time_series_plot_"
        + to_string(folderNumber) + ".plt"));
  TEST(remove("time_series_" + to_string(folderNumber)));


  TestOutputGenerator generatorVerticalTwo;

  TEST(generatorVerticalTwo.getOutputFolderNumber() == folderNumber);
  TEST(generatorVerticalTwo.getOutputFileName() == ("./time_series_"
        + to_string(folderNumber) + "/time_series"));
  TEST(generatorVerticalTwo.getOutputFileMotifSetsName() == ("./time_series_"
        + to_string(folderNumber) + "/time_series_meta"));
  TEST(generatorVerticalTwo.getOutputFileGNUPlotScriptName() ==
      ("./time_series_" + to_string(folderNumber) + "/time_series_plot"));
  TEST(generatorVerticalTwo.getTimeSeriesName() == ("Time Series"));

  generatorVerticalTwo.printTimeSeriesVertical(testTimeSeries,
      { topMotifSetRange, topMotifPairSimilarity }, { windowSize, windowSize },
      { topMotifSetPos, topMotifPairPos });

  //test if time series file has correct content
  correctStream.str("");
  correctStream.clear();
  correctStream << fixed;

  iTestTimeSeries.open("time_series_" + to_string(folderNumber)
      + "/time_series_" + to_string(folderNumber) + ".csv");

  for (double item : testTimeSeries) {

    correctStream << item;

    if (TEST_IF(getline(iTestTimeSeries, line) ? true : false))
      TEST(line == correctStream.str());

    correctStream.str("");
    correctStream.clear();
  }

  TEST(!getline(iTestTimeSeries, line));

  iTestTimeSeries.close();

  //test if motif positions file has correct content
  iTestPositions.open("time_series_" + to_string(folderNumber)
      + "/time_series_meta_" + to_string(folderNumber) + ".csv");

  correctStream << defaultfloat;
  correctStream.str("");
  correctStream.clear();

  correctStream << " , \"Top Latent Motif Locations\"," <<
    " \"Top Pair Motif Locations\"";
  if(TEST_IF(getline(iTestPositions, line) ? true : false))
    TEST(line == correctStream.str());

  correctStream.str("");
  correctStream.clear();

  correctStream << fixed << "\"range/similarity\", " <<
    0.240684 << ", " << 0.000000 << defaultfloat;

  if(TEST_IF(getline(iTestPositions, line) ? true : false))
    TEST(line == correctStream.str());

  correctStream.str("");
  correctStream.clear();

  correctStream << "\"window size\", " << windowSize << ", " << windowSize;

  if(TEST_IF(getline(iTestPositions, line) ? true : false))
    TEST(line == correctStream.str());

  correctStream.str("");
  correctStream.clear();

  correctStream << "\"position 0\", " << topMotifSetPos[0] << ", " <<
    topMotifPairPos[0];

  if(TEST_IF(getline(iTestPositions, line) ? true : false))
    TEST(line == correctStream.str());

  correctStream.str("");
  correctStream.clear();

  correctStream << "\"position 1\", " << topMotifSetPos[1] << ", " <<
    topMotifPairPos[1];

  if(TEST_IF(getline(iTestPositions, line) ? true : false))
    TEST(line == correctStream.str());

  correctStream.str("");
  correctStream.clear();

  correctStream << "\"position 2\", " << topMotifSetPos[2] << ", ";

  if(TEST_IF(getline(iTestPositions, line) ? true : false))
    TEST(line == correctStream.str());

  TEST(!getline(iTestPositions, line));

  iTestPositions.close();

  iTestPlot.open("time_series_" + to_string(folderNumber)
      + "/time_series_plot_" + to_string(folderNumber) + ".plt");

  if(TEST_IF(getline(iTestPlot, line) ? true : false))
#ifdef _WIN32
    TEST(line == "set terminal wxt size 1600, 300");
#elif __APPLE__
    TEST(line == "set terminal aqua size 1600, 300");
#else
    TEST(line == "set terminal qt size 1600, 300 persist");
#endif
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set xrange[0:299]");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set yrange[-67.793538:110.480620]");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set grid");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set obj rect from 35, graph 0 to 54,"
        " graph 1 back fc rgb \"#c5c5c5\" fs border rgb \"#c5c5c5\"");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set obj rect from 150, graph 0 to 169,"
        " graph 1 back fc rgb \"#c5c5c5\" fs border rgb \"#c5c5c5\"");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set obj rect from 250, graph 0 to 269,"
        " graph 1 back fc rgb \"#c5c5c5\" fs border rgb \"#c5c5c5\"");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set obj rect from 150, graph 0 to 169,"
        " graph 1 back fc rgb \"#656565\" "
        "fs pattern 1 border rgb \"#656565\"");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "set obj rect from 250, graph 0 to 269,"
        " graph 1 back fc rgb \"#656565\" "
        "fs pattern 1 border rgb \"#656565\"");
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == ("plot \"time_series_" + to_string(folderNumber) + ".csv\""
        " title \"Time Series " + to_string(folderNumber) + "\" with lines"));
#ifdef _WIN32
#elif __APPLE__
#else
  TEST_IF(getline(iTestPlot, line) ? true : false);
  if(TEST_IF(getline(iTestPlot, line) ? true : false))
    TEST(line == "pause mouse close");
#endif

  TEST(!getline(iTestPlot, line));

  iTestPlot.close();

  //remove files and folder
  TEST(remove("time_series_" + to_string(folderNumber) + "/time_series_"
        + to_string(folderNumber) + ".csv"));
  TEST(remove("time_series_" + to_string(folderNumber) + "/time_series_meta_"
        + to_string(folderNumber) + ".csv"));
  TEST(remove("time_series_" + to_string(folderNumber) + "/time_series_plot_"
        + to_string(folderNumber) + ".plt"));
  TEST(remove("time_series_" + to_string(folderNumber)));


  TestOutputGenerator generatorHorizontalThree;

  TEST(generatorHorizontalThree.getOutputFolderNumber() == folderNumber);
  TEST(generatorHorizontalThree.getOutputFileName() == ("./time_series_"
        + to_string(folderNumber) + "/time_series"));
  TEST(generatorHorizontalThree.getOutputFileMotifSetsName() ==
      ("./time_series_" + to_string(folderNumber) + "/time_series_meta"));
  TEST(generatorHorizontalThree.getOutputFileGNUPlotScriptName() ==
      ("./time_series_" + to_string(folderNumber) + "/time_series_plot"));
  TEST(generatorHorizontalThree.getTimeSeriesName() == ("Time Series"));

  //test if output generator throw exceptions correctly
  try {

    generatorHorizontalThree.printTimeSeriesHorizontal({},
        { topMotifPairSimilarity }, { 8 }, { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorHorizontalThree.printTimeSeriesHorizontal(testTimeSeries, {},
      { 8 }, { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorHorizontalThree.printTimeSeriesHorizontal(testTimeSeries,
        { topMotifPairSimilarity, topMotifPairSimilarity }, { 8 },
        { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorHorizontalThree.printTimeSeriesHorizontal(testTimeSeries,
        { topMotifPairSimilarity }, {}, { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorHorizontalThree.printTimeSeriesHorizontal(testTimeSeries,
        { topMotifPairSimilarity }, { 8, 8 }, { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorHorizontalThree.printTimeSeriesHorizontal(testTimeSeries,
        { topMotifPairSimilarity }, { 8 }, { {} });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  //remove folder
  TEST(remove("time_series_" + to_string(folderNumber)));


  TestOutputGenerator generatorVerticalThree;

  TEST(generatorVerticalThree.getOutputFolderNumber() == folderNumber);
  TEST(generatorVerticalThree.getOutputFileName() == ("./time_series_"
        + to_string(folderNumber) + "/time_series"));
  TEST(generatorVerticalThree.getOutputFileMotifSetsName() == ("./time_series_"
        + to_string(folderNumber) + "/time_series_meta"));
  TEST(generatorVerticalThree.getOutputFileGNUPlotScriptName() ==
      ("./time_series_" + to_string(folderNumber) + "/time_series_plot"));
  TEST(generatorVerticalThree.getTimeSeriesName() == ("Time Series"));

  //test if output generator throw exceptions correctly
  try {

    generatorVerticalThree.printTimeSeriesVertical({}, { topMotifPairSimilarity
        }, { 8 }, { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorVerticalThree.printTimeSeriesVertical(testTimeSeries, {},
      { 8 }, { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorVerticalThree.printTimeSeriesVertical(testTimeSeries,
        { topMotifPairSimilarity, topMotifPairSimilarity }, { 8 },
        { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorVerticalThree.printTimeSeriesVertical(testTimeSeries,
        { topMotifPairSimilarity }, {}, { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorVerticalThree.printTimeSeriesVertical(testTimeSeries,
        { topMotifPairSimilarity }, { 8, 8 }, { topMotifPairPos });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  try {

    generatorVerticalThree.printTimeSeriesVertical(testTimeSeries,
        { topMotifPairSimilarity }, { 8 }, { {} });
    TEST(!"Has to throw an error");
  }
  catch (int e) {

    TEST(e == EXIT_FAILURE);
  }

  //remove folder
  TEST(remove("time_series_" + to_string(folderNumber)));


  //reset cerr
  cerr.rdbuf(cerr_buff);
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

    TEST_SECTION("global units");
    test_check_if_float();
    test_parseArgs();
    test_checkArg();
    test_print_version();
    test_print_help();

    TEST_SECTION("tsgenerator units");
    test_motifsetcollection();
    test_outputgenerator();
    test_freepositions();
    test_basets();
    test_tsgenerator();

    TEST_SUMMARY;
  }
  catch (...) {
    TEST_EXCEPTION;
  }

  return EXIT_SUCCESS;
}
