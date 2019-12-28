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
tsg::par argTokens;
tsg::par payload;

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
  char arg7[] = "tsg::word1 tsg::word2";
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
  TEST(argTokens[7] == "tsg::word1 tsg::word2");
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
  TEST(payload[1] == "tsg::word1 tsg::word2");
  TEST(payload[2] == "end");
}

void test_print_version() {

  TEST_GROUP_FUNCTION;

  std::ostringstream local;
  std::ostringstream testStream;

  auto cout_buff = std::cout.rdbuf(local.rdbuf());

  print_version();

  testStream << PROGNAME << " " << VERSION << std::endl;

  testStream << "Copyright \302\251 2018 Rafael Moczalla" << std::endl;
  testStream << "Lizenz: Creative Commons - Attribution" <<
    " - Non-Commercial - Share Alike" << std::endl;
  testStream << "There are no guarantees as far as the law permits." <<
    std::endl << std::endl;
  testStream << "Written by Rafael Moczalla" << std::endl;

  std::cout.rdbuf(cout_buff);

  TEST(testStream.str() == local.str());
}

void test_print_help() {

  TEST_GROUP_FUNCTION;

  std::ostringstream local;
  std::ostringstream testStream;

  auto cout_buff = std::cout.rdbuf(local.rdbuf());

  print_help();

  testStream << "Call: " << std::endl;
  testStream << "    " << PROGNAME << " [Options]" << std::endl << std::endl;
  testStream << "Options:" << std::endl;
  testStream << "    -g,    --generator tsg::word  " <<
    "         Sets the generator, the method for injecting sequences" <<
    std::endl;
  testStream << "                                  " <<
    "       into the time series matching the synthetic motif. The" <<
    std::endl;
  testStream << "                                  " <<
    "       available methods are pair motif, set motif and latent" <<
    std::endl;
  testStream << "                                  " <<
    "       motif." << std::endl;
  testStream << "    -ty,   --type tsg::word       " <<
    "         Sets the motif type, the shape of the injected motif and" <<
    std::endl;
  testStream << "                                  " <<
    "       the inserted sequences. The available methods are box," <<
    std::endl;
  testStream << "                                  " <<
    "       triangle, semicircle, trapezoid, positiveflank," << std::endl;
  testStream << "                                  " <<
    "       negativeflank, sine and cosine." << std::endl;
  testStream << "    -me,   --method tsg::word     " <<
    "         Sets the motif type, the shape of the injected motif and" <<
    std::endl;
  testStream << "                                  " <<
    "       the inserted sequences. The available methods are" << std::endl;
  testStream << "                                  " <<
    "       simpleRandomWalk, realRandomWalk, normalRandomWalk," << std::endl;
  testStream << "                                  " <<
    "       linearRandomWalk, boundedSimpleRandomWalk," << std::endl;
  testStream << "                                  " <<
    "       boundedRealRandomWalk, boundedNormalRandomWalk," << std::endl;
  testStream << "                                  " <<
    "       boundedLinearRandomWalk, uniformRandom, normalRandom," <<
    std::endl;
  testStream << "                                  " <<
    "       piecewiseLinearRandom and splineRepeated." << std::endl;
  testStream << "    -l,    --length INTEGER       " <<
    "         Sets the length of the time series." << std::endl;
  testStream << "    -w,    --windowSize INTEGER   " <<
    "         Sets the window size of the time series." << std::endl;
  testStream << "    -si,   --size INTEGER         " <<
    "         Sets the size of the motif, the number of inserted" <<
    std::endl;
  testStream << "                                  " <<
    "       sequences non-self matched by the motif." << std::endl;
  testStream << "    -no,   --noise FLOAT          " <<
    "         Sets the noise value. A random value in the range from" <<
    std::endl;
  testStream << "                                  " <<
    "       -FLOAT to FLOAT is added to the base time series and the" <<
    std::endl;
  testStream << "                                  " <<
    "       inserted sequences." << std::endl;
  testStream << "    -d,    --delta FLOAT          " <<
    "         Sets the maximum absolute difference between two" << std::endl;
  testStream << "                                  " <<
    "       consecutive values in the time series." << std::endl;
  testStream << "    -he,   --delta FLOAT          " <<
    "         Sets the maximum absolute difference between two values" <<
    std::endl;
  testStream << "                                  " <<
    "       of the base motif." << std::endl;
  testStream << "    -st,   --step FLOAT           " <<
    "         Sets the maximum step size in x direction from two" << std::endl;
  testStream << "                                  " <<
    "       consecutive values when creating a splined base times" <<
    std::endl;
  testStream << "                                  " <<
    "       series." << std::endl;
  testStream << "    -ti,   --times INTEGER        " <<
    "         Sets the number of values computed to generate a" << std::endl;
  testStream << "                                  " <<
    "       repeating" << std::endl;
  testStream << "                                  " <<
    "       pattern when generating a splined base time series." << std::endl;
  testStream << "    -ma,   --maxi FLOAT           " <<
    "         Sets the maximum absolute value in the base times series." <<
    std::endl;
  testStream << "    -o,    --out NAME             " <<
    "         Sets the output file name." << std::endl;
  testStream << "    -tsn,  --timeSeriesName NAME  " <<
    "         Sets the time series name." << std::endl;
  testStream << "    -ho,   --horizontalOutput     " <<
    "         Sets the output mode to horizontal." << std::endl;
  testStream << "    -r,    --range FLOAT FLOAT    " <<
    "         Sets the range of the time series values." << std::endl;
  testStream << "    -h,    --help                 " <<
    "         Prints these help text." << std::endl;
  testStream << "    -v,    --version              " <<
    "         Prints the version information." << std::endl;

  std::cout.rdbuf(cout_buff);

  TEST(testStream.str() == local.str());
}

void test_outputgenerator() {

  //redirect std::cerr to keep test output clean
  std::ostringstream local;

  auto cerr_buff = std::cerr.rdbuf(local.rdbuf());

  TEST_GROUP_FUNCTION;
  int folderNumber = 0;
  int secondFolderNumber;

  while(std::filesystem::exists("time_series_" + std::to_string(folderNumber)))
    folderNumber++;

  std::filesystem::create_directory("time_series_"
      + std::to_string(folderNumber));

  secondFolderNumber = folderNumber + 1;

  while(std::filesystem::exists("time_series_"
        + std::to_string(secondFolderNumber)))
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
        + std::to_string(secondFolderNumber) + "/time_series"));
  TEST(generator.getOutputFileMotifSetsName() == ("./time_series_"
        + std::to_string(secondFolderNumber) + "/time_series_meta"));
  TEST(generator.getOutputFileGNUPlotScriptName() == ("./time_series_"
        + std::to_string(secondFolderNumber) + "/time_series_plot"));
  TEST(generator.getTimeSeriesName() == ("Test Name"));
  TEST(std::filesystem::remove("time_series_"
        + std::to_string(secondFolderNumber) + "/time_series_"
        + std::to_string(secondFolderNumber) + ".csv"));
  TEST(std::filesystem::remove("time_series_"
        + std::to_string(secondFolderNumber)
        + "/time_series_meta_" + std::to_string(secondFolderNumber) + ".csv"));
  TEST(std::filesystem::remove("time_series_"
        + std::to_string(secondFolderNumber)
        + "/time_series_plot_" + std::to_string(secondFolderNumber) + ".plt"));
  TEST(std::filesystem::remove("time_series_"
        + std::to_string(secondFolderNumber)));

  std::filesystem::remove("time_series_" + std::to_string(folderNumber));


  TestOutputGenerator generatorHorizontal;

  TEST(generatorHorizontal.getOutputFolderNumber() == folderNumber);
  TEST(generatorHorizontal.getOutputFileName() == ("./time_series_"
        + std::to_string(folderNumber) + "/time_series"));
  TEST(generatorHorizontal.getOutputFileMotifSetsName() == ("./time_series_"
        + std::to_string(folderNumber) + "/time_series_meta"));
  TEST(generatorHorizontal.getOutputFileGNUPlotScriptName() ==
      ("./time_series_" + std::to_string(folderNumber) + "/time_series_plot"));
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
  TEST(std::filesystem::remove("time_series_" + std::to_string(folderNumber)));


  TestOutputGenerator generatorVertical;

  TEST(generatorVertical.getOutputFolderNumber() == folderNumber);
  TEST(generatorVertical.getOutputFileName() == ("./time_series_"
        + std::to_string(folderNumber) + "/time_series"));
  TEST(generatorVertical.getOutputFileMotifSetsName() == ("./time_series_"
        + std::to_string(folderNumber) + "/time_series_meta"));
  TEST(generatorVertical.getOutputFileGNUPlotScriptName() == ("./time_series_"
        + std::to_string(folderNumber) + "/time_series_plot"));
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
  TEST(std::filesystem::remove("time_series_" + std::to_string(folderNumber)));


  TestOutputGenerator generatorHorizontalTwo;

  TEST(generatorHorizontalTwo.getOutputFolderNumber() == folderNumber);
  TEST(generatorHorizontalTwo.getOutputFileName() == ("./time_series_"
        + std::to_string(folderNumber) + "/time_series"));
  TEST(generatorHorizontalTwo.getOutputFileMotifSetsName() == ("./time_series_"
        + std::to_string(folderNumber) + "/time_series_meta"));
  TEST(generatorHorizontalTwo.getOutputFileGNUPlotScriptName() ==
      ("./time_series_" + std::to_string(folderNumber) + "/time_series_plot"));
  TEST(generatorHorizontalTwo.getTimeSeriesName() == ("Time Series"));

  generatorHorizontalTwo.printTimeSeriesHorizontal(testTimeSeries,
      { topMotifSetRange, topMotifPairSimilarity }, { windowSize, windowSize },
      { topMotifSetPos, topMotifPairPos });

  //test if time series file has correct content
  std::ostringstream correctStream;
  correctStream << std::fixed;
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

  std::ifstream iTestTimeSeries("time_series_" + std::to_string(folderNumber)
      + "/time_series_" + std::to_string(folderNumber) + ".csv");

  tsg::word line;

  if (TEST_IF(getline(iTestTimeSeries, line) ? true : false))
    TEST(line == correctStream.str());

  TEST(!getline(iTestTimeSeries, line));

  iTestTimeSeries.close();

  //test if motif positions file has correct content
  std::ifstream iTestPositions("time_series_" + std::to_string(folderNumber)
      + "/time_series_meta_" + std::to_string(folderNumber) + ".csv");

  correctStream << std::defaultfloat;
  correctStream.str("");
  correctStream.clear();

  correctStream << " , \"range/similarity\", \"window size\"," <<
    " \"position 0\", \"position 1\", \"position 2\"";

  if(TEST_IF(getline(iTestPositions, line) ? true : false))
    TEST(line == correctStream.str());

  correctStream.str("");
  correctStream.clear();

  correctStream << std::fixed << "\"Top Latent Motif Locations\", " <<
    0.240684 << ", " << windowSize << ", ";
  correctStream << std::defaultfloat << topMotifSetPos[0] << ", " <<
    topMotifSetPos[1] << ", " << topMotifSetPos[2];

  if(TEST_IF(getline(iTestPositions, line) ? true : false))
    TEST(line == correctStream.str());

  correctStream.str("");
  correctStream.clear();

  correctStream << std::fixed << "\"Top Pair Motif Locations\", " <<
    0.000000 << ", " << windowSize << ", ";
  correctStream << std::defaultfloat << topMotifPairPos[0] << ", " <<
    topMotifPairPos[1];

  if(TEST_IF(getline(iTestPositions, line) ? true : false))
    TEST(line == correctStream.str());

  TEST(!getline(iTestPositions, line));

  iTestPositions.close();

  std::ifstream iTestPlot("time_series_" + std::to_string(folderNumber)
      + "/time_series_plot_" + std::to_string(folderNumber) + ".plt");

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
    TEST(line == ("plot \"time_series_" + std::to_string(folderNumber)
          + ".csv\" matrix title \"Time Series " + std::to_string(folderNumber)
          + "\" with lines"));
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
  TEST(std::filesystem::remove("time_series_" + std::to_string(folderNumber)
        + "/time_series_"
        + std::to_string(folderNumber) + ".csv"));
  TEST(std::filesystem::remove("time_series_" + std::to_string(folderNumber)
        + "/time_series_meta_"
        + std::to_string(folderNumber) + ".csv"));
  TEST(std::filesystem::remove("time_series_" + std::to_string(folderNumber)
        + "/time_series_plot_"
        + std::to_string(folderNumber) + ".plt"));
  TEST(std::filesystem::remove("time_series_" + std::to_string(folderNumber)));


  TestOutputGenerator generatorVerticalTwo;

  TEST(generatorVerticalTwo.getOutputFolderNumber() == folderNumber);
  TEST(generatorVerticalTwo.getOutputFileName() == ("./time_series_"
        + std::to_string(folderNumber) + "/time_series"));
  TEST(generatorVerticalTwo.getOutputFileMotifSetsName() == ("./time_series_"
        + std::to_string(folderNumber) + "/time_series_meta"));
  TEST(generatorVerticalTwo.getOutputFileGNUPlotScriptName() ==
      ("./time_series_" + std::to_string(folderNumber) + "/time_series_plot"));
  TEST(generatorVerticalTwo.getTimeSeriesName() == ("Time Series"));

  generatorVerticalTwo.printTimeSeriesVertical(testTimeSeries,
      { topMotifSetRange, topMotifPairSimilarity }, { windowSize, windowSize },
      { topMotifSetPos, topMotifPairPos });

  //test if time series file has correct content
  correctStream.str("");
  correctStream.clear();
  correctStream << std::fixed;

  iTestTimeSeries.open("time_series_" + std::to_string(folderNumber)
      + "/time_series_" + std::to_string(folderNumber) + ".csv");

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
  iTestPositions.open("time_series_" + std::to_string(folderNumber)
      + "/time_series_meta_" + std::to_string(folderNumber) + ".csv");

  correctStream << std::defaultfloat;
  correctStream.str("");
  correctStream.clear();

  correctStream << " , \"Top Latent Motif Locations\"," <<
    " \"Top Pair Motif Locations\"";
  if(TEST_IF(getline(iTestPositions, line) ? true : false))
    TEST(line == correctStream.str());

  correctStream.str("");
  correctStream.clear();

  correctStream << std::fixed << "\"range/similarity\", " <<
    0.240684 << ", " << 0.000000 << std::defaultfloat;

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

  iTestPlot.open("time_series_" + std::to_string(folderNumber)
      + "/time_series_plot_" + std::to_string(folderNumber) + ".plt");

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
    TEST(line == ("plot \"time_series_" + std::to_string(folderNumber)
          + ".csv\" title \"Time Series " + std::to_string(folderNumber)
          + "\" with lines"));
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
  TEST(std::filesystem::remove("time_series_" + std::to_string(folderNumber)
        + "/time_series_"
        + std::to_string(folderNumber) + ".csv"));
  TEST(std::filesystem::remove("time_series_" + std::to_string(folderNumber)
        + "/time_series_meta_"
        + std::to_string(folderNumber) + ".csv"));
  TEST(std::filesystem::remove("time_series_" + std::to_string(folderNumber)
        + "/time_series_plot_"
        + std::to_string(folderNumber) + ".plt"));
  TEST(std::filesystem::remove("time_series_" + std::to_string(folderNumber)));


  TestOutputGenerator generatorHorizontalThree;

  TEST(generatorHorizontalThree.getOutputFolderNumber() == folderNumber);
  TEST(generatorHorizontalThree.getOutputFileName() == ("./time_series_"
        + std::to_string(folderNumber) + "/time_series"));
  TEST(generatorHorizontalThree.getOutputFileMotifSetsName() ==
      ("./time_series_" + std::to_string(folderNumber) + "/time_series_meta"));
  TEST(generatorHorizontalThree.getOutputFileGNUPlotScriptName() ==
      ("./time_series_" + std::to_string(folderNumber) + "/time_series_plot"));
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
  TEST(std::filesystem::remove("time_series_" + std::to_string(folderNumber)));


  TestOutputGenerator generatorVerticalThree;

  TEST(generatorVerticalThree.getOutputFolderNumber() == folderNumber);
  TEST(generatorVerticalThree.getOutputFileName() == ("./time_series_"
        + std::to_string(folderNumber) + "/time_series"));
  TEST(generatorVerticalThree.getOutputFileMotifSetsName() == ("./time_series_"
        + std::to_string(folderNumber) + "/time_series_meta"));
  TEST(generatorVerticalThree.getOutputFileGNUPlotScriptName() ==
      ("./time_series_" + std::to_string(folderNumber) + "/time_series_plot"));
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
  TEST(std::filesystem::remove("time_series_" + std::to_string(folderNumber)));


  //reset cerr
  std::cerr.rdbuf(cerr_buff);
}

int main() {

  try {
    TEST_WAIT(false);

    TEST_SECTION("check if float function");
    test_check_if_float();
    TEST_SECTION("parse arguments function");
    test_parseArgs();
    TEST_SECTION("check arguments function");
    test_checkArg();
    TEST_SECTION("print version function");
    test_print_version();
    TEST_SECTION("print help function");
    test_print_help();

    TEST_SUMMARY;
  }
  catch (...) {
    TEST_EXCEPTION;
  }

  return EXIT_SUCCESS;
}
