///\file unittests.cpp
///
///\brief File contains all unit test.
///
///This is the source file to run all unit tests with stest.

#include <unittests.hpp>


void test_check_if_float() {

  TEST_GROUP_FUNCTION;
  TEST_R(check_if_float("123456789"));
  TEST_R(check_if_float("-123456789"));
  TEST_R(check_if_float("123456789.123456789"));
  TEST_R(check_if_float("-123456789.123456789"));
  TEST_R(!check_if_float("abcd"));
  TEST_R(!check_if_float(""));
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

  TEST_R(argTokens[0] == "--arg0");
  TEST_R(argTokens[1] == "payload");
  TEST_R(argTokens[2] == "-arg2");
  TEST_R(argTokens[3] == "-1000.1");
  TEST_R(argTokens[4] == "-empty");
  TEST_R(argTokens[5] == "--multiplePayload");
  TEST_R(argTokens[6] == "-100.000");
  TEST_R(argTokens[7] == "tsg::word1 tsg::word2");
  TEST_R(argTokens[8] == "end");
}

void test_checkArg() {

  TEST_GROUP_FUNCTION;

  TEST_R(!checkArg(argTokens, "", payload));
  TEST_R(!checkArg(argTokens, "-1000.1", payload));
  TEST_R(checkArg(argTokens, "-empty", payload));
  TEST_R(payload.empty());
  TEST_R(checkArg(argTokens, "--arg0", payload));
  TEST_R(payload.size() == 1);
  TEST_R(payload[0] == "payload");
  TEST_R(checkArg(argTokens, "-arg2", payload));
  TEST_R(payload.size() == 1);
  TEST_R(payload[0] == "-1000.1");
  TEST_R(checkArg(argTokens, "--multiplePayload", payload));
  TEST_R(payload.size() == 3);
  TEST_R(payload[0] == "-100.000");
  TEST_R(payload[1] == "tsg::word1 tsg::word2");
  TEST_R(payload[2] == "end");
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

  TEST_R(testStream.str() == local.str());
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
  testStream << "    -ho,   --horizontalOutput     " <<
    "         Sets the output mode to horizontal." << std::endl;
  testStream << "    -r,    --range FLOAT FLOAT    " <<
    "         Sets the range of the time series values." << std::endl;
  testStream << "    -h,    --help                 " <<
    "         Prints these help text." << std::endl;
  testStream << "    -v,    --version              " <<
    "         Prints the version information." << std::endl;

  std::cout.rdbuf(cout_buff);

  TEST_R(testStream.str() == local.str());
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

  generator.open();
  generator.close();

  TEST_R(generator.getFolderNumber() == secondFolderNumber);
  TEST_R(generator.getFileName() == ("./time_series_"
        + std::to_string(secondFolderNumber) + "/time_series"));
  TEST_R(generator.getMetaFileName() == ("./time_series_"
        + std::to_string(secondFolderNumber) + "/time_series_meta"));
  TEST_R(std::filesystem::remove("time_series_"
        + std::to_string(secondFolderNumber) + "/time_series_"
        + std::to_string(secondFolderNumber) + ".csv"));
  TEST_R(std::filesystem::remove("time_series_"
        + std::to_string(secondFolderNumber)
        + "/time_series_meta_" + std::to_string(secondFolderNumber) + ".csv"));
  TEST_R(std::filesystem::remove("time_series_"
        + std::to_string(secondFolderNumber)));

  generator.open();
  generator.printTimeSeriesHorizontal(testTimeSeries);
  generator.close();

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
    TEST_R(line == correctStream.str());

  TEST_R(!getline(iTestTimeSeries, line));

  //remove folder
  TEST_R(std::filesystem::remove("time_series_" + std::to_string(folderNumber)));

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
