///\file global.cpp
///
///\brief File contains the global function definitions.
///
///This is the header file of the global function definitions. It contains
///several functions to parse data.

#include<global.hpp>

bool check_if_float(const tsg::word string_in) {

  bool isFloat = false;
  int stringItr = 0;

  if (string_in[stringItr] == '-')
    stringItr++;

  if (stringItr < (int)string_in.length() && '0' <= string_in[stringItr] &&
      string_in[stringItr] <= '9') {

    while (stringItr < (int)string_in.length() && '0' <= string_in[stringItr]
        && string_in[stringItr] <= '9')
      stringItr++;

    if (stringItr == (int)string_in.length())
      isFloat = true;
    else {

      if (string_in[stringItr] == '.') {

        stringItr++;

        if (stringItr < (int)string_in.length() && '0' <= string_in[stringItr]
            && string_in[stringItr] <= '9') {

          while (stringItr < (int)string_in.length() && '0' <=
              string_in[stringItr] && string_in[stringItr] <= '9')
            stringItr++;

          if (stringItr == (int)string_in.length())
            isFloat = true;
        }
      }
    }
  }

  return isFloat;
}

void parseArgs(int argc, char **argv, tsg::par &argTokens) {

  for (int i = 0; i < argc; i++) {

    tsg::word token(argv[i]);

    if (token != "")
      argTokens.push_back(token);
  }
}

bool checkArg(tsg::par &argTokens, const tsg::word arg_in, tsg::par
    &payload_out) {

  if (!payload_out.empty()) {

    payload_out.clear();
    payload_out.resize(0);
  }

  tsg::par::iterator payload = std::find(argTokens.begin(),
      argTokens.end(), arg_in);

  if (payload == argTokens.end() || check_if_float(*payload))
    return false;

  payload++;

  while (payload != argTokens.end() && ((*payload).front() != '-' ||
        check_if_float(*payload))) {
    if (*payload != "")
      payload_out.push_back(*payload);

    payload++;
  }

  return true;
}

void print_version() {

  std::cout << PROGNAME << " " << VERSION << std::endl;

  std::cout << "Copyright \302\251 2018 Rafael Moczalla" << std::endl;
  std::cout << "Lizenz: Creative Commons - Attribution - Non-Commercial" <<
    " - Share Alike" << std::endl;
  std::cout << "There are no guarantees as far as the law permits." <<
    std::endl << std::endl;
  std::cout << "Written by Rafael Moczalla" << std::endl;
}

void print_help() {

  std::cout << "Call: " << std::endl;
  std::cout << "    " << PROGNAME << " [Options]" << std::endl << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "    -g,    --generator tsg::word  " <<
    "         Sets the generator, the method for injecting sequences" <<
    std::endl;
  std::cout << "                                  " <<
    "       into the time series matching the synthetic motif. The" <<
    std::endl;
  std::cout << "                                  " <<
    "       available methods are pair motif, set motif and latent" <<
    std::endl;
  std::cout << "                                  " <<
    "       motif." << std::endl;
  std::cout << "    -ty,   --type tsg::word       " <<
    "         Sets the motif type, the shape of the injected motif and" <<
    std::endl;
  std::cout << "                                  " <<
    "       the inserted sequences. The available methods are box," <<
    std::endl;
  std::cout << "                                  " <<
    "       triangle, semicircle, trapezoid, positiveflank," << std::endl;
  std::cout << "                                  " <<
    "       negativeflank, sine and cosine." << std::endl;
  std::cout << "    -me,   --method tsg::word     " <<
    "         Sets the motif type, the shape of the injected motif and" <<
    std::endl;
  std::cout << "                                  " <<
    "       the inserted sequences. The available methods are" << std::endl;
  std::cout << "                                  " <<
    "       simpleRandomWalk, realRandomWalk, normalRandomWalk," << std::endl;
  std::cout << "                                  " <<
    "       linearRandomWalk, boundedSimpleRandomWalk," << std::endl;
  std::cout << "                                  " <<
    "       boundedRealRandomWalk, boundedNormalRandomWalk," << std::endl;
  std::cout << "                                  " <<
    "       boundedLinearRandomWalk, uniformRandom, normalRandom," <<
    std::endl;
  std::cout << "                                  " <<
    "       piecewiseLinearRandom and splineRepeated." << std::endl;
  std::cout << "    -l,    --length INTEGER       " <<
    "         Sets the length of the time series." << std::endl;
  std::cout << "    -w,    --windowSize INTEGER   " <<
    "         Sets the window size of the time series." << std::endl;
  std::cout << "    -si,   --size INTEGER         " <<
    "         Sets the size of the motif, the number of inserted" <<
    std::endl;
  std::cout << "                                  " <<
    "       sequences non-self matched by the motif." << std::endl;
  std::cout << "    -no,   --noise FLOAT          " <<
    "         Sets the noise value. A random value in the range from" <<
    std::endl;
  std::cout << "                                  " <<
    "       -FLOAT to FLOAT is added to the base time series and the" <<
    std::endl;
  std::cout << "                                  " <<
    "       inserted sequences." << std::endl;
  std::cout << "    -d,    --delta FLOAT          " <<
    "         Sets the maximum absolute difference between two" << std::endl;
  std::cout << "                                  " <<
    "       consecutive values in the time series." << std::endl;
  std::cout << "    -he,   --delta FLOAT          " <<
    "         Sets the maximum absolute difference between two values" <<
    std::endl;
  std::cout << "                                  " <<
    "       of the base motif." << std::endl;
  std::cout << "    -st,   --step FLOAT           " <<
    "         Sets the maximum step size in x direction from two" << std::endl;
  std::cout << "                                  " <<
    "       consecutive values when creating a splined base times" <<
    std::endl;
  std::cout << "                                  " <<
    "       series." << std::endl;
  std::cout << "    -ti,   --times INTEGER        " <<
    "         Sets the number of values computed to generate a" << std::endl;
  std::cout << "                                  " <<
    "       repeating" << std::endl;
  std::cout << "                                  " <<
    "       pattern when generating a splined base time series." << std::endl;
  std::cout << "    -ma,   --maxi FLOAT           " <<
    "         Sets the maximum absolute value in the base times series." <<
    std::endl;
  std::cout << "    -o,    --out NAME             " <<
    "         Sets the output file name." << std::endl;
  std::cout << "    -tsn,  --timeSeriesName NAME  " <<
    "         Sets the time series name." << std::endl;
  std::cout << "    -ho,   --horizontalOutput     " <<
    "         Sets the output mode to horizontal." << std::endl;
  std::cout << "    -r,    --range FLOAT FLOAT    " <<
    "         Sets the range of the time series values." << std::endl;
  std::cout << "    -h,    --help                 " <<
    "         Prints these help text." << std::endl;
  std::cout << "    -v,    --version              " <<
    "         Prints the version information." << std::endl;
}
