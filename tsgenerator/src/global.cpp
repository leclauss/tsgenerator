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
  std::cout << "    " << PROGNAME << " -l INTEGER -w INTEGER -rd FLOAT" <<
    " [Options]" << std::endl;
  std::cout << std::endl;
  std::cout << "    -l,    --length INTEGER                           " <<
    "       Sets the length of the time series." << std::endl;
  std::cout << "    -w,    --windowSize INTEGER                       " <<
    "       Sets the window size of the time series." << std::endl;
  std::cout << "    -rd,   --randomness FLOAT                         " <<
    "       Sets the randomness of the time series." << std::endl;
  std::cout << std::endl;
  std::cout << "Options:" << std::endl;
  std::cout << "    -lm,   --latentMotif (STRING INTEGER FLOAT)+      " <<
    "       Sets the operating mode to latent motif and" << std::endl;
  std::cout << "                                                      " <<
    "       the parameters for the latent motif namely" << std::endl;
  std::cout << "                                                      " <<
    "       type, count and height." << std::endl;
  std::cout << "    -o,    --out NAME                                 " <<
    "       Sets the output file name." << std::endl;
  std::cout << "    -tsn,  --timeSeriesName NAME                      " <<
    "       Sets the time series name." << std::endl;
  std::cout << "    -ho,   --horizontalOutput                         " <<
    "       Sets the output mode to horizontal." << std::endl;
  std::cout << "    -r,    --range FLOAT FLOAT                        " <<
    "       Sets the range of the time series values." << std::endl;
  std::cout << "    -h,    --help                                     " <<
    "       Prints these help text." << std::endl;
  std::cout << "    -v,    --version                                  " <<
    "       Prints the version information." << std::endl;
}
