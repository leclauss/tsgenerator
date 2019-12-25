///\file global.cpp
///
///\brief File contains the global function definitions.
///
///This is the header file of the global function definitions. It contains
///several functions to parse data.

#include<global.hpp>

bool check_if_float(const string string_in) {

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

void parseArgs(int argc, char **argv, vector<string> &argTokens) {

  for (int i = 0; i < argc; i++) {

    string token(argv[i]);

    if (token != "")
      argTokens.push_back(token);
  }
}

bool checkArg(vector<string> &argTokens, const string arg_in, vector<string>
    &payload_out) {

  if (!payload_out.empty()) {

    payload_out.clear();
    payload_out.resize(0);
  }

  vector<string>::iterator payload = find(argTokens.begin(), argTokens.end(),
      arg_in);

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

  cout << PROGNAME << " " << VERSION << endl;

  cout << "Copyright \302\251 2018 Rafael Moczalla" << endl;
  cout << "Lizenz: Creative Commons - Attribution - Non-Commercial" <<
    " - Share Alike" << endl;
  cout << "There are no guarantees as far as the law permits." << endl
    << endl;
  cout << "Written by Rafael Moczalla" << endl;
}

void print_help() {

  cout << "Call: " << endl;
  cout << "    " << PROGNAME << " -l INTEGER -w INTEGER -rd FLOAT" <<
    " [Options]" << endl;
  cout << endl;
  cout << "    -l,    --length INTEGER                           " <<
    "       Sets the length of the time series." << endl;
  cout << "    -w,    --windowSize INTEGER                       " <<
    "       Sets the window size of the time series." << endl;
  cout << "    -rd,   --randomness FLOAT                         " <<
    "       Sets the randomness of the time series." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "    -lm,   --latentMotif (STRING INTEGER FLOAT)+      " <<
    "       Sets the operating mode to latent motif and" << endl;
  cout << "                                                      " <<
    "       the parameters for the latent motif namely" << endl;
  cout << "                                                      " <<
    "       type, count and height." << endl;
  cout << "    -o,    --out NAME                                 " <<
    "       Sets the output file name." << endl;
  cout << "    -tsn,  --timeSeriesName NAME                      " <<
    "       Sets the time series name." << endl;
  cout << "    -ho,   --horizontalOutput                         " <<
    "       Sets the output mode to horizontal." << endl;
  cout << "    -r,    --range FLOAT FLOAT                        " <<
    "       Sets the range of the time series values." << endl;
  cout << "    -h,    --help                                     " <<
    "       Prints these help text." << endl;
  cout << "    -v,    --version                                  " <<
    "       Prints the version information." << endl;
}
