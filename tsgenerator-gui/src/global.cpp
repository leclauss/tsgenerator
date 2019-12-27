///\file global.cpp
///
///\brief File contains the global function definitions.
///
///This is the header file of the global function definitions. It contains
///several functions to parse data.

#include<global.hpp>

bool check_if_float(const std::string string_in) {

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

