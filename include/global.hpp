///\file main.hpp
///
///\brief File contains the global declarations.
///
///This is the source file of the global declarations. It contains several
///functions to parse data and macros.

#ifndef MAIN_HPP
#define MAIN_HPP

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;


// ----------------------------------------------------------------------------
///\brief Defines the program name.
///
///Makro that defines the name of the program.
// ----------------------------------------------------------------------------
#define PROGNAME "TSGenerator"

// ----------------------------------------------------------------------------
///\brief Defines the program version.
///
///Makro that defines the version of the program.
// ----------------------------------------------------------------------------
#define VERSION  "3.0.0.8"

// ----------------------------------------------------------------------------
///\brief Checks if string represents float.
///
///\param [in] string_in The string to be checked.
///
///\return true if string is float otherwise false.
///
///This function checks if the string is a float. Float representation of form
///[-](123456789)+[.(123456789)+] are accepted.
// ----------------------------------------------------------------------------
bool check_if_float(const string string_in);

// ----------------------------------------------------------------------------
///\brief Parses arguments into a vector.
///
///\param [in] argc Hands over the number of arguments.
///\param [in] **argv Hands over the arguments.
///\param [out] *argTokens Returns the parsed arguments vector.
///
///This is a argument parser function that parses arguments into a string
///vector.
// ----------------------------------------------------------------------------
void parseArgs(int argc, char **argv, vector<string> &argTokens);

// ----------------------------------------------------------------------------
///\brief Checks whether argument is contained in argument vector or not.
///
///\param [in] *argTokens Hands over the parsed arguments vector.
///\param [in] arg_in Hands over the argument token.
///\param [out] *payload_out Returns the parsed payload.
///
///\return true if argument was found otherwise false.
///
///Function that checks whether argument is contained in the parsed argument
///vector or not. Additional it returns the payload, if available. Returns
///true, if the argument was found and sets the payload string. Otherwise the
///function returns false.
// ----------------------------------------------------------------------------
bool checkArg(vector<string> &argTokens, const string arg_in, vector<string>
    &payload_out);

// ----------------------------------------------------------------------------
///\brief Prints out the version text.
///
///A text containing the program name and version as well as additional
///informations like the authors name.
// ----------------------------------------------------------------------------
void print_version();

// ----------------------------------------------------------------------------
///\brief Prints out the help text.
///
///A text containing the information about the use and the arguments of the
///TSGenerator.
// ----------------------------------------------------------------------------
void print_help();

#endif
