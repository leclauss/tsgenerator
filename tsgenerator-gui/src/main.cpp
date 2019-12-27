///\file main.cpp
///
///\brief Here starts the TSGenerator.
///
///This is the main.cpp file of the TSGenerator. The program uses a
///configuration file for the time series configuration.

#ifdef _WIN32
#include <Windows.h>
#endif

#include <moc_tsggui.cpp>


///\brief Function called at program start.
///
///\param [in] argc Hands over the number of arguments to main.
///\param [in] *argv[] Hands over the arguments to main.
///
///This is the typical main function of C++. The main function is called at the
///program start. First all non-local objects with static storage duration are
///initialized. We use the main function to parse the argument list to set up
///and to start the TSGenerator. Also a configuration file is created, if
///a configuration file do not exist.
int main(int argc, char *argv[]) {

#ifdef _WIN32
  FreeConsole();
#endif

  int ret;

  {
    //run gui
    TsgGui tsggui(argc, argv);

    ret = tsggui.exec();
  }

  return ret;
}
