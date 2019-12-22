///\file motifsetcollection.cpp
///
///\brief File contains the motif set collection definitions.
///
///This is the source file of the motif set collection. The motif set collection consists of motif set class definitions.

#include <motifsetcollection.hpp>


void generateBoxMotif(vector<double> &subsequence_out, const double &warp_in,
    const double &scale_in, const int &length_in, const double &height_in) {

  if (!subsequence_out.empty()) {

    subsequence_out.clear();
    subsequence_out.resize(0);
  }

  double warp = warp_in;
  double scale = scale_in;
  int length = length_in;
  double height = height_in;

  int realLength = (int)(warp * length);

  subsequence_out.push_back(0.0);

  if (realLength > 2) {

    for (int itr = 0; itr < realLength - 2; itr++)
      subsequence_out.push_back(scale * height);
  }

  subsequence_out.push_back(0.0);
}

void generateTriangleMotif(vector<double> &subsequence_out, const double
    &warp_in, const double &scale_in, const int &length_in, const double
    &height_in) {

  if (!subsequence_out.empty()) {

    subsequence_out.clear();
    subsequence_out.resize(0);
  }

  double warp = warp_in;
  double scale = scale_in;
  int length = length_in;
  double height = height_in;

  int realLength = (int)(warp * length);

  subsequence_out.push_back(0.0);

  if (realLength > 2) {

    double increase = scale * height / ((realLength - 1.0) / 2.0);

    for (int itr = 2; itr < realLength - 1; itr += 2)
      subsequence_out.push_back(increase * (itr / 2.0));

    for (int itr = realLength - (realLength % 2 == 0 ? 2 : 1);
        itr > 1;
        itr -= 2)
      subsequence_out.push_back(increase * (itr / 2.0));
  }

  subsequence_out.push_back(0.0);
}

void generateSemicircleMotif(vector<double> &subsequence_out, const double
    &warp_in, const double &scale_in, const int &length_in, const double
    &height_in) {

  if (!subsequence_out.empty()) {

    subsequence_out.clear();
    subsequence_out.resize(0);
  }

  double warp = warp_in;
  double scale = scale_in;
  int length = length_in;
  double height = height_in;

  if (height < 0.0) {

    height = -height;
    scale = -scale;
  }

  int realLength = (int)(warp * length);

  subsequence_out.push_back(0.0);

  if (realLength > 2) {

    int diameter = realLength - 2;

    double scaleFactor = 2.0 * scale * height / diameter;
    int itr = diameter - 1;

    for (; itr > 0; itr -= 2)
      subsequence_out.push_back(scaleFactor * sqrt((diameter / 2.0) * (diameter
              / 2.0) - (itr / 2.0) * (itr / 2.0)));

    if (itr < 0)
      itr += 2;

    for (; itr <= diameter; itr += 2)
      subsequence_out.push_back(scaleFactor * sqrt((diameter / 2.0) * (diameter
              / 2.0) - (itr / 2.0) * (itr / 2.0)));
  }

  subsequence_out.push_back(0.0);
}

void generateTrapezoidMotif(vector<double> &subsequence_out, const double
    &warp_in, const double &scale_in, const int &length_in, const double
    &height_in) {

  if (!subsequence_out.empty()) {

    subsequence_out.clear();
    subsequence_out.resize(0);
  }

  double warp = warp_in;
  double scale = scale_in;
  int length = length_in;
  double height = height_in;

  int realLength = (int)(warp * length);

  subsequence_out.push_back(0.0);

  if (realLength > 2) {
    int quarter = (realLength - 2) / 4;
    double increase = scale * height / (quarter + 1);

    for (int itr = 1; itr <= quarter; itr++)
      subsequence_out.push_back(increase * itr);

    for (int itr = 0; itr < (realLength - 2) - 2 * quarter; itr++)
      subsequence_out.push_back(scale * height);

    for (int itr = quarter; itr > 0; itr--)
      subsequence_out.push_back(increase * itr);
  }

  subsequence_out.push_back(0.0);
}

void generatePositiveFlankMotif(vector<double> &subsequence_out, const double
    &warp_in, const double &scale_in, const int &length_in, const double
    &height_in) {

  if (!subsequence_out.empty()) {

    subsequence_out.clear();
    subsequence_out.resize(0);
  }

  double warp = warp_in;
  double scale = scale_in;
  int length = length_in;
  double height = height_in;

  int realLength = (int)(warp * length);

  subsequence_out.push_back(0.0);

  if (realLength > 2) {

    double increase = scale * height / (realLength - 2.0);

    for (int itr = 1; itr < realLength - 1; itr++)
      subsequence_out.push_back(increase * itr);
  }

  subsequence_out.push_back(0.0);
}

void generateNegativeFlankMotif(vector<double> &subsequence_out, const double
    &warp_in, const double &scale_in, const int &length_in, const double
    &height_in) {

  if (!subsequence_out.empty()) {

    subsequence_out.clear();
    subsequence_out.resize(0);
  }

  double warp = warp_in;
  double scale = scale_in;
  int length = length_in;
  double height = height_in;

  int realLength = (int)(warp * length);

  subsequence_out.push_back(0.0);

  if (realLength > 2) {

    double increase = scale * height / (realLength - 2.0);

    for (int itr = realLength - 2; itr > 0; itr--)
      subsequence_out.push_back(increase * itr);
  }

  subsequence_out.push_back(0.0);
}

void generateSineMotif(vector<double> &subsequence_out, const double &warp_in,
    const double &scale_in, const int &length_in, const double &height_in) {

  if (!subsequence_out.empty()) {

    subsequence_out.clear();
    subsequence_out.resize(0);
  }

  double warp = warp_in;
  double scale = scale_in;
  int length = length_in;
  double height = height_in / 2.0;

  int realLength = (int)(warp * length);
  int middle = (realLength - 1) / 2;

  if (realLength > 2) {

    for (int itr = 0; itr < realLength - 1; itr++)
      if (realLength % 2 == 1 && itr == middle)
        subsequence_out.push_back(0.0);
      else
        subsequence_out.push_back(scale * height * sin((2 * M_PI) / (realLength
                - 1.0) * itr));
  }

  subsequence_out.push_back(0.0);
}

void generateCosineMotif(vector<double> &subsequence_out, const double
    &warp_in, const double &scale_in, const int &length_in, const double
    &height_in) {

  if (!subsequence_out.empty()) {

    subsequence_out.clear();
    subsequence_out.resize(0);
  }

  double warp = warp_in;
  double scale = scale_in;
  int length = length_in;
  double height = height_in / 2.0;

  int realLength = (int)(warp * length);

  if (realLength > 2)
    for (int itr = 0; itr < realLength; itr++)
      subsequence_out.push_back((scale * height * sin(M_PI / 2.0 + (2 * M_PI)
              / (realLength - 1.0) * itr)) - scale * height);
}
