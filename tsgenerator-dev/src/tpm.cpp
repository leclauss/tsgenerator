///\file tpm.cpp
///
///\brief File contains a top pair motif discovery algorithm.
///
///This algorithm computes the top pair motif of a time series by evaluating
///iteratively the diagonals of the distance matrix of the time series.

#include <tpm.hpp>


namespace tsg {

  double tpm(const rseq &timeSeries_in, const rseq &sums_in, const rseq
      &sumSquares_in, int &pos0_out, int &pos1_out, const int window_in) {

    int window = window_in;
    double rWindow = 1.0 / window;

    int length = (int)(timeSeries_in.size());

    //precompute running mean and standard deviation
    rseq mean;
    for (int i = 0; i < length - window + 1; i++)
      mean.push_back(sums_in[i] * rWindow);
    rseq var;
    for (int i = 0; i < length - window + 1; i++)
      var.push_back(sumSquares_in[i] * rWindow - mean[i] * mean[i]);
    rseq sigma;
    for (int i = 0; i < length - window + 1; i++)
      sigma.push_back(var[i] > 1.0 ? sqrt(var[i]) : 1.0);

    double q = 0.0;
    double distance;
    double bsf = std::numeric_limits<double>::infinity();

    //for each diagonal of non overlapping subsequences
    for (int k = window - 1; k < length - window + 1; k++) {

      //the dot product
      q = 0.0;

      //iterate throw the diagonal
      for (int i = 0; i < length - window + 1 - k; i++) {

        if (i == 0)
          //compute the initial dot product
          for (int j = 0; j < window; j++)
            q += timeSeries_in[i + j] * timeSeries_in[k + j];
        else
          //compute dot product iteratively
          q += timeSeries_in[i + window - 1] * timeSeries_in[i + k + window
            - 1] - timeSeries_in[i - 1] * timeSeries_in[i + k - 1];

        //compute distance with dot product
        distance = 2 * (window - (q - window * mean[i] * mean[i + k])
            / (sigma[i] * sigma[i + k]));

        //save the best so far positions and distance
        if (distance < bsf) {

          bsf = distance;
          pos0_out = i;
          pos1_out = i + k;
        }
      }
    }

    return bsf;
  }
}
