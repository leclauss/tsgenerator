///\file tsgutils.cpp
///
///\brief File contains additional utilities.
///
///This is a set of additional utilities for time series modifications and
///analysis.

#include <tsgutils.hpp>


namespace tsg {

  void minMaxDist(const rseq &timeSeries_in, const rseq &sums_in, const rseq
      &sumSquares_in, const int window_in, double &min_out, double &max_out) {

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
    min_out = std::numeric_limits<double>::infinity();
    max_out = 0.0;

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

        //save the best so far distances
        if (distance < min_out)
          min_out = distance;

        if (distance > max_out)
          max_out = distance;
      }
    }

    min_out = sqrt(min_out);
    max_out = sqrt(max_out);
  }

  void calcRunnings(const rseq &sequence_in, const int window_in, rseq
      &sums_out, rseq &sumSquares_out) {

    // obtain running sum and sum of square
    if (!sums_out.empty() || !sumSquares_out.empty()) {

      sums_out.clear();
      sumSquares_out.clear();
    }

    sums_out.resize(sequence_in.size() - window_in + 1);
    sumSquares_out.resize(sequence_in.size() - window_in + 1);

    sums_out[0] = sequence_in[0];
    sumSquares_out[0] = sequence_in[0] * sequence_in[0];

    for (int i = 1; i < window_in; i++) {

      sums_out[0] += sequence_in[i];
      sumSquares_out[0] += sequence_in[i] * sequence_in[i];
    }

    for (int i = 0; i < (int)sequence_in.size() - window_in; i++) {

      sums_out[i + 1] = sums_out[i] - sequence_in[i] +
        sequence_in[i + window_in];

      sumSquares_out[i + 1] = sumSquares_out[i] - sequence_in[i] *
        sequence_in[i] + sequence_in[i + window_in] *
        sequence_in[i + window_in];
    }
  }

  void generateBaseTimeSeries(rseq &timeSeries_out, const word method_in, const
      int length_in, const double delta_in, const double step_in, const int
      times_in, const double noise_in) {

    if (!timeSeries_out.empty()) {

      timeSeries_out.clear();
      timeSeries_out.resize(0);
    }

    baseTS BaseTS;

    // get the method
    int method = (int)(std::distance(TSGenerator.methods.begin(),
          std::find(TSGenerator.methods.begin(),
          TSGenerator.methods.end(), method_in)));

    // check if method exists
    if (method >= (int) baseTS.methods.size()) {

      std::cerr << "ERROR: Unknown method: " << method_in << std::endl;
      throw(EXIT_FAILURE);
    }

    switch (method) {

      case 0:
        baseTS.simpleRandomWalk(timeSeries_out, length_in, delta_in, noise_in);
        break;
      case 1:
        baseTS.realRandomWalk(timeSeries_out, length_in, delta_in, noise_in);
        break;
      case 2:
        baseTS.normalRandomWalk(timeSeries_out, length_in, delta_in, noise_in);
        break;
      case 3:
        baseTS.linearRandomWalk(timeSeries_out, length_in, delta_in, step_in,
            noise_in);
        break;
      case 4:
        baseTS.simpleRandomWalk(timeSeries_out, length_in, delta_in, maxi_in,
            noise_in);
        break;
      case 5:
        baseTS.realRandomWalk(timeSeries_out, length_in, delta_in, maxi_in,
            noise_in);
        break;
      case 6:
        baseTS.normalRandomWalk(timeSeries_out, length_in, delta_in, maxi_in,
            noise_in);
        break;
      case 7:
        baseTS.linearRandomWalk(timeSeries_out, length_in, delta_in, step_in, maxi_in,
            noise_in);
        break;
      case 8:
        baseTS.uniformRandom(timeSeries_out, length_in, delta_in, noise_in);
        break;
      case 9:
        baseTS.normalRandom(timeSeries_out, length_in, delta_in, noise_in);
        break;
      case 10:
        baseTS.piecewiseLinearRandom(timeSeries_out, length_in, delta_in,
            noise_in);
        break;
      case 11:
        baseTS.splineRepeated(timeSeries_out, length_in, delta_in, (int)step_in, times_in,
            noise_in);
        break;
        default:
        std::cerr << "ERROR: Unknown method: " << method << std::endl;
        throw(EXIT_FAILURE);
    }
  }
}
