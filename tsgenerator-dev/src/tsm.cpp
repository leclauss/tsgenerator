///\file tsm.cpp
///
///\brief File contains a top set motif discovery algorithm.
///
///This algorithm computes the top set motif of a time series by evaluating
///iteratively the diagonals of the distance matrix of the time series.

#include <tsm.hpp>


namespace tsg {

  void zNormPAA(const rseq &timeSeries_in, const rseq &sums_in, const rseq
      &sumSquares_in, rseq &paa_out, const int pos_in, const int window_in,
      const int paa_in) {

    if (!paa_out.empty()) {

      paa_out.clear();
      paa_out.resize(0);
    }

    double rW = 1.0 / (double)window_in;
    double c = (double)timeSeries_in.size() / (double)paa_in;
    double rC = (double)paa_in / (double)timeSeries_in.size();
    double s = 0.0;
    double den = 0.0;

    for (int i = 0; i < paa_in; i++) {

      s = 0;
      for (int j = c * i; j < c * i + c; j++) {

        den = (rW * (sumSquares_in[pos_in + i] - rW * (sums_in[pos_in + i]
            * sums_in[pos_in + i])));

        if (den < 1.0)
          s += (timeSeries_in[pos_in + j] - rW * sums_in[pos_in + i]) / 1.0;
        else
          s += (timeSeries_in[pos_in + j] - rW * sums_in[pos_in + i]) / den;
      }

      paa_out.push_back(rC * s);
    }
  }

  double invNormalCDF(const double p_in, const double prec_in) {

    if (p_in <= 0.0)
      return -std::numeric_limits<double>::infinity();

    if (p_in >= 1.0)
      return std::numeric_limits<double>::infinity();

    double p = 0.0;
    double err = 1.0;
    double conq = 0.25;
    double inv = 0.5;
    double sqrtH = sqrt(0.5);
    double last = -1.0;

    while (fabs(err) > prec_in) {

      p = 0.5 * erfc(-inv * sqrtH);
      err = p_in - p;

      if (p_in == p)
        return inv;
      else if (p_in < p) {

        inv -= conq;

        if (last > 0.0) {

          conq *= 0.5;
          last *= -1.0;
        }
      }
      else {

        inv += conq;

        if (last < 0.0) {

          conq *= 0.5;
          last *= -1.0;
        }
      }
    }

    return inv;
  }

  void invNormalCDF(const int alphabet_in, rseq &p_out, const double prec_in) {

    if (!p_out.empty()) {

      p_out.clear();
      p_out.resize(0);
    }

    double n = alphabet_in;
    double rN = 1.0 / n;

    for (int i = 0; i < (int)n + 1; i++)
      p_out.push_back(invNormalCDF(rN * i, prec_in));
  }

  void zNormSAX(const rseq &timeSeries_in, const rseq &sums_in, const rseq
      &sumSquares_in, word &sax_out, const int pos_in, const int window_in,
      const int paa_in, const rseq &breakpoints_in) {

    if (!sax_out.empty()) {

      sax_out.clear();
      sax_out.resize(0);
    }

    rseq paa;

    zNormPAA(timeSeries_in, sums_in, sumSquares_in, paa, pos_in, window_in,
        paa_in);

    for (int i = 0; i < (int)paa.size(); i++)
      for (int j = 0; j < (int)breakpoints_in.size(); j++)
        if (paa[j] < breakpoints_in[j]) {

          sax_out.push_back(j);
          j = (int)breakpoints_in.size();
        }
  }

  void tsm(const rseq &timeSeries_in, const rseq &sums_in, const rseq
      &sumSquares_in, iseqs &motif_out, const int window_in, const int paa_in,
      const rseq &breakpoints_in) {

    if (!motif_out.empty()) {

      motif_out.clear();
      motif_out.resize(0);
    }

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

    std::map<int, iseq> buckets;
    std::map<int, iseq>::iterator itB;
    word sax = "";
    int buf[length - window];

    for (int i = 0; i < length - window + 1; i++) {

      zNormSAX(timeSeries_in, sums_in, sumSquares_in, sax, i, window, paa_in,
          breakpoints_in);

      buf[i] = 0;

      for (int j = (int)sax.size(); j > 0; j--)
        buf[i] += sax[j] * j * (int)breakpoints_in.size();

      itB = buckets.find(buf[i]);

      if (itB == buckets.end())
        buckets.insert(std::pair<int, iseq>(buf[i], { i }));
      else
        itB->second.push_back(i);
    }




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
        }
      }
    }
  }
}
