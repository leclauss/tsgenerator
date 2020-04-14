///\file tsgtypes.hpp
///
///\brief File contains the TSGenerator type definitions.
///
///This is file contains the type defintions of the tsg library.

#ifndef TSGTYPES_HPP
#define TSGTYPES_HPP

#include <vector>
#include <set>
#include <string>


namespace tsg {

  ///\brief This is a vector of reals.
  ///
  ///The vector of reals is used as a mathimatical vector of reals or
  ///a sequence of reals.
  typedef std::vector<double> rseq;

  ///\brief This is a vector of vectors of reals.
  ///
  ///The vector of vectors of reals is used as a mathimatical matrix of reals
  ///or a sequence of sequences of reals.
  typedef std::vector<std::vector<double>> rseqs;

  ///\brief This is a vector of indices.
  ///
  ///The vector of indices is used as a mathimatical vector of integers or as
  ///a sequence of integers.
  typedef std::vector<int> iseq;

  ///\brief This is a vector of vectors of integers.
  ///
  ///The vector of vectors of integers is used as a mathimatical matrix of
  ///integers or a sequence of sequences of integers.
  typedef std::vector<std::vector<int>> iseqs;

  ///\brief This is a set of unique integers.
  ///
  ///The set of unique integers contains integers without duplicates.
  typedef std::set<int> iset;

  ///\brief This is a set of unique reals.
  ///
  ///The set of unique reals contains reals without duplicates.
  typedef std::set<double> rset;

  ///\brief This is a word.
  ///
  ///The word is just a regular std string.
  typedef std::string word;

  ///\brief This is a paragraph.
  ///
  ///The paragraph is a sequence of words.
  typedef std::vector<std::string> par;

  ///\brief This struct represents a time series subsequence.
  ///
  ///The time series subsequence consists of a postion and length value.
  struct subsequence {

    int position = 0;
    int length = 0;
  };

  ///\brief This is a time series.
  ///
  ///The time series is a seuqence of subsequences.
  typedef std::vector<subsequence> subsequences;

  ///\brief This struct an interval of the time series.
  ///
  ///The time series interval consists of a postion and length.
  struct interval {

    int start = 0;
    int end = 0;
  };

  ///\brief This is a interval sequence.
  ///
  ///The interval sequence is a seuqence of intervals.
  typedef std::vector<interval> intervals;

  ///\brief This is the default time series length default.
  ///
  ///The default time series length sets the default number of values of the
  ///synthetic times series.
  const int defaultLength = 4000;

  ///\brief This is the default window size value.
  ///
  ///The default windows size sets the default number of values of
  ///a subsequnece of the synthetic times series.
  const int defaultWindow = 30;

  ///\brief This is the default delta value.
  ///
  ///The default delta sets the default maximal absolute of difference of two
  ///values in the base time series.
  const double defaultDelta = 1.0;

  ///\brief This is the default noise value.
  ///
  ///The default noise sets the default maximal absolute of difference to
  ///a value of the synthetic time series or a subsequence.
  const double defaultNoise = 2.0;

  ///\brief This is the default type.
  ///
  ///The default type sets the default shape of the top motif injected into the
  ///base time series.
  const word defaultType = "box";

  ///\brief This is the default motif size value.
  ///
  ///The default motif size sets the default number of subsequnces non-self
  ///matched by the motif. Ignored for pair motifs.
  const int defaultMotifSize = 3;

  ///\brief This is the default motif height value.
  ///
  ///The default motif height sets the default absolute base difference in
  ///a motif sequence.
  const double defaultHeight = 10.0;

  ///\brief This is the default step value.
  ///
  ///The default step sets the default maximal absolut difference in
  ///x direction between two values when generating a base time series some
  ///kind of approximations between the values.
  const double defaultStep = 1.0;

  ///\brief This is the default times value.
  ///
  ///The default times sets the default number of values generated for
  ///a repeating pattern when generating the base time series.
  const int defaultTimes = 3;

  ///\brief This is the default maxi value.
  ///
  ///The default maxi sets the default maximal absolut difference between any
  ///pair of values in the synthetic times series.
  const double defaultMaxi = 20.0;

  ///\brief This is the default smaller value.
  ///
  ///The default smaller sets the default number of smaller motifs injected
  ///into the time series after injecting the top pair motif.
  const int defaultSmaller = 1;

  ///\brief This is the default method value.
  ///
  ///The default method sets the default base time series generation procedure.
  const word defaultMethod = "boundedNormalRandomWalk";

  ///\brief This is the default generator value.
  ///
  ///The default generator sets the default injection procedure like pair
  ///motif, set motif or latent motif injection.
  const word defaultGen = "latent motif";
}

#endif
