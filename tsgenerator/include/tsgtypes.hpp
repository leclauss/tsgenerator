///\file tsgtypes.hpp
///
///\brief File contains the TSGenerator type definitions.
///
///This is file contains the type defintions of the tsg library.

#ifndef TSGTYPES_HPP
#define TSGTYPES_HPP

#include <vector>
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
}

#endif
