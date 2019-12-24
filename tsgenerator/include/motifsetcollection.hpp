///\file motifsetcollection.hpp
///
/// \brief File contains the motif set collection declarations.
///
/// This is the header file of the motif set collection. The motif set
///collection consists of motif set class declarations.

#ifndef MOTIFCOLLECTION_HPP
#define MOTIFCOLLECTION_HPP

#define _USE_MATH_DEFINES

#include <cmath>
#include <list>
#include <tsgtypes.hpp>


namespace tsg {

  ///\brief Generates a box motif.
  ///
  ///\param [out] &subsequence_out Hands over the subsequence of values.
  ///\param [in] &warp_in Hands over the warp factor for height scaling.
  ///\param [in] &scale_in Hands over the scale factor for length scaling.
  ///\param [in] &length_in Hands over the base length.
  ///\param [in] &height_in Hands over the base heigth.
  ///
  /// The box motif seqeuence looks similiar to:
  /// \code
  /// ^
  /// |
  /// |      .........
  /// |
  /// |
  /// |
  /// | .....         .....
  /// O---------------------->
  /// \endcode
  ///
  ///This is the box motif function. The function generates a box motif. The
  ///box motif starts and ends with 0.0.
  void generateBoxMotif(rseq &subsequence_out, const double warp_in = 1.0,
      const double scale_in = 1.0, const int length_in = 100, const double
      height_in = 50.0);

  ///\brief Generates a triangle motif.
  ///
  ///\param [out] &subsequence_out Hands over the subsequence of values.
  ///\param [in] &warp_in Hands over the warp factor for height scaling.
  ///\param [in] &scale_in Hands over the scale factor for length scaling.
  ///\param [in] &length_in Hands over the base length.
  ///\param [in] &height_in Hands over the base heigth.
  ///
  /// The TriangleMotifSet looks similiar to:
  /// \code
  /// ^
  /// |
  /// |          .
  /// |        .   .
  /// |      .       .
  /// |    .           .
  /// | ..               ..
  /// O---------------------->
  /// \endcode
  ///
  ///This is the triangle motif function. The function generates a triangle
  ///motif. The triangle motif starts and ends with 0.0.
  void generateTriangleMotif(rseq &subsequence_out, const double warp_in = 1.0,
      const double scale_in = 1.0, const int length_in = 100, const double
      height_in = 50.0);

  ///\brief Generates a semicircle motif.
  ///
  ///\param [out] &subsequence_out Hands over the subsequence of values.
  ///\param [in] &warp_in Hands over the warp factor for height scaling.
  ///\param [in] &scale_in Hands over the scale factor for length scaling.
  ///\param [in] &length_in Hands over the base length.
  ///\param [in] &height_in Hands over the base heigth.
  ///
  /// The SemicircleMotifSet looks similiar to:
  /// \code
  /// ^
  /// |
  /// |       .......
  /// |    .           .
  /// |   .             .
  /// |
  /// | ..               ..
  /// O---------------------->
  /// \endcode
  ///
  ///This is the semicircle motif function. The function generates a semicircle
  ///motif. The semicircle motif starts and ends with 0.0.
  void generateSemicircleMotif(rseq &subsequence_out, const double warp_in
      = 1.0, const double scale_in = 1.0, const int length_in = 100, const
      double height_in = 50.0);

  ///\brief Generates a trapezoid motif.
  ///
  ///\param [out] &subsequence_out Hands over the subsequence of values.
  ///\param [in] &warp_in Hands over the warp factor for height scaling.
  ///\param [in] &scale_in Hands over the scale factor for length scaling.
  ///\param [in] &length_in Hands over the base length.
  ///\param [in] &height_in Hands over the base heigth.
  ///
  /// The TrapezoidMotifSet looks similiar to:
  /// \code
  /// ^
  /// |
  /// |      .........
  /// |     .         .
  /// |    .           .
  /// |   .             .
  /// | ..               ..
  /// O---------------------->
  /// \endcode
  ///
  ///This is the trapezoid motif function. The function generates a trapezoid
  ///motif. The trapezoid motif starts and ends with 0.0.
  void generateTrapezoidMotif(rseq &subsequence_out, const double warp_in
      = 1.0, const double scale_in = 1.0, const int length_in = 100, const
      double height_in = 50.0);

  ///\brief Generates a positive flank motif.
  ///
  ///\param [out] &subsequence_out Hands over the subsequence of values.
  ///\param [in] &warp_in Hands over the warp factor for height scaling.
  ///\param [in] &scale_in Hands over the scale factor for length scaling.
  ///\param [in] &length_in Hands over the base length.
  ///\param [in] &height_in Hands over the base heigth.
  ///
  /// The PositiveFlankMotifSet looks similiar to:
  /// \code
  /// ^
  /// |
  /// |                 .
  /// |             .
  /// |         .
  /// |     .
  /// | ..               ..
  /// O---------------------->
  /// \endcode
  ///
  ///This is the positive flank motif function. The function generates
  ///a positive flank motif. The positive flank motif starts and ends with 0.0.
  void generatePositiveFlankMotif(rseq &subsequence_out, const double warp_in
      = 1.0, const double scale_in = 1.0, const int length_in = 100, const
      double height_in = 50.0);

  ///\brief Generates a negative flank motif.
  ///
  ///\param [out] &subsequence_out Hands over the subsequence of values.
  ///\param [in] &warp_in Hands over the warp factor for height scaling.
  ///\param [in] &scale_in Hands over the scale factor for length scaling.
  ///\param [in] &length_in Hands over the base length.
  ///\param [in] &height_in Hands over the base heigth.
  ///
  /// The NegetiveFlankMotifSet looks similiar to:
  /// \code
  /// ^
  /// |
  /// |   .
  /// |       .
  /// |           .
  /// |               .
  /// | ..               ..
  /// O---------------------->
  /// \endcode
  ///
  ///This is the negative flank motif function. The function generates
  ///a negative flank motif. The negative flank motif starts and ends with 0.0.
  void generateNegativeFlankMotif(rseq &subsequence_out, const double warp_in
      = 1.0, const double scale_in = 1.0, const int length_in = 100, const
      double height_in = 50.0);

  ///\brief Generates a sine motif.
  ///
  ///\param [out] &subsequence_out Hands over the subsequence of values.
  ///\param [in] &warp_in Hands over the warp factor for height scaling.
  ///\param [in] &scale_in Hands over the scale factor for length scaling.
  ///\param [in] &length_in Hands over the base length.
  ///\param [in] &height_in Hands over the base heigth.
  ///
  /// The SineMotifSet looks similiar to:
  /// \code
  /// ^
  /// |     . .
  /// |   .     .
  /// |  .
  /// O-.--------.--------.-->
  /// |                  .
  /// |           .     .
  /// |             . .
  /// \endcode
  ///
  ///This is the sine motif function. The function generates a sine motif. The
  ///sine motif starts and ends with 0.0.
  void generateSineMotif(rseq &subsequence_out, const double warp_in = 1.0,
      const double scale_in = 1.0, const int length_in = 100, const double
      height_in = 50.0);

  ///\brief Generates a cosine motif.
  ///
  ///\param [out] &subsequence_out Hands over the subsequence of values.
  ///\param [in] &warp_in Hands over the warp factor for height scaling.
  ///\param [in] &scale_in Hands over the scale factor for length scaling.
  ///\param [in] &length_in Hands over the base length.
  ///\param [in] &height_in Hands over the base heigth.
  ///
  ///The CosineMotifSet looks similiar to:
  /// \code
  /// ^
  /// O--.---------------.--->
  /// |    .           .
  /// |
  /// |     .         .
  /// |
  /// |       .     .
  /// |         . .
  /// \endcode
  ///
  ///This is the cosine motif function. The function generates a cosine motif.
  ///The cosine motif starts and ends with 0.0.
  void generateCosineMotif(rseq &subsequence_out, const double warp_in = 1.0,
      const double scale_in = 1.0, const int length_in = 100, const double
      height_in = 50.0);
}

#endif
