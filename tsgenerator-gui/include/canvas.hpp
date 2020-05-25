#ifndef CANVAS_HPP
#define CANVAS_HPP

#include <tsgtypes.hpp>
#include <QLabel>
#include <QPainter>
#include <QPen>
#include <QImage>
#include <QPixmap>
#include <QPoint>
#include <QColor>
#include <QMouseEvent>

class Canvas : public QLabel
{

private:
  tsg::rseq shape;
  QPoint lastPoint;
  int penWidth = 1;
  int drawing = false;
  int res = 300;

  QImage img;
  QPixmap pix;
  QPainter pen;

public:
  Canvas();
  ~Canvas();

  ///\brief This function returns the drawn shape.
  ///
  ///\param [in] &shape_out Hands over the shape values.
  ///
  ///This function returns the drawn shape as a list.
  void getShape(tsg::rseq &shape_out);

protected:
  ///\brief This function handles the mouse press action from custom label.
  ///
  ///\param [in] *event_in Hands over the event.
  ///
  ///This function saves the first point when the mouse is clicked on the
  ///custom label.
  void mousePressEvent(QMouseEvent *event_in) override;

  ///\brief This function handles the mouse move action from custom label.
  ///
  ///\param [in] *event_in Hands over the event.
  ///
  ///This function draws the a line segment when the mouse is moved on the
  ///custom label.
  void mouseMoveEvent(QMouseEvent *event_in) override;

  ///\brief This function handles the mouse released action from custom label.
  ///
  ///\param [in] *event_in Hands over the event.
  ///
  ///This function draws the last line when the mouse is released from the
  ///custom label.
  void mouseReleaseEvent(QMouseEvent *event_in) override;

private:
  ///\brief This function draws a line on the canvas.
  ///
  ///\param [in] endPoint_in Hands over the end point of the line.
  ///
  ///This function draws a line on the canvas depending on the mouse positions.
  void drawLineTo(const QPoint &endPoint_in);

  ///\brief This function converts the image into a list.
  ///
  ///This function converts the drawn image into a list of values.
  void imageToPoints();
};

#endif // CANVAS_HPP
