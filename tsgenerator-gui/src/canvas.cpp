///\file canvas.cpp
///
///\brief Source file of a canvas.
///
///This is the source file of a canvas with drawing functionality.

#include <canvas.hpp>

Canvas::Canvas() : QLabel() {

  pix = QPixmap(res, res);
  pix.fill(qRgb(255, 255, 255));

  this->setPixmap(pix);
}

Canvas::~Canvas() { }

void Canvas::getShape(tsg::rseq &shape_out) {

  this->imageToPoints();

  shape_out = shape;
}

void Canvas::mousePressEvent(QMouseEvent *event_in) {

  if (event_in->button() == Qt::LeftButton) {

    lastPoint = event_in->pos();
    drawing = true;

    pix.fill(qRgb(255, 255, 255));

    this->setPixmap(pix);
  }
}

void Canvas::mouseMoveEvent(QMouseEvent *event_in) {

  if ((event_in->buttons() & Qt::LeftButton) && drawing) {

    drawLineTo(event_in->pos());
  }
}

void Canvas::mouseReleaseEvent(QMouseEvent *event_in) {

  if (event_in->button() == Qt::LeftButton && drawing) {

    drawLineTo(event_in->pos());
    drawing = false;
  }

}

void Canvas::drawLineTo(const QPoint &endPoint_in) {

  pen.begin(&pix);
  pen.setPen(QPen(Qt::black, penWidth, Qt::SolidLine, Qt::RoundCap,
      Qt::RoundJoin));
  pen.drawLine(lastPoint, endPoint_in);
  pen.end();

  this->setPixmap(pix);

  lastPoint = endPoint_in;
}

void Canvas::imageToPoints() {

  img = pix.toImage();

  if (!shape.empty()) {

    shape.clear();
    shape.resize(0);
  }

  double mean;
  int start = -1;
  int count;
  bool end;

  for (int x = 0; x < res; x++) {

    count = 0;
    end = true;
    mean = 0.0;

    for ( int y = 0; y < res; y++) {

      if (img.pixel(x, y) == qRgb(0, 0, 0)) {

        if (start < 0)
          start = x;

        count++;
        end = false;
        mean += (double) res - y;
      }
    }

    if (start < 0)
      continue;

    if (start >= 0 && end)
      break;

    mean /= (double) count;

    //push back y coordinate
    shape.push_back(mean);
  }

  //level the y values such that the first y is 0.0
  for (size_t i = 1; i < shape.size(); i++)
    shape[i] -= shape[0];

  shape[0] = 0.0;
}
