#pragma once

class ROI {
 public:
  virtual bool isInside(double x, double y) const = 0;
};
