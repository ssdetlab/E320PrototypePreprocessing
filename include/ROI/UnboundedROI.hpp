#pragma once

#include "ROI/ROI.hpp"

class UnboundedROI : public ROI {
 public:
  UnboundedROI() = default;

  bool isInside(double x, double y) const override { return true; }
};
