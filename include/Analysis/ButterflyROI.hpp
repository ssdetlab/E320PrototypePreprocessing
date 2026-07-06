#pragma once

#include <cmath>

#include "Analysis/ROI.hpp"

class ButterflyROI : public ROI {
 public:
  ButterflyROI(double x0, double y0, double a, double b, double theta,
               double curv)
      : m_x0(x0), m_y0(y0), m_a(a), m_b(b), m_theta(theta), m_curv(curv) {}

  bool isInside(double x, double y) const override {
    // Translate point to ellipse center
    double dx = x - m_x0;
    double dy = y - m_y0;

    // Precompute trig
    double cosT = std::cos(m_theta);
    double sinT = std::sin(m_theta);

    // Rotate into ellipse frame (inverse rotation)
    double xPrime = cosT * dx + sinT * dy;
    double yPrime = -sinT * dx + cosT * dy;

    double boundY = m_b * std::sqrt(1 + (xPrime / m_curv) * (xPrime / m_curv));

    return (std::abs(yPrime) < boundY);
  }

 private:
  double m_x0;     // center x
  double m_y0;     // center y
  double m_a;      // semi-major axis
  double m_b;      // semi-minor axis
  double m_theta;  // rotation angle (radians)
  double m_curv;   // curvature
};
