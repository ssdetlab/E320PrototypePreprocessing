#include <cmath>
#include <iostream>

#include "DetectorEvent.h"

namespace EPICS {

/// @brief convert EPICS dipole PV to T
///
/// @param val PV value -- beam energy in eV
/// @param theta deflection angle in rad
/// @param bz dipole length in m
///
/// @return field strength in T
double convertEpicsDipole(double val, double theta, double bz) {
  const double c_light = 2.99792458E8;  // speed of light in m/s
  const double me = 0.5109989461e6;     // electron mass in eV
  return std::sqrt(val * val - me * me) * std::sin(theta) / (c_light * bz);
}

/// @brief convert EPICS x-corrector PV to T
///
/// @param val PV value -- kG * m
/// @param bz dipole length in m
///
/// @return field strength in T
double convertEpicsXCorrector(double val, double bz) {
  return val / bz * 0.1;
}

/// @brief convert EPICS quadrupole PV to T
///
/// @param val PV value -- kG
///
/// @return field gradient in T/m
double convertEpicsQuad(double val) {
  return val * 0.1;
}

/// @brief extract EPICS PV value as a double
///
/// @param detEvent detector event pointer
/// @param pvname name of the PV to extract
///
/// @return value of the PV

double getPv(const detector_event *detEvent, const std::string &pvname) {
  try {
    return std::stod(detEvent->ev_epics_frame.pv_map.at(pvname).value);
  } catch (std::out_of_range &err) {
    std::cout << "MISSING PV " << pvname << "\n";
    throw err;
  }
}

}  // namespace EPICS
