#pragma once

#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>

#include <TObject.h>

namespace E320Io {

/// @brief ALPIDE readout frame
class ChipEvent {
 public:
  /// @brief Pixel cluster struct
  struct Cluster {
    /// Cluster geometric center in
    /// local chip coordinates
    double xCenter;
    double yCenter;

    /// Cluster length in local chip
    /// coordinates
    std::size_t dx = 0;
    std::size_t dy = 0;

    /// Number of pixels in the cluster
    std::size_t clSize = 0;

    /// Cluster pixels list
    std::vector<std::pair<std::uint16_t, std::uint16_t>> pixels;
  };

  /// Hit map of the chip in the event
  std::vector<Cluster> hits;

  /// ALPIDE readout flags
  bool isBusyViolation = false;
  bool isFlushedIncomplete = false;
  bool isStrobeExtended = false;
  bool isBusyTransition = false;

  /// MOSAIC readout flags
  bool endOfRun = false;
  bool overflow = false;
  bool timeout = false;
  bool headerError = false;
  bool decoder10b8bError = false;
  bool eventOversizeError = false;

  /// Chip identificators
  std::uint8_t chipId = 0;
  std::uint8_t channel = 0;

  ClassDef(ChipEvent, 5);
};

class StaveEvent {
 public:
  /// Chip event storage
  std::vector<ChipEvent> chipEventBuffer;

  /// Stave identificator
  std::uint8_t staveId = 0;

  ClassDef(StaveEvent, 4);
};

class DetectorEvent {
 public:
  /// Stave event storage
  std::vector<StaveEvent> staveEventBuffer;

  /// Event MOSAIC trigger ID
  std::uint64_t trgN = 0;

  /// UNIX timestamp in ms
  std::uint64_t tsBegin = 0;
  std::uint64_t tsEnd = 0;

  /// EPICS parity
  std::size_t epicsParity = 0;
  /// EPICS PID
  std::size_t epicsPID = 0;
  /// EPICS DAQ number
  std::size_t epicsDAQNumber = 0;

  /// Run number
  std::uint32_t runNumber = 0;

  ClassDef(DetectorEvent, 8);
};

class ChipRunMeta {
 public:
  /// Chip identificator
  std::uint8_t chipId = 0;

  /// Number of received triggers
  std::uint16_t nTrigs = 0;

  /// Number of generated strobes
  std::uint16_t nStrobes = 0;

  /// Number of matrix readouts
  std::uint16_t nMatrixReadouts = 0;

  /// Number of sent frames
  std::uint16_t nFrames = 0;

  ClassDef(ChipRunMeta, 3);
};

class StaveRunMeta {
 public:
  /// Stave identificator
  std::uint8_t staveId;

  /// Frames with data
  std::uint32_t nonEmptyFrames;

  /// Number of triggers received
  /// by MOSAIC
  std::uint32_t nTrigsMosaic;

  /// ALPIDE meta data
  std::vector<ChipRunMeta> chipsRmd;

  ClassDef(StaveRunMeta, 3);
};

class RunMeta {
 public:
  using TimeServerMd = std::tuple<std::uint64_t, std::uint64_t, std::uint64_t>;

  /// Current run number
  std::uint32_t runNumber;

  /// UNIX timestamp in ms
  std::uint64_t runStart;
  std::uint64_t runEnd;

  /// Stave meta data
  std::vector<StaveRunMeta> staveRmd;

  /// Time server meta data
  std::vector<TimeServerMd> tsMd;

  ClassDef(RunMeta, 3);
};

}  // namespace E320Io
