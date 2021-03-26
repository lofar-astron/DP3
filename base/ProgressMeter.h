// ProgressMeter.h: Visual indication of a tasks progress.
// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LOFAR_COMMON_PROGRESSMETER_H
#define LOFAR_COMMON_PROGRESSMETER_H

#include <string>

namespace dp3 {
namespace base {

/// @brief Visual indication of a tasks progress.

/// This class shows the progress of a task.
/// It is copied from casacore.
///
/// It shows the progress on stdout using a line with dots and percentages.
///
/// It is possible to attach the progressmeter to function in, say, a GUI
/// that can show the progress in a more visual way.
///
/// The progress meter will usually be removed from the screen once the maximum
/// value is set, so you should not reuse the ProgressMeter after that has
/// happened. It is harmless, but it will not result in any visual feedback for
/// the user.
///
/// While the "min" is usually less than "max", if in fact it is greater than
/// "max" the progress meter will count down correctly.
///
/// For example:
/// \code
/// void calculate(unsigned int n) {
///   int skip = n / 200;
///   ProgressMeter meter(0, n, "Title", "Subtitle", "", "", true, skip);
///   for (unsigned int i=0; i<n; i++) {
///       ... calculate ...
///       meter.update(i);
///   }
/// }
/// \endcode

class ProgressMeter {
 public:
  /// Makes a null progress meter, i.e. no updates to the screen are
  /// generated.
  ProgressMeter();

  /// Create a progress meter with the given min and max values and labels.
  /// if <tt>estimateTime=true</tt>, an estimate of the
  /// time remaining will be made for the user. This estimate assumes that
  /// the remaining portion will compute at the same rate as the portion
  /// completed so far, so the time should not be estimated for processes
  /// which are non-linear.
  ///
  /// Any labels which are set to the empty string will have sensible defaults
  /// supplied. In particular, <tt>minlabel</tt> and <tt>maxlabel</tt>
  /// will be set to the display the minimum and maximum values.
  ///
  /// Normally the progress bar will be updated with every call to
  /// <tt>update()</tt>. If however you will be sending many events
  /// then you might want to update the GUI every <tt>updateEvery</tt>'th
  /// event for efficiency. Generally there's no point updating more than
  /// a couple of hundred times since the eye can't distinguish differences
  /// in the progress bar position at that level. If updateEvery is <=0, it
  /// is set to 1 for you.
  ProgressMeter(double min, double max, const std::string& title,
                const std::string& subtitle, const std::string& minlabel,
                const std::string& maxlabel, bool estimateTime = true,
                int updateEvery = 1);

  /// The destruction of the meter will cause an update to be sent with the
  /// maximum value. This will usually cause the GUI window to be removed
  /// from the screen. Thus the progress meter should generally live as long
  /// as the calculation it is tracking.
  ~ProgressMeter();

  void update(double value, bool force = false);

  /// Get the min and max values of the progress meter.
  ///@{
  double min() const { return min_p; }
  double max() const { return max_p; }
  ///@}

  friend class ObjectController;

 private:
  int id_p;
  double min_p, max_p;
  int update_every_p, update_count_p;

  /// These are set by ObjectController for executables that have the tasking
  /// system in them, otherwise they are null and this class just does no-ops.
  static int (*creation_function_p)(double, double, const std::string&,
                                    const std::string&, const std::string&,
                                    const std::string&, bool);
  static void (*update_function_p)(int, double);

  /// Undefined and inaccessible
  ProgressMeter(const ProgressMeter&);
  ProgressMeter& operator=(const ProgressMeter&);
};

}  // namespace base
}  // namespace dp3

#endif
