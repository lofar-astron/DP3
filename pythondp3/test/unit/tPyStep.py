# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
Script can be invoked in two ways:
- as standalone from the build/pythondp3/test/integration directory,
  using `pytest source/tPyStep.py` (extended with pytest options of your choice)
- using ctest, see DP3/pythondp3/test/integration/CMakeLists.txt
"""
# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf

sys.path.insert(0, tcf.PYTHONDIR)

try:
    "The import may fail while running pytest --collect-only"
    import dp3
except ImportError:
    pass


def test_next_step_basic():
    step1 = dp3.Step()
    step2 = dp3.Step()
    step3 = dp3.Step()
    assert step1.get_next_step() is None

    step1.set_next_step(step3)
    assert step1.get_next_step() == step3

    step1.set_next_step(step2)
    assert step1.get_next_step() == step2

    step2.set_next_step(step3)
    assert step1.get_next_step() == step2
    assert step2.get_next_step() == step3
    assert step3.get_next_step() is None

    # DP3 does not detect circular graphs.
    circle_step = dp3.Step()
    circle_step.set_next_step(circle_step)
    assert circle_step.get_next_step() == circle_step


def test_next_step_process():
    class FirstStep(dp3.Step):
        def __init__(self):
            dp3.Step.__init__(self)
            self.process_count = 0

        def process(self, dpbuffer):
            self.process_count += 1
            self.get_next_step().process(dpbuffer)

    class LastStep(dp3.Step):
        def __init__(self):
            dp3.Step.__init__(self)
            self.process_count = 0

        def process(self, dpbuffer):
            self.process_count += 1

    first_step = FirstStep()
    last_step = LastStep()
    first_step.set_next_step(last_step)
    assert first_step.process_count == 0
    assert last_step.process_count == 0

    for i in range(1, 5):
        first_step.process(dp3.DPBuffer())
        assert first_step.process_count == i
        assert last_step.process_count == i


def test_info():
    class TestInfoStep(dp3.Step):
        def __init__(self):
            dp3.Step.__init__(self)

        def _update_info(self, info):
            # Modify the info a bit showing that this method is
            # indeed called when set_info() is called.
            # An example of a step that makes actual changes to the info object
            # is the Averager step.
            first_time = info.first_time + 1.0
            last_time = info.last_time + 2.0
            time_interval = info.time_interval + 3.0
            info.set_times(first_time, last_time, time_interval)
            super()._update_info(info)

    step = TestInfoStep()

    info = dp3.DPInfo()

    first_time = 42.0
    last_time = 141.0
    time_interval = 5.0
    info.set_times(first_time, last_time, time_interval)

    # Set info, this will call _update_info()
    step.set_info(info)

    # Assert that info was set and modified by _update_info()
    # Storing step.info in local variable step_info, because
    # step.info creates a copy
    step_info = step.info
    assert step_info.first_time == first_time + 1.0
    assert step_info.last_time == last_time + 2.0
    assert step_info.time_interval == time_interval + 3.0
