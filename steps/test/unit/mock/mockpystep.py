# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# MockPyStep.py: simple python step, where the MockPyStep inherits from
# the (py)dp3.DPStep class.

import parameterset
from pydp3 import Step
import numpy as np

class MockPyStep(Step):
    """Example python DPStep that multiplies DATA and WEIGHT_SPECTRUM"""
    def __init__(self, parset, prefix):
        """
        Set up the step (constructor). Read the parset here.
        Set fetch_weights to True if the weights need to be read.
        Similarly for fetch_uvw.

        Args:
          parset: Parameter set for the entire pipeline
          prefix: Prefix for this step, e.g. "thisstepname."
        """
        super().__init__()
        self.datafactor = parset.getDouble(prefix + "datafactor")
        self.weightsfactor = parset.getDouble(prefix + "weightsfactor")

        self.fetch_weights = True
        self.fetch_uvw = False

    def update_info(self, dpinfo):
        """
        Process metadata. This will be called before any call to process.

        Args:
          dpinfo: DPInfo object with all metadata, see docs in pydp3.cc
        """
        super().update_info(dpinfo)

        # Make sure data is read
        self.info().set_need_vis_data()

        # Make sure data will be written
        self.info().set_write_data()

    def show(self):
        """Print a summary of the step and its settings"""
        print("\nMockPyStep")
        print(f"  data factor:    {self.datafactor}")
        print(f"  weights factor: {self.weightsfactor}")

    def process(self, dpbuffer):
        """
        Process one time slot of data. This function MUST call process_next_step.

        Args:
          dpbuffer: DPBuffer object which can contain data, flags and weights
                    for one time slot.
        """

        data = np.array(dpbuffer.get_data(), copy=False)
        weights = np.array(dpbuffer.get_weights(), copy=False)

        # Do the operation on data
        data *= self.datafactor
        weights *= self.weightsfactor

        # Send processed data to the next step
        self.process_next_step(dpbuffer)

    def finish(self):
        """
        If there is any remaining data, process it. This can be useful if the
        step accumulates multiple time slots.
        """
        pass
