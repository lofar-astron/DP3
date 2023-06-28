# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# MockPyStep.py: simple python step, where the MockPyStep inherits from
# the (py)dp3.DPStep class.

from dp3 import parameterset, Fields

# Import raw step from pydp3
# The extra checks in the higher level interface (dp3.Step) are thereby omitted
# The raw Step is needed because some unit tests trigger these additional checks
from dp3.pydp3 import Step

import numpy as np


class MockPyStep(Step):
    """Example python DPStep that multiplies DATA and WEIGHT_SPECTRUM"""

    def __init__(self, parset, prefix):
        """
        Set up the step (constructor). Read the parset here.

        Args:
          parset: Parameter set for the entire pipeline
          prefix: Prefix for this step, e.g. "thisstepname."
        """

        super().__init__()
        self.datafactor = parset.get_double(prefix + "datafactor")
        self.weightsfactor = parset.get_double(prefix + "weightsfactor")

    def show(self):
        """Print a summary of the step and its settings"""
        print("\nMockPyStep")
        print(f"  data factor:    {self.datafactor}")
        print(f"  weights factor: {self.weightsfactor}")

    def get_required_fields(self):
        return Fields.DATA | Fields.WEIGHTS

    def get_provided_fields(self):
        # The Flags field is added as provided fields only for testing purposes.
        return Fields.DATA | Fields.FLAGS | Fields.WEIGHTS

    def process(self, dpbuffer):
        """
        Process one time slot of data. This function MUST call self.get_next_step().process

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
        self.get_next_step().process(dpbuffer)

    def finish(self):
        """
        If there is any remaining data, process it. This can be useful if the
        step accumulates multiple time slots.
        """
        pass
