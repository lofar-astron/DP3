# MockPyStep.py: simple python step, where the MockPyStep inherits from
# the (py)dppp.DPStep class.
# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import parameterset
from pydppp import DPStep
import numpy as np


class MockPyStep(DPStep):
    def __init__(self,parset, prefix):
        super().__init__()
        self.parset = parset
        self.dpbuffers = []
        self.prefix = prefix

    def show(self) :
        print()
        print("MockPyStep")

    def process(self, dpbuffer) :
        # Accumulate buffers
        self.dpbuffers.append(dpbuffer)

        # If we have accumulated enough data, process it
        if len(self.dpbuffers) == 10:
            self.process_buffers()

            # Send processed data to the next step
            for dpbuffer in self.dpbuffers:
                self.process_next_step(dpbuffer)

            # Clear accumulated data
            self.dpbuffers = []

    def finish(self):
        # If there is any remaining data, process it
        if len(self.dpbuffers):
            self.process_buffers()
            for dpbuffer in self.dpbuffers:
                self.process_next_step(dpbuffer)
            self.dpbuffers = []

    def process_buffers(self):
        print()
        print(f"Processing a chunk of length {len(self.dpbuffers)}.")

        # Enumerate accumulated data and display just the shapes
        for i, dpbuffer in enumerate(self.dpbuffers):
            # Multiply (complex-valued) visibilities by a factor 2
            data = np.array(dpbuffer.get_data(), copy=False)
            data *= 2.
            # Divide weight by a factor 2
            weights = np.array(dpbuffer.get_weights(), copy=False)
            weights /= 2.

    def update_info(self, dpinfo) :
        super().update_info(dpinfo)
        print("MockPyStep.update_info")
        self.info().set_need_vis_data()
        self.fetch_uvw = True
        self.fetch_weigths = True
