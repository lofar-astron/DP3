# Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import dppp
import numpy as np

class DummyDPStep(dppp.DPStep):

    def __init__(self,parset, prefix):
        super().__init__()
        self.parset = parset
        self.dpbuffers = []
        self.prefix = prefix

    def show(self) :
        print()
        print("DummyDPStep")
        v = self.parset.getString(self.prefix + "somekey")
        print("    somekey: {v}")

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
        print(f"Processing a chunk of {len(self.dpbuffers)}.")

        # Enumerate accumulated data and display just the shapes
        for i, dpbuffer in enumerate(self.dpbuffers):
            data = np.array(dpbuffer.get_data(), copy=False)
            flags = np.array(dpbuffer.get_flags(), copy=False)
            uvw = np.array(dpbuffer.get_uvw(), copy=False)
            print(f"{i:4d} shape of data: ", data.shape)
            print(f"     shape of flags: ", flags.shape)
            print(f"     shape uvw: ", uvw.shape)

    def update_info(self, dpinfo) :
        super().update_info(dpinfo)
        print("DummyDPStep.update_info")
        self.info().set_need_vis_data()
        self.fetch_uvw = True
        self.fetch_weigths = True
