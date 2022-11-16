# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

# queue.py: simple python step, which puts the data it receives in a queue

from dp3 import parameterset
from dp3 import Step, Fields
import numpy as np
import queue


class QueueOutput(Step):
    """Python Step that collects data send to it in a queue"""

    def __init__(self, parset, prefix):
        """
        Set up the step (constructor). Read the parset here.

        Args:
          parset: Parameter set, not used by this step
          prefix: Prefix for this step, e.g. "thisstepname.",
        """

        super().__init__()
        self.queue = queue.Queue()

    def update_info(self, dpinfo):
        """
        Process metadata. This will be called before any call to process.

        Args:
          dpinfo: DPInfo object with all metadata, see docs in pydp3.cc
        """
        super().update_info(dpinfo)

    def show(self):
        """Print a summary of the step and its settings"""
        print("QueueOutput")

    def get_required_fields(self):
        return Fields()

    def get_provided_fields(self):
        return Fields()

    def process(self, dpbuffer):
        """
        Process one time slot of data. This step may or may not be the
        end point of a step chain. If there is a next step,
        this function will call the next step's process() function.

        Args:
          dpbuffer: DPBuffer object which can contain data, flags and weights
                    for one time slot.
        """

        # Send data to the queue
        self.queue.put(dpbuffer)
        next_step = self.get_next_step()
        if next_step is not None:
            next_step.process(dpbuffer)

    def finish(self):
        """
        If there is any remaining data, process it. This can be useful if the
        step accumulates multiple time slots.
        """
        pass
