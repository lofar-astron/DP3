# Copyright (C) 2022 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

"""
These tests check the collection of output data (DPBuffers) 
in a queue by the QueueOutput step. 

Script can be invoked in two ways:
- as standalone from the build/pythondp3/test/integration directory,
  using `pytest source/tQueueOutput.py` (extended with pytest options of your choice)
- using ctest, see pythondp3/test/integration/CMakeLists.txt
"""

import numpy as np

# Append current directory to system path in order to import testconfig
import sys

sys.path.append(".")

import testconfig as tcf

sys.path.insert(0, tcf.PYTHONDIR)
sys.path.insert(0, tcf.PYTHONMOCKDIR)

try:
    "The import may fail while running pytest --collect-only"
    import dp3
    import dp3.steps
except ImportError:
    pass


def test_queue_output():
    """
    Creates a QueueOuput step.
    Checks that data is properly collected by
    the QueueOutput step.
    """

    parset = dp3.parameterset.ParameterSet()

    queue_step = dp3.steps.QueueOutput(parset, "")

    # Newly created QueueStep should have an empty queue
    assert queue_step.queue.empty()

    # Define dimentsions of input data
    n_baselines = 10
    n_channels = 2
    n_correlations = 4

    # loop over three time steps
    n_time = 3
    for i in range(n_time):
        # Create an empty buffer
        dpbuffer_in = dp3.DPBuffer()

        # Create input data
        data_in = i * np.ones(
            (n_baselines, n_channels, n_correlations), dtype=np.complex64
        )
        weights_in = (
            10
            * i
            * np.ones(
                (n_baselines, n_channels, n_correlations), dtype=np.float32
            )
        )

        # Put input data in buffer
        dpbuffer_in.set_data(data_in)
        dpbuffer_in.set_weights(weights_in)

        # Process input data
        # The output should end up in the queue of the QueueOutput step
        queue_step.process(dpbuffer_in)

    for i in range(n_time):
        assert not queue_step.queue.empty()
        # Get data out of queue in QueueOutput step
        dpbuffer_from_queue = queue_step.queue.get()

        # Extract numpy arrays from the DPBuffer
        data_out = np.array(dpbuffer_from_queue.get_data(), copy=False)
        weights_out = np.array(dpbuffer_from_queue.get_weights(), copy=False)

        # Check whether the output values are as expected
        assert data_out.shape[0] == n_baselines
        assert data_out.shape[1] == n_channels
        assert data_out.shape[2] == n_correlations
        assert weights_out.shape[0] == n_baselines
        assert weights_out.shape[1] == n_channels
        assert weights_out.shape[2] == n_correlations
        assert np.all(data_out == i)
        assert np.all(weights_out == 10 * i)

    # All elements have been read from the queue
    # Queue should be empty again
    assert queue_step.queue.empty()


def test_chain_queues():
    """
    Creates two QueueOuput steps. Connects them and checks
    that both properly collect the data.
    """
    parset = dp3.parameterset.ParameterSet()

    queue_step1 = dp3.steps.QueueOutput(parset, "")
    queue_step2 = dp3.steps.QueueOutput(parset, "")
    queue_step1.set_next_step(queue_step2)

    # loop over three time steps
    n_time = 3
    for i in range(n_time):
        # Create an empty buffer
        dpbuffer_in = dp3.DPBuffer()

        # Create input data
        data_in = i * np.ones((1, 1, 1), dtype=np.complex64)
        weights_in = 10 * i * np.ones((1, 1, 1), dtype=np.float32)

        # Put input data in buffer
        dpbuffer_in.set_data(data_in)
        dpbuffer_in.set_weights(weights_in)

        # Process input data
        # The output should end up in the queues of the QueueOutput steps
        queue_step1.process(dpbuffer_in)

    for i in range(n_time):
        assert not queue_step1.queue.empty()
        assert not queue_step2.queue.empty()

        # Get data out of queue in QueueOutput step
        dpbuffer_from_queue1 = queue_step1.queue.get()
        dpbuffer_from_queue2 = queue_step2.queue.get()

        # Extract numpy arrays from the DPBuffer
        data_out1 = np.array(dpbuffer_from_queue1.get_data(), copy=False)
        weights_out1 = np.array(dpbuffer_from_queue1.get_weights(), copy=False)
        data_out2 = np.array(dpbuffer_from_queue2.get_data(), copy=False)
        weights_out2 = np.array(dpbuffer_from_queue2.get_weights(), copy=False)

        # Check whether the output values are as expected
        assert np.all(data_out1 == i)
        assert np.all(weights_out1 == 10 * i)
        assert np.all(data_out2 == i)
        assert np.all(weights_out2 == 10 * i)

    # All elements have been read from both queues
    # Queues should be empty again
    assert queue_step1.queue.empty()
    assert queue_step2.queue.empty()
