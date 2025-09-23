# Copyright (C) 2025 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import re

# Append current directory to system path in order to import testconfig
import sys

import pytest

sys.path.append(".")

import os
import socket
import struct
import subprocess

import numpy as np
import testconfig as tcf
from utils import spawn_dp3

# Define the socket file path
server_address = "./svpsock"
# Default MS
output_ms = "./svpout.ms"
# Default timeout
default_timeout = 40


def test_running():
    # Remove the file if it already exists
    try:
        os.unlink(server_address)
    except FileNotFoundError:
        pass
    # Create a Unix domain socket
    socket.setdefaulttimeout(default_timeout)
    server_socket = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)

    # Bind the socket to the address
    server_socket.bind(server_address)

    # Listen for incoming connections
    server_socket.listen(1)
    # Spawn the client process
    dp3_process = spawn_dp3(
        [
            f"steps=[]",
            "msout=" + output_ms,
            "msout.overwrite=True",
            "stream.socket=" + server_address,
        ]
    )

    # Wait for a connection within given timeout
    try:
        connection, _ = server_socket.accept()
    except OSError as msg:
        server_socket.close()
        dp3_process.kill()
        raise

    kTelescope = b"ALMA"
    kChan = 2
    kAnt = 3
    kPol = int(2)
    kScale = 1.0
    kStart = 4.3e5
    kInt = 1.0
    kTime = 3
    kSource = b"TESTSRC"
    kRA = 1.11
    kDec = 0.11
    kFrame = b"ICRS"
    kAntCoord = np.array([1.0, 2.0, 3.0, 4, 5, 6, 7, 8, 9], dtype=float)
    kAntDiam = np.array([70, 70, 70], dtype=float)
    kChanFreq = np.array([1e9, 1.1e9], dtype=float)
    kChanWidth = np.array([1e3, 1.1e3], dtype=float)
    # Send data
    connection.send(struct.pack("i", len(kTelescope)))
    connection.send(kTelescope)
    connection.send(struct.pack("N", kChan))
    connection.send(struct.pack("N", kAnt))
    connection.send(struct.pack("i", kPol))
    connection.send(struct.pack("d", kScale))
    connection.send(struct.pack("d", kStart))
    connection.send(struct.pack("d", kInt))
    connection.send(struct.pack("N", kTime))
    connection.send(struct.pack("i", len(kSource)))
    connection.send(kSource)
    connection.send(struct.pack("d", kRA))
    connection.send(struct.pack("d", kDec))
    connection.send(struct.pack("i", len(kFrame)))
    connection.send(kFrame)
    connection.sendall(kAntCoord.data)
    connection.sendall(kAntDiam.data)
    connection.sendall(kChanFreq.data)
    connection.sendall(kChanWidth.data)

    # Clean up the connection
    connection.close()
