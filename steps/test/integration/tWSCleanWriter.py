# Copyright (C) 2024 ASTRON (Netherlands Institute for Radio Astronomy)
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import re
import struct

# Append current directory to system path in order to import testconfig
import sys
from subprocess import check_call, check_output

import numpy as np
import pytest
from casacore.tables import table
from numpy.testing import assert_almost_equal

sys.path.append(".")

import testconfig as tcf
from utils import run_in_tmp_path, untar

MSIN = "tNDPPP-generic.MS"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")


SECS_IN_DAY = 86400

# Reordered files C struct packing information

# struct MetaHeaderBuffer {
#   double startTime;
#   uint64_t selectedRowCount;
#   uint32_t filenameLength;
# };
_meta_header_fmt = "=dQL"
_meta_header_binary_size = 20

# struct MetaRecordBuffer {
#   double u;
#   double v;
#   double w;
#   double time;
#   uint16_t antenna1;
#   uint16_t antenna2;
#   uint16_t fieldId;
# };
_meta_record_fmt = "=ddddHHH"
_meta_record_binary_size = 38

# struct PartHeaderBuffer {
#   uint64_t channelCount;
#   uint64_t channelStart;
#   uint32_t dataDescId;
#   bool hasModel;
# };
_part_header_fmt = "=QQL?"
_part_header_binary_size = 21

_float_bin_size = 4
_complex_float_bin_size = 2 * 4

# Enums to keep the extract functions meaningful

XX = 0
XY = 1
YX = 2
YY = 3

RR = 0
RL = 1
LR = 2
LL = 3


def assert_reorder_ms_meta_file(ms_filename, ms_table):
    """
    Reads the Measurement Set with python-casacore, reorder
    the metadata information and assert against the actual reordered
    output created from DP3
    """
    reorder_metaname = f"{ms_filename}-spw0-parted-meta.tmp"

    assert os.path.exists(reorder_metaname)

    nr_rows_expected = ms_table.nrows()

    uvws = ms_table.getcol("UVW")
    times = ms_table.getcol("TIME")
    antenna_1 = ms_table.getcol("ANTENNA1")
    antenna_2 = ms_table.getcol("ANTENNA2")

    start_time_expected = times[0] / SECS_IN_DAY

    with open(reorder_metaname, "rb") as f:
        data_raw = f.read(_meta_header_binary_size)
        (
            start_time_actual,
            nr_rows_actual,
            filename_length_actual,
        ) = struct.unpack_from(_meta_header_fmt, data_raw)

        data_raw = f.read(filename_length_actual)
        filename_actual = struct.unpack_from(f"{len(data_raw)}s", data_raw)[
            0
        ].decode("utf-8")

        assert_almost_equal(start_time_actual, start_time_expected, decimal=3)
        assert nr_rows_actual == nr_rows_expected
        assert (
            ms_filename in filename_actual
        )  # filename_actual will contain the tmp directory of the test, so assert with a contains

        for row in range(nr_rows_actual):
            data_raw = f.read(_meta_record_binary_size)
            (u, v, w, time, antenna1, antenna2, fieldid) = struct.unpack_from(
                _meta_record_fmt, data_raw
            )

            assert uvws[row][0] == u
            assert uvws[row][1] == v
            assert uvws[row][2] == w
            assert times[row] == time
            assert antenna_1[row] == antenna1
            assert antenna_2[row] == antenna2
            assert fieldid == 0


def _extract_pol(data, pol, corr_type="linear"):
    """
    Perform reorder of visibility data given a polarization.
    """
    instrumental = False
    if pol == "instr":
        # pol row chan -> row chan pol
        pol_data = data.transpose(1, 2, 0)
        instrumental = True
    elif pol == "diag_instr":
        pol_data = np.stack([data[0], data[3]])
        # pol row chan -> row chan pol
        pol_data = pol_data.transpose(1, 2, 0)
        instrumental = True
    elif corr_type == "linear":
        # I = (XX + YY)/2
        if pol == "I":
            pol_data = (data[XX] + data[YY]) / 2
        # Q = (XX - YY)/2
        elif pol == "Q":
            pol_data = (data[XX] - data[YY]) / 2
        # U = (XY + YX)/2
        elif pol == "U":
            pol_data = (data[XY] + data[YX]) / 2
        # V = -i(XY - YX)/2
        elif pol == "V":
            pol_data = -1j * (data[XY] - data[YX]) / 2

        # For cases without stokes pol
        if pol == "XX":
            pol_data = data[XX]
        elif pol == "XY":
            pol_data = data[XY]
        elif pol == "YX":
            pol_data = data[YX]
        elif pol == "YY":
            pol_data = data[YY]

    elif corr_type == "circular":
        # I = (LL + RR)/2
        if pol == "I":
            pol_data = (data[LL] + data[RR]) / 2
        # Q = (RL + LR)/2
        elif pol == "Q":
            pol_data = (data[RL] + data[LR]) / 2
        # U = -i (RL - LR)/2
        elif pol == "U":
            pol_data = -1j * (data[RL] - data[LR]) / 2
        # V = (RR - LL)/2
        elif pol == "V":
            pol_data = (data[RR] - data[LL]) / 2

        # For cases without stokes pol
        if pol == "LL":
            pol_data = data[LL]
        elif pol == "LR":
            pol_data = data[LR]
        elif pol == "RL":
            pol_data = data[RL]
        elif pol == "RR":
            pol_data = data[RR]

    if not instrumental:
        pol_data = pol_data.reshape((pol_data.shape[0], pol_data.shape[1], 1))

    return pol_data


def _extract_weights(weights, flags, pol, corr_type="linear"):
    """
    Perform reorder of weights data given a polarization.
    """
    instrumental = False
    if pol == "instr":
        weights_data = np.where(flags, 0.0, weights)
        weights_data = weights_data.transpose(1, 2, 0)
        instrumental = True
    elif pol == "diag_instr":
        selected_weights = np.stack([weights[0], weights[3]])
        selected_flags = np.stack([flags[0], flags[3]])
        weights_data = np.where(selected_flags, 0.0, selected_weights)
        weights_data = weights_data.transpose(1, 2, 0)
        instrumental = True
    elif corr_type == "linear":
        if pol == "I":
            corr1 = np.where(flags[XX], 0.0, weights[XX])
            corr2 = np.where(flags[YY], 0.0, weights[YY])
            weights_data = np.minimum(corr1, corr2)
        elif pol == "Q":
            corr1 = np.where(flags[XX], 0.0, weights[XX])
            corr2 = np.where(flags[YY], 0.0, weights[YY])
            weights_data = np.minimum(corr1, corr2)
        elif pol == "U":
            corr1 = np.where(flags[XY], 0.0, weights[XY])
            corr2 = np.where(flags[YX], 0.0, weights[YX])
            weights_data = np.minimum(corr1, corr2)
        elif pol == "V":
            corr1 = np.where(flags[XY], 0.0, weights[XY])
            corr2 = np.where(flags[YX], 0.0, weights[YX])
            weights_data = np.minimum(corr1, corr2)

        if pol == "XX":
            weights_data = np.where(flags[XX], 0.0, weights[XX])
        elif pol == "XY":
            weights_data = np.where(flags[XY], 0.0, weights[XY])
        elif pol == "YX":
            weights_data = np.where(flags[YX], 0.0, weights[YX])
        elif pol == "YY":
            weights_data = np.where(flags[YY], 0.0, weights[YY])
    elif corr_type == "circular":
        if pol == "I":
            corr1 = np.where(flags[LL], 0.0, weights[LL])
            corr2 = np.where(flags[RR], 0.0, weights[RR])
            weights_data = np.minimum(corr1, corr2)
        elif pol in ["Q", "U"]:
            corr1 = np.where(flags[RL], 0.0, weights[RL])
            corr2 = np.where(flags[LR], 0.0, weights[LR])
            weights_data = np.minimum(corr1, corr2)
        elif pol == "V":
            corr1 = np.where(flags[RR], 0.0, weights[RR])
            corr2 = np.where(flags[LL], 0.0, weights[LL])
            weights_data = np.minimum(corr1, corr2)

        if pol == "LL":
            weights_data = np.where(flags[LL], 0.0, weights[LL])
        elif pol == "LR":
            weights_data = np.where(flags[LR], 0.0, weights[LR])
        elif pol == "RL":
            weights_data = np.where(flags[RL], 0.0, weights[RL])
        elif pol == "RR":
            weights_data = np.where(flags[RR], 0.0, weights[RR])

    needs_scaling = pol not in ["XX", "XY", "YX", "YY", "LL", "LR", "RL", "RR"]
    if needs_scaling:
        weights_data = weights_data * 4

    if not instrumental:
        weights_data = weights_data.reshape(
            (weights_data.shape[0], weights_data.shape[1], 1)
        )

    return weights_data


def assert_reorder_ms_data_files(
    ms_filename,
    ms_table,
    pols=["I"],
    corr_type="linear",
):
    """
    Reorders the data and weights of a MS and asserts it against the
    actual results from DP3
    """
    npol_per_file = 1
    if pols == ["instr"]:
        npol_per_file = 4
    elif pols == ["diag_instr"]:
        npol_per_file = 2

    data = ms_table.getcol("DATA")
    nr_rows, nr_chan, nr_corr = data.shape

    if "WEIGHT_SPECTRUM" in ms_table.colnames():
        weight = ms_table.getcol("WEIGHT_SPECTRUM")
    else:
        weight = ms_table.getcol("WEIGHT")
    flags = ms_table.getcol("FLAG")

    if weight.shape != data.shape:
        # If the weights don't have a frequency dimension, replicate the
        # values per each channel
        weight = np.repeat(weight, nr_chan).reshape(nr_rows, nr_chan, nr_corr)

    # Transpose data for ease of use
    # "rows chan pol -> pol rows chan"
    data_transposed = data.transpose(2, 0, 1)
    weight_transposed = weight.transpose(2, 0, 1)
    flags_transposed = flags.transpose(2, 0, 1)

    for pol in pols:
        # Data array
        data_file_name = f"{ms_filename}-part0000-{pol}-b0.tmp"
        assert os.path.exists(data_file_name)

        pol_data = _extract_pol(data_transposed, pol, corr_type)

        with open(data_file_name, "rb") as f:
            data_raw = f.read(_part_header_binary_size)

            (nr_chan_actual, chan_start, ddi, has_model) = struct.unpack_from(
                _part_header_fmt, data_raw
            )

            assert nr_chan_actual == nr_chan
            assert chan_start == 0
            assert ddi == 0
            assert has_model == False

            for row in range(nr_rows):
                for chan in range(nr_chan):
                    for pol_idx in range(npol_per_file):
                        data_raw = f.read(_complex_float_bin_size)
                        (real, imag) = struct.unpack_from("ff", data_raw)
                        vis_actual = real + 1j * imag
                        assert_almost_equal(
                            vis_actual, pol_data[row][chan][pol_idx]
                        )

        weights_file_name = f"{ms_filename}-part0000-{pol}-b0-w.tmp"
        assert os.path.exists(weights_file_name)
        weights_data = _extract_weights(
            weight_transposed, flags_transposed, pol, corr_type
        )

        # Weights
        with open(weights_file_name, "rb") as f:
            for row in range(nr_rows):
                for chan in range(nr_chan):
                    for pol_idx in range(npol_per_file):
                        data_raw = f.read(_float_bin_size)
                        weight = struct.unpack_from("f", data_raw)[0]
                        assert_almost_equal(
                            weight, weights_data[row][chan][pol_idx]
                        )


def test_reorder():
    """
    Test that reordering operation yields the expected results.
    """
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"wscleanwriter.name={MSIN}",
            "steps=[wscleanwriter]",
        ]
    )

    # Assert the content of the reordered files against a reference
    with table(MSIN) as t:
        ms_table = t.query("ANTENNA1 != ANTENNA2")

        assert_reorder_ms_meta_file(MSIN, ms_table)
        assert_reorder_ms_data_files(MSIN, ms_table, pols=["I"])


def test_reorder_resolve_dot_as_name():
    """
    Test that reordering operation uses msin name when '.' is given as name.
    """
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "wscleanwriter.name=.",
            "steps=[wscleanwriter]",
        ]
    )

    # Assert the content of the reordered files against a reference
    with table(MSIN) as t:
        ms_table = t.query("ANTENNA1 != ANTENNA2")

        assert_reorder_ms_meta_file(MSIN, ms_table)
        assert_reorder_ms_data_files(MSIN, ms_table, pols=["I"])


def test_reorder_use_msout_name():
    """
    Test that reordering operation uses the name from msout when name
    is not provided.
    """
    msout_name = "test.ms"

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout={msout_name}",
            "steps=[wscleanwriter,msout]",
        ]
    )

    # Assert the content of the reordered files against a reference
    with table(MSIN) as t:
        ms_table = t.query("ANTENNA1 != ANTENNA2")

        assert_reorder_ms_meta_file(msout_name, ms_table)
        assert_reorder_ms_data_files(msout_name, ms_table, pols=["I"])


def test_reorder_handle_dot_msout_name():
    """
    Test that reordering operation uses the name from msin when msout
    name '.'
    """
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout=.",
            "steps=[wscleanwriter,msout]",
        ]
    )

    # Assert the content of the reordered files against a reference
    with table(MSIN) as t:
        ms_table = t.query("ANTENNA1 != ANTENNA2")

        assert_reorder_ms_meta_file(MSIN, ms_table)
        assert_reorder_ms_data_files(MSIN, ms_table, pols=["I"])


def test_reorder_use_msin_name():
    """
    Test that reordering operation uses the name from msin when name
    is not provided and output step is not present
    """
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "steps=[wscleanwriter]",
        ]
    )

    # Assert the content of the reordered files against a reference
    with table(MSIN) as t:
        ms_table = t.query("ANTENNA1 != ANTENNA2")

        assert_reorder_ms_meta_file(MSIN, ms_table)
        assert_reorder_ms_data_files(MSIN, ms_table, pols=["I"])


def test_reorder_change_output_name():
    """
    Test that reordering operation yields the expected results when a new
    filename is given as msout.
    """
    msout_name = "test.ms"
    msout_path = f"{msout_name}/"  # Create a path with trailing separator
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            f"msout={msout_path}",
            "steps=[wscleanwriter,msout]",
        ]
    )

    # Assert the content of the reordered files against a reference
    with table(MSIN) as t:
        ms_table = t.query("ANTENNA1 != ANTENNA2")

        assert_reorder_ms_meta_file(msout_name, ms_table)
        assert_reorder_ms_data_files(msout_name, ms_table, pols=["I"])


def test_reorder_prints_warning_when_meta_data_changes():
    """
    Test that reordering operation prints warning when meta data is
    changed in a previous step
    """
    dp3_output = check_output(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=test.ms",
            "steps=[averager, wscleanwriter, msout]",
            "averager.freqstep=2",
        ]
    )

    assert "Meta data changes detected" in dp3_output.decode("utf-8")


def test_reorder_multiple_stokes_pol():
    """
    Test reordering when multiple polarizations are given
    """
    pols = "IQUV"
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "wscleanwriter.name=.",
            f"wscleanwriter.polarization={pols}",
            "steps=[wscleanwriter]",
        ]
    )

    # Assert the content of the reordered files against a reference
    with table(MSIN) as t:
        ms_table = t.query("ANTENNA1 != ANTENNA2")

        assert_reorder_ms_meta_file(MSIN, ms_table)
        assert_reorder_ms_data_files(MSIN, ms_table, pols=list(pols))


def test_reorder_tmp_dir():
    """
    Test reordering when writing the reordered files in a temporary directory.
    """
    tmp_dir = "tmp"
    os.system(f"mkdir {tmp_dir}")
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "wscleanwriter.name=.",
            f"wscleanwriter.temporary_directory={tmp_dir}",
            "steps=[wscleanwriter]",
        ]
    )

    assert os.path.exists(f"./{tmp_dir}/{MSIN}-spw0-parted-meta.tmp")
    assert os.path.exists(f"./{tmp_dir}/{MSIN}-part0000-I-b0.tmp")
    assert os.path.exists(f"./{tmp_dir}/{MSIN}-part0000-I-b0-w.tmp")


def test_reorder_circular_pol_to_stokes():
    """
    Test reordering given a MS with circular polarization
    and reorder to be done into Stokes polarization
    """
    taqlcommand_run = f"UPDATE {MSIN}/POLARIZATION SET CORR_TYPE=[5, 6, 7, 8]"
    check_call([tcf.TAQLEXE, "-noph", taqlcommand_run])

    pols = "IQUV"

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "wscleanwriter.name=.",
            f"wscleanwriter.polarization={pols}",
            "steps=[wscleanwriter]",
        ]
    )

    with table(MSIN) as t:
        ms_table = t.query("ANTENNA1 != ANTENNA2")

        assert_reorder_ms_meta_file(MSIN, ms_table)
        assert_reorder_ms_data_files(
            MSIN,
            ms_table,
            pols=list(pols),
            corr_type="circular",
        )


def test_reorder_linear_pol():
    """
    Test reordering when output polarization is expected to be linear.
    """
    pols = "XX,XY,YX,YY"

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "wscleanwriter.name=.",
            f"wscleanwriter.polarization={pols}",
            "steps=[wscleanwriter]",
        ]
    )

    with table(MSIN) as t:
        ms_table = t.query("ANTENNA1 != ANTENNA2")

        assert_reorder_ms_meta_file(MSIN, ms_table)
        assert_reorder_ms_data_files(
            MSIN,
            ms_table,
            pols=pols.split(","),
        )


def test_reorder_circular_pol():
    """
    Test reordering when given a measurement set with
    circular polarization and output polarization is also circular polarization
    """
    taqlcommand_run = f"UPDATE {MSIN}/POLARIZATION SET CORR_TYPE=[5, 6, 7, 8]"
    check_call([tcf.TAQLEXE, "-noph", taqlcommand_run])

    pols = "LL,LR,RL,RR"

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "wscleanwriter.name=.",
            f"wscleanwriter.polarization={pols}",
            "steps=[wscleanwriter]",
        ]
    )

    with table(MSIN) as t:
        ms_table = t.query("ANTENNA1 != ANTENNA2")

        assert_reorder_ms_meta_file(MSIN, ms_table)
        assert_reorder_ms_data_files(
            MSIN, ms_table, pols=pols.split(","), corr_type="circular"
        )


def test_reorder_for_instrumental_pol():
    """
    Test reordering to output instrumental polarization.
    """
    pol = "instr"
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "wscleanwriter.name=.",
            f"wscleanwriter.polarization={pol}",
            "steps=[wscleanwriter]",
        ]
    )

    # Assert the content of the reordered files against a reference
    with table(MSIN) as t:
        ms_table = t.query("ANTENNA1 != ANTENNA2")

        assert_reorder_ms_meta_file(MSIN, ms_table)
        assert_reorder_ms_data_files(MSIN, ms_table, pols=[pol])


def test_reorder_for_diagonal_instrumental_pol():
    """
    Test reordering to output diagonal instrumental polarization.
    """
    pol = "diag_instr"
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "wscleanwriter.name=.",
            f"wscleanwriter.polarization={pol}",
            "steps=[wscleanwriter]",
        ]
    )

    # Assert the content of the reordered files against a reference
    with table(MSIN) as t:
        ms_table = t.query("ANTENNA1 != ANTENNA2")

        assert_reorder_ms_meta_file(MSIN, ms_table)
        assert_reorder_ms_data_files(MSIN, ms_table, pols=[pol])


def test_reorder_for_diagonal_instrumental_pol_for_circular_corr():
    """
    Test reordering when given a MS with circular
    polarization and is expected to output diagonal instrumental polarization
    """
    taqlcommand_run = f"UPDATE {MSIN}/POLARIZATION SET CORR_TYPE=[5, 6, 7, 8]"
    check_call([tcf.TAQLEXE, "-noph", taqlcommand_run])

    pol = "diag_instr"
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "wscleanwriter.name=.",
            f"wscleanwriter.polarization={pol}",
            "steps=[wscleanwriter]",
        ]
    )

    # Assert the content of the reordered files against a reference
    with table(MSIN) as t:
        ms_table = t.query("ANTENNA1 != ANTENNA2")

        assert_reorder_ms_meta_file(MSIN, ms_table)
        assert_reorder_ms_data_files(MSIN, ms_table, pols=[pol])


def test_reorder_split_channel():
    """
    Test reordering when chanperfile is provided and is expected to
    write into multiple reordered file parts.
    """
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "wscleanwriter.name=.",
            "wscleanwriter.chanperfile=4",
            "steps=[wscleanwriter]",
        ]
    )

    assert os.path.exists(f"./{MSIN}-spw0-parted-meta.tmp")
    assert os.path.exists(f"./{MSIN}-part0000-I-b0.tmp")
    assert os.path.exists(f"./{MSIN}-part0000-I-b0-w.tmp")
    assert os.path.exists(f"./{MSIN}-part0001-I-b0.tmp")
    assert os.path.exists(f"./{MSIN}-part0001-I-b0-w.tmp")
    assert not os.path.exists(f"./{MSIN}-part0002-I-b0.tmp")
    assert not os.path.exists(f"./{MSIN}-part0002-I-b0-w.tmp")


def test_reorder_split_channel_uneven():
    """
    Test reordering when chanperfile is not divisible by nchan is provided
    and is expected to write reorder into multiple file parts.
    """
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "wscleanwriter.name=.",
            "wscleanwriter.chanperfile=3",
            "steps=[wscleanwriter]",
        ]
    )

    assert os.path.exists(f"./{MSIN}-spw0-parted-meta.tmp")
    assert os.path.exists(f"./{MSIN}-part0000-I-b0.tmp")
    assert os.path.exists(f"./{MSIN}-part0000-I-b0-w.tmp")
    assert os.path.exists(f"./{MSIN}-part0001-I-b0.tmp")
    assert os.path.exists(f"./{MSIN}-part0001-I-b0-w.tmp")
    assert os.path.exists(f"./{MSIN}-part0002-I-b0.tmp")
    assert os.path.exists(f"./{MSIN}-part0002-I-b0-w.tmp")
    assert not os.path.exists(f"./{MSIN}-part0003-I-b0.tmp")
    assert not os.path.exists(f"./{MSIN}-part0003-I-b0-w.tmp")


def test_reorder_split_channel_nchan_per_file_larger_than_nchan():
    """
    Test reordering when chanperfile is greater than nchan, DP3 is
    expected to write a single file part.
    """
    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "wscleanwriter.name=.",
            "wscleanwriter.chanperfile=10",
            "steps=[wscleanwriter]",
        ]
    )

    assert os.path.exists(f"./{MSIN}-spw0-parted-meta.tmp")
    assert os.path.exists(f"./{MSIN}-part0000-I-b0.tmp")
    assert os.path.exists(f"./{MSIN}-part0000-I-b0-w.tmp")
    assert not os.path.exists(f"./{MSIN}-part0001-I-b0.tmp")
    assert not os.path.exists(f"./{MSIN}-part0001-I-b0-w.tmp")
