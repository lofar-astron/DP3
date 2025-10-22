import sys
from subprocess import check_call

import pytest

# Append current directory to system path in order to import testconfig
sys.path.append(".")

import testconfig as tcf
from utils import get_taql_result, run_in_tmp_path, untar

"""
Tests for predicting gaussian sources

Script can be invoked in two ways:
- as standalone from the build/steps/test/integration directory,
  using `pytest source/tPhaseShiftPredict.py` (extended with pytest options of your choice)
- using ctest, see DP3/steps/test/integration/CMakeLists.txt
"""

MSIN = "tNDPPP-generic.MS"


@pytest.fixture(autouse=True)
def source_env(run_in_tmp_path):
    untar(f"{tcf.RESOURCEDIR}/{MSIN}.tgz")


testdata = [
    ("POINT", 0, 0, 0, 2, 10),
    ("GAUSSIAN", 600, 60, 0, 2, 10),
    ("GAUSSIAN", 600, 60, 30, 0, 0),
    ("GAUSSIAN", 600, 60, 30, 2, 0),
    ("GAUSSIAN", 600, 60, 30, 0, 10),
    ("GAUSSIAN", 600, 60, 30, 2, 10),
    ("GAUSSIAN", 60, 60, 0, 2, 10),
]


@pytest.mark.parametrize(
    "source_type,major_axis,minor_axis,orientation,offset_ra_hour,offset_dec_degree",
    testdata,
)
@pytest.mark.parametrize(
    "use_fast_predict",
    [
        False,
        pytest.param(
            True,
            marks=pytest.mark.skipif(
                tcf.USE_FAST_PREDICT != "ON",
                reason="FastPredict is not available",
            ),
        ),
    ],
)
def test_phaseshift_predict(
    source_type,
    major_axis,
    minor_axis,
    orientation,
    offset_ra_hour,
    offset_dec_degree,
    use_fast_predict,
):
    """
    - Phaseshift MS to the position of a source
    - Predict the source with the new phase center
    - Phaseshift back to the original phase center
    - Predict the source again, subtract it
    Since the source model should be independent of the phase center, the result should be zero
    """
    MS_PHASECENTER_RA = "01h37m41.299"
    MS_PHASECENTER_DEC = "+033d09m35.132"
    assert -1 < offset_ra_hour < 22 and int(offset_ra_hour) == offset_ra_hour
    assert (
        -33 < offset_dec_degree < 56
        and int(offset_dec_degree) == offset_dec_degree
    )
    source_position_ra = MS_PHASECENTER_RA.replace(
        "01", f"{1 + offset_ra_hour:02d}"
    )
    source_position_dec = MS_PHASECENTER_DEC.replace(
        "+033", f"{33 + offset_dec_degree:+04d}"
    )
    with open("test.skymodel", "w") as f:
        f.write(
            f"""\
FORMAT = Name, Type, Ra, Dec, I, MajorAxis, MinorAxis, Orientation, OrientationIsAbsolute
dummysource, {source_type}, {source_position_ra}, {source_position_dec}, {major_axis}, {minor_axis}, 2, 112, True
"""
        )

    check_call(
        [
            tcf.DP3EXE,
            f"msin={MSIN}",
            "msout=out.MS",
            "steps=[shift,predict,shiftback,subtract]",
            "shift.type=phaseshift",
            f"shift.phasecenter=[{source_position_ra},{source_position_dec}]",
            "predict.type=predict",
            "predict.operation=replace",
            "predict.sourcedb=test.skymodel",
            f"predict.usefastpredict={use_fast_predict}",
            "shiftback.type=phaseshift",
            "shiftback.phasecenter=[]",
            "subtract.type=predict",
            "subtract.sourcedb=test.skymodel",
            "subtract.operation=subtract",
            f"subtract.usefastpredict={use_fast_predict}",
        ]
    )

    taql_command = (
        f"select gmean(abs(DATA)) from out.MS WHERE ANTENNA1!=ANTENNA2"
    )
    residual = float(get_taql_result(taql_command))
    print("Residual:", residual)
    assert residual < 0.001
