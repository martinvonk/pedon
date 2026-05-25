"""Test soilmodels."""

import pytest
from numpy import allclose, array, logspace

import pedon as pe
from pedon._typing import FloatArray

h = -logspace(-2, 6, num=9, dtype=float)

theta = array([0.1, 0.2, 0.3, 0.4])


def assert_close(actual: FloatArray, expected: FloatArray) -> None:
    """Assert if two arrays are close within a reasonable tolerance."""
    assert allclose(actual, expected, rtol=1e-8, atol=1e-8)


def test_theta_genuchten(genuchten: pe.SoilModel, h: FloatArray = h) -> None:
    """Test water content calculation for the van Genuchten model."""
    expected = array(
        [
            0.42999674186266773,
            0.4299590045910226,
            0.42948737007237936,
            0.4240392694915515,
            0.3884681132693863,
            0.3202485060610565,
            0.25718963002041634,
            0.20639799991988222,
            0.166007528803122,
        ]
    )
    assert_close(genuchten.theta(h=h), expected)


def test_s_genuchten(genuchten: pe.SoilModel, h: FloatArray = h) -> None:
    """Test degree of saturation calculation for the van Genuchten model."""
    expected = array(
        [
            0.9999922425301613,
            0.9999023918833871,
            0.9987794525532842,
            0.9858077845036941,
            0.9011145554033008,
            0.7386869191929917,
            0.5885467381438484,
            0.4676142855235291,
            0.3714464971502905,
        ]
    )
    assert_close(genuchten.s(h=h), expected)


def test_k_genuchten(genuchten: pe.SoilModel, h: FloatArray = h) -> None:
    """Test hydraulic conductivity calculation for the van Genuchten model."""
    expected = array(
        [
            3.2869753941674986,
            2.1425970739706632,
            1.0528919260290555,
            0.25654138503050566,
            0.011109108559820792,
            9.373293594158181e-05,
            5.47578373681202e-07,
            3.0887429051900053e-09,
            1.7373526209279636e-11,
        ]
    )
    assert_close(genuchten.k(h=h), expected)


def test_h_genuchten(genuchten: pe.SoilModel, theta: FloatArray = theta) -> None:
    """Test pressure head calculation for the van Genuchten model."""
    expected = array(
        [244927639.31052694, 139271.6694845786, 1998.5290326524432, 61.672352811394724]
    )
    assert_close(genuchten.h(theta=theta), expected)


def test_theta_brooks(brooks: pe.SoilModel, h: FloatArray = h) -> None:
    """Test water content calculation for the Brooks-Corey model."""
    expected = array(
        [
            0.43,
            0.43,
            0.43,
            0.43,
            0.0142,
            0.010042,
            0.01000042,
            0.0100000042,
            0.010000000042000001,
        ]
    )
    assert_close(brooks.theta(h=h), expected)


def test_s_brooks(brooks: pe.SoilModel, h: FloatArray = h) -> None:
    """Test degree of saturation calculation for the Brooks-Corey model."""
    expected = array([1.0, 1.0, 1.0, 1.0, 0.01, 0.0001, 1e-06, 1e-08, 1e-10])
    assert_close(brooks.s(h=h), expected)


def test_k_brooks(brooks: pe.SoilModel, h: FloatArray = h) -> None:
    """Test hydraulic conductivity calculation for the Brooks-Corey model."""
    expected = array(
        [
            10.0,
            10.0,
            10.0,
            10.0,
            1e-07,
            1.0000000000000003e-15,
            9.999999999999998e-24,
            1e-31,
            1.0000000000000001e-39,
        ]
    )
    assert_close(brooks.k(h=h), expected)


def test_h_brooks(brooks: pe.SoilModel, theta: FloatArray = theta) -> None:
    """Test pressure head calculation for the Brooks-Corey model."""
    expected = array([21.60246899, 14.86783883, 12.03443336, 10.37749043])
    assert_close(brooks.h(theta=theta), expected)


def test_theta_panday(panday: pe.SoilModel, h: FloatArray = h) -> None:
    """Test water content calculation for the Panday model."""
    expected = array(
        [
            0.4299967418626678,
            0.4299590045910227,
            0.4294873700723793,
            0.4240392694915515,
            0.3884681132693863,
            0.3202485060610564,
            0.2571896300204164,
            0.2063979999198822,
            0.16600752880312195,
        ]
    )
    assert_close(panday.theta(h=h), expected)


def test_s_panday(panday: pe.SoilModel, h: FloatArray = h) -> None:
    """Test degree of saturation calculation for the Panday model."""
    expected = array(
        [
            0.9999922425301614,
            0.9999023918833874,
            0.9987794525532842,
            0.9858077845036941,
            0.9011145554033008,
            0.7386869191929916,
            0.5885467381438485,
            0.46761428552352907,
            0.37144649715029043,
        ]
    )
    assert_close(panday.s(h=h), expected)


def test_k_panday(panday: pe.SoilModel, h: FloatArray = h) -> None:
    """Test hydraulic conductivity calculation for the Panday model."""
    expected = array(
        [
            9.999767277710188,
            9.997072042312654,
            9.96342825049768,
            9.58024751871584,
            7.317117250456762,
            4.030706962685341,
            2.0386509376131543,
            1.022499986411543,
            0.5124940191915482,
        ]
    )
    assert_close(panday.k(h=h), expected)


def test_h_panday(panday: pe.SoilModel, theta: FloatArray = theta) -> None:
    """Test pressure head calculation for the Panday model."""
    expected = array(
        [244927639.31052694, 139271.6694845786, 1998.5290326524432, 61.672352811394724]
    )
    assert_close(panday.h(theta=theta), expected)


def test_theta_gardner(gardner: pe.SoilModel, h: FloatArray = h) -> None:
    """Test water content calculation for the Gardner model."""
    expected = array(
        [
            0.42529591987340853,
            0.38520867817750715,
            0.1431345659901742,
            7.181731339805633e-06,
            7.262321084965386e-49,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )
    assert_close(gardner.theta(h=h), expected)


def test_s_gardner(gardner: pe.SoilModel, h: FloatArray = h) -> None:
    """Test degree of saturation calculation for the Gardner model."""
    expected = array(
        [
            0.9890602787753687,
            0.8958341352965282,
            0.33287108369807955,
            1.670170079024566e-05,
            1.6889118802245084e-48,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )
    assert_close(gardner.s(h=h), expected)


def test_k_gardner(gardner: pe.SoilModel, h: FloatArray = h) -> None:
    """Test hydraulic conductivity calculation for the Gardner model."""
    expected = array(
        [
            9.998000199986667,
            9.98001998667333,
            9.801986733067553,
            8.187307530779819,
            1.353352832366127,
            2.061153622438558e-08,
            1.3838965267367375e-86,
            0.0,
            0.0,
        ]
    )
    assert_close(gardner.k(h=h), expected)


def test_h_gardner(gardner: pe.SoilModel, theta: FloatArray = theta) -> None:
    """Test pressure head calculation for the Gardner model."""
    expected = array(
        [
            1.3260136569995606,
            0.6958798564905193,
            0.32727521275582455,
            0.06574605598147819,
        ]
    )
    assert_close(gardner.h(theta=theta), expected)


def test_theta_rucker(rucker: pe.SoilModel, h: FloatArray = h) -> None:
    """Test water content calculation for the Gardner-Rucker model."""
    expected = array(
        [
            0.42999999864525157,
            0.4299998646064057,
            0.4299865414790768,
            0.4287311305825537,
            0.3545658264218955,
            0.013113261891671047,
            0.01,
            0.01,
            0.01,
        ]
    )
    assert_close(rucker.theta(h=h), expected)


def test_s_rucker(rucker: pe.SoilModel, h: FloatArray = h) -> None:
    """Test degree of saturation calculation for the Gardner-Rucker model."""
    expected = array(
        [
            0.9999999967744085,
            0.9999996776342993,
            0.9999679559025638,
            0.9969788823394136,
            0.8203948248140369,
            0.007412528313502493,
            0.0,
            0.0,
            0.0,
        ]
    )
    assert_close(rucker.s(h=h), expected)


def test_k_rucker(rucker: pe.SoilModel, h: FloatArray = h) -> None:
    """Test hydraulic conductivity calculation for the Gardner-Rucker model."""
    expected = array(
        [
            9.998000199986667,
            9.98001998667333,
            9.801986733067553,
            8.187307530779819,
            1.353352832366127,
            2.061153622438558e-08,
            1.3838965267367375e-86,
            0.0,
            0.0,
        ]
    )
    assert_close(rucker.k(h=h), expected)


def test_h_rucker(rucker: pe.SoilModel, theta: FloatArray = theta) -> None:
    """Test pressure head calculation for the Gardner-Rucker model."""
    expected = array(
        [399.640933493472, 247.51431683188103, 148.39155295044985, 55.87543934217171]
    )
    assert_close(rucker.h(theta=theta), expected)


def test_theta_genuchtengardner(
    genuchtengardner: pe.SoilModel, h: FloatArray = h
) -> None:
    """Test water content calculation for the GenuchtenGardner model."""
    expected = array(
        [
            0.4299979007684259,
            0.42994728398473125,
            0.4286847798061082,
            0.4016372327065841,
            0.2421543984264464,
            0.10587778884646776,
            0.04822936278663853,
            0.025220333539544004,
            0.01605933898039554,
        ]
    )
    assert_close(genuchtengardner.theta(h=h), expected)


def test_s_genuchtengardner(genuchtengardner: pe.SoilModel, h: FloatArray = h) -> None:
    """Test degree of saturation calculation for the GenuchtenGardner model."""
    expected = array(
        [
            0.9999950018295856,
            0.9998744856779316,
            0.9968685233478767,
            0.932469601682343,
            0.5527485676820152,
            0.22828044963444707,
            0.09102229234913936,
            0.036238889379866676,
            0.01442699757237033,
        ]
    )
    assert_close(genuchtengardner.s(h=h), expected)


def test_k_genuchtengardner(genuchtengardner: pe.SoilModel, h: FloatArray = h) -> None:
    """Test hydraulic conductivity calculation for the GenuchtenGardner model."""
    expected = array(
        [
            9.998000199986667,
            9.98001998667333,
            9.801986733067553,
            8.187307530779819,
            1.353352832366127,
            2.061153622438558e-08,
            1.3838965267367375e-86,
            0.0,
            0.0,
        ]
    )
    assert_close(genuchtengardner.k(h=h), expected)


def test_h_genuchtengardner(
    genuchtengardner: pe.SoilModel, theta: FloatArray = theta
) -> None:
    """Test pressure head calculation for the GenuchtenGardner model."""
    expected = array(
        [1172.3053976502565, 173.47442816736208, 50.22660724364475, 10.481432805126213]
    )
    assert_close(genuchtengardner.h(theta=theta), expected)


def test_theta_haverkamp(haverkamp: pe.SoilModel, h: FloatArray = h) -> None:
    """Test water content calculation for the Haverkamp model."""
    expected = array(
        [
            0.3602762555045123,
            0.11108822138463235,
            0.018235294117647058,
            0.010529336192246969,
            0.010033438339917915,
            0.010002109974002438,
            0.010000133130985966,
            0.010000008399999833,
            0.01000000053000417,
        ]
    )
    assert_close(haverkamp.theta(h=h), expected)


def test_s_haverkamp(haverkamp: pe.SoilModel, h: FloatArray = h) -> None:
    """Test degree of saturation calculation for the Haverkamp model."""
    expected = array(
        [
            0.8339910845345532,
            0.24068624139198183,
            0.019607843137254898,
            0.0012603242672546866,
            7.961509504265527e-05,
            5.023747624852072e-06,
            3.169785380148769e-07,
            1.999999960205131e-08,
            1.2619146881543392e-09,
        ]
    )
    assert_close(haverkamp.s(h=h), expected)


def test_k_haverkamp(haverkamp: pe.SoilModel, h: FloatArray = h) -> None:
    """Test hydraulic conductivity calculation for the Haverkamp model."""
    expected = array(
        [
            9.92100751538024,
            8.879484773404124,
            3.333333333333333,
            0.30583037614054576,
            0.019865814910850377,
            0.0012557854962273476,
            7.924403165642413e-05,
            4.999997500001253e-06,
            3.1547866228741785e-07,
        ]
    )
    assert_close(haverkamp.k(h=h), expected)


def test_k_r_haverkamp_accepts_s_kwarg(
    haverkamp: pe.soilmodel.Haverkamp,
) -> None:
    """Test that Haverkamp k_r method accepts a pre-computed saturation argument."""
    kr_out = haverkamp.k_r(h=array([1.0, 2.0, 3.0]), s=array([0.1, 0.2, 0.3]))
    expected = array([0.7352941176470588, 0.8620689655172414, 0.9146341463414634])
    assert_close(kr_out, expected)


def test_h_haverkamp_inverse(
    haverkamp: pe.soilmodel.Haverkamp, theta: FloatArray = theta
) -> None:
    """Test that Haverkamp h is the inverse of theta."""
    expected = array(
        [
            0.11334904230389142,
            0.04501288603422439,
            0.019670381895844538,
            0.0045279908036708555,
        ]
    )
    h_out = haverkamp.h(theta=theta)
    assert_close(h_out, expected)
    assert_close(haverkamp.theta(h=h_out), theta)


@pytest.mark.parametrize(
    "fixture_name",
    [
        "genuchten",
        "brooks",
        "panday",
        "gardner",
        "rucker",
        "genuchtengardner",
        "haverkamp",
    ],
)
def test_h_theta_roundtrip(fixture_name: str, request: pytest.FixtureRequest) -> None:
    """Test that h and theta are inverse pairs for models with closed-form inverses."""
    model = request.getfixturevalue(fixture_name)
    h_out = model.h(theta=theta)
    theta_out = model.theta(h=h_out)
    assert_close(theta, theta_out)
