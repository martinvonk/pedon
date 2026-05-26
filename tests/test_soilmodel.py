"""Test soilmodels."""

from typing import get_args

import pytest
from numpy import allclose, array, logspace

import pedon as pe
from pedon._typing import FloatArray, SoilModelNames

h = -logspace(-2, 6, num=9, dtype=float)

theta = array([0.1, 0.2, 0.3, 0.4])


def assert_close(actual: FloatArray, expected: FloatArray) -> None:
    """Assert if two arrays are close within a reasonable tolerance."""
    assert allclose(actual, expected, rtol=1e-8, atol=1e-8)


@pytest.mark.parametrize("soilmodel_name", get_args(SoilModelNames))
def test_get_soilmodel(
    soilmodel_name: SoilModelNames,
) -> None:
    """Test that get_soilmodel maps each supported name to the right class."""
    smt = pe.soilmodel.get_soilmodel(soilmodel_name)
    assert smt.__name__ == soilmodel_name


def test_theta_genuchten(genuchten: pe.SoilModel) -> None:
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


def test_s_genuchten(genuchten: pe.SoilModel) -> None:
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


def test_k_genuchten(genuchten: pe.SoilModel) -> None:
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


def test_theta_brooks(brooks: pe.SoilModel) -> None:
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


def test_s_brooks(brooks: pe.SoilModel) -> None:
    """Test degree of saturation calculation for the Brooks-Corey model."""
    expected = array([1.0, 1.0, 1.0, 1.0, 0.01, 0.0001, 1e-06, 1e-08, 1e-10])
    assert_close(brooks.s(h=h), expected)


def test_k_brooks(brooks: pe.SoilModel) -> None:
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


def test_theta_panday(panday: pe.SoilModel) -> None:
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


def test_s_panday(panday: pe.SoilModel) -> None:
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


def test_k_panday(panday: pe.SoilModel) -> None:
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


def test_theta_gardner(gardner: pe.SoilModel) -> None:
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


def test_s_gardner(gardner: pe.SoilModel) -> None:
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


def test_k_gardner(gardner: pe.SoilModel) -> None:
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


def test_theta_rucker(rucker: pe.SoilModel) -> None:
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


def test_s_rucker(rucker: pe.SoilModel) -> None:
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


def test_k_rucker(rucker: pe.SoilModel) -> None:
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


def test_theta_genuchtengardner(genuchtengardner: pe.SoilModel) -> None:
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


def test_s_genuchtengardner(genuchtengardner: pe.SoilModel) -> None:
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


def test_k_genuchtengardner(genuchtengardner: pe.SoilModel) -> None:
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


def test_theta_fredlund(fredlund: pe.SoilModel) -> None:
    """Test water content calculation for the Fredlund model."""
    expected = array(
        [
            0.4292739372053212,
            0.4127503403050247,
            0.2432576190859083,
            0.08790598851970158,
            0.0474319313818822,
            0.03189870288483172,
            0.02384298091305578,
            0.01894497637562149,
            0.01566526178712822,
        ]
    )
    assert_close(fredlund.theta(h=h), expected)


def test_s_fredlund(fredlund: pe.SoilModel) -> None:
    """Test degree of saturation calculation for the Fredlund model."""
    expected = array(
        [
            0.99831148187284,
            0.9598845123372668,
            0.5657153932230425,
            0.20443253144116647,
            0.1103068171671679,
            0.07418302996472494,
            0.05544879282105995,
            0.04405808459446858,
            0.03643084136541447,
        ]
    )
    assert_close(fredlund.s(h=h), expected)


def test_k_fredlund(fredlund: pe.SoilModel) -> None:
    """Test hydraulic conductivity calculation for the Fredlund model."""
    expected = array(
        [
            10.154616004536054,
            5.218778175589122,
            0.12213707500185565,
            8.003798775846895e-05,
            1.0251971820534869e-07,
            2.5630136211214665e-10,
            9.079007831684122e-13,
            3.816854013113528e-15,
            0.0,
        ]
    )
    assert_close(fredlund.k(h=h), expected)


def test_h_fredlund(fredlund: pe.SoilModel, theta: FloatArray = theta) -> None:
    """Test pressure head calculation for the Fredlund model."""
    expected = array(
        [7.031778582692749, 1.5128170925180793, 0.5981696228905503, 0.15333981041674402]
    )
    h_out = fredlund.h(theta=theta)
    assert_close(h_out, expected)


def test_k_r_fredlund_rejects_s_kwarg(fredlund: pe.SoilModel) -> None:
    """Test that Fredlund k_r rejects saturation input."""
    with pytest.raises(NotImplementedError, match="using the pressure head"):
        fredlund.k_r(h=array([1.0, 2.0, 3.0]), s=array([0.1, 0.2, 0.3]))


def test_theta_kosugi(kosugi: pe.SoilModel) -> None:
    """Test water content calculation for the Kosugi model."""
    h_abs = abs(h)
    expected = array(
        [
            0.43,
            0.42999999819621027,
            0.42997391351309045,
            0.4184485298896685,
            0.22,
            0.021551470110331506,
            0.010026086486909508,
            0.010000001803789757,
            0.010000000000003467,
        ]
    )
    assert_close(kosugi.theta(h=h_abs), expected)


def test_s_kosugi(kosugi: pe.SoilModel) -> None:
    """Test degree of saturation calculation for the Kosugi model."""
    h_abs = abs(h)
    expected = array(
        [
            1.0,
            0.9999999957052626,
            0.999937889316882,
            0.972496499737306,
            0.5,
            0.02750350026269398,
            6.211068311787629e-05,
            4.294737515113387e-09,
            8.254664573734071e-15,
        ]
    )
    assert_close(kosugi.s(h=h_abs), expected)


def test_k_kosugi(kosugi: pe.SoilModel) -> None:
    """Test hydraulic conductivity calculation for the Kosugi model."""
    h_abs = abs(h)
    expected = array(
        [
            10.0,
            9.99994795699993,
            9.916381557070964,
            5.754238967213946,
            0.09362821398232724,
            1.3669488922238534e-06,
            4.376455940449844e-15,
            1.9946478748160554e-27,
            1.1051085183358813e-43,
        ]
    )
    assert_close(kosugi.k(h=h_abs), expected)


def test_h_kosugi(kosugi: pe.SoilModel, theta: FloatArray = theta) -> None:
    """Test pressure head calculation for the Kosugi model."""
    expected = array(
        [258.5622598127402, 115.43965450658366, 55.06583675346047, 17.234094056734545]
    )
    h_out = kosugi.h(theta=theta)
    assert_close(h_out, expected)
    assert_close(kosugi.theta(h=h_out), theta)


def test_k_r_kosugi_accepts_s_kwarg(kosugi: pe.SoilModel) -> None:
    """Test that Kosugi k_r method accepts a pre-computed saturation argument."""
    kr_out = kosugi.k_r(h=array([1.0, 2.0, 3.0]), s=array([0.1, 0.2, 0.3]))
    expected = array([1.3528004311960868e-05, 0.0001896793173053034, 0.0009808586637367542])
    assert_close(kr_out, expected)


def test_alpha_w_genuchtenkool(genuchtenkool: pe.soilmodel.GenuchtenKool) -> None:
    """Test scaled alpha parameter for the GenuchtenKool model."""
    assert genuchtenkool.alpha_w == 0.05


def test_theta_genuchtenkool(genuchtenkool: pe.SoilModel) -> None:
    """Test water content calculation for the GenuchtenKool model."""
    expected = array(
        [
            0.4299910737668094,
            0.42988779112842795,
            0.4286129877445311,
            0.41563208021830816,
            0.3624881596094997,
            0.2936755393053078,
            0.23558483318061144,
            0.18920450408643255,
            0.1523482139298137,
        ]
    )
    assert_close(genuchtenkool.theta(h=h), expected)


def test_s_genuchtenkool(genuchtenkool: pe.SoilModel) -> None:
    """Test degree of saturation calculation for the GenuchtenKool model."""
    expected = array(
        [
            0.9999787470638319,
            0.9997328360200666,
            0.9966975898679312,
            0.965790667186448,
            0.8392575228797612,
            0.6754179507269233,
            0.5371067456681224,
            0.42667739068198224,
            0.33892431888050883,
        ]
    )
    assert_close(genuchtenkool.s(h=h), expected)


def test_k_genuchtenkool(genuchtenkool: pe.SoilModel) -> None:
    """Test hydraulic conductivity calculation for the GenuchtenKool model."""
    expected = array(
        [
            2.8343153705439206,
            1.6927117058050811,
            0.6817168120309531,
            0.09609620636559924,
            0.0018452194138960125,
            1.2243219408738212e-05,
            6.982358918415397e-08,
            3.9308819964993456e-10,
            2.2106947229727662e-12,
        ]
    )
    assert_close(genuchtenkool.k(h=h), expected)


def test_h_genuchtenkool(
    genuchtenkool: pe.SoilModel, theta: FloatArray = theta
) -> None:
    """Test pressure head calculation for the GenuchtenKool model."""
    expected = array(
        [97971055.72421077, 55708.66779383144, 799.4116130609773, 24.66894112455789]
    )
    assert_close(genuchtenkool.h(theta=theta), expected)


def test_theta_haverkamp(haverkamp: pe.SoilModel) -> None:
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


def test_s_haverkamp(haverkamp: pe.SoilModel) -> None:
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


def test_k_haverkamp(haverkamp: pe.SoilModel) -> None:
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
        "fredlund",
        "kosugi",
        "genuchtenkool",
        "haverkamp",
    ],
)
def test_h_theta_roundtrip(fixture_name: str, request: pytest.FixtureRequest) -> None:
    """Test that h and theta are inverse pairs for models with closed-form inverses."""
    model = request.getfixturevalue(fixture_name)
    h_out = model.h(theta=theta)
    theta_out = model.theta(h=h_out)
    assert_close(theta, theta_out)
