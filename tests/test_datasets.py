"""Tests for pedon datasets."""

import pedon as pe

pe.soil.Soil("Del Monte Sand").from_name(pe.Brooks)


def test_sample_staring_2018() -> None:
    """Test loading soil sample from Staring 2018 dataset."""
    ss = pe.soil.SoilSample().from_staring("B01", year="2018")
    assert ss.silt_p == 5.0
    assert ss.clay_p == 0.0
    assert ss.om_p == 7.5
    assert ss.m50 == 0.01575


def test_sample_staring_2001() -> None:
    """Test loading soil sample from Staring 2001 dataset."""
    ss = pe.soil.SoilSample().from_staring("B01", year="2001")
    assert ss.silt_p == 7.0
    assert ss.clay_p == 0.0
    assert ss.rho == 1.55
    assert ss.om_p == 2.5
    assert ss.m50 == 0.0155


def test_soil_from_name() -> None:
    """Test loading soil from name."""
    s = pe.soil.Soil("Del Monte Sand").from_name(pe.Brooks)
    assert s.source == "VS2D"
    assert s.description == "20 mesh"
    assert isinstance(s.model, pe.Brooks)
    assert s.model.k_s == 7e5
    assert s.model.theta_r == 0.011
    assert s.model.theta_s == 0.36
    assert s.model.h_b == 0.112
    assert s.model.l == 2.5


def test_soil_from_staring() -> None:
    """Test loading soil from Staring dataset."""
    s = pe.soil.Soil("O01").from_staring()
    assert s.source == "Staring_2018"
    assert s.description == "leemarm zeer fijn tot matig fijn zand"
    assert isinstance(s.model, pe.Genuchten)
    assert s.model.k_s == 22.32
    assert s.model.theta_r == 0.01
    assert s.model.theta_s == 0.366
    assert s.model.alpha == 0.016
    assert s.model.n == 2.16
    assert s.model.l == 2.868


def test_soil_list_genuchten() -> None:
    """Test listing soil names for Genuchten model."""
    gen_soils = pe.soil.Soil.list_names(pe.Genuchten)
    gen_soils_current = [
        "Sand",
        "Loamy Sand",
        "Sandy Loam",
        "Loam",
        "Silt",
        "Silt Loam",
        "Sandy Clay Loam",
        "Clay Loam",
        "Silty Clay Loam",
        "Sandy Clay",
        "Silty Clay",
        "Clay",
        "B01",
        "B02",
        "B03",
        "B04",
        "B05",
        "B06",
        "B07",
        "B08",
        "B09",
        "B10",
        "B11",
        "B12",
        "B13",
        "B14",
        "B15",
        "B16",
        "B17",
        "B18",
        "O01",
        "O02",
        "O03",
        "O04",
        "O05",
        "O06",
        "O07",
        "O08",
        "O09",
        "O10",
        "O11",
        "O12",
        "O13",
        "O14",
        "O15",
        "O16",
        "O17",
        "O18",
        "Medium Sand",
        "Del Monte Sand",
        "Fresno Medium Sand",
        "Unconsolidated Sand",
        "Fine Sand",
        "Columbia Sandy Loam",
        "Touchet Silt Loam",
        "Hygiene Sandstone",
        "Adelanto Loam",
        "Limon Silt",
        "Yolo Light Clay",
    ]
    assert all([s in gen_soils for s in gen_soils_current])
