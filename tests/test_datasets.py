"""Tests for pedon datasets."""

import pytest

import pedon as pe


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


def test_soil_from_clapp() -> None:
    """Test loading soil from the Clapp database."""
    s = pe.soil.Soil("Sand").from_name(pe.Campbell, source="Clapp")
    assert s.source == "Clapp"
    assert s.description == "Sand"
    assert isinstance(s.model, pe.Campbell)
    assert s.model.k_s == 1520.64
    assert s.model.theta_s == 0.395
    assert s.model.h_b == 3.5
    assert s.model.b == 3.0


def test_soil_from_rawls() -> None:
    """Test loading soil from the Rawls PTF."""
    s = pe.soil.Soil("Sand").from_name(pe.Brooks, source="Rawls")
    assert s.source == "Rawls"
    assert s.description == "Sand"
    assert isinstance(s.model, pe.Brooks)
    assert s.model.k_s == 504.0
    assert s.model.theta_r == 0.02
    assert s.model.theta_s == 0.437
    assert s.model.h_b == 7.26
    assert s.model.l == 0.592


def test_soil_from_staring_with_int_year() -> None:
    """Test integer year input for from_staring."""
    s = pe.soil.Soil("O01").from_staring(2018)
    assert s.source == "Staring_2018"


def test_soil_from_staring_invalid_year() -> None:
    """Invalid Soil.from_staring year should raise a clear ValueError."""
    with pytest.raises(ValueError, match="Year must either be '2001' or '2018'"):
        pe.soil.Soil("O01").from_staring(1999)


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


def test_soilsample_from_staring_invalid_year() -> None:
    """Invalid Staring year should raise a clear ValueError."""
    with pytest.raises(ValueError, match="No Staring series available"):
        pe.soil.SoilSample().from_staring("B01", year="1999")


def test_soil_from_name_requires_source_when_multiple() -> None:
    """Names with multiple sources should require an explicit source."""
    with pytest.raises(ValueError, match="Please provide the source"):
        pe.soil.Soil("O01").from_name(pe.Genuchten)


def test_soil_from_name_hydrus_prefix_sets_source(caplog) -> None:
    """HYDRUS_ prefix should be stripped and source set automatically."""
    with caplog.at_level("WARNING"):
        s = pe.soil.Soil("HYDRUS_Del Monte Sand").from_name(pe.Brooks)

    assert s.name == "Del Monte Sand"
    assert s.source is not None
    assert "Removed 'HYDRUS_' from soil name" in caplog.text


def test_soil_list_names_invalid_type_raises_typeerror() -> None:
    """Invalid soil model argument types should raise TypeError."""
    with pytest.raises(TypeError, match=r"Type\[SoilModel\] \| SoilModel \| str"):
        pe.soil.Soil.list_names(123)  # type: ignore[arg-type]
