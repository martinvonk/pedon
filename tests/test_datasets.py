import pedon as pe


def test_sample_staring_2018() -> None:
    pe.soil.SoilSample().from_staring("B01", year="2018")


def test_sample_staring_2001() -> None:
    pe.soil.SoilSample().from_staring("B02", year="2001")


def test_soil_from_name() -> None:
    pe.soil.Soil("Del Monte Sand").from_name(pe.Brooks)


def test_soil_from_staring() -> None:
    pe.soil.Soil("O01").from_staring()


def test_soil_list_genuchten() -> None:
    pe.soil.Soil.list_names(pe.Genuchten)
