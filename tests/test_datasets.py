import pedon as pe


def test_sample_staring_2018():
    pe.soil.SoilSample().from_staring("B01", year="2018")


def test_sample_staring_2001():
    pe.soil.SoilSample().from_staring("B02", year="2001")


def test_soil_from_name():
    pe.soil.Soil("VS2D_Del Monte Sand").from_name(pe.Brooks)


def test_soil_from_staring():
    pe.soil.Soil("O01").from_staring()
