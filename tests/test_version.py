import re

from pedon._version import __version__, get_versions, show_versions


def test_get_versions_returns_dict() -> None:
    """Test that get_versions returns a dictionary."""
    versions = get_versions()
    assert isinstance(versions, dict)


def test_get_versions_contains_python() -> None:
    """Test that get_versions includes python version."""
    versions = get_versions()
    assert "python" in versions
    assert isinstance(versions["python"], str)
    # Check it's a valid version string like "3.13.5"
    assert re.match(r"\d+\.\d+\.\d+", versions["python"])


def test_get_versions_contains_pedon() -> None:
    """Test that get_versions includes pedon version."""
    versions = get_versions()
    assert "pedon" in versions
    assert versions["pedon"] == __version__
    assert isinstance(versions["pedon"], str)


def test_get_versions_contains_dependencies() -> None:
    """Test that get_versions includes dependency versions."""
    versions = get_versions()
    # Should contain at least numpy, pandas, scipy (core dependencies)
    # Check at least one common dependency is present
    has_deps = any(dep in versions for dep in ["numpy", "pandas", "scipy"])
    assert has_deps, f"Expected dependencies in versions, got: {versions.keys()}"


def test_get_versions_all_values_are_strings() -> None:
    """Test that all version values are strings."""
    versions = get_versions()
    for name, version in versions.items():
        assert isinstance(version, str), (
            f"Version for {name} is not a string: {type(version)}"
        )


def test_show_versions_runs_without_error(capsys) -> None:
    """Test that show_versions runs and prints output."""
    show_versions()
    captured = capsys.readouterr()
    assert len(captured.out) > 0
    assert "python" in captured.out
    assert "pedon" in captured.out


def test_show_versions_output_format(capsys) -> None:
    """Test that show_versions output is properly formatted."""
    show_versions()
    captured = capsys.readouterr()
    lines = captured.out.strip().split("\n")

    # Should have at least a few lines (python, pedon, dependencies)
    assert len(lines) >= 2

    # Each line should contain a colon separator
    for line in lines:
        assert ":" in line, f"Expected ':' separator in line: {line}"

    # Check version numbers are present
    assert re.search(r"\d+\.\d+", captured.out), "Expected version numbers in output"


def test_show_versions_alignment(capsys) -> None:
    """Test that show_versions output is properly aligned."""
    show_versions()
    captured = capsys.readouterr()
    lines = captured.out.strip().split("\n")

    # Extract the position of colons (should be aligned)
    colon_positions = [line.index(":") for line in lines if ":" in line]

    # All colons should be at the same position
    assert len(set(colon_positions)) == 1, "Colons should be aligned"
