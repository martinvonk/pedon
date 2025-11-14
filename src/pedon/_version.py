from importlib import metadata
from platform import python_version

__version__ = "0.0.10"


def get_versions() -> dict[str, str]:
    """Get version information.

    Returns:
        A dictionary with version information.
    """
    version_d = {
        "python": python_version(),
        "pedon": __version__,
    }
    requirements = metadata.requires("pedon")
    if requirements:
        deps = [x for x in requirements if "extra" not in x]
        for dep in deps:
            version_d[dep] = metadata.version(dep)
    return version_d


def show_versions() -> None:
    """Print version information."""
    versions = get_versions()
    max_name_len = max(len(name) for name in versions.keys())
    msg = ""
    for name, version in versions.items():
        msg += f"{name:<{max_name_len}} : {version}\n"

    print(msg.strip())
