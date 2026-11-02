"""
Local data locations, read from paths.toml. See paths.example.toml.

The file is found through the LU_PATHS environment variable, or next to the kits.
A location that is not set gives a placeholder path, so importing always works
and the error appears only when the data are opened.
"""
import os
import tomllib

_DEFAULT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "paths.toml")


def get(name):
    r"""
    get: the location called `name` in paths.toml.

    Data, meshes and outputs live in different places on each machine. Keeping their locations in one
    file, outside the code, lets the same code run anywhere.

    The file is the one named by the LU_PATHS environment variable, or py-kits/paths.toml. An unknown name
    returns a placeholder that names it, rather than failing.

        Parameters:
            name (str): Key in paths.toml, e.g. 'glorys04'

        Returns:
            path (str): The location, with ~ expanded, or '<name: not set in paths.toml>' if absent

        Formula:
            \text{path} = \text{paths.toml}[\text{name}],\ \text{with } \sim \text{ expanded}

        See also:
            datakit.source
    """
    path = os.environ.get("LU_PATHS", _DEFAULT)
    table = {}
    if os.path.exists(path):
        with open(path, "rb") as f:
            table = tomllib.load(f)
    return os.path.expanduser(table.get(name, f"<{name}: not set in paths.toml>"))
