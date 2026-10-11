"""Package version; readable by the build backend without runtime imports."""

VERSION = "1.9.7"
__version__ = tuple(
    int(part) if part.isdigit() else part for part in VERSION.split(".")
)
