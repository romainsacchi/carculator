import os
from pathlib import Path

from setuptools import setup

packages = []
root_dir = os.path.dirname(__file__)
if root_dir:
    os.chdir(root_dir)

# read the contents of your README file
this_directory = Path(__file__).parent
README = (this_directory / "README.md").read_text()

# Probably should be changed, __init__.py is no longer required for Python 3
for dirpath, dirnames, filenames in os.walk("carculator_utils"):
    # Ignore dirnames that start with '.'
    if "__init__.py" in filenames:
        pkg = dirpath.replace(os.path.sep, ".")
        if os.path.altsep:
            pkg = pkg.replace(os.path.altsep, ".")
        packages.append(pkg)


def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join("..", path, filename))
    return paths


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="carculator_utils",
    version="1.0.0",
    packages=packages,
    author="Romain Sacchi <romain.sacchi@psi.ch>",
    python_requires=">=3.9",
    license=open("LICENSE").read(),
    package_data={
        "carculator_utils": package_files(os.path.join("carculator_utils", "data"))
    },
    install_requires=[
        "pandas",
        "xarray",
        "numpy<1.24.0",
        "klausen",
        "xlrd",
        "numexpr",
        "bw2io",
        "pycountry",
        "wurst",
        "pyyaml",
    ],
    url="https://github.com/romainsacchi/carculator_utils",
    description="Provides convenience functions for carculator_utils",
    long_description_content_type="text/markdown",
    long_description=README,
    classifiers=[
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
)
