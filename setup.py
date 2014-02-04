# databrowser's setup.py
from distutils.core import setup

files = ["favicon.ico"]

setup(
    name = "offspec_databrowser",
    packages = ["offspec_databrowser"],
    package_data = {"offspec_databrowser" : files },
    scripts=['postinstall.py'],
    version = "0.9.7.1",
    description = "Offspecular neutron reflectometry reduction (ICP-format files)",
    author = "Brian Ben Maranville",
    author_email = "brian.maranville@nist.gov",
    url = "http://code.google.com/p/offspec-databrowser/",
    download_url = "http://code.google.com/p/offspec-databrowser/downloads/databrowser-0.9.4.tgz",
    keywords = ["neutron", "offspecular", "reflectometry"],
    classifiers = [
        "Programming Language :: Python",
        "Development Status :: 3 - Alpha",
        "Environment :: Other Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization",
        ],
    long_description = """\
Offspecular neutron scattering data reduction tool
--------------------------------------------------

Imports: 2-d ICP-style data files (must be prebinned in one axis if using PSD)

Features:
- graphical editor for He3 cell parameters
- manual markup of data files with extra information:
  - polarization state ('++', '+-', '-+', '--' or None)
  - total elapsed time for measurement (needed to calculate time-dependent polarization)
  - extra comments per file
- function to auto-calculate elapsed time from "last-modified" timestamp on
    unrebinned files in source directory
- subtraction of background from files: specify angular range to treat as "background",
    and then average that angular region and subtract from whole range (normalized to monit.)
- polarization correction, with automatic adjustment for incomplete input states:
  - if all 4 cross-sections present, does proper inversion of 4x4 NT matrix 
      outputs all 4 Reflectivity cross-sections
  - if one spin-flip cross-section is missing, sets R(missing) to R(present) for SF and
      outputs only 3 cross sections, corresponding to input cross-sections
  - if both spin-flip cross sections missing, assume R(SF) = 0 and reduce NT matrix to 2x2
      outputs only 2 cross sections (non-spin-flip), as the SF graphs would be all zeros!
- advanced plotting:
  - right click on 2d colormaps to copy color range or slice region from another plot
  - right click on line graphs to pop-out to separate graph window
  - click and drag legends around on graph

This version requires Python 2.6 or later; a forward-compatible (Python 3) is underway.
"""
)