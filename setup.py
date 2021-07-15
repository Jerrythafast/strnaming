#!/usr/bin/env python3

#
# Copyright (C) 2021 Jerry Hoogenboom
#
# This file is part of STRNaming, an algorithm for generating simple,
# informative names for sequenced STR alleles in a standardised and
# automated manner.
#
# STRNaming is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# STRNaming is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with STRNaming.  If not, see <http://www.gnu.org/licenses/>.
#

import setuptools

with open("README.md", "rt", encoding="UTF-8") as fh:
    long_description = fh.read()

version = {}
with open("strnaming/__init__.py", "rt", encoding="UTF-8") as fh:
    exec(fh.read(), version)

setuptools.setup(
    name="strnaming",
    version=version['__version__'],
    description="STRNaming STR Sequence Nomenclature",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://fdstools.nl",
    project_urls={
        "Bug Tracker": "https://github.com/Jerrythafast/STRNaming/issues",
        "Source Code": "https://github.com/Jerrythafast/STRNaming",
    },
    author="Jerry Hoogenboom",
    author_email="jerryhoogenboom@outlook.com",
    license="LGPLv3+",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Legal Industry",
        "Topic :: Scientific/Engineering :: Bio-Informatics"],
    keywords='bioinformatics forensics NGS MPS DNA sequencing STR nomenclature',
    packages=setuptools.find_packages(),
    install_requires=[],
    python_requires='>=3.5',
    package_data={
        "strnaming": [
            "data/reference-structures-*.bin.gz",
            "data/refseq-*-*.txt.gz",
            "data/ranges_*.txt"]
    },
    entry_points={
        'console_scripts': [
            "strnaming=strnaming.strnaming:main"
        ]
    }
)
