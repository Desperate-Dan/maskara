from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from maskara import __version__, _program

setup(
    name='maskara',
    version=__version__,    
    description='maskara',
    url='https://github.com/Desperate-Dan/maskara',
    author='Daniel Maloney',
    author_email='dmaloney@ed.ac.uk',
    packages=find_packages(),
    scripts=['maskara/maskara.py'],
    install_requires=['biopython',
                        "pysam",
                        "numpy"
                      ],
    entry_points="""
    [console_scripts]
    {program} = maskara:main
    """.format(program = _program)
)
