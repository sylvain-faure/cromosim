from os import path
from setuptools import setup

CLASSIFIERS = [
    "Development Status :: 5 - Production/Stable",
    "Intended Audience :: Education",
    "Intended Audience :: Science/Research",
    "Intended Audience :: Developers",
    "License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)",
    "Programming Language :: Python :: 3.5",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Topic :: Software Development",
    "Topic :: Scientific/Engineering",
    "Operating System :: Microsoft :: Windows",
    "Operating System :: Unix",
    "Operating System :: MacOS",
    "Topic :: Scientific/Engineering :: Mathematics",
    "Topic :: Software Development :: Libraries :: Python Modules"]

MAJOR = "2"
MINOR = "0"
PATCH = "2"
VERSION = "{0}.{1}.{2}".format(MAJOR, MINOR, PATCH)


def write_version_py(filename='cromosim/version.py'):
    a = open(filename, 'w')
    try:
        a.write("version = '{}'".format(VERSION))
    finally:
        a.close()


# read the contents of your README file
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

# README = open("README.rst").readlines()

write_version_py()

setup(
    name="cromosim",
    version=VERSION,
    # description=README[0],
    # long_description_content_type='text/x-rst',
    # long_description="".join(README[1:]),
    long_description=long_description,
    long_description_content_type='text/markdown',
    author="Sylvain Faure, Bertrand Maury",
    author_email="sylvain.faure@math.u-psud.fr, bertrand.maury@math.u-psud.fr",
    url="http://www.cromosim.fr",
    license="GPL",
    keywords="Crowd Motion Simulator",
    classifiers=CLASSIFIERS,
    # packages=find_packages(exclude=['doc','examples']),
    packages=["cromosim"],
    include_package_data=True,
    install_requires=[
                      'numpy',
                      'scipy',
                      'Pillow',
                      'matplotlib',
                      'numpydoc',
                      'sphinx',
                      'scikit-fmm',
                      'cvxopt',
                      'imageio'
                      ],
)
