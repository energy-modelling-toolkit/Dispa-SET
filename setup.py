#!/usr/bin/env python

from setuptools import setup, find_packages
import codecs
import os
HERE = os.path.abspath(os.path.dirname(__file__))

def read(*parts):
    """
    Build an absolute path from *parts* and and return the contents of the
    resulting file.  Assume UTF-8 encoding.
    """
    with codecs.open(os.path.join(HERE, *parts), "rb", "utf-8") as f:
        return f.read()

# Sets the __version__ variable
__version__ = None
exec(open('DispaSET/_version.py').read())

setup(
    name='dispaset',
    version=__version__,
    author='Sylvain Quoilin, Konstantinos Kavvadias',
    author_email='squoilin@uliege.be',
    description='An open-source unit commitment and optimal dispatch model ',
    license='EUPL v1.1.',
    url='http://www.dispaset.eu',
    download_url='http://www.dispaset.eu/en/latest/releases.html',
    packages=find_packages(),
    long_description=read('README.md'),
    include_package_data=True,
    install_requires=[
        "future >= 0.15",
        "click >= 3.3",
        "numpy >= 1.12",
        "pandas >= 0.19",
        "xlrd >= 0.9",
        "matplotlib >= 1.5.1",
        "gdxcc >= 7",
        "gamsxcc",
        "optcc"
    ],
    extras_require={'pyomo': ['pyomo>=5.2'],
                    },
    entry_points={
        'console_scripts': [
            'dispaset = DispaSET.cli:cli'
        ]},
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: EUPL Software License',
        'Programming Language :: Python'
    ],
    keywords=['energy systems', 'optimization', 'mathematical programming']
)