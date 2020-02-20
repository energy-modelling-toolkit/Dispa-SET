#!/usr/bin/env python

from setuptools import setup, find_packages
import codecs
import os
HERE = os.path.abspath(os.path.dirname(__file__))

# FINAL_RELEASE is the last stable version of Dispa-SET
# A more precisely version try to be automatically determined from the git repository using setuptools_scm.
# If it's not possible (using git archive tarballs for example), FINAL_RELEASE will be used as fallback version.
# edited manually when a new release is out (git tag -a)
FINAL_RELEASE = open(os.path.join(HERE, 'VERSION')).read().strip()

def read(*parts):
    """
    Build an absolute path from *parts* and and return the contents of the
    resulting file.  Assume UTF-8 encoding.
    """
    with codecs.open(os.path.join(HERE, *parts), "rb", "utf-8") as f:
        return f.read()

setup(
    name='dispaset',
    author='Sylvain Quoilin, Konstantinos Kavvadias',
    author_email='squoilin@uliege.be',
    description='An open-source unit commitment and optimal dispatch model ',
    license='EUPL v1.2.',
    url='http://www.dispaset.eu',
    download_url='http://www.dispaset.eu/en/latest/releases.html',
    packages=find_packages(),
    long_description=read('README.md'),
    include_package_data=True,
    use_scm_version={
        'version_scheme': 'post-release',
        'local_scheme': lambda version: version.format_choice("" if version.exact else "+{node}", "+dirty"),
        'fallback_version': FINAL_RELEASE,
    },
    setup_requires=["setuptools_scm"],
    install_requires=[
        "future >= 0.15",
        "click >= 3.3",
        "numpy >= 1.12",
        "pandas >= 0.19",
        "xlrd >= 0.9",
        "matplotlib >= 1.5.1",
        "gdxcc >= 7",
        "gamsxcc",
        "optcc",
        "setuptools_scm",
    ],
    entry_points={
        'console_scripts': [
            'dispaset = dispaset.cli:cli'
        ]},
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: EUPL Software License',
        'Programming Language :: Python'
    ],
    keywords=['energy systems', 'optimization', 'mathematical programming']
)
