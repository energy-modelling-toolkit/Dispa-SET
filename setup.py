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

# Retrieve the release tag as fallback version
exec(open(os.path.join(HERE, 'dispaset/_release.py')).read().strip())

def local_version(version):
    return version.format_choice("+{node}", "+dirty")

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
        'local_scheme': local_version,
        'write_to': os.path.join(HERE, 'dispaset/_version.py'),
        'fallback_version': release,
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
        "optcc"
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
