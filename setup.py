#!/usr/bin/env python

from pathlib import Path
from setuptools import setup, find_packages


# Sets the __version__ variable
__version__ = None
exec(open('DispaSET/_version.py').read())


def get_subdirs(path, glob_string):
    return [
        i.relative_to(str(path) + '/')
        for i in path.glob(glob_string) if i.is_dir()
    ]


def find_package_data():
    """Returns a list of found directories with package data files"""
    path = Path('./DispaSET')
    package_data = ['config/*.yaml','GAMS/*']
    for test_case_dir in get_subdirs(path, 'test/*'):
        package_data.append(str(test_case_dir) + '/*.csv')
    print(package_data)
    return package_data


setup(
    name='Dispa-SET',
    version=__version__,
    author='Sylvain Quoilin, Konstantino Kavvadias',
    author_email='squoilin@uliege.be',
    description='An open-source unit commitment and optimal dispatch model ',
    license='EUPL v1.1.',
    url='http://www.dispaset.eu',
    download_url='http://www.dispaset.eu/en/latest/releases.html',
    packages=find_packages(),
    package_data={'DispaSET': find_package_data()},
    install_requires=[
        "click >= 3.3",
        "numpy >= 1.12",
        "pandas >= 0.19, < 0.20"
    ],
    entry_points={
        'console_scripts': [
            'DispaSET = dispacli:cli'
        ]
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: EUPL Software License',
        'Programming Language :: Python'
    ],
    keywords=['energy systems', 'optimization', 'mathematical programming']
)