#!/usr/bin/env python

from setuptools import setup, find_packages


# Sets the __version__ variable
__version__ = None
exec(open('DispaSET/_version.py').read())

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
    package_data={'DispaSET': ['config/*.yaml','GAMS/*']},
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