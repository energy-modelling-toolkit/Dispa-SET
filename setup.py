#!/usr/bin/env python

from setuptools import setup, find_packages
import os


# Sets the __version__ variable
__version__ = None
exec(open('DispaSET/_version.py').read())

def recursive_listdir(folder):
    '''
    Generates a list of all files in a a specified directory
    '''
    return [(dp,[os.path.join(dp, f)]) for dp, dn, fn in os.walk(os.path.expanduser(folder)) for f in fn]

def define_data_files():
    list_files = []
    for folder in ['ConfigFile','Database','Externals']:
        list_files += recursive_listdir(folder)
    return list_files
        

setup(
    name='Dispa-SET',
    version=__version__,
    author='Sylvain Quoilin, Konstantinos Kavvadias',
    author_email='squoilin@uliege.be',
    description='An open-source unit commitment and optimal dispatch model ',
    license='EUPL v1.1.',
    url='http://www.dispaset.eu',
    download_url='http://www.dispaset.eu/en/latest/releases.html',
    packages=find_packages(),
    package_data={'DispaSET': ['config/*.yaml','GAMS/*']},
    data_files=define_data_files(),
    install_requires=[
        "click >= 3.3",
        "numpy >= 1.12",
        "pandas >= 0.19"
    ],
    scripts=['dispaset'],
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: EUPL Software License',
        'Programming Language :: Python'
    ],
    keywords=['energy systems', 'optimization', 'mathematical programming']
)