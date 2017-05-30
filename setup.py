#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys

try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup

setup(
    name='pcmdpy',
    version='0.0.1',
    author='Ben Cook',
    author_email='bcook@cfa.harvard.edu',
    packages=['pcmdpy'],
    url='',
    license='LICENSE',
    description='Tools for modelling crowded-field photometry using the Pixel Color-Magnitude Diagram technique',
    package_data={'pcmdpy':['isoc_csv/*csv']},
    include_package_data=True,
    #install_requires=['numpy','scipy','pandas','matplotlib','emcee'],
)
