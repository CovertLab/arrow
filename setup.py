import os
import glob
import setuptools
from distutils.core import setup, Extension
import numpy.distutils.misc_util

headers = glob.glob('arrow/*.h')
sources = glob.glob('arrow/*.c')
obsidian = Extension('obsidian', sources=sources)

with open("README.md", 'r') as readme:
	long_description = readme.read()

current_dir = os.getcwd()
arrow_dir = os.path.join(current_dir, 'arrow')
include = [arrow_dir] + numpy.distutils.misc_util.get_numpy_include_dirs()

setup(
	name='stochastic-arrow',
	version='0.1.7',
	packages=['arrow'],
	author='Ryan Spangler',
	author_email='spanglry@stanford.edu',
	url='https://github.com/CovertLab/arrow',
	license='MIT',
	include_dirs=include,
	ext_modules=[obsidian],
	long_description=long_description,
	long_description_content_type='text/markdown')
