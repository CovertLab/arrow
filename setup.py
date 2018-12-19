import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
	long_description = readme.read()

setup(
	name='stochastic-arrow',
	version='0.0.3',
	packages=['arrow'],
	author='Ryan Spangler',
	author_email='spanglry@stanford.edu',
	url='https://github.com/CovertLab/arrow',
	license='MIT',
	long_description=long_description,
	long_description_content_type='text/markdown')
