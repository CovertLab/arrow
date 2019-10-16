import os
# from glob import glob
import setuptools  # used indirectly for bdist_wheel cmd and long_description_content_type
from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext
import numpy.distutils.misc_util

with open("README.md", 'r') as readme:
	long_description = readme.read()

current_dir = os.getcwd()
arrow_dir = os.path.join(current_dir, 'arrow')
include = [arrow_dir] + numpy.distutils.misc_util.get_numpy_include_dirs()

arrowhead = cythonize([
	Extension('arrow.arrowhead',
			  sources=['arrow/arrowhead.pyx', 'arrow/mersenne.c', 'arrow/obsidian.c',],
			  include_dirs=['arrow'],
			  define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
			  )],
	include_path=['arrow'],
	# annotate=True,  # to get an HTML code listing
	)

setup(
	name='stochastic-arrow',
	version='0.2.0',
	packages=['arrow'],
	author='Ryan Spangler, John Mason, Jerry Morrison',
	author_email='spanglry@stanford.edu',
	url='https://github.com/CovertLab/arrow',
	license='MIT',
	include_dirs=include,
	ext_modules=arrowhead,
	long_description=long_description,
	long_description_content_type='text/markdown',
	cmdclass={'build_ext': build_ext},
	requires=['numpy (>=1.14)', 'six'],
	classifiers=[
		'Development Status :: 3 - Alpha',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python',
		'Programming Language :: Python :: 2.7',
		'Programming Language :: Python :: 3',
		'Topic :: Scientific/Engineering',
	])
