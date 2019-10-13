import os
# from glob import glob
from distutils.core import setup, Extension
from Cython.Build import cythonize
# from Cython.Distutils import build_ext
import numpy.distutils.misc_util

obsidian = [Extension(
	'obsidian',
	sources=['arrow/mersenne.c', 'arrow/obsidian.c', 'arrow/obsidianmodule.c'])]

with open("README.md", 'r') as readme:
	long_description = readme.read()

current_dir = os.getcwd()
arrow_dir = os.path.join(current_dir, 'arrow')
include = [arrow_dir] + numpy.distutils.misc_util.get_numpy_include_dirs()

arrowhead = cythonize([
	Extension('arrow.arrowhead',
			  sources=['arrow/arrowhead.pyx'],
			  include_dirs=['arrow'],
			  )],
	# annotate=True,  # to get an HTML code listing
	)

arrow = [Extension('arrow.arrowhead',
		   sources=['arrow/mersenne.c', 'arrow/obsidian.c', 'arrow/arrowhead.c'],
		   include_dirs=['arrow'])]

setup(
	name='stochastic-arrow',
	version='0.2.0',
	packages=['arrow'],
	author='Ryan Spangler',
	author_email='spanglry@stanford.edu',
	url='https://github.com/CovertLab/arrow',
	license='MIT',
	include_dirs=include,
	ext_modules=arrow + obsidian,
	long_description=long_description,
	long_description_content_type='text/markdown',
	# cmdclass={'build_ext': build_ext},
	classifiers=[
		'Development Status :: 3 - Alpha',
		'License :: OSI Approved :: MIT License',
		'Programming Language :: Python',
		'Programming Language :: Python :: 2.7',
		# 'Programming Language :: Python :: 3',
		'Topic :: Scientific/Engineering',
	])
