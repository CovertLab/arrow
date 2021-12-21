import os
import setuptools  # used indirectly for bdist_wheel cmd and long_description_content_type
from distutils.core import setup
from distutils.extension import Extension
import numpy.distutils.misc_util

with open("README.md", 'r') as readme:
    long_description = readme.read()

current_dir = os.getcwd()
arrow_dir = os.path.join(current_dir, 'arrow')
include = [arrow_dir] + numpy.distutils.misc_util.get_numpy_include_dirs()

# Compile the Cython code to C for development builds:
#    USE_CYTHON=1 python setup.py build_ext --inplace
# and for building source distribution packages:
#    USE_CYTHON=1 python setup.py sdist
# and *not* when installing a distribution package.
# See http://docs.cython.org/en/latest/src/userguide/source_files_and_compilation.html#distributing-cython-modules
USE_CYTHON = 'USE_CYTHON' in os.environ

ext = '.pyx' if USE_CYTHON else '.c'

cython_extensions = [
    Extension('arrow.arrowhead',
              sources=['arrow/mersenne.c', 'arrow/obsidian.c', 'arrow/arrowhead'+ext,],
              include_dirs=['arrow'] + numpy.distutils.misc_util.get_numpy_include_dirs(),
              define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
              )]

if USE_CYTHON:
    from Cython.Build import cythonize
    cython_extensions = cythonize(
        cython_extensions,
        include_path=['arrow'],
        annotate=True,  # to get an HTML code listing
    )

setup(
    name='stochastic-arrow',
    version='0.5.0',
    packages=['arrow'],
    author='Ryan Spangler, John Mason, Jerry Morrison, Chris Skalnik, Travis Ahn-Horst',
    author_email='ryan.spangler@gmail.com',
    url='https://github.com/CovertLab/arrow',
    license='MIT',
    include_dirs=include,
    ext_modules=cython_extensions,
    long_description=long_description,
    long_description_content_type='text/markdown',
    requires=['numpy (>=1.14)', 'six'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
    ])
